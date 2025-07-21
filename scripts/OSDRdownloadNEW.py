#!/usr/bin/env python3

# USAGE EXAMPLES 
# 
# python OSDRdownload.py \
#   --ids GLDS-168 GLDS-245 OSD-871 \
#   --dest data/nasa_studies
#
# python OSDRdownload.py \
#   --filter project_type=spaceflight \
#   --output results/spaceflight_studies.csv \
#   --dest data/spaceflight
# 
# python OSDRdownload.py \
#   --ids GLDS-214 OSD-883 \
#   --dest data/redo \
#   --retry-failed

"""
NASA GLDS + OSDR Downloader (refactored)
- DRY download logic via a generic download_file()
- Robust retries via urllib3 Retry
- Pathlib for filesystem operations
- Type hints and PEP8 compliance
"""

import argparse
import csv
import os
import logging
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional

import requests
from requests.adapters import HTTPAdapter
from tqdm import tqdm
from urllib3.util.retry import Retry
from urllib.parse import urlencode

summary_lock = threading.Lock()

# Constants
CHUNK_SIZE = 8_192
TIMEOUT = 60
MAX_RETRIES = 3
BACKOFF_FACTOR = 1

# Configure root logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
log = logging.getLogger(__name__)


def setup_session() -> requests.Session:
    """
    Create a requests.Session with a Retry-backed HTTPAdapter.
    """
    session = requests.Session()
    retry = Retry(
        total=MAX_RETRIES,
        backoff_factor=BACKOFF_FACTOR,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"],
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=100, pool_maxsize=100)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def download_file(
    session: requests.Session,
    url: str,
    dest: Path,
    summary: Dict[str, int],
    key: str,
    summary_field: str,
    desc: Optional[str] = None,
    retry_failed: bool = False
) -> None:
    """
    Download a file with retries, size-check, and progress bar.

    Parameters
    ----------
    session
        Configured requests.Session
    url
        Full URL to download
    dest
        Local Path where file will be saved
    summary
        Shared summary dict to increment upon success/skip
    key
        Top‐level key in `summary` for this study/accession
    summary_field
        Field name under summary[key] to increment
    desc
        Description for progress bar (e.g., "GLDS-123 metadata")
    retry_failed
        If True, re-download incomplete files
    """
    dest.parent.mkdir(parents=True, exist_ok=True)

    # 1. HEAD to get remote size
    try:
        head = session.head(url, timeout=TIMEOUT)
        head.raise_for_status()
        remote_size = int(head.headers.get("Content-Length", 0))
    except Exception:
        remote_size = None

    # 2. Skip if already exists and sizes match
    if dest.exists() and remote_size is not None:
        local_size = dest.stat().st_size
        if local_size == remote_size:
            log.info("⏩ %s already exists (%s)", dest.name, desc or "")
            summary[key][summary_field] += 1
            return
        if not retry_failed:
            log.warning("⚠️ %s incomplete (%d/%d bytes), skipping", dest.name, local_size, remote_size)
            return
        log.info("🔁 Retrying %s due to --retry-failed", dest.name)
        dest.unlink(missing_ok=True)

    # 3. Stream download with tqdm
    for attempt in range(MAX_RETRIES):
        try:
            with session.get(url, stream=True, timeout=TIMEOUT) as r:
                r.raise_for_status()
                total = remote_size or 0
                with tqdm(
                    total=total,
                    unit="B",
                    unit_scale=True,
                    desc=desc or dest.name,
                ) as pbar, dest.open("wb") as f:
                    for chunk in r.iter_content(CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                            pbar.update(len(chunk))
            with summary_lock:
                summary[key][summary_field] += 1
            log.info("✔ Downloaded %s", dest.name)
            return
        except Exception as e:
            log.warning("✘ %s failed (attempt %d): %s", dest.name, attempt + 1, e)
            time.sleep(2 ** attempt)
    log.error("❌ Gave up on %s after %d attempts", dest.name, MAX_RETRIES)


def build_search_url(base_url: str, params: Dict[str, str]) -> str:
    """Construct a URL with encoded query parameters."""
    return f"{base_url}?{urlencode(params, doseq=True)}"


def fetch_all_results(
    session: requests.Session,
    base_url: str,
    params: Dict[str, str],
    page_size: int = 100
) -> List[Dict]:
    """
    Fetch all paginated JSON results from an API endpoint.
    """
    params = params.copy()
    params.update({"size": page_size})
    page = 0
    all_results: List[Dict] = []

    while True:
        params["page"] = page
        url = build_search_url(base_url, params)
        r = session.get(url, timeout=TIMEOUT)
        r.raise_for_status()
        data = r.json().get("content", [])
        if not data:
            break
        all_results.extend(data)
        page += 1

    return all_results


def extract_metadata(results: List[Dict]) -> List[Dict[str, str]]:
    """
    Extract accession, title, and description from API results.
    """
    out: List[Dict[str, str]] = []
    for entry in results:
        acc = entry.get("study_accession")
        if acc:
            out.append({
                "accession": acc,
                "title": entry.get("title", "").strip(),
                "description": entry.get("description", "").strip(),
            })
    return out


def save_metadata(metadata: List[Dict[str, str]], output_path: Path) -> None:
    """
    Save study metadata to a CSV file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["accession", "title", "description"])
        writer.writeheader()
        writer.writerows(metadata)


def download_glds_studies(
    accessions: List[str],
    output_dir: Path,
    session: requests.Session,
    summary: Dict[str, Dict[str, int]],
    retry_failed: bool = False
) -> None:
    base_url = "https://genelab-data.ndc.nasa.gov/genelab/static/data/"

    def task(acc: str):
        summary[acc] = {"metadata": 0, "processed": 0}
        study_dir = output_dir / acc
        study_url = f"{base_url}{acc}/"

        # Metadata ZIP
        download_file(
            session,
            f"{study_url}{acc}_metadata.zip",
            study_dir / f"{acc}_metadata.zip",
            summary, acc, "metadata",
            desc=f"{acc} metadata",
            retry_failed=retry_failed
        )

        # Processed files
        for ext in [".txt", ".csv", ".tsv"]:
            download_file(
                session,
                f"{study_url}processed/{acc}_processed{ext}",
                study_dir / f"{acc}_processed{ext}",
                summary, acc, "processed",
                desc=f"{acc} processed{ext}",
                retry_failed=retry_failed
            )

    with ThreadPoolExecutor(max_workers=8) as pool:
        for _ in tqdm(pool.map(lambda a: task(a), accessions),
                      total=len(accessions),
                      desc="GLDS downloads"):
            pass


def download_osd_studies(
    accessions: List[str],
    output_dir: Path,
    session: requests.Session,
    summary: Dict[str, Dict[str, int]],
    retry_failed: bool = False
) -> None:
    base_api = "https://osdr.nasa.gov/osdr/data/osd/files/"

    def task(osd_id: str):
        info = session.get(f"{base_api}{osd_id}", timeout=TIMEOUT)
        info.raise_for_status()
        data = info.json()
        study_key = f"OSD-{osd_id}"
        files = data.get("studies", {}).get(study_key, {}).get("study_files", [])
        summary[study_key] = {"files": 0}

        for f in files:
            url = f"https://osdr.nasa.gov{f['remote_url']}"
            dest = output_dir / study_key / f["file_name"]
            download_file(
                session, url, dest, summary,
                study_key, "files",
                desc=f"{study_key} → {dest.name}",
                retry_failed=retry_failed
            )

    with ThreadPoolExecutor(max_workers=8) as pool:
        for _ in tqdm(pool.map(lambda o: task(o), accessions),
                      total=len(accessions),
                      desc="OSDR downloads"):
            pass

from pathlib import Path
from typing import List, Dict

def choose_studies(
    metadata: List[Dict[str, str]],
    csv_path: Path,
    dest_dir: Path
) -> List[str]:
    """
    Prompt user to inspect metadata CSV and select studies to download.

    Returns a list of selected accession strings.
    """
    # Offer to open the CSV
    ans = input(f"\nOpen '{csv_path}' now for inspection? [y/N]: ").strip().lower()
    if ans == "y":
        import webbrowser
        webbrowser.open(str(csv_path))
        input("Press Enter to continue...")

    # List studies
    print("\nStudy List:")
    for i, entry in enumerate(metadata):
        title = entry["title"]
        print(f"[{i}] {entry['accession']}: {title[:60]}{'...' if len(title)>60 else ''}")

    # Menu
    print("\nOptions:\n[A] Download all\n[S] Select specific\n[N] Skip download")
    choice = input("Choose [A/S/N]: ").strip().lower()

    if choice == "n":
        print("❌ Skipping download.")
        return []
    if choice == "a":
        selected = metadata
    elif choice == "s":
        idxs = [
            int(i) for i in input("Enter indices (e.g., 0,2): ")
            .split(",") if i.strip().isdigit()
        ]
        selected = [metadata[i] for i in idxs if 0 <= i < len(metadata)]
    else:
        print("⚠ Invalid choice. Aborting.")
        return []

    confirm = input("Proceed to download these datasets? [y/N]: ").strip().lower()
    if confirm != "y":
        print("❌ Aborted download.")
        return []

    # Return just the accession strings
    return [entry["accession"] for entry in selected]


def main():
    parser = argparse.ArgumentParser(description="Download NASA GeneLab and OSDR datasets")
    parser.add_argument("-i", "--ids", nargs="+", help="GLDS-### or OSD-### IDs")
    parser.add_argument("-f", "--filter", action="append", help="KEY=VALUE filters")
    parser.add_argument("-o", "--output", help="CSV path for metadata")
    parser.add_argument("-d", "--dest", default="data/OSDR", help="Destination dir")
    parser.add_argument("--retry-failed", action="store_true", help="Re-download incomplete")

    args = parser.parse_args()
    if not args.ids and not args.filter:
        parser.error("Specify --ids or --filter")
    if args.filter and not args.output:
        parser.error("--output required with --filter")

    session = setup_session()
    summary: Dict[str, Dict[str, int]] = {}
    dest = Path(args.dest)

    if args.ids:
        glds = [x.upper().split("GLDS-")[-1] for x in args.ids if x.upper().startswith("GLDS-")]
        osd  = [x.upper().split("OSD-")[-1]  for x in args.ids if x.upper().startswith("OSD-")]
        dest.mkdir(parents=True, exist_ok=True)

        if glds:
            download_glds_studies(glds, dest, session, summary, retry_failed=args.retry_failed)
        if osd:
            download_osd_studies(osd, dest, session, summary, retry_failed=args.retry_failed)

    if args.filter:
        filters = dict(kv.split("=", 1) for kv in args.filter)
        filters["api_key"] = os.getenv("NASA_API_KEY", "DEMO_KEY")  # get API key from environment, falling back to DEMO_KEY
        base_search = "https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/"
        results = fetch_all_results(session, base_search, filters)
        metadata = extract_metadata(results)
        save_metadata(metadata, Path(args.output))
        print(f"✔ Saved {len(metadata)} studies to {args.output}")

        # selection tool
        selected_ids = choose_studies(metadata, Path(args.output), dest)
        if selected_ids:
            download_glds_studies(selected_ids, dest, session, summary)

    # Summary
    print("\nDownload Summary:")
    for sid, stats in summary.items():
        counts = ", ".join(f"{k}: {v}" for k, v in stats.items())
        print(f"{sid} → {counts}")


if __name__ == "__main__":
    main()

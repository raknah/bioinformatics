#!/usr/bin/env python3

"""
NASA GLDS + OSDR Downloader

Supports:
- Direct download via --ids [GLDS-###, OSD-###]
- OSDR metadata discovery via --filter (future use)
"""

import os, argparse, logging, urllib.parse, requests, csv, time
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger(__name__)

# -------------------- API Utilities --------------------

# 1. Build a search query based on provided filters

def build_search_url(base_url: str, params: dict) -> str:
    """
    Construct a search URL with query parameters.

    base_url: str
        BASE OSDR endpoint URL
    params: dict
        Dictionary of filter key-values


    Returns
    -------
    full_url: str
        Encoded search URL with parameters

    """
    qs = urllib.parse.urlencode(params, doseq=True)
    return f"{base_url}?{qs}"

# 2. Fetch all results available given search query

def fetch_all_results(session: requests.Session, base_url: str, params: dict,
                      page_size: int=100) -> list:
    """
    Fetch all paginated JSON results from the OSDR study search API.

    Parameters
    ----------
    session : requests.Session
        Reusable HTTP session with retry settings.
    base_url : str
        Base URL of the OSDR API (e.g., https://osdr.nasa.gov/osdr/data/osd/query/studies/)
    params : dict
        Dictionary of filter key-value pairs (e.g., {"project_type": "spaceflight"}).
        May also include `api_key`.
    page_size : int, optional
        Number of results per page (default = 100).

    Returns
    -------
    results : list of dict
        Combined list of all study entries returned by the API.
    """
    params.update({"size": page_size, "api_key": params.get("api_key")})
    page = 0
    results = []
    while True:
        params["page"] = page
        url = build_search_url(base_url, params)
        r = session.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()
        content = data.get("content", [])
        if not content:
            break
        results.extend(content)
        page += 1
    return results

# 3. Extract accessions and metadata that correspond to the search criteria

def extract_metadata(results: list) -> list:
    """
    Extract GLDS study accession strings from OSDR API results
    
    Parameters
    ----------
    results: list
        List of study metadata dictionaries
        
    Returns
    -------
    metadata: list of dict
        Each item: {"accession": str, "title": str, "description": str}
    """

    metadata = []

    for entry in results:
        accession = entry.get("study_accession")

        if accession:
            metadata.append({
                "accession": accession,
                "title": entry.get("title", "").strip(),
                "description": entry.get("description", "").strip()
            })

    return metadata

# 4. Save extracted accessions to a file

def save_metadata(metadata: list, output_path: str) -> None:
    """
    Save study metadata to a CSV file.

    Parameters
    ----------
    metadata : list of dict
        Each dict contains keys: "accession", "title", "description".
    output_path : str
        Destination CSV file path (e.g. studies.csv)
    """

    fieldnames = ["accession", "title", "description"]

    with open(output_path, "w", newline = "", encoding = "utf-8") as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()
        writer.writerows(metadata)


# -------------------- GLDS Downloader --------------------

def download_glds_studies(accessions: list[str], output_dir: str, session: requests.Session, summary: dict, retry_failed: bool = False) -> None:
    """
    Download metadata and processed files for all GLDS studies in parallel,
    showing one overall progress bar.
    """
    base_url = "https://genelab-data.ndc.nasa.gov/genelab/static/data/"

    def task(accession: str):
        summary[accession] = {"metadata": False, "processed": 0}
        study_url = f"{base_url}{accession}/"
        local_dir = os.path.join(output_dir, accession)
        os.makedirs(local_dir, exist_ok=True)

        # ── Metadata ─────────────────────────────────────────────
        meta_url = f"{study_url}{accession}_metadata.zip"
        meta_path = os.path.join(local_dir, f"{accession}_metadata.zip")

        for attempt in range(3):
            try:
                r = session.get(meta_url, timeout=60)
                if r.status_code == 200:
                    remote_size = int(r.headers.get("Content-Length", 0))
                    if os.path.exists(meta_path):
                        local_size = os.path.getsize(meta_path)
                        if local_size == remote_size:
                            log.info("⏩ %s: metadata already downloaded", accession)
                            summary[accession]["metadata"] = True
                            break
                        elif not retry_failed:
                            log.warning("⏭ %s: metadata incomplete (%d/%d bytes), skipping (use --retry-failed)", accession, local_size, remote_size)
                            break
                        else:
                            log.info("🔁 %s: metadata incomplete, retrying due to --retry-failed", accession)
                            try:
                                os.remove(meta_path)
                            except Exception as e:
                                log.warning("⚠️ Could not delete %s: %s", meta_path, e)


                    with open(meta_path, "wb") as f:
                        f.write(r.content)
                    summary[accession]["metadata"] = True
                    log.info("✔ %s: metadata downloaded", accession)
                else:
                    log.warning("⚠ %s: metadata not found", accession)
                break
            except requests.RequestException as e:
                log.warning("⏳ %s: metadata error (attempt %d) → %s", accession, attempt + 1, e)
                time.sleep(2 ** attempt)
        else:
            log.error("❌ %s: failed to download metadata after 3 attempts", accession)

        # ── Processed ─────────────────────────────────────────────
        for ext in [".txt", ".csv", ".tsv"]:
            proc_url = f"{study_url}processed/{accession}_processed{ext}"
            proc_path = os.path.join(local_dir, f"{accession}_processed{ext}")

            for attempt in range(3):
                try:
                    r = session.get(proc_url, timeout=60)
                    if r.status_code == 200:
                        remote_size = int(r.headers.get("Content-Length", 0))
                        if os.path.exists(proc_path):
                            local_size = os.path.getsize(proc_path)
                            if local_size == remote_size:
                                log.info("⏩ %s: processed%s already downloaded", accession, ext)
                                summary[accession]["processed"] += 1
                                break
                            elif not retry_failed:
                                log.warning("⏭ %s: processed%s incomplete (%d/%d bytes), skipping (use --retry-failed)", accession, ext, local_size, remote_size)
                                break
                            else:
                                log.info("🔁 %s: processed%s incomplete, retrying due to --retry-failed", accession, ext)
                                try:
                                    os.remove(proc_path)
                                except Exception as e:
                                    log.warning("⚠️ Could not delete %s: %s", proc_path, e)

                        with open(proc_path, "wb") as f:
                            f.write(r.content)
                        summary[accession]["processed"] += 1
                        log.info("✔ %s: processed%s downloaded", accession, ext)
                        break
                except requests.RequestException as e:
                    log.warning("⏳ %s: processed%s error (attempt %d) → %s", accession, ext, attempt + 1, e)
                    time.sleep(2 ** attempt)
            else:
                log.error("❌ %s: failed to download processed%s after 3 attempts", accession, ext)

# -------------------- OSDR Downloader --------------------

def download_osd_studies(accessions: list, output_dir: str, session: requests.Session, summary: dict, retry_failed=False):
    """
    Download all files for each OSDR study using the OSDR Data File API in parallel.
    """
    base_api = "https://osdr.nasa.gov/osdr/data/osd/files/"

    def task(osd_id: str):
        try:
            url = f"{base_api}{osd_id}"
            r = session.get(url, timeout=30)
            r.raise_for_status()
            info = r.json()
            study_key = f"OSD-{osd_id}"
            files = info.get("studies", {}).get(study_key, {}).get("study_files", [])
            osd_dir = os.path.join(output_dir, study_key)
            os.makedirs(osd_dir, exist_ok=True)
            summary[study_key] = {"files": 0}

            for f in files:
                fname = f["file_name"]
                remote = f["remote_url"]
                full_url = f"https://osdr.nasa.gov{remote}"
                local_path = os.path.join(osd_dir, fname)

                if os.path.exists(local_path):
                    local_size = os.path.getsize(local_path)
                    remote_size = int(rr.headers.get("Content-Length", 0))
                    
                    if local_size == remote_size:
                        log.info("⏩ Skipping %s (already complete)", fname)
                        summary[study_key]["files"] += 1
                        return  # skip this file entirely
                    
                    elif not retry_failed:
                        log.warning("⏭ %s exists but incomplete (%.2f%%); skipping (use --retry-failed to re-download)",
                                    fname, 100 * local_size / remote_size)
                        return
                    else:
                        log.info("🔁 %s exists but incomplete; re-downloading due to --retry-failed", fname)
                        try:
                            os.remove(local_path)
                        except Exception as e:
                            log.warning("⚠ Could not delete %s: %s", local_path, e)

                for attempt in range(3):
                    try:
                        rr = session.get(full_url, stream=True, timeout=60)
                        rr.raise_for_status()
                        content_len = int(rr.headers.get("Content-Length", 0))

                        with tqdm(total=content_len, unit='B', unit_scale=True, desc=f"{study_key} → {fname}") as progress:
                            with open(local_path, "wb") as out:
                                for chunk in rr.iter_content(chunk_size=8192):
                                    if chunk:
                                        out.write(chunk)
                                        progress.update(len(chunk))

                        summary[study_key]["files"] += 1
                        log.info("✔ OSD-%s → %s", osd_id, fname)
                        break  # success, exit retry loop

                    except requests.RequestException as e:
                        log.warning("✘ OSD-%s file %s → %s (attempt %d)", osd_id, fname, str(e), attempt + 1)
                        time.sleep(2 ** attempt)  # backoff: 1s → 2s → 4s

                else:
                    log.error("❌ Gave up on OSD-%s file %s after 3 attempts", osd_id, fname)

        except Exception as e:
            log.error("❌ Failed to fetch OSD-%s: %s", osd_id, str(e))

    # Run all tasks in parallel with one overall bar
    osd_ids = accessions
    with ThreadPoolExecutor(max_workers=16) as pool:
        futures = [pool.submit(task, oid) for oid in osd_ids]
        for _ in tqdm(as_completed(futures),
                      total=len(futures),
                      desc="Downloading OSD-studies"):
            pass

# MAIN

def main():
    """
    Entry point for NASA data downloader.
    Supports:
    - --ids: Direct GLDS/OSD download mode (recommended)
    - --filter + --output: OSDR metadata discovery (API currently down)
    """
    # ── 1. CLI ─────────────────────────────────────────────────────────────
    parser = argparse.ArgumentParser(description="Download NASA GeneLab and OSDR datasets")
    parser.add_argument(
        "--ids", "-i", nargs="+", metavar="ID",
        help="List of GLDS or OSD study IDs, e.g. GLDS-168 OSD-871"
    )
    parser.add_argument(
        "--filter", "-f", action="append", metavar="KEY=VALUE",
        help="Query filters for OSDR metadata (e.g. project_type=spaceflight)"
    )
    parser.add_argument(
        "--output", "-o",
        help="Output CSV path for filtered metadata (required with --filter)"
    )
    parser.add_argument(
        "--dest", "-d", default="data/OSDR",
        help="Download destination directory"
    )
    parser.add_argument(
    "--retry-failed",
    action="store_true",
    help="Force re-download of incomplete or partially downloaded files"
    )

    args = parser.parse_args()

    if not args.ids and not args.filter:
        parser.error("You must specify either --ids or --filter.")

    if args.filter and not args.output:
        parser.error("--output is required when using --filter.")

    # ── 2. HTTP session ────────────────────────────────────────────────────
    session = requests.Session()
    session.mount("https://", requests.adapters.HTTPAdapter(pool_connections=100, pool_maxsize=100, max_retries=3))

    summary = {}

    # ── 3. Direct ID Mode ──────────────────────────────────────────────────
    if args.ids:
        glds_ids = [x.upper() for x in args.ids if x.upper().startswith("GLDS-")]
        osd_ids = [x.replace("OSD-", "") for x in args.ids if x.upper().startswith("OSD-")]

        os.makedirs(args.dest, exist_ok=True)

        if glds_ids:
            download_glds_studies(glds_ids, args.dest, session, summary, retry_failed=args.retry_failed)

        if osd_ids:
            download_osd_studies(osd_ids, args.dest, session, summary, retry_failed=args.retry_failed)


    # ── 4. OSDR Metadata Search Mode ───────────────────────────────────────
    if args.filter:
        filters = {}
        for kv in args.filter:
            key, val = kv.split("=", 1)
            filters[key.strip()] = val.strip()
        filters["api_key"] = os.getenv("NASA_API_KEY", "DEMO_KEY")

        base_search = "https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/"
        results = fetch_all_results(session, base_search, filters)
        metadata = extract_metadata(results)

        os.makedirs(os.path.dirname(args.output), exist_ok=True)
        save_metadata(metadata, args.output)
        print(f"✔ Saved {len(metadata)} studies to {args.output}")

        # Interactive selection
        open_now = input(f"\nOpen '{args.output}' now for inspection? [y/N]: ").strip().lower()
        if open_now == "y":
            import webbrowser
            webbrowser.open(args.output)
            input("Press Enter to continue...")

        print("\nStudy List:")
        for i, entry in enumerate(metadata):
            print(f"[{i}] {entry['accession']}: {entry['title'][:60]}...")

        print("\nOptions:")
        print("[A] Download all")
        print("[S] Select specific studies")
        print("[N] Skip download")
        choice = input("Choose [A/S/N]: ").strip().lower()

        if choice == "n":
            print("❌ Skipping download.")
            return
        elif choice == "a":
            selected = metadata
        elif choice == "s":
            idxs = input("Enter comma-separated indices (e.g., 0,2,5): ")
            idxs = [int(i.strip()) for i in idxs.split(",") if i.strip().isdigit()]
            selected = [metadata[i] for i in idxs if 0 <= i < len(metadata)]
        else:
            print("⚠ Invalid choice. Aborting.")
            return

        confirm = input("Proceed to download these datasets? [y/n]: ").strip().lower()
        if confirm != "y":
            print("❌ Aborted download.")
            return

        os.makedirs(args.dest, exist_ok=True)
        accessions = [entry["accession"] for entry in selected]
        download_glds_studies(accessions, args.dest, session, summary)

    # ── 5. Summary ─────────────────────────────────────────────────────────
    print("\n────────────── Download Summary ──────────────")
    for study_id, stats in summary.items():
        if study_id.startswith("GLDS"):
            meta = "✔" if stats.get("metadata") else "✘"
            proc = stats.get("processed", 0)
            print(f"{study_id}: metadata {meta}, processed files: {proc}")
        elif study_id.startswith("OSD"):
            count = stats.get("files", 0)
            print(f"{study_id}: {count} files downloaded")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# EXAMPLE USES
#
# python GEOdownload.py \
#   --csv GSE_shortlist.csv \
#   --project space \
#   --threads 3

# python GEOdownload.py \
#   --list shortlist.txt \
#   --project space

# python GEOdownload.py \
#   --gse GSE63952 GSE64375 GSE12345 \
#   --project space \
#   --threads 3

"""
GEOdownload.py
Download GEO Series datasets (GSE IDs) in batch using GEOparse,
with progress bars for GSE-level and GSM-level completion.
"""

import os
import sys
import argparse
import pandas as pd
from GEOparse import get_GEO
from tqdm import tqdm
from multiprocessing import Pool

def download_gse(args):
    gse_id, project_name, base = args
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    out_dir = os.path.join(root, base, project_name, "data", "GEO", gse_id)
    os.makedirs(out_dir, exist_ok=True)

    try:
        gse = get_GEO(geo=gse_id, destdir=out_dir, annotate_gpl=True)
    except Exception as e:
        print(f"[ERROR] {gse_id}: {e}", file=sys.stderr)
        return

    # save expression files
    for gsm, gsmobj in tqdm(
        gse.gsms.items(),
        desc=f"[{gse_id}] expr",
        leave=False,
        position=1
    ):
        try:
            out_file = os.path.join(out_dir, f"{gsm}_expression.csv")
            gsmobj.table.to_csv(out_file, index=False)
        except Exception as e:
            print(f"[ERROR] {gsm}: {e}", file=sys.stderr)

    # save metadata
    rows = []
    for gsm, gsmobj in tqdm(
        gse.gsms.items(),
        desc=f"[{gse_id}] meta",
        leave=False,
        position=2
    ):
        rows.append({
            "GSM_ID": gsm,
            "Title": gsmobj.metadata.get("title", [""])[0],
            "Source": gsmobj.metadata.get("source_name_ch1", [""])[0],
            "Characteristics": "; ".join(gsmobj.metadata.get("characteristics_ch1", []))
        })
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, "metadata.csv"), index=False)
    print(f"[DONE] {gse_id}")

def parse_args():
    p = argparse.ArgumentParser(description="Download GEO Series in batch.")
    p.add_argument("-p", "--project", required=True,
                   help="Project folder name (e.g., 'space')")
    p.add_argument("-b", "--base", default="projects",
                   help="Base path where projects/ folder lives")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("-g", "--gse", nargs="+",
                       help="One or more GSE IDs to download")
    group.add_argument("-l", "--list",
                       help="Path to .txt file with GSE IDs (one per line)")
    group.add_argument("-c", "--csv",
                       help="Path to CSV (e.g. from GEOselector) with GSE IDs")
    p.add_argument("--id-col", default="GSE_ID",
                   help="Column name in CSV that contains GSE IDs")
    p.add_argument("-t", "--threads", type=int, default=3,
                   help="Number of concurrent downloads (default: 3)")
    return p.parse_args()

def main():
    args = parse_args()

    # collect GSE IDs
    if args.csv:
        df = pd.read_csv(args.csv)
        if args.id_col not in df.columns:
            print(f"[ERROR] Column '{args.id_col}' not found in CSV.", file=sys.stderr)
            sys.exit(1)
        ids = df[args.id_col].dropna().astype(str).tolist()
    elif args.list:
        if args.list == "-":
            ids = [line.strip() for line in sys.stdin if line.strip()]
        else:
            if not os.path.isfile(args.list):
                print(f"[ERROR] List file not found: {args.list}", file=sys.stderr)
                sys.exit(1)
            with open(args.list) as f:
                ids = [line.strip() for line in f if line.strip()]
    else:
        ids = args.gse

    if not ids:
        print("[ERROR] No GSE IDs provided.", file=sys.stderr)
        sys.exit(1)

    tasks = [(gse, args.project, args.base) for gse in ids]

    print(f"📦 Downloading {len(tasks)} datasets into: {args.base}/{args.project}/data/GEO")

    if args.threads > 1:
        with Pool(processes=args.threads) as pool:
            for _ in tqdm(pool.imap_unordered(download_gse, tasks),
                          total=len(tasks),
                          desc="Overall GSEs"):
                pass
    else:
        for t in tqdm(tasks, desc="Overall GSEs"):
            download_gse(t)

if __name__ == "__main__":
    main()


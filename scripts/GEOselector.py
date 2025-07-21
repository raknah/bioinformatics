#!/usr/bin/env python3

# EXAMPLE USAGE 1
#
# python GEOselector.py \
#   --email fomo@cbigspace.org \
#   --query '("space radiation" OR "ionizing radiation" OR "Gy") AND "Homo sapiens"[Organism] AND gse[Entry Type]' \
#   --min-total 6 \
#   --max-ratio 3.0 \
#   --control-keywords control sham unexposed \
#   --treated-keywords irradiated radiation gy cosmic \
#   --output GSE_shortlist.csv



#!/usr/bin/env python3
"""
GEOselector.py
Shortlist GEO Series (GSE) datasets for human space-radiation transcriptomic analysis,
using group keyword detection, sample balance, and matrix availability filters.
"""

import argparse
import time
import pandas as pd
from Bio import Entrez

# ──────── Helper Functions ─────────────────────────────────────────────────────

def fetch_gse_ids(query: str, retmax: int) -> list[str]:
    handle = Entrez.esearch(
        db="gds",
        term=query,
        retmax=retmax,
        retmode="xml"
    )
    result = Entrez.read(handle)
    handle.close()
    return result["IdList"]

def fetch_gse_summary(gse_id: str) -> dict:
    handle = Entrez.esummary(db="gds", id=gse_id, retmode="xml")
    recs = Entrez.read(handle)
    handle.close()
    return recs[0]

# ──────── Main CLI ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Shortlist GEO Series datasets for DE analysis using keyword group detection and balance filters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--email", "-e", required=True, help="Your email for Entrez API")
    parser.add_argument("--query", required=True, help="Entrez search query for GSEs")
    parser.add_argument("--min-total", type=int, default=6, help="Minimum total sample count")
    parser.add_argument("--max-ratio", type=float, default=3.0, help="Maximum allowed group imbalance (e.g. 3.0 = 3:1)")
    parser.add_argument("--retmax", type=int, default=2000, help="Maximum number of GSE IDs to retrieve")
    parser.add_argument("--output", "-o", default="GSE_shortlist.csv", help="Output path for filtered CSV")

    parser.add_argument("--control-keywords", nargs="+", default=["control", "sham"],
                        help="Terms identifying control group samples")
    parser.add_argument("--treated-keywords", nargs="+", default=["irradiated", "radiation", "gy"],
                        help="Terms identifying treated/irradiated samples")

    args = parser.parse_args()
    Entrez.email = args.email

    print(f"🔍 Querying GEO: {args.query}")
    gse_ids = fetch_gse_ids(args.query, args.retmax)
    print(f"⚡ Retrieved {len(gse_ids)} GSE IDs")

    records = []

    for gse in gse_ids:
        try:
            summ = fetch_gse_summary(gse)
            n_samples = int(summ.get("n_samples", 0))
            taxon = summ.get("taxon", "")
            title = summ.get("title", "")
            summary = summ.get("summary", "")
            supp = summ.get("supplementaryfiles", "")

            if taxon != "Homo sapiens":
                continue
            if n_samples < args.min_total:
                continue

            # Combine text
            txt = (title + " " + summary).lower()

            # Group detection
            n_control = sum(txt.count(word) for word in args.control_keywords)
            n_treated = sum(txt.count(word) for word in args.treated_keywords)

            if min(n_control, n_treated) == 0:
                continue

            ratio = max(n_control, n_treated) / min(n_control, n_treated)
            if ratio > args.max_ratio:
                continue

            if not any(ext in supp.lower() for ext in [".csv", ".txt", ".soft.gz"]):
                continue

            records.append({
                "GSE_ID": gse,
                "Samples": n_samples,
                "Control_Est": n_control,
                "Treated_Est": n_treated,
                "Ratio": round(ratio, 2),
                "Title": title,
                "Summary": summary,
                "SuppFiles": supp
            })

            time.sleep(0.34)  # Stay under NCBI limit

        except Exception as err:
            print(f"⚠️  Error on {gse}: {err}")

    df = pd.DataFrame(records)
    df = df.sort_values("Samples", ascending=False)
    df.to_csv(args.output, index=False)
    print(f"\n✔︎ Saved {len(df)} datasets to {args.output}")

if __name__ == "__main__":
    main()

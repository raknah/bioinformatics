#!/usr/bin/env python3
"""
Extract top 10 DEGs from GEO2R *.top.table.tsv files and create a single
Enrichr-compatible gene list. Grouped by GSE ID and dose.

Usage:
    python top10_enrichr.py --gse GSE63952 --dir "/path/to/data"
"""

import os
import glob
import re
import argparse
import pandas as pd

def get_sort_column(df):
    for col in ("adj.P.Val", "P.Value", "pvalue", "p_value"):
        if col in df.columns:
            return col
    return None

def main():
    parser = argparse.ArgumentParser(description="Extract top 10 DEGs and make Enrichr gene list.")
    parser.add_argument("--gse", required=True, help="GSE ID (e.g., GSE63952)")
    parser.add_argument("--dir", required=True, help="Directory containing *.top.table.tsv files")
    args = parser.parse_args()

    gse_id = args.gse
    base_dir = args.dir

    input_dir = os.path.join(base_dir, gse_id)
    output_dir = os.path.join(base_dir, gse_id)
    enrichr_file = os.path.join(base_dir, gse_id, is f"{gse_id}_all_top10_genes.txt")

    os.makedirs(output_dir, exist_ok=True)
    print(f"Looking in: {os.path.abspath(input_dir)}")

    pattern = os.path.join(input_dir, "*.top.table.tsv")
    files = glob.glob(pattern)
    print(f"Found {len(files)} files.\n")

    gene_set = set()

    for filepath in files:
        basename = os.path.basename(filepath)
        m = re.match(rf"^{gse_id}(low|med|high)\.top\.table\.tsv$", basename, re.IGNORECASE)
        if not m:
            print(f"Skipping: {basename} (filename doesn't match expected pattern)")
            continue
        dose = m.group(1)
        print(f"Processing {gse_id} [{dose}]")

        df = pd.read_csv(filepath, sep="\t")
        sort_col = get_sort_column(df)
        if not sort_col:
            print(f"  ✘ No p-value column found in {basename}")
            continue

        top10 = df.sort_values(by=sort_col, ascending=True).head(10)
        out_name = f"{gse_id}{dose}.top10.table.tsv"
        top_path = os.path.join(output_dir, out_name)
        top10.to_csv(top_path, sep="\t", index=False)
        print(f"  ✔ Saved top 10 to {top_path}")

        symbol_col = next((col for col in df.columns if col.lower().replace(".", "") == "genesymbol"), None)
        if symbol_col:
            genes = top10[symbol_col].dropna().astype(str).str.strip()
            gene_set.update(genes)

    if gene_set:
        with open(enrichr_file, "w") as f:
            for gene in sorted(gene_set):
                f.write(f"{gene}\n")
        print(f"\n✔ Saved Enrichr gene list: {enrichr_file}")
    else:
        print("\n⚠️ No gene symbols extracted. Nothing written.")

if __name__ == "__main__":
    main()

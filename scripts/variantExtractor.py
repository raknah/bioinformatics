#!/usr/bin/env python3

import os
import pandas as pd
import requests
import subprocess

# === Defaults ===
DISGENET_FILE = "genomes/all_gene_disease_associations.tsv.gz"
DISGENET_URL = "https://www.disgenet.org/static/disgenet_ap1/files/downloads/all_gene_disease_associations.tsv.gz"
DEFAULT_VCF = "genomes/variants/clinvar.vcf.gz"

# === Download DisGeNET if needed ===
if not os.path.exists(DISGENET_FILE):
    print(f"⚠️ DisGeNET file not found.")
    print(f"Please download it manually from https://www.disgenet.org/downloads")
    print(f"And place it at: {DISGENET_FILE}")
    exit(1)

# === 1. Prompt for disease ===
query = input("🧠 Enter disease/condition name (e.g. autism): ").strip().lower()
filename_stub = query.replace(" ", "_")
DEFAULT_BED_DIR = f"projects/{query}"

# === 2. Prompt for BED output location ===
print(f"\n📂 Default BED folder: {DEFAULT_BED_DIR}")
bed_dir = input("   Enter custom folder (or press Enter to use default): ").strip() or DEFAULT_BED_DIR
bedfile = os.path.join(bed_dir, f"{filename_stub}-genes.bed")

# === 3. Read and filter DisGeNET ===
print(f"\n📦 Loading DisGeNET from: {DISGENET_FILE}")
df = pd.read_csv(DISGENET_FILE, sep="\t", compression="gzip")
filtered = df[df["diseaseName"].str.contains(query, case=False)]
genes = sorted(set(filtered["geneSymbol"]))
print(f"✅ Found {len(genes)} genes linked to '{query}'\n")

# === 4. Ensembl coordinate lookup ===
def get_coords(symbol):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
    r = requests.get(url)
    if r.status_code == 200:
        j = r.json()
        return (j["display_name"], j["seq_region_name"], j["start"], j["end"])
    return None

# === 5. Write BED ===
print(f"💾 Writing BED to: {bedfile}")
os.makedirs(os.path.dirname(bedfile), exist_ok=True)

with open(bedfile, "w") as out:
    for gene in genes:
        coords = get_coords(gene)
        if coords:
            name, chrom, start, end = coords
            out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

print("✅ BED file complete.\n")

# === 6. Optional VCF Extraction ===
extract = input("🧬 Extract matching variants from a VCF? [y/n]: ").strip().lower()
if extract in ("y", ""):
    print(f"\n📦 Default VCF path: {DEFAULT_VCF}")
    vcf_path = input("   Enter custom VCF path (or press Enter to use default): ").strip() or DEFAULT_VCF
    out_vcf = os.path.splitext(bedfile)[0] + "-variants.vcf.gz"

    print(f"\n🚀 Extracting variants to: {out_vcf}")
    subprocess.run(["bcftools", "view", "-R", bedfile, vcf_path, "-Oz", "-o", out_vcf])
    subprocess.run(["bcftools", "index", out_vcf])
    print(f"✅ Variant file written and indexed: {out_vcf}")
else:
    print("✔️ BED file ready. You can extract variants later using bcftools.\n")

import os
import sys

from GEOparse import get_GEO



# step 1: parse CLI input
if len(sys.argv) < 3:
    print("Incorrect Usage. Expected: python geoDownload.py <GSE_ID> <project_name>")
    sys.exit(1)

gse_id = sys.argv[1]
project_name = sys.argv[2]

print(f"Arguments received: GSE_ID={gse_id}, project_name={project_name}")

# step 2: destination path
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(project_root, "projects", project_name, "data", "geo", gse_id)
os.makedirs(output_dir, exist_ok = True)

# step 3: download GEO data
print(f"Downloading {gse_id} into {output_dir}...")
gse = get_GEO(geo = gse_id, destdir=output_dir, annotate_gpl=True)

# step 4: save expression data
for gsm_name, gsm in gse.gsms.items():
    gsm.table.to_csv(os.path.join(output_dir, f"{gsm_name}_expression.csv"), index = False)

# step 5: save metadata
meta_path = os.path.join(output_dir, "metadata.tsv")
with open(meta_path, 'w') as f:
    f.write("GSM_ID\tTitle\tSource\tCharacteristics\n")
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [""])[0]
        source = gsm.metadata.get("source_name_ch1", [""])[0]
        characteristics = "; ".join(gsm.metadata.get("characteristics_ch1", []))
        f.write(f"{gsm_name}\t{title}\t{source}\t{characteristics}\n")

print("job done big man")


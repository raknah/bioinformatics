#!/usr/bin/env python3

"""
GEOselector.py
Shortlist GEO Series (GSE) datasets for human space-radiation transcriptomic analysis
by querying Entrez then loading each GSE with GEOparse for reliable sample metadata.

Dependencies:
    pip install biopython GEOparse pandas tqdm
"""

import argparse, time
import pandas as pd
from Bio import Entrez
import GEOparse
from tqdm import tqdm

def fetch_gse_uids(query: str, retmax: int):
    handle = Entrez.esearch(db="gds", term=query, retmax=retmax, retmode="xml")
    rec = Entrez.read(handle); handle.close()
    return rec["IdList"]

def fetch_gse_accession(uid: str) -> str:
    rec = Entrez.esummary(db="gds", id=uid, retmode="xml")
    acc = Entrez.read(rec)[0]["Accession"]
    rec.close()
    return acc

def classify_samples(gse, ctrl_keys, trt_keys, audit):
    ctrl = trt = amb = uncls = 0
    for gsm, sample in gse.gsms.items():
        desc = " ".join(sample.metadata.get("characteristics_ch1", [])).lower()
        ic = any(k.lower() in desc for k in ctrl_keys)
        it = any(k.lower() in desc for k in trt_keys)
        if ic and not it:
            ctrl += 1
        elif it and not ic:
            trt += 1
        elif ic and it:
            amb += 1
        else:
            uncls += 1
        if audit:
            print(f"  {gsm}: ctrl={ic} trt={it} » {desc[:60]}…")
    return ctrl, trt, amb, uncls

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--email","-e", required=True)
    p.add_argument("--query",  required=True,
                   help='Entrez query, e.g. \'("space radiation" OR "ionizing radiation") AND "Homo sapiens"[Organism] AND gse[Entry Type]\'')
    p.add_argument("--min-total",   type=int,   default=6)
    p.add_argument("--max-ratio",   type=float, default=3.0)
    p.add_argument("--retmax",      type=int,   default=2000)
    p.add_argument("--control-keywords", nargs="+", required=True)
    p.add_argument("--treated-keywords", nargs="+", required=True)
    p.add_argument("--verbose",     action="store_true")
    p.add_argument("--audit",       action="store_true")
    p.add_argument("--output","-o", default="GSE_shortlist.csv")
    args = p.parse_args()

    Entrez.email = args.email
    print(f"🔍 Searching GEO DataSets: {args.query}")
    uids = fetch_gse_uids(args.query, args.retmax)
    print(f"⚡ Retrieved {len(uids)} series\n")

    stats = {"passed":0,"missing_group":0,"bad_ratio":0,"too_few":0,"load_fail":0}
    out = []

    bar = tqdm(uids, desc="🔬 Filtering", unit="UID", dynamic_ncols=True)
    for uid in bar:
        try:
            # fetch accession (e.g. "GSE63952")
            acc = fetch_gse_accession(uid)

            # load via GEOparse
            gse = GEOparse.get_GEO(geo=acc, silent=True)

            total = len(gse.gsms)
            if total < args.min_total:
                stats["too_few"] += 1
                reason = f"{total} samples"
                raise StopIteration

            # classify
            c,t,a,u = classify_samples(gse, args.control_keywords, args.treated_keywords, args.audit)

            if min(c,t) == 0:
                stats["missing_group"] += 1
                reason = f"ctrl={c},trt={t}"
                raise StopIteration

            ratio = max(c,t)/min(c,t)
            if ratio > args.max_ratio:
                stats["bad_ratio"] += 1
                reason = f"ratio={ratio:.2f}"
                raise StopIteration

            # passed
            stats["passed"] += 1
            reason = "✅ passed"
            out.append({
                "GSE":    acc,
                "UID":    uid,
                "Samples":total,
                "Control":c,
                "Treated":t,
                "Ratio": round(ratio,2),
                "Title":  gse.metadata.get("title",[""])[0]
            })

        except StopIteration:
            pass
        except Exception as e:
            stats["load_fail"] += 1
            reason = f"ERROR:{e}"

        bar.set_postfix_str(reason)
        if args.verbose:
            passed = stats["passed"]
            skipped = sum(v for k,v in stats.items() if k!="passed")
            print(f"📊 passed={passed} skipped={skipped}: {reason} ({acc})")

        time.sleep(0.5)

    bar.close()

    print("\n📊 Final counts:")
    for k,v in stats.items():
        print(f"  - {k}: {v}")

    if out:
        df = pd.DataFrame(out).sort_values("Samples",ascending=False)
        df.to_csv(args.output,index=False)
        print(f"\n✔︎ Saved {len(out)} series to {args.output}")
    else:
        print("\n❌ No series passed filters.")

if __name__=="__main__":
    main()

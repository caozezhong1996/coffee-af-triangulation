"""Extract the 40 coffee-instrument SNPs from 9 caffeine-family metabolite
GWAS (Chen et al. 2023, European-ancestry metabolomics, N~8,299) on the EBI
GWAS Catalog FTP, using the remote tabix client. Resumable per metabolite.

Exposure alignment: beta is flipped so it is per coffee-INTAKE-INCREASING
allele (ea from spectrum_harmonized_snps.csv, phenocode I9_AF).
"""
import csv
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

from remote_tabix import RemoteTabix, get_header

PANEL = {
    "GCST90200436": ("GCST90200001-GCST90201000", "caffeine"),
    "GCST90199654": ("GCST90199001-GCST90200000", "paraxanthine"),
    "GCST90199644": ("GCST90199001-GCST90200000", "theobromine"),
    "GCST90199647": ("GCST90199001-GCST90200000", "theophylline"),
    "GCST90200285": ("GCST90200001-GCST90201000", "xanthine"),
    "GCST90199754": ("GCST90199001-GCST90200000", "7-methylxanthine"),
    "GCST90199699": ("GCST90199001-GCST90200000", "3-methylxanthine"),
    "GCST90199756": ("GCST90199001-GCST90200000", "1,7-dimethyluric acid"),
    "GCST90199760": ("GCST90199001-GCST90200000", "AAMU"),
}

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{bucket}/{g}/harmonised/{g}.h.tsv.gz"
SRC = "../spectrum_20260901/spectrum_harmonized_snps.csv"
TBI_DIR = "metab_tbi"
OUT_DIR = "metab_raw"

OUT_COLS = ["query_rsid", "chrom", "pos", "ea", "oa", "bx", "bxse",
            "metab_rsid", "effect_allele", "other_allele", "beta_raw",
            "se_raw", "eaf", "p_value", "flip", "allele_match"]


def load_snps():
    snps = {}
    for r in csv.DictReader(open(SRC)):
        if r["phenocode"] != "I9_AF":
            continue
        snps[r["SNP"]] = dict(chrom=r["chr"], pos=int(r["pos"]),
                              ea=r["ea"].upper(), oa=r["oa"].upper(),
                              bx=r["bx"], bxse=r["bxse"])
    return snps


def extract_one(g, bucket, label, snps):
    out_path = os.path.join(OUT_DIR, f"{g}.csv")
    if os.path.exists(out_path):
        print(f"[skip] {label} {g} already done", flush=True)
        return
    url = BASE.format(bucket=bucket, g=g)
    hdr = get_header(url)
    idx = {name: i for i, name in enumerate(hdr)}
    need = ["chromosome", "base_pair_location", "effect_allele",
            "other_allele", "beta", "standard_error",
            "effect_allele_frequency", "p_value", "rsid"]
    missing_cols = [c for c in need if c not in idx]
    if missing_cols:
        print(f"[hdr-fail] {label}: missing {missing_cols}; hdr={hdr}", flush=True)
        return
    rt = RemoteTabix(url, os.path.join(TBI_DIR, f"{g}.tbi"))

    def query_snp(item):
        rsid, s = item
        rec = None
        for refname in (s["chrom"], f"chr{s['chrom']}"):
            if refname not in rt.refs:
                continue
            got = rt.query(refname, s["pos"] - 1, s["pos"])
            if got:
                for f in got:
                    if f[idx["rsid"]] == rsid or int(f[idx["base_pair_location"]]) == s["pos"]:
                        rec = f
                        break
            if rec:
                break
        row = dict(query_rsid=rsid, chrom=s["chrom"], pos=s["pos"],
                   ea=s["ea"], oa=s["oa"], bx=s["bx"], bxse=s["bxse"])
        if rec:
            eff, oth = rec[idx["effect_allele"]].upper(), rec[idx["other_allele"]].upper()
            beta = float(rec[idx["beta"]])
            if eff == s["ea"] and oth == s["oa"]:
                flip, match = 0, "exact"
            elif eff == s["oa"] and oth == s["ea"]:
                flip, match = 1, "swapped"
                beta = -beta
            else:
                flip, match = 0, "mismatch"
            row.update(metab_rsid=rec[idx["rsid"]], effect_allele=eff,
                       other_allele=oth, beta_raw=rec[idx["beta"]],
                       se_raw=rec[idx["standard_error"]],
                       eaf=rec[idx["effect_allele_frequency"]],
                       p_value=rec[idx["p_value"]], flip=flip,
                       allele_match=match)
            # beta_raw column stores ALIGNED beta after flip for direct use
            row["beta_raw"] = beta
        return rsid, row

    rows, missing = [], []
    with ThreadPoolExecutor(max_workers=6) as ex:
        for rsid, row in ex.map(query_snp, snps.items()):
            if "metab_rsid" in row:
                rows.append(row)
            else:
                missing.append(rsid)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"[done] {label} {g}: {len(rows)}/40 found; missing={missing}", flush=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    snps = load_snps()
    print(f"{len(snps)} instrument SNPs", flush=True)
    for g, (bucket, label) in PANEL.items():
        t0 = time.time()
        extract_one(g, bucket, label, snps)
        print(f"  ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    sys.exit(main())

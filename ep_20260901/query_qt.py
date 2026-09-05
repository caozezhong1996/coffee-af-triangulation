"""Query QT-interval summary stats (GCST90165290, harmonised) for the 40
coffee-associated SNPs using the remote tabix client. Output: CSV rows.
"""
import csv
import sys
import time

from remote_tabix import RemoteTabix

BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90165001-GCST90166000/GCST90165290/harmonised/GCST90165290.h.tsv.gz"
TBI = "GCST90165290.h.tsv.gz.tbi"
SRC = "../spectrum_20260901/spectrum_harmonized_snps.csv"
OUT = "qt_raw_40snps.csv"

COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele",
        "beta", "standard_error", "effect_allele_frequency", "p_value",
        "variant_id", "variant_id_when_present", "n",
        "hm_coordinate_conversion", "hm_code", "rsid"]


def main():
    snps = {}
    for r in csv.DictReader(open(SRC)):
        if r["phenocode"] != "I9_AF":
            continue
        snps[r["SNP"]] = (r["chr"], int(r["pos"]))
    print(f"{len(snps)} SNPs to query", flush=True)

    rt = RemoteTabix(BASE, TBI)
    print("index refs:", rt.ref_names[:30], flush=True)
    print("col_seq,col_beg,format:", rt.col_seq, rt.col_beg, rt.format, flush=True)

    rows_out, missing = [], []
    for i, (rsid, (chrom, pos)) in enumerate(snps.items(), 1):
        hit = None
        got = None
        for refname in (chrom, f"chr{chrom}"):
            if refname not in rt.refs:
                continue
            got = rt.query(refname, pos - 1, pos)
            if got:
                break
        if got:
            for f in got:
                rec = dict(zip(COLS, f))
                if rec.get("rsid") == rsid or int(rec["base_pair_location"]) == pos:
                    hit = rec
                    break
        if hit:
            hit["query_rsid"] = rsid
            rows_out.append(hit)
        else:
            missing.append(rsid)
        print(f"[{i}/40] {rsid} chr{chrom}:{pos} -> {'OK' if hit else 'MISS'}", flush=True)
        time.sleep(0.2)

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["query_rsid"] + COLS)
        w.writeheader()
        w.writerows(rows_out)
    print(f"\n{len(rows_out)}/40 found; missing: {missing}")
    print("saved:", OUT)


if __name__ == "__main__":
    sys.exit(main())

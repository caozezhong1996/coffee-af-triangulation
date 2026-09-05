"""Molecular-fingerprint analysis: coffee-consumption instrument -> caffeine-
family metabolite levels (Chen 2023 metabolomics, N~8,299).

Per metabolite: multiplicative random-effects IVW over the extracted SNPs,
overall (39) and by mechanism cluster (M=caffeine pharmacokinetics,
A=adiposity/intake propensity, O=other). BH FDR across the 9 metabolites.
"""
import csv
import glob
import math
import os
import json

RAW = "metab_raw"
CLU = "cluster_assignments.csv"
LABELS = {
    "GCST90200436": "caffeine",
    "GCST90199654": "paraxanthine",
    "GCST90199644": "theobromine",
    "GCST90199647": "theophylline",
    "GCST90200285": "xanthine",
    "GCST90199754": "7-methylxanthine",
    "GCST90199699": "3-methylxanthine",
    "GCST90199756": "1,7-dimethyluric acid",
    "GCST90199760": "AAMU",
}


def norm_cdf(z):
    return 0.5 * (1 + math.erf(z / math.sqrt(2)))


def ivw(rows):
    """rows: list of (beta, se). Multiplicative random-effects IVW."""
    n = len(rows)
    if n < 3:
        return None
    w = [1.0 / (se * se) for _, se in rows]
    b = [x for x, _ in rows]
    sw = sum(w)
    biv = sum(wi * bi for wi, bi in zip(w, b)) / sw
    q = sum(wi * (bi - biv) ** 2 for wi, bi in zip(w, b))
    phi = max(1.0, q / (n - 1))
    se = math.sqrt(phi / sw)
    z = biv / se
    p = 2 * (1 - norm_cdf(abs(z)))
    return dict(n=n, beta=biv, se=se, p=p, Q=q, Q_df=n - 1,
                Q_p=2 * (1 - norm_cdf(abs(math.sqrt(q)))) if False else None,
                I2=max(0.0, (q - (n - 1)) / q) * 100 if q > 0 else 0.0)


def bh_fdr(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [1.0] * m
    prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        val = min(prev, pvals[i] * m / rank)
        q[i] = val
        prev = val
    return q


def main():
    clusters = {}
    for r in csv.DictReader(open(CLU)):
        clusters[r["SNP"]] = r["cluster"]

    results = {}
    for path in sorted(glob.glob(os.path.join(RAW, "*.csv"))):
        g = os.path.basename(path).replace(".csv", "")
        label = LABELS.get(g, g)
        rows = []
        for r in csv.DictReader(open(path)):
            if r["allele_match"] == "mismatch":
                continue
            rows.append((r["query_rsid"], float(r["beta_raw"]), float(r["se_raw"])))
        subsets = {
            "all": rows,
            "M": [x for x in rows if clusters.get(x[0]) == "M"],
            "A": [x for x in rows if clusters.get(x[0]) == "A"],
            "O": [x for x in rows if x[0] not in clusters],
        }
        res = {}
        for k, sub in subsets.items():
            r = ivw([(b, se) for _, b, se in sub])
            if r:
                res[k] = r
        results[label] = res
        print(f"{label}: all n={res.get('all',{}).get('n')} "
              f"beta={res.get('all',{}).get('beta', float('nan')):+.4f} "
              f"p={res.get('all',{}).get('p', 1):.2e}", flush=True)

    # FDR across the 9 'all' tests
    labels = list(results)
    p_all = [results[l]["all"]["p"] for l in labels]
    q_all = bh_fdr(p_all)
    for l, qv in zip(labels, q_all):
        results[l]["all"]["q"] = qv

    with open("metab_results.json", "w") as fh:
        json.dump(results, fh, indent=2)

    with open("metab_results.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metabolite", "subset", "n", "beta", "se", "p", "q_FDR",
                    "Q", "Q_df", "I2_pct"])
        for l in labels:
            for k in ("all", "M", "A", "O"):
                r = results[l].get(k)
                if not r:
                    continue
                w.writerow([l, k, r["n"], f"{r['beta']:.5f}", f"{r['se']:.5f}",
                            f"{r['p']:.3e}",
                            f"{r.get('q', ''):.3e}" if r.get("q") else "",
                            f"{r['Q']:.1f}", r["Q_df"], f"{r['I2']:.0f}"])

    print("\n=== Overall (all SNPs), FDR-adjusted ===")
    for l in labels:
        r = results[l]["all"]
        print(f"{l:22s} beta={r['beta']:+.4f} ({r['beta']-1.96*r['se']:+.4f}.."
              f"{r['beta']+1.96*r['se']:+.4f}) p={r['p']:.2e} q={r['q']:.2e} I2={r['I2']:.0f}%")
    print("\n=== M cluster (caffeine pharmacokinetics) ===")
    for l in labels:
        r = results[l].get("M")
        if r:
            print(f"{l:22s} beta={r['beta']:+.4f} p={r['p']:.2e}")
    print("\nsaved metab_results.csv / metab_results.json")


if __name__ == "__main__":
    main()

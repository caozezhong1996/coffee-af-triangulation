"""Mechanism-stratified (cluster) MR, 2026-09-02.
Cluster M: caffeine pharmacokinetic loci (CYP1A1/CYP1A2, CYP2A6, AHR x2, POR).
Cluster A: adiposity/intake-propensity loci (FTO, MC4R, TMEM18, SEC16B, GDF15,
GCKR, MLXIPL). Cluster O: all remaining instruments.
Outcomes: FinnGen R12 AF (40 SNP), BBJ AF (34 SNP), PR/QT/HR (ep_harmonized).
Statistics: IVW multiplicative RE per cluster + Cochran Q-between (subgroup
difference) + mean F per cluster.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

CL_M = {"rs2472297": "CYP1A1/CYP1A2", "rs56113850": "CYP2A6",
        "rs4410790": "AHR", "rs73075167": "AHR", "rs1057868": "POR"}
CL_A = {"rs1421085": "FTO", "rs476828": "MC4R", "rs13387939": "TMEM18",
        "rs516636": "SEC16B", "rs75347775": "GDF15", "rs780093": "GCKR",
        "rs34060476": "MLXIPL"}


def ivw(bx, by, sy):
    w = 1 / sy**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    se_fe = 1 / np.sqrt(np.sum(w * bx**2))
    Q = np.sum(w * (by - b * bx) ** 2)
    phi = max(Q / (len(bx) - 1), 1.0) if len(bx) > 1 else 1.0
    se = se_fe * np.sqrt(phi)
    return b, se, Q


def cluster_test(df, scheme):
    """scheme: dict name -> set(snp). Returns rows + Q-between."""
    rows = []
    Qs, ks, bs, ses = [], [], {}, {}
    for name, members in scheme.items():
        sub = df[df.SNP.isin(members)]
        if len(sub) == 0:
            continue
        b, se, Q = ivw(sub.bx.values, sub.by.values, sub.sy.values)
        F = float(np.mean((sub.bx / sub.bxse) ** 2))
        p = 2 * stats.norm.sf(abs(b / se))
        rows.append(dict(cluster=name, nsnp=len(sub), F=round(F, 1),
                         b=round(b, 4), se=round(se, 4),
                         lo=round(b - 1.96 * se, 4), hi=round(b + 1.96 * se, 4),
                         p=f"{p:.2e}", Q=round(Q, 1)))
        Qs.append(Q); ks.append(len(sub)); bs[name] = b; ses[name] = se
    bt, st, Qt = ivw(df.bx.values, df.by.values, df.sy.values)
    Q_between = Qt - sum(Qs)
    df_between = len(scheme) - 1
    p_between = stats.chi2.sf(Q_between, df_between)
    return rows, Q_between, df_between, p_between


def load_fng():
    h = pd.read_csv("../spectrum_20260901/spectrum_harmonized_snps.csv")
    return h[h.phenocode == "I9_AF"][["SNP", "bx", "bxse", "by", "sy"]].copy()


def load_bbj():
    b = pd.read_csv("../key_outputs/bbj_harmonized.csv")
    b = b.rename(columns={"sx": "bxse"})
    return b[["SNP", "bx", "bxse", "by", "sy"]].copy()


def load_ep(trait_key):
    e = pd.read_csv("ep_harmonized.csv")
    e = e[e.trait.str.contains(trait_key)].copy()
    return e[["SNP", "bx", "bxse", "by", "sy"]].copy()


def schemes(df):
    allsnps = set(df.SNP)
    m = set(CL_M) & allsnps
    a = set(CL_A) & allsnps
    rest = allsnps - m
    three = {"M (caffeine metabolism)": m,
             "A (adiposity/intake propensity)": a,
             "Other": allsnps - m - a}
    return {"metabolism_vs_rest": {"M (caffeine metabolism)": m, "Other": rest},
            "three_cluster": three}


out_rows = []
for trait, df in [("AF (FinnGen R12)", load_fng()), ("AF (BBJ)", load_bbj()),
                  ("PR interval", load_ep("PR")), ("QT interval", load_ep("QT")),
                  ("Resting heart rate", load_ep("heart rate"))]:
    for sname, scheme in schemes(df).items():
        rows, Qb, dfb, pb = cluster_test(df, scheme)
        for r in rows:
            r.update(trait=trait, scheme=sname,
                     Q_between=round(Qb, 2), df_between=dfb,
                     p_between=f"{pb:.3f}")
            out_rows.append(r)
        print(f"\n== {trait} [{sname}] Q_between={Qb:.2f} df={dfb} P={pb:.3g}")
        for r in rows:
            print(f"  {r['cluster']:35s} n={r['nsnp']:2d} F={r['F']:6.1f} "
                  f"b={r['b']:+.4f} ({r['lo']:+.4f},{r['hi']:+.4f}) P={r['p']}")

res = pd.DataFrame(out_rows)
res.to_csv("cluster_mr_results.csv", index=False)
pd.DataFrame([{"SNP": k, "cluster": "M", "gene": v} for k, v in CL_M.items()] +
             [{"SNP": k, "cluster": "A", "gene": v} for k, v in CL_A.items()]
             ).to_csv("cluster_assignments.csv", index=False)
print("\nsaved cluster_mr_results.csv, cluster_assignments.csv")

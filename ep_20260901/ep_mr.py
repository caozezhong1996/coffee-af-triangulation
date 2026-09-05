"""Electrophysiology intermediate-phenotype MR (module B, 2026-09-02).
Exposure: 40-SNP coffee-consumption instrument (key_outputs/instruments.csv).
Outcomes: QT (GCST90165290 harmonised, GRCh38), resting HR (Zhu 2019,
GCST007609 BOLT raw, b37), PR (GCST010321 METAL raw, b37).
Methods: IVW (multiplicative random effects), weighted median (parametric
bootstrap), MR-Egger, Cochran Q. Outputs ep_mr_results.csv + ep_harmonized.csv.
"""
import json

import numpy as np
import pandas as pd
from scipy import stats

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

ins = pd.read_csv("../key_outputs/instruments.csv")[
    ["SNP", "effect_allele.exposure", "other_allele.exposure",
     "eaf.exposure", "beta.exposure", "se.exposure"]]
ins.columns = ["SNP", "ea", "oa", "eaf", "bx", "bxse"]
ins["ea"] = ins.ea.str.upper()
ins["oa"] = ins.oa.str.upper()


def ivw(bx, by, sy):
    w = 1 / sy**2
    b = np.sum(w * bx * by) / np.sum(w * bx**2)
    se_fe = 1 / np.sqrt(np.sum(w * bx**2))
    Q = np.sum(w * (by - b * bx) ** 2)
    k = len(bx)
    phi = max(Q / (k - 1), 1.0)
    return b, se_fe * np.sqrt(phi), Q, stats.chi2.sf(Q, k - 1)


def wmedian(bx, by, sy, nboot=5000, seed=42):
    theta = by / bx
    w = 1 / sy**2

    def med(t, ww):
        o = np.argsort(t)
        t, ww = t[o], ww[o]
        cw = (np.cumsum(ww) - 0.5 * ww) / np.sum(ww)
        return np.interp(0.5, cw, t)

    rng = np.random.default_rng(seed)
    ests = [med(rng.normal(by, sy) / bx, w) for _ in range(nboot)]
    return med(theta, w), np.std(ests, ddof=1)


def egger(bx, by, sy):
    s = np.sign(bx)
    byo, bxo = by * s, np.abs(bx)
    w = 1 / sy**2
    X = np.column_stack([np.ones_like(bxo), bxo])
    XtW = X.T * w
    beta = np.linalg.solve(XtW @ X, XtW @ byo)
    resid = byo - X @ beta
    phi = max(np.sum(w * resid**2) / (len(bx) - 2), 1.0)
    cov = phi * np.linalg.inv(XtW @ X)
    return beta[0], np.sqrt(cov[0, 0]), beta[1], np.sqrt(cov[1, 1])


def harmonise(x, o_ea, o_oa, o_beta, o_eaf, trait):
    """Return (by, keep, note) aligning outcome beta to exposure effect allele."""
    xa, xo = x.ea, x.oa
    a1, a2 = o_ea.upper(), o_oa.upper()
    pair_x = {xa, xo}
    pal = pair_x in ({"A", "T"}, {"C", "G"})
    if {a1, a2} == pair_x:
        note = "direct"
        by = o_beta if a1 == xa else -o_beta
    elif {COMP[a1], COMP[a2]} == pair_x:
        note = "strand"
        by = o_beta if COMP[a1] == xa else -o_beta
    else:
        return np.nan, False, "allele_mismatch"
    if pal and o_eaf is not None and np.isfinite(o_eaf):
        ea_f = o_eaf if (a1 == xa or COMP.get(a1) == xa) else 1 - o_eaf
        if abs(ea_f - x.eaf) > 0.2:
            return np.nan, False, "pal_freq_mismatch"
        if 0.42 < ea_f < 0.58 and 0.42 < x.eaf < 0.58:
            return np.nan, False, "pal_ambiguous"
    return by, True, note


def run(trait, df, ea_col, oa_col, b_col, se_col, eaf_col):
    m = ins.merge(df, left_on="SNP", right_on="q", how="inner")
    res = [harmonise(r, r[ea_col], r[oa_col], r[b_col],
                     r[eaf_col] if eaf_col else np.nan, trait)
           for _, r in m.iterrows()]
    m["by"] = [x[0] for x in res]
    m["keep"] = [x[1] for x in res]
    m["hnote"] = [x[2] for x in res]
    m["sy"] = m[se_col].astype(float)
    kept = m[m.keep].copy()
    kept["trait"] = trait
    bx, by, sy = kept.bx.values, kept.by.values, kept.sy.values
    b, se, Q, Qp = ivw(bx, by, sy)
    p = 2 * stats.norm.sf(abs(b / se))
    bm, sem = wmedian(bx, by, sy)
    pm = 2 * stats.norm.sf(abs(bm / sem))
    ei, eis, es, ess = egger(bx, by, sy)
    pe = 2 * stats.norm.sf(abs(es / ess))
    pi = 2 * stats.norm.sf(abs(ei / eis))
    F = float(np.mean((bx / ins.loc[kept.index, "bxse"]) ** 2))
    row = dict(trait=trait, n_instruments=len(kept), dropped=len(m) - len(kept),
               F_mean=round(F, 1),
               ivw_b=round(b, 4), ivw_se=round(se, 4), ivw_p=f"{p:.2e}",
               ivw_lo=round(b - 1.96 * se, 4), ivw_hi=round(b + 1.96 * se, 4),
               Q=round(Q, 1), Q_p=f"{Qp:.3f}",
               wm_b=round(bm, 4), wm_se=round(sem, 4), wm_p=f"{pm:.2e}",
               egger_b=round(es, 4), egger_se=round(ess, 4), egger_p=f"{pe:.2e}",
               egger_int=round(ei, 4), egger_int_se=round(eis, 4),
               egger_int_p=f"{pi:.3f}")
    return row, kept


qt = pd.read_csv("qt_raw_40snps.csv").rename(columns={"query_rsid": "q"})
hr = pd.read_csv("hr_raw_40snps.csv").rename(columns={"query_rsid": "q"})
pr = pd.read_csv("pr_raw_40snps.csv").rename(columns={"query_rsid": "q"})
for c in ["beta", "standard_error", "effect_allele_frequency"]:
    qt[c] = pd.to_numeric(qt[c], errors="coerce")
for c in ["BETA", "SE"]:
    hr[c] = pd.to_numeric(hr[c], errors="coerce")
for c in ["Effect", "StdErr", "EAF"]:
    pr[c] = pd.to_numeric(pr[c], errors="coerce")

rows, kept_all = [], []
r, k = run("QT interval (GCST90165290)", qt, "effect_allele", "other_allele",
           "beta", "standard_error", "effect_allele_frequency")
rows.append(r); kept_all.append(k)
r, k = run("Resting heart rate (Zhu 2019)", hr, "A1", "A0", "BETA", "SE", None)
rows.append(r); kept_all.append(k)
r, k = run("PR interval (GCST010321)", pr, "Allele1", "Allele2",
           "Effect", "StdErr", "EAF")
rows.append(r); kept_all.append(k)

out = pd.DataFrame(rows)
out.to_csv("ep_mr_results.csv", index=False)
pd.concat(kept_all)[["SNP", "trait", "ea", "oa", "bx", "bxse", "by", "sy",
                     "hnote"]].to_csv("ep_harmonized.csv", index=False)
print(out.to_string(index=False))

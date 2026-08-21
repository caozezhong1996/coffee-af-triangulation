# -*- coding: utf-8 -*-
# Plasma caffeine biomarker MR (post hoc), run 2026-08-20.
# Instruments: independent GW-significant loci from Cornelis et al. 2016
# (Hum Mol Genet 25:5472-5482; Table 1 Z statistics), outcomes extracted from
# local Nielsen 2018 (coloc/nielsen_gwas.tbl.gz) and FinnGen R11 I9_AF files.
# Exposure betas derived from Z under standardized scaling at conservative
# minimum N (9054 plasma caffeine; 5323 paraxanthine/caffeine ratio).
import json
from math import exp, sqrt, erfc

outc = json.load(open(r"C:\Users\64915\Documents\kimi\workspace\mr_data\biomarker_outcome_snps.json"))
cornelis = {
  "rs4410790":  ("T", 0.36,  7.36, -8.07),
  "rs6968554":  ("A", 0.37,  7.22, -8.44),
  "rs2472297":  ("T", 0.27, -9.34,  9.58),
  "rs56113850": ("T", 0.43,  None,  9.59),
}
def exp_beta(z, eaf, n):
    den = sqrt(2*eaf*(1-eaf)*n); return z/den, 1/den
def out_beta(rs, ea):
    n, g = outc["nielsen"][rs], outc["finngen"][rs]
    bn, sn = float(n["beta"]), float(n["se"])
    if n["a2"] != ea: bn = -bn
    bg, sg = float(g["beta"]), float(g["se"])
    if g["alt"] != ea: bg = -bg
    return (bn, sn), (bg, sg)
def ivw(pairs):
    ws = [1/se**2 for _, se in pairs]
    b = sum(w*bi for w, (bi, _) in zip(ws, pairs)) / sum(ws)
    se = sqrt(1/sum(ws))
    q = sum(w*(bi-b)**2 for w, (bi, _) in zip(ws, pairs))
    return b, se, q
def report(tag, snps, N, zidx):
    print(f"== {tag} ==")
    per_ds = []
    for ds in ["nielsen", "finngen"]:
        pairs = []
        for rs in snps:
            ea, eaf, z137, zr = cornelis[rs]
            be, see = exp_beta((z137, zr)[zidx], eaf, N)
            bo, so = out_beta(rs, ea)[0 if ds == "nielsen" else 1]
            pairs.append((bo/be, so/abs(be)))
            print(f"  {ds} {rs}: WR={pairs[-1][0]:+.4f} ({pairs[-1][1]:.4f})")
        b, se, q = ivw(pairs)
        per_ds.append((b, se))
        print(f"  {ds} IVW: OR {exp(b):.3f} ({exp(b-1.96*se):.3f}-{exp(b+1.96*se):.3f}), Q={q:.2f}")
    b, se, _ = ivw(per_ds)
    z = b/se
    print(f"  POOLED: OR {exp(b):.3f} ({exp(b-1.96*se):.3f}-{exp(b+1.96*se):.3f}), P={erfc(abs(z)/sqrt(2)):.3f}")

report("Plasma caffeine per SD -> AF", ["rs4410790", "rs2472297"], 9054, 0)
report("Paraxanthine/caffeine ratio per SD -> AF", ["rs6968554", "rs2472297", "rs56113850"], 5323, 1)

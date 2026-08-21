# -*- coding: utf-8 -*-
# Replication MR using Miyazawa 2023 Nat Genet cross-ancestry AF meta-GWAS
# (GCST90204201; 77,690 cases / 1,167,040 controls; combines BBJ Japan + Nielsen 2018 + FinnGen).
# NOTE: outcome overlaps with our Nielsen/FinnGen/BBJ datasets -> superset, not fully independent.
# Run: 2026-08-20
import gzip, csv, math, json

WS = r"C:\Users\64915\Documents\kimi\workspace"
GZ = WS + r"\mr_data\miyazawa_af_2023.h.tsv.gz"
INST = WS + r"\coffee_af_code\key_outputs\instruments.csv"

# ---- instrument sets ----
# 1) coffee intake (cups/day, UKB ukb-b-5237), 40 SNPs, from instruments.csv
coffee = {}
with open(INST, newline='', encoding='utf-8') as f:
    for r in csv.DictReader(f):
        if r['mr_keep.exposure'] == 'TRUE':
            coffee[r['SNP']] = (r['effect_allele.exposure'].upper(),
                                r['other_allele.exposure'].upper(),
                                float(r['beta.exposure']), float(r['se.exposure']))

# 2) plasma caffeine biomarker instruments (exposure betas per SD, from biomarker_mr_plasma_caffeine.py)
def exp_beta(z, eaf, n):
    se = 1.0 / math.sqrt(2 * eaf * (1 - eaf) * n)
    return abs(z) * se, se

biomarker = {
  # rsid: (EA, OEA, exp_beta, exp_se)  plasma caffeine, N=9054
  "rs4410790": dict(EA="T", eaf=0.36, z=7.36, N=9054, trait="plasma"),
  "rs2472297": dict(EA="T", eaf=0.27, z=-9.34, N=9054, trait="plasma"),
  "rs6968554": dict(EA="A", eaf=0.37, z=7.22, N=5323, trait="ratio"),
  "rs56113850": dict(EA="T", eaf=0.43, z=9.59, N=5323, trait="ratio"),
}
# note: rs2472297 also serves ratio (z=9.58 per original script); use same SNP row
ratio_extra_z = {"rs2472297": 9.58}

want = set(coffee) | set(biomarker)
rows = {}
with gzip.open(GZ, 'rt', encoding='utf-8', errors='replace') as f:
    header = f.readline().rstrip('\n').split('\t')
    idx = {c: i for i, c in enumerate(header)}
    for line in f:
        p = line.rstrip('\n').split('\t')
        rs = p[idx['rsid']]
        if rs in want and rs not in rows:
            rows[rs] = dict(ea=p[idx['effect_allele']].upper(),
                            oa=p[idx['other_allele']].upper(),
                            beta=float(p[idx['beta']]), se=float(p[idx['standard_error']]),
                            eaf=p[idx['effect_allele_frequency']], p=p[idx['p_value']])

print("outcome rows found:", len(rows), "of", len(want))
missing = sorted(want - set(rows))
print("missing:", missing)

COMP = {'A':'T','T':'A','C':'G','G':'C'}
def flip(a): return ''.join(COMP.get(x, x) for x in a)

def harmonize_out(rs, exp_ea, exp_oa):
    """return (beta_out, se_out) aligned to exposure EA, or None if mismatch."""
    o = rows[rs]
    if o['ea'] == exp_ea and o['oa'] == exp_oa:
        return o['beta'], o['se']
    if o['ea'] == exp_oa and o['oa'] == exp_ea:
        return -o['beta'], o['se']
    # strand flip
    if flip(o['ea']) == exp_ea and flip(o['oa']) == exp_oa:
        return o['beta'], o['se']
    if flip(o['ea']) == exp_oa and flip(o['oa']) == exp_ea:
        return -o['beta'], o['se']
    return None

def ivw(pairs):
    # pairs: list of (bx, sx, by, sy)
    num = sum((bx * by) / sy**2 for bx, sx, by, sy in pairs)
    den = sum((bx**2) / sy**2 for bx, sx, by, sy in pairs)
    b = num / den
    se_fe = 1.0 / math.sqrt(den)
    # Cochran Q
    q = sum(((by - b * bx) / sy)**2 for bx, sx, by, sy in pairs)
    k = len(pairs)
    c = max(1.0, q / (k - 1)) if k > 1 else 1.0
    se = se_fe * math.sqrt(c)
    z = b / se
    pval = math.erfc(abs(z) / math.sqrt(2))
    return b, se, pval, q, k

def wald(bx, sx, by, sy):
    b = by / bx
    se = sy / abs(bx)  # delta method, no measurement error term (NO ME assumption)
    z = b / se
    return b, se, math.erfc(abs(z) / math.sqrt(2))

# ---- coffee intake ----
pairs, skipped = [], []
for rs, (ea, oa, bx, sx) in coffee.items():
    if rs not in rows:
        skipped.append(rs); continue
    h = harmonize_out(rs, ea, oa)
    if h is None:
        skipped.append(rs + "(allele mismatch)"); continue
    pairs.append((bx, sx, h[0], h[1]))
b, se, p, q, k = ivw(pairs)
print(f"\nCoffee intake cups/day -> AF (Miyazawa meta, n={k} SNPs):")
print(f"  IVW(re) beta={b:.4f} OR={math.exp(b):.3f} ({math.exp(b-1.96*se):.3f}-{math.exp(b+1.96*se):.3f}) P={p:.2e}  Q={q:.1f}")
print("  skipped:", skipped)

# weighted median (simple bootstrap-free implementation)
def wmed(pairs):
    ws = sorted(((by/bx, (abs(bx)/sy)) for bx, sx, by, sy in pairs), key=lambda t: t[0])
    w = [x[1] for x in ws]; tot = sum(w); c = 0
    for (r, wi) in ws:
        c += wi
        if c >= tot / 2:
            return r
    return ws[-1][0]

def wmed_boot(pairs, B=10000, seed=42):
    import random
    rnd = random.Random(seed)
    ests = []
    k = len(pairs)
    for _ in range(B):
        samp = [pairs[rnd.randrange(k)] for _ in range(k)]
        ests.append(wmed(samp))
    ests.sort()
    return wmed(pairs), ests[int(0.025*B)], ests[int(0.975*B)]

def egger(pairs):
    # MR-Egger: regress by on bx, weights 1/sy^2
    k = len(pairs)
    w = [1/sy**2 for _, _, _, sy in pairs]
    x = [bx for bx, _, _, _ in pairs]; y = [by for _, _, by, _ in pairs]
    W = sum(w); Wx = sum(wi*xi for wi, xi in zip(w, x)); Wy = sum(wi*yi for wi, yi in zip(w, y))
    Wxx = sum(wi*xi*xi for wi, xi in zip(w, x)); Wxy = sum(wi*xi*yi for wi, xi, yi in zip(w, x, y))
    slope = (Wxy - Wx*Wy/W) / (Wxx - Wx*Wx/W)
    inter = (Wy - slope*Wx) / W
    # residual variance (multiplicative random effects)
    res = sum(wi*(yi - inter - slope*xi)**2 for wi, xi, yi in zip(w, x, y))
    sig2 = max(1.0, res/(k-2))
    se_slope = (sig2/(Wxx - Wx*Wx/W))**0.5
    se_inter = (sig2*(1/W + (Wx/W)**2/(Wxx - Wx*Wx/W)))**0.5
    z_i = inter/se_inter
    return slope, se_slope, math.erfc(abs(slope/se_slope)/math.sqrt(2)), inter, se_inter, math.erfc(abs(z_i)/math.sqrt(2))
print(f"  Weighted median OR={math.exp(wmed(pairs)):.3f}")
wm, wmlo, wmhi = wmed_boot(pairs)
print(f"  Weighted median (bootstrap) OR={math.exp(wm):.3f} ({math.exp(wmlo):.3f}-{math.exp(wmhi):.3f})")
es, ese, esp, ei, eie, eip = egger(pairs)
print(f"  MR-Egger slope OR={math.exp(es):.3f} ({math.exp(es-1.96*ese):.3f}-{math.exp(es+1.96*ese):.3f}) P={esp:.3f}; intercept={ei:.4f} (P={eip:.3f})")

# ---- biomarkers ----
def biomarker_run(name, rsids, trait):
    pairs = []
    for rs in rsids:
        d = biomarker[rs]
        z = d['z'] if (trait == 'plasma' or rs != 'rs2472297') else ratio_extra_z[rs]
        bx, sx = exp_beta(z, d['eaf'], d['N'])
        if rs not in rows:
            print(f"  {rs} missing in outcome"); continue
        h = harmonize_out(rs, d['EA'], None) if False else None
        o = rows[rs]
        # harmonize: outcome EA vs exposure EA (other allele may differ in file; check both)
        if o['ea'] == d['EA']:
            by, sy = o['beta'], o['se']
        elif o['oa'] == d['EA']:
            by, sy = -o['beta'], o['se']
        elif flip(o['ea']) == d['EA']:
            by, sy = o['beta'], o['se']
        elif flip(o['oa']) == d['EA']:
            by, sy = -o['beta'], o['se']
        else:
            print(f"  {rs} allele mismatch out {o['ea']}/{o['oa']} vs exp {d['EA']}"); continue
        pairs.append((bx, sx, by, sy))
    if len(pairs) == 1:
        b, se, p = wald(*pairs[0]); k = 1; q = float('nan')
    else:
        b, se, p, q, k = ivw(pairs)
    print(f"{name} (Miyazawa meta, n={k} SNPs): OR={math.exp(b):.3f} ({math.exp(b-1.96*se):.3f}-{math.exp(b+1.96*se):.3f}) P={p:.3f}")

print()
biomarker_run("Plasma caffeine per SD -> AF", ["rs4410790", "rs2472297"], "plasma")
biomarker_run("Paraxanthine/caffeine ratio per SD -> AF", ["rs6968554", "rs2472297", "rs56113850"], "ratio")

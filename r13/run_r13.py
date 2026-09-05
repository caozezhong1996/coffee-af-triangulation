# -*- coding: utf-8 -*-
"""FinnGen R13 I9_AF sensitivity: harmonise 40 instruments, IVW + Egger."""
import pandas as pd, numpy as np, math, gzip

def norm_sf(z):
    return math.erfc(abs(z) / math.sqrt(2))

comp = {'A':'T','T':'A','C':'G','G':'C'}
def flip(a): return ''.join(comp[x] for x in a)

inst = pd.read_csv(r'mr_results\instruments.csv')
targets = {r['SNP']: (r['effect_allele.exposure'], r['other_allele.exposure'],
                      r['beta.exposure'], r['se.exposure'])
           for _, r in inst.iterrows()}

found = {}
with gzip.open(r'r13\finngen_R13_I9_AF.gz', 'rt') as g:
    g.readline()
    for line in g:
        p = line.split('\t')
        rs = p[4]
        if rs in targets and rs not in found:
            found[rs] = p

print(f'matched {len(found)} / {len(targets)} instruments by rsID')

rows = []
for rsid, (ea, oa, bx, sx) in targets.items():
    if rsid not in found:
        rows.append(dict(SNP=rsid, ok=False)); continue
    p = found[rsid]
    ref, alt = p[2], p[3]
    beta, seb, af = float(p[8]), float(p[9]), float(p[10])
    if alt == ea and ref == oa:
        by = beta
    elif alt == oa and ref == ea:
        by = -beta
    elif flip(alt) == ea and flip(ref) == oa:
        by = beta
    elif flip(alt) == oa and flip(ref) == ea:
        by = -beta
    else:
        rows.append(dict(SNP=rsid, ok=False, note=f'allele mismatch {ref}/{alt} vs {oa}/{ea}')); continue
    rows.append(dict(SNP=rsid, ok=True, bx=bx, sx=sx, by=by, sy=seb, af_alt=af))

h = pd.DataFrame([r for r in rows if r.get('ok')])
miss = [r['SNP'] for r in rows if not r.get('ok')]
print('harmonized:', len(h), ' missing:', miss)

theta = (h.by / h.bx).values
w = ((h.bx / h.sy) ** 2).values
ivw = np.sum(w * theta) / np.sum(w)
se = math.sqrt(1 / np.sum(w))
Q = float(np.sum(w * (theta - ivw) ** 2)); df = len(h) - 1
I2 = max(0.0, (Q - df) / Q) * 100
phi = max(1.0, Q / df); se_re = se * math.sqrt(phi)
print(f'R13 IVW(RE): beta={ivw:.3f} se={se_re:.3f} OR={math.exp(ivw):.2f} '
      f'({math.exp(ivw-1.96*se_re):.2f}-{math.exp(ivw+1.96*se_re):.2f}) '
      f'P={norm_sf(ivw/se_re):.4f} Q={Q:.1f} df={df} I2={I2:.0f}%')

bxv, byv = h.bx.values, h.by.values
xm = np.sum(w*bxv)/np.sum(w); ym = np.sum(w*byv)/np.sum(w)
b1 = np.sum(w*(bxv-xm)*(byv-ym))/np.sum(w*(bxv-xm)**2)
b0 = ym - b1*xm
resid = byv - b0 - b1*bxv
s2 = np.sum(w*resid**2)/(len(h)-2)
se_b0 = math.sqrt(s2*(1/np.sum(w) + xm**2/np.sum(w*(bxv-xm)**2)))
print(f'R13 Egger intercept: {b0:.4f} (se {se_b0:.4f}) P={norm_sf(b0/se_b0):.2f} slope={b1:.3f} OR_slope={math.exp(b1):.2f}')

h.to_csv(r'r13\r13_harmonized.csv', index=False)
print('saved r13\r13_harmonized.csv')

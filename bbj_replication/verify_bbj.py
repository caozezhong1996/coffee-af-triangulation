# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, zipfile, gzip, math

def norm_sf(z):  # two-sided p for |z|
    return math.erfc(abs(z) / math.sqrt(2))

# 1) check rs762551 / rs4410790 in BBJ raw file
z = zipfile.ZipFile(r'coloc\bbj_AF.zip')
name = [n for n in z.namelist() if n.endswith('.txt.gz')][0]
hits = {}
with z.open(name) as f:
    with gzip.open(f, 'rt') as g:
        g.readline()
        for line in g:
            p = line.split('\t')
            if p[3] in ('rs762551', 'rs4410790'):
                hits[p[3]] = (p[9], p[10], p[13].strip())
print('BBJ lookup:', hits if hits else 'neither rs762551 nor rs4410790 found')

# 2) recompute IVW / Egger from harmonized data
h = pd.read_csv(r'coloc\bbj_harmonized.csv')
theta, w = h.theta.values, h.w.values
ivw = np.sum(w * theta) / np.sum(w)
ivw_se = math.sqrt(1 / np.sum(w))
Q = float(np.sum(w * (theta - ivw) ** 2)); df = len(h) - 1
I2 = max(0.0, (Q - df) / Q) * 100
phi = max(1.0, Q / df)
ivw_se_re = ivw_se * math.sqrt(phi)
p_ivw = norm_sf(ivw / ivw_se_re)
print(f'BBJ IVW(RE): beta={ivw:.3f} se={ivw_se_re:.3f} OR={math.exp(ivw):.2f} '
      f'({math.exp(ivw-1.96*ivw_se_re):.2f}-{math.exp(ivw+1.96*ivw_se_re):.2f}) '
      f'P={p_ivw:.3f} Q={Q:.1f} df={df} I2={I2:.0f}%')

# Egger intercept (weighted)
X = h.bx.values; y = h.by.values
xm = np.sum(w * X) / np.sum(w); ym = np.sum(w * y) / np.sum(w)
b1 = np.sum(w * (X - xm) * (y - ym)) / np.sum(w * (X - xm) ** 2)
b0 = ym - b1 * xm
resid = y - b0 - b1 * X
s2 = np.sum(w * resid ** 2) / (len(h) - 2)
se_b0 = math.sqrt(s2 * (1 / np.sum(w) + xm ** 2 / np.sum(w * (X - xm) ** 2)))
p_b0 = norm_sf(b0 / se_b0)
print(f'BBJ Egger intercept: {b0:.4f} (se {se_b0:.4f}) P={p_b0:.2f}')

# 3) three-dataset pooled fixed-effect (Nielsen, FinnGen, BBJ)
est = [(math.log(1.22), (math.log(1.53) - math.log(0.97)) / 3.9199),
       (math.log(1.63), (math.log(2.36) - math.log(1.13)) / 3.9199),
       (ivw, ivw_se_re)]
ww = np.array([1 / e[1] ** 2 for e in est]); bb = np.array([e[0] for e in est])
pool = np.sum(ww * bb) / np.sum(ww); pool_se = math.sqrt(1 / np.sum(ww))
Qp = float(np.sum(ww * (bb - pool) ** 2))
p_pool = norm_sf(pool / pool_se)
p_Q = math.exp(-Qp / 2)  # chi2 df=2 SF
print(f'Pooled: OR={math.exp(pool):.2f} ({math.exp(pool-1.96*pool_se):.2f}-{math.exp(pool+1.96*pool_se):.2f}) '
      f'P={p_pool:.3f} Q={Qp:.2f} df=2 P={p_Q:.2f}')
print(f'BBJ weight in pooled: {ww[2]/np.sum(ww)*100:.1f}%')

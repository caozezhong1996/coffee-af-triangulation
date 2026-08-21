# -*- coding: utf-8 -*-
import pandas as pd, numpy as np, math

def norm_sf(z):
    return math.erfc(abs(z) / math.sqrt(2))

h = pd.read_csv(r'coloc\bbj_harmonized.csv')
bx, by, sy = h.bx.values, h.by.values, h.sy.values
theta = by / bx
w = (bx / sy) ** 2  # correct IVW weight for Wald ratio

ivw = np.sum(w * theta) / np.sum(w)
se = math.sqrt(1 / np.sum(w))
Q = float(np.sum(w * (theta - ivw) ** 2)); df = len(h) - 1
I2 = max(0.0, (Q - df) / Q) * 100
phi = max(1.0, Q / df)
se_re = se * math.sqrt(phi)
p = norm_sf(ivw / se_re)
print(f'BBJ IVW: beta={ivw:.3f} SE={se:.3f} (RE {se_re:.3f}) OR={math.exp(ivw):.2f} '
      f'({math.exp(ivw-1.96*se_re):.2f}-{math.exp(ivw+1.96*se_re):.2f}) P={p:.3f} '
      f'Q={Q:.1f} df={df} I2={I2:.0f}%')

# Egger with correct weights
xm = np.sum(w * bx) / np.sum(w); ym = np.sum(w * by) / np.sum(w)
b1 = np.sum(w * (bx - xm) * (by - ym)) / np.sum(w * (bx - xm) ** 2)
b0 = ym - b1 * xm
resid = by - b0 - b1 * bx
s2 = np.sum(w * resid ** 2) / (len(h) - 2)
se_b0 = math.sqrt(s2 * (1 / np.sum(w) + xm ** 2 / np.sum(w * (bx - xm) ** 2)))
print(f'BBJ Egger intercept: {b0:.4f} (se {se_b0:.4f}) P={norm_sf(b0/se_b0):.2f}, slope={b1:.3f}')

# three-dataset pooled with corrected BBJ
est = [(math.log(1.22), (math.log(1.53) - math.log(0.97)) / 3.9199),
       (math.log(1.63), (math.log(2.36) - math.log(1.13)) / 3.9199),
       (ivw, se_re)]
ww = np.array([1 / e[1] ** 2 for e in est]); bb = np.array([e[0] for e in est])
pool = np.sum(ww * bb) / np.sum(ww); pool_se = math.sqrt(1 / np.sum(ww))
Qp = float(np.sum(ww * (bb - pool) ** 2))
print(f'Pooled: OR={math.exp(pool):.2f} ({math.exp(pool-1.96*pool_se):.2f}-{math.exp(pool+1.96*pool_se):.2f}) '
      f'P={norm_sf(pool/pool_se):.3f} Q={Qp:.2f} df=2 P={math.exp(-Qp/2):.2f} BBJwt={ww[2]/np.sum(ww)*100:.1f}%')

# fix the csv w column
h['w'] = w
h['theta'] = theta
h.to_csv(r'coloc\bbj_harmonized.csv', index=False)
print('csv w/theta columns corrected')

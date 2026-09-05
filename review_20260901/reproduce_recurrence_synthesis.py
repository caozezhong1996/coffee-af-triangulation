# -*- coding: utf-8 -*-
"""Reproduce the post-rhythm-control recurrence synthesis (Supplementary Table S26).

Two published ratio measures are pooled by fixed-effects inverse-variance
weighting on the log scale:

  1. DECAF (multicenter RCT, after electrical cardioversion, N=200):
     continued caffeinated coffee (>=1 cup/day) vs abstinence,
     Cox HR for AF/flutter recurrence over 6 months = 0.61 (0.42-0.89).
  2. Hagele et al. (retrospective cohort, first-time HPSD PVI, N=199):
     >=2 vs 0-1 cups/day, multivariable-adjusted OR for 12-month atrial
     arrhythmia recurrence (after 90-day blanking) = 0.47 (0.23-0.95).

Expected output (recurrence_synthesis.json):
  pooled 0.5762 (0.4134-0.8029), P = .0011, I2 = 0%, Q P = .524.
"""
import json
import math

STUDIES = [
    # name, point estimate, CI lower, CI upper
    ("DECAF (HR, 6-month AF/flutter recurrence)", 0.61, 0.42, 0.89),
    ("Hagele et al. (adj. OR, 12-month recurrence after PVI)", 0.47, 0.23, 0.95),
]

Z = 1.959963984540054  # 97.5th percentile of the standard normal


def norm_sf(z):
    return math.erfc(abs(z) / math.sqrt(2))


def main():
    rows = []
    for name, est, lo, hi in STUDIES:
        log_est = math.log(est)
        se = (math.log(hi) - math.log(lo)) / (2 * Z)
        w = 1.0 / (se * se)
        rows.append(dict(name=name, est=est, lo=lo, hi=hi,
                         log_est=log_est, se=se, w=w))
        print(f"{name}: log={log_est:+.4f}  se={se:.4f}  w={w:.2f}")

    wsum = sum(r["w"] for r in rows)
    pooled = sum(r["w"] * r["log_est"] for r in rows) / wsum
    se_pool = math.sqrt(1.0 / wsum)
    lo = math.exp(pooled - Z * se_pool)
    hi = math.exp(pooled + Z * se_pool)
    p = norm_sf(pooled / se_pool)

    Q = sum(r["w"] * (r["log_est"] - pooled) ** 2 for r in rows)
    df = len(rows) - 1
    # chi-square with df=1 survival function
    Q_P = norm_sf(math.sqrt(Q)) if df == 1 else None
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0

    out = dict(pooled=math.exp(pooled), lo=lo, hi=hi, P=p, Q_P=Q_P, I2=I2)
    print("\nPooled (fixed effect):")
    print(f"  {out['pooled']:.4f} ({out['lo']:.4f}-{out['hi']:.4f})  "
          f"P={out['P']:.4f}  I2={out['I2']:.0f}%  Q_P={out['Q_P']:.3f}")
    for r in rows:
        print(f"  weight {r['name']}: {r['w'] / wsum * 100:.0f}%")

    with open("recurrence_synthesis.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("saved recurrence_synthesis.json")


if __name__ == "__main__":
    main()

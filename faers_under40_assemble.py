# -*- coding: utf-8 -*-
"""Assemble faers_under40_counts.json and faers_under40_summary.csv from the
counts collected on 2026-09-06 (avoids re-querying the slow endpoint)."""
import csv
import json
from math import erfc, exp, log, sqrt

RAW = {
    "lt40": dict(caf_af=11, caf_tot=3693, bg_af=1742, bg_tot=1493260),
    "m40":  dict(caf_af=28, caf_tot=7028, bg_af=17063, bg_tot=2938327),
    "g65":  dict(caf_af=63, caf_tot=3086, bg_af=43744, bg_tot=2867996),
}


def ror(a, b, c, d):
    r = a * d / (b * c)
    se = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = log(r) / se
    return r, exp(log(r) - 1.96 * se), exp(log(r) + 1.96 * se), erfc(abs(z) / sqrt(2))


counts = {"query_date": "2026-09-06",
          "note": "serious reports; age restricted to reports with onset age in years (unit=801)"}
for k, d in RAW.items():
    a, b = d["caf_af"], d["caf_tot"] - d["caf_af"]
    c, e = d["bg_af"] - a, d["bg_tot"] - d["bg_af"] - b
    r, lo, hi, p = ror(a, b, c, e)
    counts[k] = dict(d, ROR=r, lo=lo, hi=hi, P=p)
counts["lt40_nomask"] = 11
counts["lt40_ed_overlap"] = 0
counts["ed_af_lt40"] = 0

with open("faers_under40_counts.json", "w") as fh:
    json.dump(counts, fh, indent=1)

with open("faers_under40_summary.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["stratum", "caffeine_AF", "caffeine_total", "background_AF",
                "background_total", "ROR", "CI_lo", "CI_hi", "P"])
    for name in ("lt40", "m40", "g65"):
        d = counts[name]
        w.writerow([name, d["caf_af"], d["caf_tot"], d["bg_af"], d["bg_tot"],
                    f"{d['ROR']:.2f}", f"{d['lo']:.2f}", f"{d['hi']:.2f}",
                    f"{d['P']:.4f}"])
    w.writerow([])
    w.writerow(["energy-drink overlap in <40 caffeine+AF", counts["lt40_ed_overlap"]])
    w.writerow(["energy-drink AF reports <40 (any)", counts["ed_af_lt40"]])
    w.writerow(["<40 caffeine+AF surviving masking correction", counts["lt40_nomask"]])

for name in ("lt40", "m40", "g65"):
    d = counts[name]
    print(f"{name}: {d['caf_af']}/{d['caf_tot']} vs {d['bg_af']}/{d['bg_tot']} "
          f"-> ROR {d['ROR']:.2f} ({d['lo']:.2f}-{d['hi']:.2f}) P={d['P']:.4f}")
print("saved faers_under40_counts.json / faers_under40_summary.csv")

# -*- coding: utf-8 -*-
"""Final clean pipeline for the FAERS <40 case-level characterization (M6).

Stage A (counts): stratum-specific counts and RORs as of query date.
Stage B (facets): aggregate composition of the <40 caffeine+AF reports
(products, indications, substances, years, outcomes) via small count facets.
Stage C (optional): full report download is bandwidth-limited; the companion
script faers_under40_salvage.py tabulates per-report detail from partial
downloads when available.

Outputs: faers_under40_counts.json, faers_under40_facets.json,
faers_under40_summary.csv
"""
import csv
import json
import time
import urllib.parse
import urllib.request
from math import erfc, exp, log, sqrt

BASE = "https://api.fda.gov/drug/event.json"
AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
CAF = 'patient.drug.medicinalproduct:"caffeine"'
SER = "serious:1"
YRS = "patient.patientonsetageunit:801"  # age expressed in years
LT40 = "patient.patientonsetage:[0 TO 39]"
M40 = "patient.patientonsetage:[40 TO 64]"
G65 = "patient.patientonsetage:[65 TO 120]"
ED = ('patient.drug.medicinalproduct:"5-hour energy" OR '
      'patient.drug.medicinalproduct:"monster energy" OR '
      'patient.drug.medicinalproduct:"red bull" OR '
      'patient.drug.medicinalproduct:"rockstar" OR '
      'patient.drug.medicinalproduct:"bang energy" OR '
      'patient.drug.medicinalproduct:"celsius"')

FACETS = {
    "products": "patient.drug.medicinalproduct.exact",
    "indications": "patient.drug.drugindication.exact",
    "substances": "patient.drug.openfda.substance_name.exact",
    "sex": "patient.patientsex",
    "death": "seriousnessdeath",
    "reporter": "primarysource.qualification",
    "receivedate": "receivedate",
}


def get(url, tries=6):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "coffee-af-research/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if i == tries - 1:
                return {"error": str(e)}
            time.sleep(8 * (i + 1))
        except Exception as e:
            if i == tries - 1:
                return {"error": str(e)}
            time.sleep(8 * (i + 1))


def count(search):
    d = get(BASE + "?search=" + urllib.parse.quote(search) + "&limit=1")
    if d is None:
        return 0
    if "error" in d:
        return d["error"]
    return d["meta"]["results"]["total"]


def ror(a, b, c, d):
    r = a * d / (b * c)
    se = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = log(r) / se
    return r, exp(log(r) - 1.96 * se), exp(log(r) + 1.96 * se), erfc(abs(z) / sqrt(2))


def main():
    counts = {"query_date": "2026-09-06"}
    for name, s in (("lt40", LT40), ("m40", M40), ("g65", G65)):
        d = dict(
            caf_af=count(f"{SER} AND {YRS} AND {s} AND {CAF} AND {AF}"),
            caf_tot=count(f"{SER} AND {YRS} AND {s} AND {CAF}"),
            bg_af=count(f"{SER} AND {YRS} AND {s} AND {AF}"),
            bg_tot=count(f"{SER} AND {YRS} AND {s}"),
        )
        a, b = d["caf_af"], d["caf_tot"] - d["caf_af"]
        c, e = d["bg_af"] - a, d["bg_tot"] - d["bg_af"] - b
        d["ROR"], d["lo"], d["hi"], d["P"] = ror(a, b, c, e)
        counts[name] = d
        time.sleep(3)
    counts["lt40_ed_overlap"] = count(
        f"{SER} AND {YRS} AND {LT40} AND {CAF} AND {AF} AND ({ED})")
    counts["ed_af_lt40"] = count(f"{SER} AND {YRS} AND {LT40} AND ({ED}) AND {AF}")

    with open("faers_under40_counts.json", "w") as fh:
        json.dump(counts, fh, indent=1)
    print(json.dumps(counts, indent=1), flush=True)

    sel = f"{SER} AND {YRS} AND {LT40} AND {CAF} AND {AF}"
    facets = {"query": sel, "query_date": "2026-09-06"}
    for name, field in FACETS.items():
        d = get(BASE + "?search=" + urllib.parse.quote(sel) + "&count=" + field)
        facets[name] = ("no results" if d is None else d).get("results", d) \
            if isinstance(d, dict) else d
        print(f"[facet {name}] done", flush=True)
        time.sleep(4)
    with open("faers_under40_facets.json", "w", encoding="utf-8") as fh:
        json.dump(facets, fh, indent=1, ensure_ascii=False)

    # compact summary table
    with open("faers_under40_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["stratum", "caffeine_AF", "caffeine_total",
                    "background_AF", "background_total",
                    "ROR", "CI_lo", "CI_hi", "P"])
        for name in ("lt40", "m40", "g65"):
            d = counts[name]
            w.writerow([name, d["caf_af"], d["caf_tot"], d["bg_af"],
                        d["bg_tot"], f"{d['ROR']:.2f}", f"{d['lo']:.2f}",
                        f"{d['hi']:.2f}", f"{d['P']:.4f}"])
        w.writerow([])
        w.writerow(["energy-drink overlap in <40 caffeine+AF", counts["lt40_ed_overlap"]])
        w.writerow(["energy-drink AF reports <40 (any)", counts["ed_af_lt40"]])
    print("saved faers_under40_counts.json / faers_under40_facets.json / "
          "faers_under40_summary.csv")


if __name__ == "__main__":
    main()

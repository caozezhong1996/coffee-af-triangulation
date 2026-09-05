# -*- coding: utf-8 -*-
"""Salvage complete report objects from the truncated openFDA download
(faers_under40_raw.json) by incremental raw_decode; tabulate per-report
age, year, caffeine role/product/indication, co-drugs, death flag."""
import json
import re

raw = open("faers_under40_raw.json", encoding="utf-8", errors="replace").read()
i = raw.index('"results"')
i = raw.index("[", i)
dec = json.JSONDecoder()
pos = i + 1
reports = []
while pos < len(raw):
    while pos < len(raw) and raw[pos] in ", \n\r\t":
        pos += 1
    try:
        obj, end = dec.raw_decode(raw, pos)
        reports.append(obj)
        pos = end
    except json.JSONDecodeError:
        break

print("complete reports salvaged:", len(reports))
out = []
for rep in reports:
    pat = rep.get("patient", {})
    drugs = pat.get("drug", [])
    caf = [d for d in drugs if "caffeine" in d.get("medicinalproduct", "").lower()]
    out.append(dict(
        id=rep.get("safetyreportid"),
        date=rep.get("receivedate"),
        death=rep.get("seriousnessdeath"),
        age=pat.get("patientonsetage"),
        sex=pat.get("patientsex"),
        n_drugs=len(drugs),
        caf_products=[d.get("medicinalproduct") for d in caf],
        caf_roles=[d.get("drugcharacterization") for d in caf],
        caf_indications=[d.get("drugindication") for d in caf],
        all_products=[d.get("medicinalproduct") for d in drugs],
    ))
json.dump(out, open("faers_under40_salvaged.json", "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
for r in out:
    print(f"\n--- {r['id']} | {r['date']} | death={r['death']} | age={r['age']} sex={r['sex']} | n_drugs={r['n_drugs']}")
    print("  caffeine entries:", list(zip(r["caf_products"], r["caf_roles"], r["caf_indications"])))
    print("  all:", "; ".join(str(p) for p in r["all_products"])[:400])

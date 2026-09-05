# -*- coding: utf-8 -*-
"""Facet-based characterization of the <40 caffeine+AF serious FAERS reports
(small responses; robust to slow links)."""
import json
import time
import urllib.parse
import urllib.request

BASE = "https://api.fda.gov/drug/event.json"
AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
CAF = 'patient.drug.medicinalproduct:"caffeine"'
SER = "serious:1"
YRS = "patient.patientonsetageunit:801"
LT40 = "patient.patientonsetage:[0 TO 39]"
SEL = f"{SER} AND {YRS} AND {LT40} AND {CAF} AND {AF}"

FACETS = {
    "products": "patient.drug.medicinalproduct.exact",
    "caffeine_role": "patient.drug.drugcharacterization",
    "indications": "patient.drug.drugindication.exact",
    "substances": "patient.drug.openfda.substance_name.exact",
    "ages": "patient.patientonsetage.exact",
    "sex": "patient.patientsex",
    "year": "receivedate",
    "country": "occurcountry",
    "reporter": "primarysource.qualification",
    "outcomes": "seriousnessdeath",
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


def facet(field):
    url = BASE + "?search=" + urllib.parse.quote(SEL) + "&count=" + field
    d = get(url)
    if d is None:
        return {"note": "no results (404)"}
    if "error" in d:
        return d
    return d.get("results", [])


def main():
    out = {"query": SEL, "query_date": "2026-09-06"}
    for name, field in FACETS.items():
        out[name] = facet(field)
        print(f"[{name}] done", flush=True)
        time.sleep(4)
    with open("faers_under40_facets.json", "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, ensure_ascii=False)
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()

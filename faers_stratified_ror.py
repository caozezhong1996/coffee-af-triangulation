# -*- coding: utf-8 -*-
# FAERS caffeine-AF disproportionality by sex and age strata (post hoc, 2026-08-20).
# Queries openFDA drug/event counts and computes stratum-specific RORs.
import json, urllib.request, urllib.parse, time
from math import log, exp, sqrt, erfc

BASE = "https://api.fda.gov/drug/event.json"
def count(search, tries=4):
    url = BASE + "?search=" + urllib.parse.quote(search) + "&limit=1"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.load(r)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404: return 0
            if i == tries-1: raise
            time.sleep(6*(i+1))

AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
CAF = 'patient.drug.medicinalproduct:"caffeine"'
SER = 'serious:1'
strata = {
  "male":     'patient.patientsex:1',
  "female":   'patient.patientsex:2',
  "age_lt40": 'patient.patientonsetage:[0 TO 39]',
  "age_40_64":'patient.patientonsetage:[40 TO 64]',
  "age_ge65": 'patient.patientonsetage:[65 TO 120]',
}
res = {"overall": dict(caf_af=count(f"{SER} AND {CAF} AND {AF}"),
                       caf_tot=count(f"{SER} AND {CAF}"),
                       bg_af=count(f"{SER} AND {AF}"),
                       bg_tot=count(SER))}
for name, s in strata.items():
    time.sleep(2)
    res[name] = dict(caf_af=count(f"{SER} AND {s} AND {CAF} AND {AF}"),
                     caf_tot=count(f"{SER} AND {s} AND {CAF}"),
                     bg_af=count(f"{SER} AND {s} AND {AF}"),
                     bg_tot=count(f"{SER} AND {s}"))
for k, d in res.items():
    a = d["caf_af"]; b = d["caf_tot"]-a; c = d["bg_af"]-a; e = d["bg_tot"]-d["bg_af"]-b
    ror = a*e/(b*c); se = sqrt(1/a+1/b+1/c+1/e)
    z = log(ror)/se
    print(k, f"ROR {ror:.2f} ({exp(log(ror)-1.96*se):.2f}-{exp(log(ror)+1.96*se):.2f}) P={erfc(abs(z)/sqrt(2)):.3f}")

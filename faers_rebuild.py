import urllib.parse, urllib.request, json, time

BASE = "https://api.fda.gov/drug/event.json"
def count(search, tries=4):
    url = BASE + "?search=" + urllib.parse.quote(search) + "&limit=1"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.load(r)["meta"]["results"]["total"]
        except urllib.error.HTTPError as e:
            if e.code == 404: return 0
            if i == tries-1: return f"ERR:{e}"
            time.sleep(6*(i+1))
        except Exception as e:
            if i == tries-1: return f"ERR:{e}"
            time.sleep(6*(i+1))

AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
SER = 'serious:1'
V2  = 'safetyreportversion:[2 TO 100]'
MASK_DRUGS = '"amiodarone" OR "dronedarone" OR "adenosine" OR "dofetilide" OR "ibutilide" OR "flecainide" OR "propafenone" OR "sotalol" OR "digoxin" OR "dobutamine" OR "dopamine" OR "epinephrine" OR "norepinephrine" OR "milrinone" OR "theophylline"'
NOTMASK = f'NOT patient.drug.medicinalproduct:({MASK_DRUGS})'
exposures = {
  "caffeine": 'patient.drug.medicinalproduct:"caffeine"',
  "caffeine citrate": 'patient.drug.medicinalproduct:"caffeine citrate"',
  "taurine": 'patient.drug.medicinalproduct:"taurine"',
  "energy drink (brand names)": '(patient.drug.medicinalproduct:"5-hour energy" OR patient.drug.medicinalproduct:"monster energy" OR patient.drug.medicinalproduct:"red bull" OR patient.drug.medicinalproduct:"rockstar" OR patient.drug.medicinalproduct:"bang energy" OR patient.drug.medicinalproduct:"celsius")',
}
out = {}
def q(key, s):
    out[key] = count(s)
    print(key, out[key], flush=True)
    time.sleep(3)

q("bg_all", SER); q("bg_af", f"{SER} AND {AF}")
q("bg_all_v2", f"{SER} AND {V2}"); q("bg_af_v2", f"{SER} AND {AF} AND {V2}")
q("bg_all_nomask", f"{SER} AND {NOTMASK}"); q("bg_af_nomask", f"{SER} AND {AF} AND {NOTMASK}")
for name, e in exposures.items():
    k = name.split(" ")[0].replace("-","")
    q(f"{k}_total", f"{SER} AND {e}")
    q(f"{k}_af", f"{SER} AND {e} AND {AF}")
    q(f"{k}_total_v2", f"{SER} AND {e} AND {V2}")
    q(f"{k}_af_v2", f"{SER} AND {e} AND {AF} AND {V2}")
    q(f"{k}_total_nomask", f"{SER} AND {e} AND {NOTMASK}")
    q(f"{k}_af_nomask", f"{SER} AND {e} AND {AF} AND {NOTMASK}")
json.dump(out, open("faers_rebuild.json","w"), indent=1)
print("SAVED")

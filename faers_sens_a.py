import urllib.parse, urllib.request, json, time

BASE = "https://api.fda.gov/drug/event.json"
def count(search, tries=4):
    url = BASE + "?search=" + urllib.parse.quote(search) + "&limit=1"
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=40) as r:
                return json.load(r)["meta"]["results"]["total"]
        except Exception as e:
            if i == tries-1: return f"ERR:{e}"
            time.sleep(6*(i+1))

AF  = 'patient.reaction.reactionmeddrapt:"atrial fibrillation"'
SER = 'serious:1'
V2  = 'safetyreportversion:[2 TO 100]'
exposures = {
  "caffeine": 'patient.drug.medicinalproduct:"caffeine"',
  "caffeine citrate": 'patient.drug.medicinalproduct:"caffeine citrate"',
  "taurine": 'patient.drug.medicinalproduct:"taurine"',
  "5-hour ENERGY": 'patient.drug.medicinalproduct:"5-hour energy"',
  "MONSTER ENERGY": 'patient.drug.medicinalproduct:"monster energy"',
  "RED BULL": 'patient.drug.medicinalproduct:"red bull"',
}
out = {}
def q(key, s):
    out[key] = count(s)
    print(key, out[key], flush=True)
    time.sleep(4)

q("bg_all", SER); q("bg_af", f"{SER} AND {AF}")
for name, e in exposures.items():
    q(f"{name}_total", f"{SER} AND {e}")
    q(f"{name}_af", f"{SER} AND {e} AND {AF}")
q("bg_all_v2", f"{SER} AND {V2}"); q("bg_af_v2", f"{SER} AND {AF} AND {V2}")
for name, e in exposures.items():
    q(f"{name}_total_v2", f"{SER} AND {e} AND {V2}")
    q(f"{name}_af_v2", f"{SER} AND {e} AND {AF} AND {V2}")
json.dump(out, open("faers_sens_a.json","w"), indent=1)
print("SAVED")

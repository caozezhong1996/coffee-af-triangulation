import json, subprocess, time, sys
import pandas as pd

ins = pd.read_csv("../key_outputs/instruments.csv")
targets = {"I9_PAROXTAC","I9_OTHARR","I9_AVBLOCK","I9_CONDUCTIO","I9_AF"}

def curl(url):
    for _ in range(3):
        try:
            out = subprocess.run(["curl","-s","--max-time","90",url], capture_output=True, timeout=120).stdout
            return out.decode("utf-8","ignore")
        except Exception:
            time.sleep(2)
    return ""

rows = []
fails = []
for _, r in ins.iterrows():
    rs = r["SNP"]; ea = r["effect_allele.exposure"].upper(); oa = r["other_allele.exposure"].upper()
    ac = curl(f"https://r12.finngen.fi/api/autocomplete?query={rs}")
    try:
        cand = json.loads(ac)
        var = next(c["variant"] for c in cand if c["display"].startswith(rs+" "))
    except Exception:
        fails.append((rs,"autocomplete")); continue
    vj = curl(f"https://r12.finngen.fi/api/variant/{var}")
    try:
        d = json.loads(vj)
    except Exception:
        fails.append((rs,"variant")); continue
    ref = d["variant"]["ref"].upper(); alt = d["variant"]["alt"].upper()
    if ea == alt: flip = 1
    elif ea == ref: flip = -1
    else:
        fails.append((rs, f"allele_mismatch ref={ref} alt={alt} ea={ea}")); continue
    for res in d["results"]:
        if res["phenocode"] in targets and res.get("beta") is not None:
            rows.append({
                "SNP": rs, "effect_allele": ea, "other_allele": oa,
                "eaf_exposure": r["eaf.exposure"], "beta_exposure": r["beta.exposure"], "se_exposure": r["se.exposure"],
                "phenocode": res["phenocode"], "n_case": res["n_case"], "n_control": res["n_control"],
                "beta_out_alt": res["beta"], "pval": res["pval"], "maf": res.get("maf"),
                "beta_outcome": flip*res["beta"], "flip": flip,
                "beta_sign_check_alt_beta": res["beta"],
            })
    time.sleep(0.3)

df = pd.DataFrame(rows)
df.to_csv("spectrum_raw.csv", index=False)
print("rows:", len(df), "snps:", df.SNP.nunique() if len(df) else 0)
print("endpoints:", sorted(df.phenocode.unique()) if len(df) else [])
per = df.groupby("phenocode")["SNP"].nunique()
print(per)
print("FAILS:", fails)

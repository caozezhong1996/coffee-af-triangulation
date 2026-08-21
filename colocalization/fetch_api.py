import json, os, time, threading
import pandas as pd
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

TOKEN = os.environ["OPGW"]
URL = "https://api.opengwas.io/api/associations"
OUT = r"data\api_associations.csv"
LOCK = threading.Lock()

df = pd.read_csv(r"data\windows_thinned.tsv", sep="\t")
rsids = sorted(set(df["rsids"]))
CH = 32
jobs = [(ds, rsids[i:i+CH]) for ds in ("ukb-b-5237",) for i in range(0, len(rsids), CH)]

done = pd.DataFrame()
if os.path.exists(OUT) and os.path.getsize(OUT) > 10:
    done = pd.read_csv(OUT)
have = set(zip(done["id"], done["rsid"])) if len(done) else set()
todo = [(ds, c) for ds, c in jobs if not all((ds, r) in have for r in c)]
print(f"total chunks {len(jobs)}, remaining {len(todo)}, rows already {len(done)}", flush=True)

state = {"new": [], "n": 0}

def fetch(job):
    ds, chunk = job
    body = json.dumps({"variant": chunk, "id": [ds]}).encode()
    for attempt in range(4):
        try:
            req = urllib.request.Request(URL, data=body, method="POST",
                headers={"Authorization": f"Bearer {TOKEN}", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                rows = json.loads(r.read().decode())
            with LOCK:
                state["new"].extend(rows)
                state["n"] += 1
            return
        except urllib.error.HTTPError as e:
            time.sleep(6 * (attempt + 1))
        except Exception:
            time.sleep(3 * (attempt + 1))

t0 = time.time()
BATCH = 40
for bstart in range(0, len(todo), BATCH):
    batch = todo[bstart:bstart+BATCH]
    with ThreadPoolExecutor(max_workers=4) as ex:
        list(ex.map(fetch, batch))
    if state["new"]:
        done = pd.concat([done, pd.DataFrame(state["new"])], ignore_index=True).drop_duplicates(subset=["id", "rsid"])
        done.to_csv(OUT, index=False)
        state["new"] = []
    print(f"progress {bstart+len(batch)}/{len(todo)} chunks, rows {len(done)}, elapsed {time.time()-t0:.0f}s", flush=True)

print("done. total rows:", len(done), flush=True)

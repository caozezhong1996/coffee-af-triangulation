import json, os, time, threading
import pandas as pd
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor

TOKEN = os.environ["OPGW"]
URL = "https://api.opengwas.io/api/associations"
OUT = r"data\ahr_api_associations.csv"
LOCK = threading.Lock()

rsids = [l.strip() for l in open(r"data\ahr_rsids.txt") if l.strip()]
CH = 64
jobs = [rsids[i:i+CH] for i in range(0, len(rsids), CH)]

done = pd.DataFrame()
if os.path.exists(OUT) and os.path.getsize(OUT) > 10:
    done = pd.read_csv(OUT)
have = set(done["rsid"]) if len(done) else set()
todo = [c for c in jobs if not all(r in have for r in c)]
print(f"total chunks {len(jobs)}, remaining {len(todo)}, rows already {len(done)}", flush=True)

state = {"new": [], "n": 0}
def fetch(chunk):
    body = json.dumps({"variant": chunk, "id": ["ukb-b-5237"]}).encode()
    for attempt in range(4):
        try:
            req = urllib.request.Request(URL, data=body, method="POST",
                headers={"Authorization": f"Bearer {TOKEN}", "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                rows = json.loads(r.read().decode())
            with LOCK:
                state["new"].extend(rows); state["n"] += 1
            return
        except urllib.error.HTTPError:
            time.sleep(6*(attempt+1))
        except Exception:
            time.sleep(3*(attempt+1))

t0=time.time()
for bstart in range(0, len(todo), 40):
    batch = todo[bstart:bstart+40]
    with ThreadPoolExecutor(max_workers=4) as ex:
        list(ex.map(fetch, batch))
    if state["new"]:
        done = pd.concat([done, pd.DataFrame(state["new"])], ignore_index=True).drop_duplicates(subset=["id","rsid"])
        done.to_csv(OUT, index=False); state["new"]=[]
    print(f"progress {bstart+len(batch)}/{len(todo)}, rows {len(done)}, elapsed {time.time()-t0:.0f}s", flush=True)
print("done. total rows:", len(done), flush=True)

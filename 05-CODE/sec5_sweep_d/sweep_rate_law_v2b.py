# [oid:sweepD] NOT RIGOROUS - v2b: supplement, degrees 17..20 (same engine as v2).
import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from sweep_rate_law_v2 import run_degree

DEGREES = ["17", "18", "19", "20"]
OUT = r"sweep_v2b_results.json"
LOG = r"sweep_v2b_progress.log"

def main():
    t0 = time.time()
    log = open(LOG, "a", buffering=1)
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")
    w(f"SWEEP v2b START degrees={DEGREES} workers=4")
    results = []
    with ProcessPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_degree, ds): ds for ds in DEGREES}
        for f in as_completed(futs):
            rr_ = f.result()
            results.append(rr_)
            if "error" in rr_:
                w(f"d={rr_['d']}: ERROR {rr_['error']}")
            else:
                nd = {k: v.get("median", "x") for k, v in rr_["nodes"].items()}
                w(f"d={rr_['d']}: delta={rr_['delta']:.5f} node_medians={nd} "
                  f"-> {rr_['verdict']} median={rr_['median_over_nodes']:.4f}"
                  + (f" implied lambda2={rr_['implied_lambda2']}" if rr_.get("implied_lambda2") else ""))
            json.dump(sorted(results, key=lambda z: float(z["d"])), open(OUT, "w"), indent=1)
    w(f"SWEEP v2b DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

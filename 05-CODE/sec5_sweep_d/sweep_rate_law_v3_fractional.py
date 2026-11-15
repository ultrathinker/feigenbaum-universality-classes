# [oid:sweepD] NOT RIGOROUS - v3: interesting fractional degrees, separate run (v2c engine).
# Zones: approach to d->1+ (1.3, 1.5, 1.7); left of crossover (2.2, 2.4 - delta-regime expected,
# lambda2(d) crosses zero between d=2 and d=3); crossover window fine grid (2.6, e, 2.8, 2.9);
# pi for fun (should sit smoothly between d=3 and d=3.5 on the lambda2 curve).
import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from sweep_rate_law_v2c import run_degree

DEGREES = ["1.3", "1.5", "1.7", "2.2", "2.4", "2.6",
           "2.71828182845904523536", "2.8", "2.9", "3.14159265358979323846"]
OUT = r"sweep_v3_results.json"
LOG = r"sweep_v3_progress.log"

def label(ds):
    if ds.startswith("2.71828"): return "e"
    if ds.startswith("3.14159"): return "pi"
    return ds

def main():
    t0 = time.time()
    log = open(LOG, "a", buffering=1)
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")
    w(f"SWEEP v3 START degrees={[label(x) for x in DEGREES]} workers=10")
    results = []
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(run_degree, ds): ds for ds in DEGREES}
        for f in as_completed(futs):
            rr_ = f.result()
            rr_["label"] = label(rr_["d"])
            results.append(rr_)
            if "error" in rr_:
                w(f"d={rr_['label']}: ERROR {rr_['error']}")
            else:
                nd = {k: v.get("median", "x") for k, v in rr_["nodes"].items()}
                w(f"d={rr_['label']}: delta={rr_['delta']:.5f} node_medians={nd} "
                  f"-> {rr_['verdict']} median={rr_['median_over_nodes']:.4f}"
                  + (f" implied lambda2={rr_['implied_lambda2']}" if rr_.get("implied_lambda2") else ""))
            json.dump(sorted(results, key=lambda z: float(z["d"])), open(OUT, "w"), indent=1)
    w(f"SWEEP v3 DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

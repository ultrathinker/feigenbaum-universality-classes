# [oid:sweepD] NOT RIGOROUS - v2: rate-law sweep, 5 nodes per degree, degrees up to 16.
# Per degree d (family 1 - a*|x|^d, real axis):
#   cascade a_1..a_16 (own Newton, impostor rejection, min-candidate bootstrap)
#   delta(d) via Richardson on gap ratios; a_F, L, A
#   nodes u in {-1 (calibration, exact limit), +2.5, +6, +12, +20}:
#     track superstable a*(m), m=5..12, by continuation in the scaling coordinate;
#     rate from successive differences: r(m) = (u_m - u_{m+1})/(u_{m+1} - u_{m+2})  [u_inf-free]
# Verdict per degree: median rate vs delta -> regime; implied lambda2 = 1/rate in lambda2-regime.
import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

DEGREES = ["2.5", "3", "3.5", "5", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16"]
NODES = ["-1", "2.5", "6", "12", "20"]
OUT = r"sweep_v2_results.json"
LOG = r"sweep_v2_progress.log"

def run_degree(dstr):
    from mpmath import mp, mpf, fabs, sign
    mp.dps = 80
    d = mpf(dstr)
    t0 = time.time()

    def step(x, dx, a):
        ax = fabs(x)
        if ax == 0:
            return mpf(1), dx * 0
        p = ax ** d
        q = d * ax ** (d - 1) * sign(x)
        return 1 - a * p, -p - a * q * dx

    def superstable(P, a0):
        a = mpf(a0)
        for it in range(40):
            x = mpf(0); dx = mpf(0)
            for i in range(P):
                x, dx = step(x, dx, a)
            if dx == 0:
                return None, mpf(0)
            stp = x / dx
            a = a - stp
            if a < mpf("0.5") or a > 3:
                return None, mpf(0)
            if fabs(stp) < mpf(10) ** (-65):
                break
        x = mpf(0)
        for i in range(P):
            x, _ = step(x, mpf(0), a)
        if fabs(x) > mpf(10) ** (-40):
            return None, mpf(0)
        x = mpf(0); mins = []; nxt = 2
        for i in range(1, P + 1):
            x, _ = step(x, mpf(0), a)
            if i == nxt:
                mins.append(fabs(x)); nxt *= 2
        return a, (min(mins[:-1]) if len(mins) > 1 else mpf(1))

    # ---- cascade ----
    KMAX = 16
    a = {1: mpf(1)}
    d_est = mpf(6)
    ratio_hist = []
    for k in range(2, KMAX + 1):
        if k == 2:
            seeds = [mpf(s) / 100 for s in range(102, 199, 2)]
            lo, hi = a[1] + mpf("0.01"), mpf(3)
        else:
            eg = (a[k - 1] - a[k - 2]) / d_est
            fr = [mpf(f) / 8 for f in (8, 9, 7, 10, 6, 12, 4, 14, 2, 16, 11, 5)]
            seeds = [a[k - 1] + eg * f for f in fr]
            lo, hi = a[k - 1] + eg / 10, a[k - 1] + 3 * (a[k - 1] - a[k - 2])
        cands = []
        for s in seeds:
            cand, minsep = superstable(2 ** k, s)
            if cand is not None and minsep > mpf(10) ** (-8) and lo < cand < hi:
                cands.append(cand)
                if len(cands) >= 3:
                    break
        if not cands:
            return {"d": dstr, "error": f"cascade stuck at k={k}"}
        a[k] = min(cands)
        if k >= 3:
            d_est = (a[k - 1] - a[k - 2]) / (a[k] - a[k - 1])
            ratio_hist.append(d_est)

    d_top, d_prev = ratio_hist[-1], ratio_hist[-2]
    rr = 1 / d_top
    delta = d_top + (d_top - d_prev) * rr / (1 - rr)
    r1 = 1 / delta
    aF = a[KMAX] + (a[KMAX] - a[KMAX - 1]) * r1 / (1 - r1)
    L = (a[KMAX] - aF) * delta ** (KMAX - 1)
    A = -1 / L
    fdelta = float(delta)

    # ---- nodes ----
    node_out = {}
    rates_all = []
    for ustr in NODES:
        u = mpf(ustr)
        um = {}
        u_track = u
        ok = True
        for m in range(5, 13):
            seed = aF + u_track / (A * delta ** m)
            ast, minsep = superstable(2 ** (m + 1), seed)
            if ast is None or minsep <= mpf(10) ** (-8):
                ok = False; break
            um[m] = A * (ast - aF) * delta ** m
            u_track = um[m]  # continuation in scaling coordinate
        if not ok or len(um) < 8:
            node_out[ustr] = {"status": "failed"}
            continue
        diffs = {m: um[m] - um[m + 1] for m in range(5, 12)}
        rts = []
        for m in range(5, 11):
            if diffs[m + 1] != 0:
                rts.append(float(fabs(diffs[m] / diffs[m + 1])))
        med = sorted(rts)[len(rts) // 2]
        node_out[ustr] = {"status": "ok", "rates": [round(x, 4) for x in rts], "median": round(med, 4)}
        rates_all.append(med)
        if ustr == "-1":  # bonus: exact-limit ratios
            e = {m: um[m] + 1 for m in um}
            er = [float(fabs(e[m] / e[m + 1])) for m in range(5, 11) if e[m + 1] != 0]
            node_out[ustr]["exact_limit_rates"] = [round(x, 4) for x in er]

    if not rates_all:
        return {"d": dstr, "delta": fdelta, "error": "all nodes failed"}
    gmed = sorted(rates_all)[len(rates_all) // 2]
    verdict = "delta-regime" if abs(gmed - fdelta) / fdelta < 0.05 else (
        "lambda2-regime" if gmed < 0.9 * fdelta else "unclear")
    return {"d": dstr, "delta": fdelta, "nodes": node_out, "median_over_nodes": gmed,
            "implied_lambda2": (round(1.0 / gmed, 5) if verdict == "lambda2-regime" else None),
            "verdict": verdict, "elapsed_s": round(time.time() - t0, 1)}

def main():
    t0 = time.time()
    log = open(LOG, "a", buffering=1)
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")
    w(f"SWEEP v2 START degrees={DEGREES} nodes={NODES} workers=14")
    results = []
    with ProcessPoolExecutor(max_workers=14) as ex:
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
    w(f"SWEEP v2 DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

# [oid:sweepD] NOT RIGOROUS - v2c: same engine as v2, but low-k cascade bootstrap via
# sign-scan + bisection on x_{2^k}(a) (robust for near-flat-top high degrees), Newton polish.
# Degrees: the ones v2/v2b could not bootstrap.
import json, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

DEGREES = ["13", "14", "15", "16", "17", "19", "20"]
NODES = ["-1", "2.5", "6", "12", "20"]
OUT = r"sweep_v2c_results.json"
LOG = r"sweep_v2c_progress.log"

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

    def orbit_end(P, a):
        x = mpf(0)
        for i in range(P):
            x, _ = step(x, mpf(0), a)
        return x

    def newton_polish(P, a0, iters=6):
        a = mpf(a0)
        for it in range(iters):
            x = mpf(0); dx = mpf(0)
            for i in range(P):
                x, dx = step(x, dx, a)
            if dx == 0:
                break
            a = a - x / dx
        return a

    def check(P, a):
        x = mpf(0); mins = []; nxt = 2
        for i in range(1, P + 1):
            x, _ = step(x, mpf(0), a)
            if i == nxt:
                mins.append(fabs(x)); nxt *= 2
        res = fabs(x)
        minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)
        return res, minsep

    def bisect_root(P, alo, ahi, flo):
        for it in range(140):
            amid = (alo + ahi) / 2
            fm = orbit_end(P, amid)
            if fm == 0:
                return amid
            if (fm > 0) == (flo > 0):
                alo, flo = amid, fm
            else:
                ahi = amid
            if ahi - alo < mpf(10) ** (-40):
                break
        return (alo + ahi) / 2

    def superstable_newton(P, a0):
        a = newton_polish(P, a0, 40)
        res, minsep = check(P, a)
        if res > mpf(10) ** (-40) or a < mpf("0.5") or a > 3:
            return None, mpf(0)
        return a, minsep

    # ---- cascade: k=2..6 by scan+bisection, k>=7 by Newton continuation ----
    KMAX = 16
    a = {1: mpf(1)}
    d_est = mpf(6) if d < 11 else mpf("12.34") + mpf("0.49") * (d - 10)
    ratio_hist = []
    for k in range(2, KMAX + 1):
        P = 2 ** k
        if k <= 6:
            if k == 2:
                start, stop = mpf("1.02"), mpf("2.95")
                h = mpf("0.005")
            else:
                gap_guess = (a[k - 1] - a[k - 2]) / d_est
                start = a[k - 1] + gap_guess / 20
                stop = a[k - 1] + 3 * (a[k - 1] - a[k - 2])
                h = gap_guess / 10
            aprev = start; fprev = orbit_end(P, aprev)
            root = None
            aa = start + h
            while aa < stop:
                ff = orbit_end(P, aa)
                if ff == 0:
                    root = aa; break
                if (ff > 0) != (fprev > 0):
                    root = bisect_root(P, aprev, aa, fprev); break
                aprev, fprev = aa, ff
                aa += h
            if root is None:
                return {"d": dstr, "error": f"bisection bootstrap found no root at k={k}"}
            ak = newton_polish(P, root)
            res, minsep = check(P, ak)
            if res > mpf(10) ** (-40) or minsep <= mpf(10) ** (-10):
                return {"d": dstr, "error": f"polish/period check failed at k={k}"}
            a[k] = ak
        else:
            eg = (a[k - 1] - a[k - 2]) / d_est
            cands = []
            for f in (1, mpf(9) / 8, mpf(7) / 8, mpf(5) / 4, mpf(3) / 4, mpf(3) / 2, mpf(1) / 2):
                cand, minsep = superstable_newton(P, a[k - 1] + eg * f)
                if cand is not None and minsep > mpf(10) ** (-10) and a[k - 1] + eg / 10 < cand < a[k - 1] + 3 * (a[k - 1] - a[k - 2]):
                    cands.append(cand)
                    if len(cands) >= 2:
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

    # ---- nodes (same as v2) ----
    node_out = {}
    rates_all = []
    for ustr in NODES:
        u = mpf(ustr)
        um = {}
        u_track = u
        ok = True
        for m in range(5, 13):
            seed = aF + u_track / (A * delta ** m)
            ast, minsep = superstable_newton(2 ** (m + 1), seed)
            if ast is None or minsep <= mpf(10) ** (-8):
                ok = False; break
            um[m] = A * (ast - aF) * delta ** m
            u_track = um[m]
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
    w(f"SWEEP v2c START degrees={DEGREES} workers=7")
    results = []
    with ProcessPoolExecutor(max_workers=7) as ex:
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
    w(f"SWEEP v2c DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

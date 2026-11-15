# [oid:sweepD] NOT RIGOROUS - surface-level sweep of the min(delta, 1/|lambda2|) rate law
# across critical-point orders d (family f_a(x) = 1 - a*|x|^d, real axis, calibration node only).
# Per degree: own Newton cascade a_1..a_14, delta(d) from gap ratios, Richardson a_F,
# signed e(m) = A(a_{m+1}-a_F)delta^m + 1, ratio band m=5..10 -> verdict:
#   ratio ~= delta       -> delta-regime (like d=2)
#   ratio clearly below  -> lambda2-regime, implied lambda2(d) = 1/ratio (anchors: 0.2918 @ d=4, 0.3711 @ d=6)
# Circularity guard: m>=12 excluded (a13,a14 build a_F); top-band drift reported as noise flag.
import json, sys, time
from concurrent.futures import ProcessPoolExecutor

DEGREES = ["2.5", "3", "3.5", "5", "7", "8", "9", "10"]
OUT = r"sweep_results.json"
LOG = r"sweep_progress.log"

def run_degree(dstr):
    from mpmath import mp, mpf, fabs, sign
    mp.dps = 80
    d = mpf(dstr)
    t0 = time.time()

    def step(x, dx, a):
        ax = fabs(x)
        p = ax ** d if ax > 0 else mpf(0)          # |x|^d
        q = d * ax ** (d - 1) * sign(x) if ax > 0 else mpf(0)  # d/dx |x|^d
        return 1 - a * p, -p - a * q * dx

    def superstable_k(k, a0):
        P = 2 ** k
        a = mpf(a0)
        for it in range(40):
            x = mpf(0); dx = mpf(0)
            for i in range(P):
                x, dx = step(x, dx, a)
            if dx == 0:
                return None, mpf(0)
            stp = x / dx
            a = a - stp
            if a < mpf("0.5") or a > 3:   # wandered off the family's interesting range
                return None, mpf(0)
            if fabs(stp) < mpf(10) ** (-65):
                break
        # residual + exact-period check
        x = mpf(0)
        for i in range(P):
            x, _ = step(x, mpf(0), a)
        res = fabs(x)
        x = mpf(0); mins = []; nxt = 2
        for i in range(1, P + 1):
            x, _ = step(x, mpf(0), a)
            if i == nxt:
                mins.append(fabs(x)); nxt *= 2
        minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)
        return a, (minsep if res < mpf(10) ** (-40) else mpf(0))

    KMAX = 16
    a = {1: mpf(1)}  # 0 -> 1 -> 1-a: superstable period 2 at a=1 for every d
    d_est = mpf(6)
    ratio_hist = []
    for k in range(2, KMAX + 1):
        if k == 2:
            seeds = [mpf(s) / 100 for s in range(102, 199, 2)]
            lo, hi = a[1] + mpf("0.01"), mpf("3")
        else:
            eg = (a[k - 1] - a[k - 2]) / d_est
            base = a[k - 1] + eg
            seeds = [base, base + eg / 4, base - eg / 4, base + eg / 2, base - eg / 2, base + eg, a[k - 1] + eg / 2]
            lo, hi = a[k - 1] + eg / 10, a[k - 1] + 3 * (a[k - 1] - a[k - 2])
        cands = []
        for s in seeds:
            cand, minsep = superstable_k(k, s)
            if cand is not None and minsep > mpf(10) ** (-8) and lo < cand < hi:
                cands.append(cand)
        if not cands:
            return {"d": dstr, "error": f"cascade stuck at k={k}"}
        # main cascade point = smallest valid candidate above a_{k-1}
        a[k] = min(cands)
        if k >= 3:
            d_est = (a[k - 1] - a[k - 2]) / (a[k] - a[k - 1])
            ratio_hist.append(d_est)

    # delta: one Richardson step on the gap-ratio sequence (ratios converge geometrically)
    d_top, d_prev = ratio_hist[-1], ratio_hist[-2]
    rr = 1 / d_top
    delta = d_top + (d_top - d_prev) * rr / (1 - rr)
    r = 1 / delta
    aF = a[KMAX] + (a[KMAX] - a[KMAX - 1]) * r / (1 - r)
    L = (a[KMAX] - aF) * delta ** (KMAX - 1)
    A = -1 / L
    # top ~3 levels are interpolated by construction (delta, aF, L use a14..a16) -> band ends at m=11
    e = {m: A * (a[m + 1] - aF) * delta ** m + 1 for m in range(5, 13)}
    ratios = {m: float(e[m] / e[m + 1]) for m in range(5, 12) if e[m + 1] != 0}
    band = [ratios[m] for m in sorted(ratios)]
    core = band[:6]  # m=5/6 .. 10/11
    mid = sorted(core)[len(core) // 2]
    fdelta = float(delta)
    verdict = "delta-regime" if abs(mid - fdelta) / fdelta < 0.05 else (
        "lambda2-regime" if mid < 0.9 * fdelta else "unclear")
    return {"d": dstr, "delta": fdelta, "delta_conv": abs(float(d_top - d_prev)), "aF": mp.nstr(aF, 25),
            "ratios_m5_11": band, "ratio_median": mid,
            "implied_lambda2": (1.0 / mid if verdict == "lambda2-regime" else None),
            "verdict": verdict, "elapsed_s": round(time.time() - t0, 1)}

def main():
    t0 = time.time()
    log = open(LOG, "a", buffering=1)
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")
    w(f"SWEEP START degrees={DEGREES} workers=8")
    results = []
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = {ex.submit(run_degree, ds): ds for ds in DEGREES}
        for f in [x for x in futs]:
            pass
        from concurrent.futures import as_completed
        for f in as_completed(futs):
            rr = f.result()
            results.append(rr)
            w(f"d={rr['d']}: " + (f"ERROR {rr['error']}" if "error" in rr else
              f"delta={rr['delta']:.5f} ratios={['%.3f' % x for x in rr['ratios_m5_11']]} "
              f"median={rr['ratio_median']:.4f} -> {rr['verdict']}"
              + (f" implied lambda2={rr['implied_lambda2']:.5f}" if rr.get("implied_lambda2") else "")))
            json.dump(sorted(results, key=lambda z: float(z["d"])), open(OUT, "w"), indent=1)
    w(f"SWEEP DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()

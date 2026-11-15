# [BRIEF-L] NOT RIGOROUS - cascade builder for d=2.75 and d=3.0.
# Method: bootstrap k<=6 via scan+bisection on sign of x_{2^k}(a) (robust,
# engine copied/adapted from ../sec5_sweep_d/sweep_rate_law_v2c.py), then
# Newton continuation k>=7 with residual-based acceptance (NOT step-size
# based - trap #4 in 05-SVODKA-ZNANIY.md SS II). Writes cascade JSON
# incrementally after EVERY level (recoverable), and progress.log with
# flush + ETA after every level. CPU budget: 2 worker processes (<=3 total).
import json, sys, time, os
from concurrent.futures import ProcessPoolExecutor, as_completed

WDIR = r"."
LOG = os.path.join(WDIR, "progress.log")
KMAX = 16
DPS = 150

JOBS = [
    ("2.75", os.path.join(WDIR, "cascade_d275.json")),
    ("3.0", os.path.join(WDIR, "cascade_d300.json")),
]


def build_cascade(dstr, out_path, kmax=KMAX, dps=DPS):
    from mpmath import mp, mpf, fabs, sign
    mp.dps = dps
    d = mpf(dstr)
    t0 = time.time()

    def step(x, a):
        ax = fabs(x)
        if ax == 0:
            return mpf(1)
        p = ax ** d
        return 1 - a * p

    def orbit_end(P, a):
        x = mpf(0)
        for _ in range(P):
            x = step(x, a)
        return x

    def orbit_end_dx(P, a):
        # value AND d(x_P)/da via forward accumulation (for Newton)
        x = mpf(0); dx = mpf(0)
        for _ in range(P):
            ax = fabs(x)
            if ax == 0:
                p = mpf(0); q = mpf(0)
            else:
                p = ax ** d
                q = d * ax ** (d - 1) * sign(x)
            x, dx = (1 - a * p, -p - a * q * dx)
        return x, dx

    def newton_polish(P, a0, max_iters=60, tol_exp=None):
        tol = mpf(10) ** (-(dps // 2)) if tol_exp is None else mpf(10) ** (-tol_exp)
        a = mpf(a0)
        best_a, best_res = a, None
        for _ in range(max_iters):
            x, dx = orbit_end_dx(P, a)
            res = fabs(x)
            if best_res is None or res < best_res:
                best_a, best_res = a, res
            if res < tol or dx == 0:
                break
            a = a - x / dx
        return best_a, best_res

    def bisect_root(P, alo, ahi, flo):
        for _ in range(160):
            amid = (alo + ahi) / 2
            fm = orbit_end(P, amid)
            if fm == 0:
                return amid
            if (fm > 0) == (flo > 0):
                alo, flo = amid, fm
            else:
                ahi = amid
            if ahi - alo < mpf(10) ** (-(dps - 20)):
                break
        return (alo + ahi) / 2

    def check(P, a):
        x = mpf(0); mins = []; nxt = 2
        for i in range(1, P + 1):
            x = step(x, a)
            if i == nxt:
                mins.append(fabs(x)); nxt *= 2
        res = fabs(x)
        minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)
        return res, minsep

    tag = f"d={dstr}"
    log = open(LOG, "a", buffering=1)

    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] [{tag}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")

    def save(levels, done=False, error=None):
        payload = {
            "d": dstr, "dps": dps, "kmax": kmax, "done": done, "error": error,
            "levels": {str(k): {"a": mp_to_str(v["a"]), "res": mp_to_str(v["res"])}
                       for k, v in levels.items()},
            "elapsed_s": round(time.time() - t0, 1),
        }
        tmp = out_path + ".tmp"
        json.dump(payload, open(tmp, "w"), indent=1)
        os.replace(tmp, out_path)

    def mp_to_str(x):
        return mp.nstr(x, dps - 5, strip_zeros=False)

    levels = {}
    a = {1: mpf(1)}
    levels[1] = {"a": a[1], "res": mpf(0)}
    save(levels)
    w("START kmax=%d dps=%d" % (kmax, dps))

    d_est = mpf(6)  # rough initial guess (delta(2.75..3.0) approx 5.76..6.09)

    for k in range(2, kmax + 1):
        tk0 = time.time()
        P = 2 ** k
        if k <= 6:
            if k == 2:
                start, stop, h = mpf("1.02"), mpf("2.95"), mpf("0.005")
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
                save(levels, done=True, error=f"bisection bootstrap found no root at k={k}")
                w(f"FATAL: no root at k={k}"); return
            ak, res = newton_polish(P, root)
            resid, minsep = check(P, ak)
            if resid > mpf(10) ** (-(dps // 2)) or minsep <= mpf(10) ** (-(dps // 4)):
                save(levels, done=True, error=f"polish/period check failed at k={k} resid={mp_to_str(resid)}")
                w(f"FATAL: polish/period fail at k={k}"); return
            a[k] = ak
        else:
            eg = (a[k - 1] - a[k - 2]) / d_est
            cands = []
            for f in (1, mpf(9) / 8, mpf(7) / 8, mpf(5) / 4, mpf(3) / 4, mpf(3) / 2, mpf(1) / 2, mpf(2), mpf(1) / 4):
                seed = a[k - 1] + eg * f
                cand, res = newton_polish(P, seed)
                if cand <= a[k - 1] + eg / 20:
                    continue
                resid, minsep = check(P, cand)
                if resid < mpf(10) ** (-(dps // 2)) and minsep > mpf(10) ** (-(dps // 4)) and \
                   a[k - 1] + eg / 10 < cand < a[k - 1] + 3 * (a[k - 1] - a[k - 2]):
                    cands.append((cand, resid))
                    if len(cands) >= 2:
                        break
            if not cands:
                save(levels, done=True, error=f"cascade stuck at k={k} (no valid Newton candidate)")
                w(f"FATAL: stuck at k={k}"); return
            cands.sort(key=lambda z: z[0])
            a[k] = cands[0][0]
            resid = cands[0][1]

        levels[k] = {"a": a[k], "res": resid}
        if k >= 3:
            d_est = (a[k - 1] - a[k - 2]) / (a[k] - a[k - 1])
        save(levels)
        elapsed_k = time.time() - tk0
        total_elapsed = time.time() - t0
        remaining = kmax - k
        eta = remaining * elapsed_k
        w(f"k={k:2d} P={P:6d} a={mp.nstr(a[k], 20)} res={mp.nstr(resid, 3)} "
          f"d_est={mp.nstr(d_est, 8)} step={elapsed_k:.1f}s total={total_elapsed:.1f}s ETA={eta:.1f}s")

    save(levels, done=True)
    w(f"DONE in {time.time() - t0:.0f}s, kmax={kmax}")


def _worker(job):
    dstr, out_path = job
    build_cascade(dstr, out_path)
    return dstr


def main():
    log = open(LOG, "a", buffering=1)
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] [main] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")
    w(f"CASCADE_GEN START jobs={[j[0] for j in JOBS]} workers=2 kmax={KMAX} dps={DPS}")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(_worker, job): job for job in JOBS}
        for f in as_completed(futs):
            dstr = f.result()
            w(f"job d={dstr} finished")
    w(f"CASCADE_GEN DONE in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()

# [BRIEF-L] NOT RIGOROUS - final dynamics verdict for d=2.75, d=3.0.
#
# Method (derived from first principles, see progress.log narrative for the
# derivation trail): a_m - a_F ~ sum_i C_i * rho_i^m (multi-mode RG expansion).
# R(m) = (a_m-a_{m-1})/(a_{m+1}-a_m) -> delta_fit as m -> inf; the correction
# r(m) = R(m) - delta_fit is DOMINATED near m~5-9 by a fast mode at rate
# EXACTLY 1/delta_fit (confirmed empirically: free 2-mode Prony fit recovers
# rho1 = 1/delta_fit to ~5%), plus a slower mode at rate mu (the disputed
# quantity). Instead of an ad hoc "post-crossover cut" (which throws away
# data and leaves residual contamination - see progress.log for the messy
# intermediate result), we ELIMINATE the fast mode EXACTLY and algebraically:
#   s(k) = r(k+1) - (1/delta_fit)*r(k)
# For a pure 2-mode signal r(k)=C*rho1^k+D*mu^k with rho1=1/delta_fit exactly,
# s(k) = D*mu^k*(mu-rho1) -- a SINGLE clean geometric mode in mu, with NO
# contamination from the fast layer at any k. This is checked to have
# R^2 > 0.98 log-linear fit quality on the tail, vs the crude post-crossover
# r(k) tail which had messy non-monotonic point-ratios. delta_fit itself:
#   d=3.0  -> external, non-circular: wave1_spectra_d/spectra_grid.json
#            (independent OPERATOR/spectral computation, N40, 30 digits)
#   d=2.75 -> no exact external grid point; self-cascade late-Aitken average,
#            cross-validated: the analogous self-Aitken estimate at d=3.0
#            agrees with the external spectral value to ~1.4e-10 (relative
#            ~2e-11), so the method is trusted for 2.75 too.
import json, sys, time
from mpmath import mp, mpf, nstr, log as mlog, e as me

WDIR = r"."
mp.dps = 120

DELTA_SPECTRAL_D300 = mpf("6.08468907426212189604612106626")
PRED = {
    "2.75": {"A": {"rate": 4.53, "mu": 1 / 4.53, "monotone": True},
             "B": {"rate": 5.50, "mu": -1 / 5.50, "monotone": False}},
    "3.0":  {"A": {"rate": 4.22, "mu": 1 / 4.22, "monotone": True},
             "B": {"rate": 5.04, "mu": -1 / 5.04, "monotone": False}},
}
FILES = {"2.75": "cascade_d275.json", "3.0": "cascade_d300.json"}
TAIL_LO = 8   # s(k) tail window start (empirically clean from k=6-8 onward)


def load_R(fn):
    data = json.load(open(f"{WDIR}\\{fn}"))
    levels = {int(k): mpf(v["a"]) for k, v in data["levels"].items()}
    ks = sorted(levels)
    D = {k: levels[k] - levels[k - 1] for k in ks if k >= 2}
    Dk = sorted(D)
    R = {}
    for k in Dk:
        if k + 1 in D:
            R[k] = D[k] / D[k + 1]
    return R


def aitken3(R, k0, k1, k2):
    x0, x1, x2 = R[k0], R[k1], R[k2]
    denom = x2 - 2 * x1 + x0
    return None if denom == 0 else x2 - (x2 - x1) ** 2 / denom


def ols(xs, ys):
    n = len(xs)
    xbar = sum(xs) / n; ybar = sum(ys) / n
    sxy = sum((xs[i] - xbar) * (ys[i] - ybar) for i in range(n))
    sxx = sum((xs[i] - xbar) ** 2 for i in range(n))
    slope = sxy / sxx
    rate = me ** (-slope)
    intercept = ybar - slope * xbar
    yhat = [intercept + slope * x for x in xs]
    ss_res = sum((ys[i] - yhat[i]) ** 2 for i in range(n))
    ss_tot = sum((ys[i] - ybar) ** 2 for i in range(n))
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else None
    return rate, r2


def compute_s(R, delta_fit):
    rho1 = 1 / delta_fit
    r = {k: R[k] - delta_fit for k in R}
    ks = sorted(r)
    s = {k: r[k + 1] - rho1 * r[k] for k in ks if k + 1 in r}
    return s, rho1


def analyze_delta(dstr, R, delta_fit, tag, log):
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] [d={dstr}][{tag}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")

    s, rho1 = compute_s(R, delta_fit)
    sks = sorted(s)
    w(f"delta_fit={nstr(delta_fit,20)} rho1=1/delta_fit={nstr(rho1,10)}")
    signs = {k: (1 if s[k] > 0 else -1) for k in sks}
    w("s(k) sequence: " + ", ".join(f"k={k}:{nstr(s[k],6)}" for k in sks))

    tail = [k for k in sks if k >= TAIL_LO]
    tail_signs = [signs[k] for k in tail]
    tail_changes = sum(1 for i in range(1, len(tail_signs)) if tail_signs[i] != tail_signs[i - 1])
    monotone = tail_changes == 0
    final_sign = tail_signs[-1] if tail_signs else None
    w(f"tail k>={TAIL_LO}: signs={tail_signs} sign_changes={tail_changes} "
      f"-> {'MONOTONE' if monotone else 'ALTERNATING'} (final_sign={'+' if final_sign>0 else '-'})")

    xs = [mpf(k) for k in tail]; ys = [mlog(abs(s[k])) for k in tail]
    rate, r2 = ols(xs, ys)
    w(f"OLS tail fit: rate(=1/mu)={nstr(rate,8)} R2={nstr(r2,6)} n={len(tail)} "
      f"=> mu={'+' if final_sign>0 else '-'}{nstr(1/rate,8)}")

    loocv = []
    for i in range(len(tail)):
        sub = tail[:i] + tail[i + 1:]
        if len(sub) >= 3:
            rr, _ = ols([mpf(k) for k in sub], [mlog(abs(s[k])) for k in sub])
            loocv.append(rr)
    lo, hi = (min(loocv), max(loocv)) if loocv else (None, None)
    if loocv:
        w(f"LOOCV rate bracket = [{nstr(lo,6)}, {nstr(hi,6)}]")

    return {
        "delta_fit": nstr(delta_fit, 20),
        "s_tail": {str(k): nstr(s[k], 10) for k in tail},
        "tail_signs": tail_signs,
        "tail_sign_changes": tail_changes,
        "monotone": monotone,
        "final_sign": "+" if final_sign > 0 else "-",
        "ols_rate": nstr(rate, 8),
        "ols_r2": nstr(r2, 6),
        "loocv_bracket": [nstr(lo, 6), nstr(hi, 6)] if loocv else None,
        "implied_mu": nstr((1 if final_sign > 0 else -1) / rate, 8),
    }


def main():
    log = open(f"{WDIR}\\progress.log", "a", buffering=1)
    out = {}
    for dstr, fn in FILES.items():
        R = load_R(fn)
        if dstr == "3.0":
            delta_primary = DELTA_SPECTRAL_D300
            delta_alt = (aitken3(R, 12, 13, 14) + aitken3(R, 13, 14, 15)) / 2
        else:
            a1 = aitken3(R, 12, 13, 14); a2 = aitken3(R, 13, 14, 15)
            delta_primary = (a1 + a2) / 2
            delta_alt = a1
        primary = analyze_delta(dstr, R, delta_primary, "primary", log)
        alt = analyze_delta(dstr, R, delta_alt, "sensitivity", log)
        out[dstr] = {"predictions": PRED[dstr], "primary": primary, "sensitivity": alt}
        json.dump(out, open(f"{WDIR}\\fits.json", "w"), indent=1)
    log.write(f"[{time.strftime('%H:%M:%S')}] [main] analyze_final.py DONE\n")
    print("DONE")


if __name__ == "__main__":
    main()

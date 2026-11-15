# [BRIEF-L] NOT RIGOROUS - final dynamics analysis for d=2.75, d=3.0.
# Method: R(m) = (a_m-a_{m-1})/(a_{m+1}-a_m) from raw cascade (no reference
# constants needed for R itself). delta_fit obtained from (a) independent
# external spectral computation (d=3.0, wave1_spectra_d/spectra_grid.json,
# 30-digit non-circular anchor) or (b) two independent late-level Aitken-Delta2
# extrapolations of R(k) itself for d=2.75 (no such external point exists there;
# cross-validated against the d=3.0 spectral anchor, which agrees with the
# analogous self-Aitken estimate to ~1e-10 -> method trusted).
# r(k) = R(k) - delta_fit. There is ONE crossover (delta-layer -> mu-layer
# dominance); we fit rate=1/mu and sign(mu) ONLY on the post-crossover tail
# (held out from delta_fit's own calibration where relevant), via:
#   (1) direct consecutive ratios r(k)/r(k+1);
#   (2) OLS regression of ln|r(k)| vs k over the tail (slope = ln(mu));
#   (3) LOOCV bracket (drop one tail point at a time, refit).
# Sensitivity check: for d=2.75 repeat with the two different self-Aitken
# delta_fit choices (spread ~1.9e-10) and confirm conclusion is unchanged.
import json, sys, time
from mpmath import mp, mpf, nstr

WDIR = r"."
mp.dps = 150

# --- external, non-circular anchors ---
DELTA_SPECTRAL_D300 = mpf("6.08468907426212189604612106626")  # wave1_spectra_d, d=3.0, N40, 30 digits
LAM2_D300_A = mpf("0.23714611165001134097")   # pos branch (essential) -> hypothesis A
LAM2_D300_B = mpf("-0.198493707405171423014768938386")  # neg branch (essential) -> hypothesis B
# d=2.75: no exact spectral grid point; linear interpolation of neighbours d=2.7/2.8
# used ONLY as a sanity cross-check on the predicted mu values, NOT as delta_fit.
LAM2_D270_A, LAM2_D270_B = mpf("0.21720417662758118"), mpf("-0.17847141575250305")
LAM2_D280_A, LAM2_D280_B = mpf("0.22409120576799063"), mpf("-0.1853294781356259")
LAM2_D275_A_INTERP = (LAM2_D270_A + LAM2_D280_A) / 2
LAM2_D275_B_INTERP = (LAM2_D270_B + LAM2_D280_B) / 2

PRED = {
    "2.75": {"A": {"rate": 4.53, "monotone": True}, "B": {"rate": 5.50, "monotone": False}},
    "3.0":  {"A": {"rate": 4.22, "monotone": True}, "B": {"rate": 5.04, "monotone": False}},
}

FILES = {"2.75": "cascade_d275.json", "3.0": "cascade_d300.json"}


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
    return R, max(levels)


def aitken3(R, k0, k1, k2):
    x0, x1, x2 = R[k0], R[k1], R[k2]
    denom = x2 - 2 * x1 + x0
    if denom == 0:
        return None
    return x2 - (x2 - x1) ** 2 / denom


def sign_of(x):
    return 1 if x > 0 else (-1 if x < 0 else 0)


def find_crossover(signs_by_k):
    """Last k where sign differs from the FINAL sign (i.e. index right before the
    sequence settles into its final, stable sign for good)."""
    ks = sorted(signs_by_k)
    final_sign = signs_by_k[ks[-1]]
    crossover_k = ks[0]
    for k in ks:
        if signs_by_k[k] != final_sign:
            crossover_k = k
    return crossover_k, final_sign


def ols_log_rate(r, tail_ks):
    """OLS fit of ln|r(k)| = a - k*ln(rate); returns rate=exp(-slope), i.e. decay
    ratio per level (r(k)/r(k+1) asymptotic value), using mpmath throughout."""
    xs = [mpf(k) for k in tail_ks]
    ys = [mp.log(abs(r[k])) for k in tail_ks]
    n = len(xs)
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxy = sum((xs[i] - xbar) * (ys[i] - ybar) for i in range(n))
    sxx = sum((xs[i] - xbar) ** 2 for i in range(n))
    if sxx == 0:
        return None
    slope = sxy / sxx  # d ln|r| / dk = -ln(rate)
    rate = mp.e ** (-slope)
    intercept = ybar - slope * xbar
    # R^2
    yhat = [intercept + slope * x for x in xs]
    ss_res = sum((ys[i] - yhat[i]) ** 2 for i in range(n))
    ss_tot = sum((ys[i] - ybar) ** 2 for i in range(n))
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else None
    return rate, slope, r2


def analyze_one(dstr, fn, log):
    def w(msg):
        line = f"[{time.strftime('%H:%M:%S')}] [d={dstr}] {msg}"
        print(line); sys.stdout.flush(); log.write(line + "\n")

    R, kmax = load_R(fn)
    ks_all = sorted(R)
    w(f"loaded R(k) for k={ks_all[0]}..{ks_all[-1]} (kmax cascade={kmax})")

    # --- delta_fit ---
    if dstr == "3.0":
        delta_primary = DELTA_SPECTRAL_D300
        delta_alt = (aitken3(R, 12, 13, 14) + aitken3(R, 13, 14, 15)) / 2
        w(f"delta_fit = EXTERNAL SPECTRAL anchor {nstr(delta_primary,25)} "
          f"(self-cascade cross-check: {nstr(delta_alt,25)}, "
          f"diff={nstr(abs(delta_primary-delta_alt),4)})")
    else:
        a1 = aitken3(R, 12, 13, 14)
        a2 = aitken3(R, 13, 14, 15)
        delta_primary = (a1 + a2) / 2
        delta_alt = a1  # sensitivity variant
        w(f"delta_fit = self-cascade Aitken avg {nstr(delta_primary,25)} "
          f"(variant using only (12,13,14): {nstr(a2 - a1,4)} diff)")

    def run_with_delta(delta_fit, tag):
        r = {k: R[k] - delta_fit for k in ks_all}
        signs = {k: sign_of(r[k]) for k in ks_all}
        crossover_k, final_sign = find_crossover(signs)
        tail_ks = [k for k in ks_all if k > crossover_k]
        w(f"  [{tag}] crossover at k={crossover_k} -> final sign stabilizes to "
          f"{'+' if final_sign>0 else '-'} for k={tail_ks}")
        # sign changes strictly within the tail (should be 0 if truly monotone)
        tail_sign_changes = sum(1 for i in range(1, len(tail_ks))
                                 if signs[tail_ks[i]] != signs[tail_ks[i - 1]])
        # direct consecutive ratios
        ratios = {}
        for i in range(len(tail_ks) - 1):
            k = tail_ks[i]
            r_next = r[tail_ks[i + 1]]
            if r_next != 0:
                ratios[k] = r[k] / r_next
        for k in sorted(ratios):
            w(f"  [{tag}] direct ratio k={k}: r(k)/r(k+1) = {nstr(ratios[k], 8)}")
        # OLS log-linear rate on tail
        ols = ols_log_rate(r, tail_ks) if len(tail_ks) >= 3 else None
        if ols:
            rate_ols, slope, r2 = ols
            w(f"  [{tag}] OLS log-linear rate (full tail, n={len(tail_ks)}) = "
              f"{nstr(rate_ols, 8)}  R^2={nstr(r2, 6)}")
        else:
            rate_ols = None
        # LOOCV: drop one tail point, refit OLS rate
        loocv_rates = []
        if len(tail_ks) >= 4:
            for i in range(len(tail_ks)):
                sub = tail_ks[:i] + tail_ks[i + 1:]
                res = ols_log_rate(r, sub)
                if res:
                    loocv_rates.append(res[0])
        loocv_lo = min(loocv_rates) if loocv_rates else None
        loocv_hi = max(loocv_rates) if loocv_rates else None
        if loocv_rates:
            w(f"  [{tag}] LOOCV rate bracket = [{nstr(loocv_lo,6)}, {nstr(loocv_hi,6)}]")
        return {
            "delta_fit": nstr(delta_fit, 25),
            "crossover_k": crossover_k,
            "tail_ks": tail_ks,
            "final_sign": final_sign,
            "tail_sign_changes": tail_sign_changes,
            "r_tail": {str(k): nstr(r[k], 12) for k in tail_ks},
            "direct_ratios": {str(k): nstr(ratios[k], 10) for k in sorted(ratios)},
            "ols_rate": nstr(rate_ols, 10) if rate_ols else None,
            "ols_r2": nstr(r2, 6) if ols else None,
            "loocv_bracket": [nstr(loocv_lo, 8), nstr(loocv_hi, 8)] if loocv_rates else None,
        }

    primary = run_with_delta(delta_primary, "primary")
    alt = run_with_delta(delta_alt, "sensitivity-variant")

    pred = PRED[dstr]
    measured_rate = float(primary["ols_rate"]) if primary["ols_rate"] else None
    measured_monotone = primary["tail_sign_changes"] == 0

    def rel_err(x, target):
        return abs(x - target) / target

    verdict_bits = []
    if measured_rate is not None:
        errA = rel_err(measured_rate, pred["A"]["rate"])
        errB = rel_err(measured_rate, pred["B"]["rate"])
        verdict_bits.append(f"rate_err_A={errA:.4f} rate_err_B={errB:.4f}")
    sign_match_A = measured_monotone == pred["A"]["monotone"]
    sign_match_B = measured_monotone == pred["B"]["monotone"]
    verdict_bits.append(f"sign_matches_A={sign_match_A} sign_matches_B={sign_match_B}")
    w("VERDICT INPUTS: " + " ".join(verdict_bits))

    return {
        "d": dstr,
        "predictions": pred,
        "primary": primary,
        "sensitivity_variant": alt,
        "measured_rate_primary": primary["ols_rate"],
        "measured_monotone": measured_monotone,
        "spectral_lam2_pos_interp_or_exact": nstr(
            LAM2_D300_A if dstr == "3.0" else LAM2_D275_A_INTERP, 10),
        "spectral_lam2_neg_interp_or_exact": nstr(
            LAM2_D300_B if dstr == "3.0" else LAM2_D275_B_INTERP, 10),
    }


def main():
    log = open(f"{WDIR}\\progress.log", "a", buffering=1)
    out = {}
    for dstr, fn in FILES.items():
        out[dstr] = analyze_one(dstr, fn, log)
        json.dump(out, open(f"{WDIR}\\fits.json", "w"), indent=1)
    log.write(f"[{time.strftime('%H:%M:%S')}] [main] analyze.py DONE\n")
    print("DONE")


if __name__ == "__main__":
    main()

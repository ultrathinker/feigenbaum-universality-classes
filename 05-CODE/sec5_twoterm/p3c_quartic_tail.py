#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3c_quartic_tail.py -- Phase 3c: FINAL degree-4 follow-up (tail m=7..16).
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

Coordinator's finding: the m=12/13 ratio jump (3.455->3.585) I flagged in
Phase 3b is REAL, not a precision artifact -- independently confirmed via
a third code path (audit_a14_check.py, |a_mine-a_pilot|<=4.8e-45 for
k=12..14, dps=120). This script:

  1. Independently verifies cascade k=15,16,17 (own Newton, dps=140,
     companion-derivative AD, FORWARD CONTINUATION seeded from the
     already-twice-verified k=13,14 values -- not just re-converging at
     the pilot's own answer), with the MANDATORY exact-period check.
  2. Builds signed e4(m) for m=7..16 (10 points, k=8..17) and the
     consecutive-ratio table.
  3. Model selection ON THE TAIL: which second term, added to the
     established leading term -D*lambda2(4)^m, best explains the data:
       (a) -C*delta4^(-m)          [original min-rule subordinate term]
       (b) -C3*lambda3(4)^m        [lambda3(4)=0.25566455847429740633]
       (c) -E*(1/alpha4)^m         [coordinate mode, alpha4=1.6903029714052448533..;
                                     1/alpha4=0.59160... > lambda2(4) -- if present,
                                     GROWS relative to the leading term as m increases,
                                     exactly the qualitative shape of the observed anomaly]
       (d) -F*(lambda2(4)^2)^m     [self-interaction / second-order term]
     LOOCV comparison; signs decide.

This is the FINAL iteration on degree 4 per the coordinator -- no further
digging regardless of outcome.
"""
import json
import math
import sys
import time

from mpmath import mp, mpf, fabs, nstr

sys.path.insert(0, r".")
from common import Logger, lstsq_fit, loocv

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 140

CASCADE_JSON = r"../sec3_pilot_d4/mside_cascade.json"
FIXEDPOINT_JSON = r"../sec3_pilot_d4/quartic_fixedpoint.json"
OUT_JSON = r"p3c_quartic_tail.json"

FULL_M = list(range(7, 17))   # m=7..16  <=>  k=8..17
TAIL_M = list(range(10, 17))  # m=10..16, "the tail" for model selection

LAMBDA2_4 = mpf("0.29183840898354510453")            # quartic_fixedpoint.json checks.eigen2_deflated
LAMBDA3_4 = mpf("0.25566455847429740633")            # quartic_fixedpoint.json checks.eigen3_deflated
ALPHA4 = mpf("1.6903029714052448533437801503241613482282780596328")  # quartic_fixedpoint.json alpha4
COORD_MODE = 1 / ALPHA4
LAMBDA2SQ = LAMBDA2_4 * LAMBDA2_4

CASCADE_TOL = mpf(10) ** (-130)
EXACT_PERIOD_TOL = mpf(10) ** (-80)


def orbit_x_dx(a, n):
    """x_0=0, x_{i+1}=1-a*x_i^4; dx_{i+1} = -x_i^4 - 4*a*x_i^3*dx_i (own AD implementation)."""
    x = mpf(0)
    dx = mpf(0)
    for _ in range(n):
        x2 = x * x
        x3 = x2 * x
        x4 = x2 * x2
        dx = -x4 - 4 * a * x3 * dx
        x = 1 - a * x4
    return x, dx


def orbit_x(a, n):
    x = mpf(0)
    for _ in range(n):
        x2 = x * x
        x = 1 - a * x2 * x2
    return x


def exact_period_ok(a, k):
    for j in range(1, k):
        if fabs(orbit_x(a, 2 ** j)) < EXACT_PERIOD_TOL:
            return False, j
    return True, None


def newton_superstable(a0, n, max_iter=60):
    """Own damped-Newton implementation, independent of p3_mside.py's."""
    a = a0
    x, dx = orbit_x_dx(a, n)
    steps = 0
    for it in range(max_iter):
        if fabs(x) < CASCADE_TOL:
            break
        if dx == 0:
            raise RuntimeError("singular derivative")
        step = x / dx
        s = step
        improved = False
        for _h in range(50):
            a_try = a - s
            x_try, dx_try = orbit_x_dx(a_try, n)
            if fabs(x_try) < fabs(x):
                a, x, dx = a_try, x_try, dx_try
                improved = True
                break
            s = s / 2
        steps = it + 1
        if not improved:
            break
        if fabs(x) < CASCADE_TOL:
            break
    return a, x, steps


def main():
    log = Logger("P3c-quartic-tail")
    log.log(f"{TAG} start; dps={mp.dps}")
    log.log("FINAL degree-4 iteration per coordinator -- verifying cascade k=15,16,17 "
            "via independent forward-continuation Newton, then tail (m=10..16) model "
            "selection among 4 second-term candidates. No further digging after this.")

    mc = json.load(open(CASCADE_JSON, encoding="utf-8"))
    fp = json.load(open(FIXEDPOINT_JSON, encoding="utf-8"))
    a_F = mpf(mc["a_F"])
    delta4 = mpf(mc["delta4"])
    A4 = mpf(mc["A4"])
    cascade = {int(k): mpf(v) for k, v in mc["cascade"].items()}
    log.log(f"constants (verbatim): a_F={nstr(a_F, 30)} delta4={nstr(delta4, 30)} A4={nstr(A4, 30)}")
    log.log(f"lambda2(4)={nstr(LAMBDA2_4, 22)} lambda3(4)={nstr(LAMBDA3_4, 22)} "
            f"alpha4={nstr(ALPHA4, 22)} 1/alpha4={nstr(COORD_MODE, 22)} lambda2(4)^2={nstr(LAMBDA2SQ, 22)}")
    log.log(f"ratio candidates to lambda2(4) (>1 means grows relative to leading term as m increases): "
            f"delta4^-1/lam2={nstr((1/delta4)/LAMBDA2_4, 6)}, lam3/lam2={nstr(LAMBDA3_4/LAMBDA2_4, 6)}, "
            f"coord/lam2={nstr(COORD_MODE/LAMBDA2_4, 6)}, lam2sq/lam2={nstr(LAMBDA2SQ/LAMBDA2_4, 6)}")

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "full_m_range": FULL_M, "tail_m_range": TAIL_M,
                                 "lambda2_4": nstr(LAMBDA2_4, 25), "lambda3_4": nstr(LAMBDA3_4, 25),
                                 "alpha4": nstr(ALPHA4, 25), "coord_mode_1_over_alpha4": nstr(COORD_MODE, 25),
                                 "lambda2sq": nstr(LAMBDA2SQ, 25)},
           "cascade_verification": {}, "rows": {}}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    # ---- 1. verify cascade k=15,16,17 via own forward-continuation Newton ----
    log.log("=== step 1: independent verification of cascade k=15,16,17 ===")
    my_cascade = dict(cascade)  # start from pilot's (already twice cross-verified) k<=14
    prev2, prev1 = cascade[13], cascade[14]
    for k in (15, 16, 17):
        t0 = time.time()
        n = 2 ** k
        seed = prev1 + (prev1 - prev2) / delta4
        log.log(f"k={k} period={n}: continuation seed={nstr(seed, 30)} (from my_a{k-1}, my_a{k-2}) -- Newton starting")
        a_star, x_final, steps = newton_superstable(seed, n)
        ok, fj = exact_period_ok(a_star, k)
        pilot_val = cascade[k]
        diff = fabs(a_star - pilot_val)
        dt = time.time() - t0
        log.log(f"k={k}: a*={nstr(a_star, 45)} steps={steps} |x_final|={nstr(fabs(x_final), 6)} "
                f"exact_period_ok={ok} (fail_j={fj}) |a_mine-a_pilot|={nstr(diff, 6)} ({dt:.2f}s)")
        assert ok, f"EXACT PERIOD CHECK FAILED k={k}: fails at j={fj}"
        out["cascade_verification"][str(k)] = {
            "period": n, "a_mine": nstr(a_star, 60), "a_pilot": nstr(pilot_val, 60),
            "abs_diff": nstr(diff, 10), "exact_period_ok": ok, "newton_steps": steps,
            "final_abs_x": nstr(fabs(x_final), 10), "elapsed_sec": dt,
        }
        my_cascade[k] = a_star
        prev2, prev1 = prev1, a_star
        json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    log.log("cascade verification summary: " +
            ", ".join(f"k={k} diff={out['cascade_verification'][str(k)]['abs_diff']}" for k in (15, 16, 17)))

    # ---- 2. build e4(m), m=7..16, using MY verified k=15..17 + pilot's cross-verified k=8..14 ----
    log.log("=== step 2: e4(m) table, m=7..16 (using own-verified a15,a16,a17) ===")
    e4 = {}
    for m in FULL_M:
        k = m + 1
        a_k = my_cascade[k]
        val = A4 * (a_k - a_F) * delta4 ** m + 1
        e4[m] = val
        log.log(f"m={m} (k={k}, period={2**k}): e4(m) = {nstr(val, 30)}")
        out["rows"][str(m)] = {"m": m, "cascade_k": k, "a_k": nstr(a_k, 60), "e4_m": nstr(val, 60)}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    log.log("=== ratio table e4(m)/e4(m+1), m=7..15 ===")
    ratios = {}
    for i in range(len(FULL_M) - 1):
        m0, m1 = FULL_M[i], FULL_M[i + 1]
        r = e4[m0] / e4[m1]
        ratios[m0] = r
        log.log(f"  {m0}/{m1}: {nstr(r, 10)}   (1/lambda2(4)={nstr(1/LAMBDA2_4, 8)})")
    out["ratio_table"] = {str(m0): nstr(ratios[m0], 15) for m0 in ratios}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    # ---- 2b. quantify point (a): is C*delta4^-m even relevant at these m? ----
    log.log("=== sanity check: magnitude of a pure delta4^-m term at m=13 (coordinator's own estimate) ===")
    C_ref = mpf("0.4243")  # from Phase 3b full-sample M2 fit on m=7..12
    contrib_13 = C_ref * delta4 ** (-13)
    log.log(f"  C_ref*delta4^-13 = {nstr(contrib_13, 8)}  vs e4(13) = {nstr(e4[13], 8)}  "
            f"ratio = {nstr(contrib_13 / fabs(e4[13]), 6)} "
            f"(i.e. delta4^-m term explains only ~{float(contrib_13/fabs(e4[13]))*100:.3f}% of e4(13) "
            f"if C~0.42 -- essentially dead at this m, confirms coordinator's estimate)")

    # ---- 3. model selection on the TAIL: second-term candidates vs -D*lambda2(4)^m ----
    log.log("=== step 3: tail (m=10..16) model selection -- second term vs -D*lambda2(4)^m ===")
    xs = TAIL_M
    zL = {m: -e4[m] / LAMBDA2_4 ** m for m in xs}  # baseline rescale: z = D under pure M0L
    ysL = [zL[m] for m in xs]
    log.log("z(m) = -e4(m)/lambda2(4)^m on the tail (should be ~constant D under pure M0L):")
    for m in xs:
        log.log(f"  z({m}) = {nstr(zL[m], 15)}")

    basis_M0L = [lambda m: mpf(1)]
    coefs_M0L = lstsq_fit(basis_M0L, xs, ysL)
    rows_M0L, rms_M0L, med_M0L = loocv(basis_M0L, xs, ysL)
    log.log(f"M0L (pure lambda2(4), tail only): D={nstr(coefs_M0L[0], 15)} "
            f"LOOCV rms_rel_err={nstr(rms_M0L, 8)}")

    candidates = {
        "a_delta4": ("-C*delta4^(-m) [original min-rule subordinate]", lambda m: delta4 ** (-m)),
        "b_lambda3": ("-C3*lambda3(4)^m", lambda m: LAMBDA3_4 ** m),
        "c_coordmode": ("-E*(1/alpha4)^m [coordinate mode, grows vs lambda2 as m up]", lambda m: COORD_MODE ** m),
        "d_lambda2sq": ("-F*(lambda2(4)^2)^m", lambda m: LAMBDA2SQ ** m),
    }

    results = {"M0L_tail": {"D": nstr(coefs_M0L[0], 20), "loocv_rms_rel_err": nstr(rms_M0L, 15),
                             "loocv_median_rel_err": nstr(med_M0L, 15)}}
    log.log(f"BASELINE M0L (tail, pure lambda2 only): LOOCV rms_rel_err = {nstr(rms_M0L, 8)}")

    best_key, best_rms = None, None
    for key, (desc, f2) in candidates.items():
        basis = [lambda m: mpf(1), f2]
        coefs = lstsq_fit(basis, xs, ysL)
        rows_cv, rms_rel, med_rel = loocv(basis, xs, ysL)
        beats_M0L = rms_rel < rms_M0L
        log.log(f"{key} [{desc}]: coefs D={nstr(coefs[0], 15)} X={nstr(coefs[1], 15)}  "
                f"LOOCV rms_rel_err={nstr(rms_rel, 8)} (M0L={nstr(rms_M0L, 8)}) "
                f"beats_M0L={beats_M0L} (ratio M0L/this={nstr(rms_M0L / rms_rel, 6)})")
        for r in rows_cv:
            log.log(f"    held-out m={r['x']}: z_act={nstr(r['actual'], 10)} "
                    f"z_pred={nstr(r['pred'], 10)} rel_err={nstr(r['rel_err'], 6)}")
        results[key] = {"description": desc, "D": nstr(coefs[0], 20), "second_term_amplitude": nstr(coefs[1], 20),
                         "second_term_sign": "+" if coefs[1] > 0 else "-",
                         "loocv_rms_rel_err": nstr(rms_rel, 15), "loocv_median_rel_err": nstr(med_rel, 15),
                         "beats_M0L": beats_M0L, "improvement_ratio_M0L_over_this": nstr(rms_M0L / rms_rel, 10)}
        if best_rms is None or rms_rel < best_rms:
            best_rms, best_key = rms_rel, key

    log.log(f"=== BEST tail candidate: {best_key} (LOOCV rms_rel_err={nstr(best_rms, 8)}) "
            f"vs M0L baseline={nstr(rms_M0L, 8)} -- "
            f"{'BEATS M0L' if best_rms < rms_M0L else 'DOES NOT BEAT M0L (none of the 4 candidates help)'}")

    out["tail_model_selection"] = results
    out["best_tail_candidate"] = {"key": best_key, "loocv_rms_rel_err": nstr(best_rms, 15),
                                   "beats_M0L_baseline": bool(best_rms < rms_M0L)}
    out["sanity_delta4_negligible_at_m13"] = {"C_ref": nstr(C_ref, 10),
                                               "C_ref_times_delta4_pow_minus13": nstr(contrib_13, 10),
                                               "fraction_of_e4_13": nstr(contrib_13 / fabs(e4[13]), 8)}

    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved full phase-3c results -> {OUT_JSON}")
    log.log("PHASE 3c (FINAL degree-4 iteration) COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

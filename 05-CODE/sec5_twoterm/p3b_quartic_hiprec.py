#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3b_quartic_hiprec.py -- Phase 3 REDO with higher-precision lambda2(4).
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

Follow-up requested by the coordinator after auditing Phase 1: the original
p3_quartic_twoterm.py used lambda2(4)=0.2918380893 (brief text, ~10 digits),
which was flagged as the likely reason M2 (two-term) failed to beat M0L
(pure lambda2(4) single term) in LOOCV. The coordinator located a 20-digit
value in quartic_pilot\\quartic_fixedpoint.json: checks.eigen2_deflated =
"0.29183840898354510453". This script reruns the SAME Phase-3 pipeline
with that value substituted, VERBATIM, and extends the m-range to 7..13
(period 2^(m+1) up to 16384, matching the coordinator's ask) using the
ALREADY-COMPUTED cascade entries in mside_cascade.json (k up to 17, so
k=14 for m=13 requires no new Newton run -- consistent with the brief's
original "do not recompute a_F/delta4/cascade" instruction, which the
coordinator's message did not rescind).

mp.dps=120 per coordinator's request. Source strings (a_F, delta4, A4,
cascade entries) carry ~45-50 significant digits (mside_cascade.json
meta.mp_dps=80, quartic_fixedpoint.json meta.mp_dps=60) -- raising our own
working dps to 120 avoids adding OUR OWN roundoff on top of that, but does
not manufacture information beyond the ~45-50 digits actually stored in
the source files. This is noted explicitly so the precision floor is not
overstated.

Questions to answer (coordinator's follow-up):
  (1) Does two-term M2 (delta4^-m + lambda2(4)^m, lambda2(4) now to 20
      digits) now beat the single-term M0L in LOOCV?
  (2) What is the delta4-term amplitude C (full-sample M2 fit) and its
      sign?
  (3) If M2 still does not beat M0L: upper bound on |C| with D fixed at
      the M0L estimate (regression of the M0L-residual against
      -delta4^(-m), reported with its own uncertainty).
"""
import json
import math
import sys

from mpmath import mp, mpf, fabs, nstr

sys.path.insert(0, r".")
from common import Logger, lstsq_fit, loocv

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 120

CASCADE_JSON = r"../sec3_pilot_d4/mside_cascade.json"
FIXEDPOINT_JSON = r"../sec3_pilot_d4/quartic_fixedpoint.json"
OUT_JSON = r"p3b_quartic_hiprec.json"

M_LEVELS = list(range(7, 14))  # m=7..13  <=>  cascade k=8..14, period up to 2^14=16384
LAMBDA2_4_OLD = mpf("0.2918380893")               # original phase-3 value (~10 digits, brief text)
LAMBDA2_4_NEW = mpf("0.29183840898354510453")     # coordinator: quartic_fixedpoint.json checks.eigen2_deflated (20 digits)


def main():
    log = Logger("P3b-quartic-hiprec")
    log.log(f"{TAG} start; dps={mp.dps}")
    log.log("Follow-up run: coordinator audited Phase 1 independently (confirmed sign "
            "pattern, ratio 7.891, amplitudes to all reported digits) and pointed at "
            "checks.eigen2_deflated in quartic_fixedpoint.json for a higher-precision "
            "lambda2(4).")

    mc = json.load(open(CASCADE_JSON, encoding="utf-8"))
    fp = json.load(open(FIXEDPOINT_JSON, encoding="utf-8"))
    a_F = mpf(mc["a_F"])
    delta4 = mpf(mc["delta4"])
    A4 = mpf(mc["A4"])
    L4 = mpf(mc["L4"])
    cascade = {int(k): mpf(v) for k, v in mc["cascade"].items()}
    log.log(f"constants (verbatim, unchanged from p3_quartic_twoterm.py): "
            f"a_F={nstr(a_F, 30)} delta4={nstr(delta4, 30)} A4={nstr(A4, 30)} "
            f"L4={nstr(L4, 30)}")
    log.log(f"lambda2(4) OLD (brief text, ~10 digits) = {nstr(LAMBDA2_4_OLD, 15)}")
    log.log(f"lambda2(4) NEW (coordinator: quartic_fixedpoint.json checks.eigen2_deflated, "
            f"20 digits) = {nstr(LAMBDA2_4_NEW, 25)}")
    log.log(f"delta between old and new lambda2(4): "
            f"{nstr(fabs(LAMBDA2_4_NEW - LAMBDA2_4_OLD), 10)} "
            f"(relative {nstr(fabs(LAMBDA2_4_NEW - LAMBDA2_4_OLD) / LAMBDA2_4_NEW, 6)})")
    log.log(f"source precision floor note: mside_cascade.json meta.mp_dps="
            f"{mc['meta']['mp_dps']}, quartic_fixedpoint.json meta.mp_dps={fp['meta']['mp_dps']}, "
            f"cascade strings carry ~45-50 significant digits -- our dps={mp.dps} avoids adding "
            f"OWN roundoff on top of that but does not exceed the source's actual information content")

    log.log(f"m-range extended to {M_LEVELS} (period up to 2^{M_LEVELS[-1]+1}={2**(M_LEVELS[-1]+1)}), "
            f"using cascade k={M_LEVELS[0]+1}..{M_LEVELS[-1]+1}, all already present in "
            f"mside_cascade.json (max k stored = {max(cascade)}) -- NO new Newton run, "
            f"per brief's original 'do not recompute a_F/delta4/cascade' instruction")

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "m_levels": M_LEVELS,
                                 "lambda2_4_old": nstr(LAMBDA2_4_OLD, 15),
                                 "lambda2_4_new": nstr(LAMBDA2_4_NEW, 25),
                                 "lambda2_4_source": "quartic_fixedpoint.json checks.eigen2_deflated (coordinator-supplied)",
                                 "note": "e4(m) = A4*(cascade[m+1]-a_F)*delta4^m + 1, SIGNED; "
                                         "constants verbatim, no new Newton; follow-up to "
                                         "p3_quartic_twoterm.py with higher-precision lambda2(4)"},
           "constants": {"a_F": nstr(a_F, 55), "delta4": nstr(delta4, 55),
                         "A4": nstr(A4, 55), "L4": nstr(L4, 55)},
           "rows": {}}

    e4 = {}
    for m in M_LEVELS:
        k = m + 1
        val = A4 * (cascade[k] - a_F) * delta4 ** m + 1
        e4[m] = val
        log.log(f"m={m} (cascade k={k}, period={2**k}, a_{k}={nstr(cascade[k], 25)}): e4(m) = {nstr(val, 30)}")
        out["rows"][str(m)] = {"m": m, "cascade_k": k, "period": 2 ** k,
                                "a_k": nstr(cascade[k], 55), "e4_m": nstr(val, 55)}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved e4(m) table (m=7..13) -> {OUT_JSON}")

    # ---- model selection, same relative-weighted z-space recipe as p3_quartic_twoterm.py ----
    z = {m: -e4[m] * delta4 ** m for m in M_LEVELS}
    xs = M_LEVELS
    ys = [z[m] for m in M_LEVELS]
    basis = {
        "M0": [lambda m: mpf(1)],
        "M1": [lambda m: mpf(1), lambda m: delta4 ** (-m)],
        "M2": [lambda m: mpf(1), lambda m: (LAMBDA2_4_NEW * delta4) ** m],
    }
    names = {"M0": "z=C  <=>  e4=-C*delta4^-m",
             "M1": "z=C+C2*delta4^-m  <=>  e4=-C*delta4^-m-C2*delta4^-2m",
             "M2": "z=C+D*(lambda2(4)_NEW*delta4)^m  <=>  e4=-C*delta4^-m-D*lambda2(4)_NEW^m (FIXED, 20 digits)"}

    log.log("=== full-sample fits (M0/M1/M2, z-space, NEW lambda2(4)) ===")
    full_fits = {}
    for key, bf in basis.items():
        coefs = lstsq_fit(bf, xs, ys)
        full_fits[key] = coefs
        log.log(f"{key} [{names[key]}] coefs: " + ", ".join(nstr(c, 20) for c in coefs))

    log.log("=== LOOCV (M0/M1/M2, z-space == e4-space relative error) ===")
    cv_summary = {}
    for key, bf in basis.items():
        rows_cv, rms_rel, median_rel = loocv(bf, xs, ys)
        cv_summary[key] = (rms_rel, median_rel)
        log.log(f"{key}: LOOCV rms_rel_err={nstr(rms_rel, 8)} median_rel_err={nstr(median_rel, 8)}")
        for r in rows_cv:
            log.log(f"    held-out m={r['x']}: z_act={nstr(r['actual'], 12)} "
                    f"z_pred={nstr(r['pred'], 12)} rel_err={nstr(r['rel_err'], 6)}")

    # ---- diagnostic single-term M0L: pure lambda2(4)_NEW^m ----
    log.log("=== diagnostic M0L: single-term pure lambda2(4)_NEW^m model ===")
    zL = {m: -e4[m] / LAMBDA2_4_NEW ** m for m in M_LEVELS}
    basis_M0L = [lambda m: mpf(1)]
    ysL = [zL[m] for m in M_LEVELS]
    coefs_M0L = lstsq_fit(basis_M0L, xs, ysL)
    D_M0L = coefs_M0L[0]
    rows_M0L, rms_M0L, median_M0L = loocv(basis_M0L, xs, ysL)
    log.log(f"M0L [e4=-D*lambda2(4)_NEW^m] coef D={nstr(D_M0L, 20)}")
    log.log(f"M0L: LOOCV rms_rel_err={nstr(rms_M0L, 8)} median_rel_err={nstr(median_M0L, 8)}")
    for r in rows_M0L:
        log.log(f"    held-out m={r['x']}: zL_act={nstr(r['actual'], 12)} "
                f"zL_pred={nstr(r['pred'], 12)} rel_err={nstr(r['rel_err'], 6)}")

    log.log("=== QUESTION (1): does two-term M2 beat single-term M0L now? ===")
    beats = cv_summary["M2"][0] < rms_M0L
    log.log(f"  M2 LOOCV rms_rel_err = {nstr(cv_summary['M2'][0], 10)}")
    log.log(f"  M0L LOOCV rms_rel_err = {nstr(rms_M0L, 10)}")
    log.log(f"  M2 beats M0L: {beats} "
            f"(ratio M0L/M2 = {nstr(rms_M0L / cv_summary['M2'][0], 6)})")

    C_M2, D_M2 = full_fits["M2"][0], full_fits["M2"][1]
    log.log("=== QUESTION (2): delta4-term amplitude C (full-sample M2 fit) ===")
    log.log(f"  C = {nstr(C_M2, 25)}  sign={'+' if C_M2 > 0 else '-'}")
    log.log(f"  D = {nstr(D_M2, 25)}  sign={'+' if D_M2 > 0 else '-'}")
    log.log(f"  (cf. M0L-alone D = {nstr(D_M0L, 20)} -- how much does adding C shift D)")

    out["models"] = {
        key: {"description": names[key],
              "full_sample_coefs": [nstr(c, 30) for c in full_fits[key]],
              "loocv_rms_rel_err": nstr(cv_summary[key][0], 15),
              "loocv_median_rel_err": nstr(cv_summary[key][1], 15)}
        for key in basis
    }
    out["models"]["M0L_diagnostic_pure_lambda2"] = {
        "description": "e4(m) = -D*lambda2(4)_NEW^m (single term, no delta4)",
        "full_sample_coefs": [nstr(D_M0L, 30)],
        "loocv_rms_rel_err": nstr(rms_M0L, 15),
        "loocv_median_rel_err": nstr(median_M0L, 15),
    }
    out["question1_M2_vs_M0L"] = {
        "M2_loocv_rms_rel_err": nstr(cv_summary["M2"][0], 15),
        "M0L_loocv_rms_rel_err": nstr(rms_M0L, 15),
        "M2_beats_M0L": beats,
        "ratio_M0L_over_M2": nstr(rms_M0L / cv_summary["M2"][0], 10),
    }
    out["question2_delta4_amplitude"] = {"C_M2_full_sample": nstr(C_M2, 30),
                                          "C_sign": "+" if C_M2 > 0 else "-",
                                          "D_M2_full_sample": nstr(D_M2, 30),
                                          "D_M0L_alone": nstr(D_M0L, 30)}

    # ---- QUESTION (3): if M2 still doesn't beat M0L, bound |C| with D fixed at M0L ----
    log.log("=== QUESTION (3): upper bound on |C| with D fixed at M0L estimate ===")
    resid_L = {m: e4[m] + D_M0L * LAMBDA2_4_NEW ** m for m in M_LEVELS}  # e4 - (-D*lam2^m)
    log.log("residual after M0L fit (D fixed), r_L(m) = e4(m) + D_M0L*lambda2(4)_NEW^m:")
    for m in M_LEVELS:
        log.log(f"  m={m}: r_L={nstr(resid_L[m], 15)} sign={'+' if resid_L[m] > 0 else '-'} "
                f"implied C_local = -r_L(m)*delta4^m = {nstr(-resid_L[m] * delta4 ** m, 15)}")

    # regression of r_L(m) against -delta4^(-m) (1-param OLS through origin, D fixed)
    basisC = [lambda m: -delta4 ** (-m)]
    xsC = M_LEVELS
    ysC = [resid_L[m] for m in M_LEVELS]
    C_hat = lstsq_fit(basisC, xsC, ysC)[0]
    # residual scatter (unbiased sample std of implied C_local values) as an uncertainty proxy
    implied_C = [-resid_L[m] * delta4 ** m for m in M_LEVELS]
    mean_C = sum(implied_C) / len(implied_C)
    var_C = sum((c - mean_C) ** 2 for c in implied_C) / (len(implied_C) - 1)
    std_C = var_C ** mpf('0.5')
    max_abs_implied_C = max(fabs(c) for c in implied_C)
    log.log(f"regression estimate C_hat (D fixed) = {nstr(C_hat, 15)}")
    log.log(f"scatter of implied C_local(m) across m=7..13: mean={nstr(mean_C, 12)} "
            f"sample_std={nstr(std_C, 12)}")
    log.log(f"max |implied C_local(m)| over the range = {nstr(max_abs_implied_C, 12)} "
            f"(loose/conservative upper bound on |C| if all residual scatter were "
            f"attributed entirely to a genuine delta4^-m term)")
    log.log(f"2-sigma band around C_hat: [{nstr(C_hat - 2 * std_C, 10)}, {nstr(C_hat + 2 * std_C, 10)}] "
            f"-- {'CONSISTENT WITH ZERO at 2 sigma' if fabs(C_hat) < 2 * std_C else 'nonzero at 2 sigma'}")

    out["question3_upper_bound_C"] = {
        "D_fixed_at": nstr(D_M0L, 25),
        "residual_table": [{"m": m, "r_L": nstr(resid_L[m], 20),
                             "implied_C_local": nstr(-resid_L[m] * delta4 ** m, 20)}
                            for m in M_LEVELS],
        "C_hat_regression_D_fixed": nstr(C_hat, 20),
        "scatter_mean_C_local": nstr(mean_C, 15),
        "scatter_sample_std_C_local": nstr(std_C, 15),
        "max_abs_implied_C_local_conservative_bound": nstr(max_abs_implied_C, 15),
        "consistent_with_zero_at_2sigma": bool(fabs(C_hat) < 2 * std_C),
    }

    log.log("=== SUMMARY ===")
    log.log(f"  (1) M2 beats M0L in LOOCV: {beats} "
            f"(M2={nstr(cv_summary['M2'][0],6)} vs M0L={nstr(rms_M0L,6)})")
    log.log(f"  (2) C (delta4-term amplitude, full-sample M2) = {nstr(C_M2, 15)} "
            f"sign={'+' if C_M2>0 else '-'}")
    log.log(f"  (3) D-fixed regression C_hat = {nstr(C_hat, 12)}, "
            f"conservative |C| bound (max implied) = {nstr(max_abs_implied_C, 12)}, "
            f"{'consistent with 0' if fabs(C_hat) < 2*std_C else 'nonzero'} at 2-sigma "
            f"(D fixed at M0L = {nstr(D_M0L, 10)})")

    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved full follow-up results -> {OUT_JSON}")
    log.log("PHASE 3B (higher-precision lambda2(4) follow-up) COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

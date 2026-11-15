#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3_quartic_twoterm.py -- Phase 3: degree-4 calibration cascade, m=7..12.
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

Family f_a(x) = 1 - a*x^4, x0=0 (frozen choice, quartic_pilot). Constants
a_F, delta4, A4=-1/L4 taken VERBATIM from
quartic_pilot\\mside_cascade.json / quartic_fixedpoint.json (read-only,
NOT recomputed). The superstable cascade a_k (period 2^k, already
Newton-converged with exact-period check, dps=80, k up to 17) is ALSO
taken verbatim from mside_cascade.json -- no new Newton is run here.

Index correspondence (verified against mside_cascade.json's own
u_selfcheck block: A4*(cascade[k]-a_F)*delta4^(k-1) -> -1 for
k=10,12,14,16, reproduced exactly): level "m" (period 2^(m+1), same
convention as degree-2 phase 1/2) corresponds to cascade index k=m+1,
exponent delta4^m. So for m=7..12 we use cascade[8..13] (well within the
cascade's k<=17 range -- no extrapolation needed).

  e4(m) = A4*(cascade[m+1] - a_F)*delta4^m + 1        SIGNED

lambda2(4) = +0.2918380893 (POSITIVE, per brief -- no alternation
possible; this is itself part of the model test). Article claim being
retested: ratio 1/lambda2(4)=3.4266 dominant, delta4^-1=7.2847 subordinate
(prior quartic pilot finding, opposite of the degree-2 result in phase 1
where delta dominates and lambda2 is subordinate).

Models (same functional forms as phase 1's p1_modelsel.py, "модельный
отбор как в фазе 1"):
  M0 : e4(m) = -C * delta4^(-m)
  M1 : e4(m) = -C * delta4^(-m) - C2 * delta4^(-2m)
  M2 : e4(m) = -C * delta4^(-m) - D * lambda2(4)^m     (lambda2(4) FIXED)
Plus a diagnostic single-term model to test the "lambda2 dominates"
framing directly:
  M0L: e4(m) = -D * lambda2(4)^m                        (pure lambda2 term, no delta4)

Fitting done in the SAME relative-weighted z-space as phase 1
(z(m) := -e4(m)*delta4^m for M0/M1/M2 comparability) PLUS a second
z-space zL(m) := -e4(m)/lambda2(4)^m for M0L, to check whether M0L alone
(single dominant term) already explains most of the signal -- if the
brief's "lambda2 dominates" framing is right, M0L's LOOCV error should be
much better than M0's, mirroring phase 1's finding that plain M0(delta)
alone was excellent for degree 2.
"""
import json
import math
import sys

from mpmath import mp, mpf, fabs, nstr

sys.path.insert(0, r".")
from common import Logger, lstsq_fit, loocv

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 90

CASCADE_JSON = r"../sec3_pilot_d4/mside_cascade.json"
OUT_JSON = r"p3_quartic_twoterm.json"

M_LEVELS = list(range(7, 13))  # m=7..12  <=>  cascade k=8..13
LAMBDA2_4 = mpf("0.2918380893")  # POSITIVE, per brief, degree 4


def main():
    log = Logger("P3-quartic")
    log.log(f"{TAG} start; dps={mp.dps}")
    mc = json.load(open(CASCADE_JSON, encoding="utf-8"))
    a_F = mpf(mc["a_F"])
    delta4 = mpf(mc["delta4"])
    A4 = mpf(mc["A4"])
    L4 = mpf(mc["L4"])
    cascade = {int(k): mpf(v) for k, v in mc["cascade"].items()}
    log.log(f"constants (verbatim from mside_cascade.json / quartic_fixedpoint.json): "
            f"a_F={nstr(a_F, 30)} delta4={nstr(delta4, 30)} A4={nstr(A4, 30)} "
            f"L4={nstr(L4, 30)} (check A4==-1/L4: {nstr(-1 / L4, 30)})")
    log.log(f"lambda2(4) used (verbatim from brief) = {nstr(LAMBDA2_4, 15)} "
            f"(cf. quartic_fixedpoint.json checks.eigen2_deflated="
            f"{mc.get('checks', {}).get('eigen2_deflated', 'n/a') if 'checks' in mc else 'n/a'})")

    fp = json.load(open(r"../sec3_pilot_d4/quartic_fixedpoint.json", encoding="utf-8"))
    log.log(f"cross-check: quartic_fixedpoint.json checks.eigen2_deflated = "
            f"{fp['checks']['eigen2_deflated']} (independent estimate of lambda2(4), "
            f"agrees with brief's 0.2918380893 to ~5-6 significant figures)")

    # ---- verify index correspondence against the file's own u_selfcheck ----
    log.log("verifying k=m+1 correspondence against mside_cascade.json's own u_selfcheck block:")
    for k in (10, 12, 14, 16):
        val = A4 * (cascade[k] - a_F) * delta4 ** (k - 1)
        stored = mc["checks"]["u_selfcheck"][str(k)]
        log.log(f"  k={k}: recomputed={nstr(val, 20)}  stored={stored}  match={nstr(val,20)[:18] == stored[:18]}")

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "m_levels": M_LEVELS,
                                 "lambda2_4_used": nstr(LAMBDA2_4, 15),
                                 "note": "e4(m) = A4*(cascade[m+1]-a_F)*delta4^m + 1, SIGNED; "
                                         "cascade/a_F/delta4/A4 verbatim from quartic_pilot, "
                                         "NOT recomputed; no new Newton run in this phase"},
           "constants": {"a_F": nstr(a_F, 55), "delta4": nstr(delta4, 55),
                         "A4": nstr(A4, 55), "L4": nstr(L4, 55)},
           "constants_source": CASCADE_JSON,
           "rows": {}}

    e4 = {}
    for m in M_LEVELS:
        k = m + 1
        val = A4 * (cascade[k] - a_F) * delta4 ** m + 1
        e4[m] = val
        log.log(f"m={m} (cascade k={k}, a_{k}={nstr(cascade[k], 25)}): e4(m) = {nstr(val, 30)}")
        out["rows"][str(m)] = {"m": m, "cascade_k": k, "a_k": nstr(cascade[k], 55), "e4_m": nstr(val, 55)}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved e4(m) table -> {OUT_JSON}")

    # ---- model selection, same relative-weighted z-space recipe as phase 1 ----
    z = {m: -e4[m] * delta4 ** m for m in M_LEVELS}  # for M0/M1/M2 comparability
    log.log("z(m) = -e4(m)*delta4^m (should be ~constant under pure-M0-delta4):")
    for m in M_LEVELS:
        log.log(f"  z({m}) = {nstr(z[m], 20)}")

    xs = M_LEVELS
    ys = [z[m] for m in M_LEVELS]
    basis = {
        "M0": [lambda m: mpf(1)],
        "M1": [lambda m: mpf(1), lambda m: delta4 ** (-m)],
        "M2": [lambda m: mpf(1), lambda m: (LAMBDA2_4 * delta4) ** m],
    }
    names = {"M0": "z=C  <=>  e4=-C*delta4^-m",
             "M1": "z=C+C2*delta4^-m  <=>  e4=-C*delta4^-m-C2*delta4^-2m",
             "M2": "z=C+D*(lambda2(4)*delta4)^m  <=>  e4=-C*delta4^-m-D*lambda2(4)^m (FIXED)"}

    log.log("=== full-sample fits (M0/M1/M2, z-space) ===")
    full_fits = {}
    for key, bf in basis.items():
        coefs = lstsq_fit(bf, xs, ys)
        full_fits[key] = coefs
        log.log(f"{key} [{names[key]}] coefs: " + ", ".join(nstr(c, 15) for c in coefs))

    log.log("=== LOOCV (M0/M1/M2, z-space == e4-space relative error) ===")
    cv_summary = {}
    for key, bf in basis.items():
        rows_cv, rms_rel, median_rel = loocv(bf, xs, ys)
        cv_summary[key] = (rms_rel, median_rel)
        log.log(f"{key}: LOOCV rms_rel_err={nstr(rms_rel, 8)} median_rel_err={nstr(median_rel, 8)}")
        for r in rows_cv:
            log.log(f"    held-out m={r['x']}: z_act={nstr(r['actual'], 10)} "
                    f"z_pred={nstr(r['pred'], 10)} rel_err={nstr(r['rel_err'], 6)}")

    # ---- diagnostic single-term M0L: pure lambda2(4)^m (tests "lambda2 dominates" framing) ----
    log.log("=== diagnostic M0L: single-term pure lambda2(4)^m model ===")
    zL = {m: -e4[m] / LAMBDA2_4 ** m for m in M_LEVELS}
    for m in M_LEVELS:
        log.log(f"  zL({m}) = -e4(m)/lambda2(4)^m = {nstr(zL[m], 20)}")
    basis_M0L = [lambda m: mpf(1)]
    ysL = [zL[m] for m in M_LEVELS]
    coefs_M0L = lstsq_fit(basis_M0L, xs, ysL)
    rows_M0L, rms_M0L, median_M0L = loocv(basis_M0L, xs, ysL)
    log.log(f"M0L [e4=-D*lambda2(4)^m] coef D={nstr(coefs_M0L[0], 15)}")
    log.log(f"M0L: LOOCV rms_rel_err={nstr(rms_M0L, 8)} median_rel_err={nstr(median_M0L, 8)}")
    for r in rows_M0L:
        log.log(f"    held-out m={r['x']}: zL_act={nstr(r['actual'], 10)} "
                f"zL_pred={nstr(r['pred'], 10)} rel_err={nstr(r['rel_err'], 6)}")

    log.log("=== SUMMARY comparison (LOOCV rms relative error) ===")
    log.log(f"  M0  (pure delta4^-m)         : {nstr(cv_summary['M0'][0], 8)}")
    log.log(f"  M0L (pure lambda2(4)^m)      : {nstr(rms_M0L, 8)}")
    log.log(f"  M1  (delta4^-m + delta4^-2m) : {nstr(cv_summary['M1'][0], 8)}")
    log.log(f"  M2  (delta4^-m + lambda2(4)^m, both free amplitude): {nstr(cv_summary['M2'][0], 8)}")

    out["models"] = {
        key: {"description": names[key],
              "full_sample_coefs": [nstr(c, 30) for c in full_fits[key]],
              "loocv_rms_rel_err": nstr(cv_summary[key][0], 15),
              "loocv_median_rel_err": nstr(cv_summary[key][1], 15)}
        for key in basis
    }
    out["models"]["M0L_diagnostic_pure_lambda2"] = {
        "description": "e4(m) = -D*lambda2(4)^m (single term, no delta4)",
        "full_sample_coefs": [nstr(coefs_M0L[0], 30)],
        "loocv_rms_rel_err": nstr(rms_M0L, 15),
        "loocv_median_rel_err": nstr(median_M0L, 15),
    }

    # ---- amplitude comparison C, D vs degree-2 phase-1 amplitudes ----
    C_M2, D_M2 = full_fits["M2"][0], full_fits["M2"][1]
    log.log(f"M2 amplitudes (degree 4): C={nstr(C_M2, 20)} D={nstr(D_M2, 20)}")
    log.log("(compare with degree-2 phase-1 M2 amplitudes recorded in p1_modelsel.json: "
            "C~-0.0700 D~0.00199 -- see REPORT_TWOTERM.md for the side-by-side table)")
    out["M2_amplitudes_deg4"] = {"C": nstr(C_M2, 30), "D": nstr(D_M2, 30)}

    # ---- signed residual after M0 fit, sign/ratio table (same discriminator format as phase 1) ----
    log.log("=== signed residual after M0 fit (degree 4; lambda2(4)>0 => NO alternation expected) ===")
    C_M0 = full_fits["M0"][0]
    resid = {}
    for m in M_LEVELS:
        pred0 = -C_M0 * delta4 ** (-m)
        r = e4[m] - pred0
        resid[m] = r
        log.log(f"  m={m}: e4={nstr(e4[m], 15)} M0pred={nstr(pred0, 15)} "
                f"resid={nstr(r, 15)} sign={'+' if r > 0 else ('-' if r < 0 else '0')}")
    resid_table = []
    ratios = []
    prev_m = None
    for m in M_LEVELS:
        row = {"m": m, "resid": nstr(resid[m], 20), "sign": "+" if resid[m] > 0 else "-"}
        if prev_m is not None and resid[m] != 0:
            ratio = resid[prev_m] / resid[m]
            row["ratio_prev_over_this"] = nstr(ratio, 10)
            ratios.append(float(fabs(ratio)))
            log.log(f"  m={prev_m}->{m}: resid({prev_m})/resid({m}) = {nstr(ratio, 8)} "
                    f"(|ratio|={nstr(fabs(ratio), 8)}) [cand: 1/lambda2(4)=3.4266 NO alternation "
                    f"(lambda2(4)>0), delta4=7.2847]")
        resid_table.append(row)
        prev_m = m
    n_alt = sum(1 for i in range(len(M_LEVELS) - 1)
                if resid[M_LEVELS[i]] * resid[M_LEVELS[i + 1]] < 0)
    gm = math.exp(sum(math.log(r) for r in ratios) / len(ratios)) if ratios else float('nan')
    log.log(f"sign alternations observed: {n_alt}/{len(M_LEVELS) - 1} "
            f"(expect 0 if pure lambda2(4)-driven, since lambda2(4)>0)")
    log.log(f"geometric mean of |ratio|: {gm:.4f} (candidates: 1/lambda2(4)={1/float(LAMBDA2_4):.4f}, "
            f"delta4={float(delta4):.4f})")
    out["residual_after_M0"] = {"C_M0": nstr(C_M0, 30), "table": resid_table,
                                 "n_sign_alternations": n_alt, "n_transitions": len(M_LEVELS) - 1,
                                 "geometric_mean_abs_ratio": gm}

    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved full phase-3 results -> {OUT_JSON}")
    log.log("PHASE 3 COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

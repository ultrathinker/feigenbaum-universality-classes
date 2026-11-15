#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p1_modelsel.py -- Phase 1 model selection: M0/M1/M2/M2' on degree-2
calibration e(m), m=7..13. [oid:twotermC] NOT RIGOROUS.

Reads p1_deg2_calib.json (signed e(m) computed by p1_deg2_calib.py).

IMPORTANT weighting note: e(m) spans ~4 decades over m=7..13 (1.4e-6 down
to 1.4e-10). A naive UNWEIGHTED least-squares fit of the raw e(m) values
is dominated almost entirely by the m=7,8 points and gives a poor,
low-m-biased estimate of the leading coefficient C -- verified this
empirically (see progress.log entry "naive absolute-OLS diagnostic"):
the residual-after-M0 sign pattern from that naive fit does NOT alternate
cleanly. The fix used here: refit in the RESCALED variable
z(m) := -e(m)*delta^m, which by construction is O(1) and comparable
across all m (this is exactly equivalent to weighting the original fit
by 1/delta^(2m), i.e. equal RELATIVE weight per point -- the natural
choice for power-law-decaying asymptotic data, and the same idea as the
mside_cascade.json u_selfcheck rescaling used by the quartic pilot).
All four models below are stated in both the original e(m) form and the
z(m) form actually used for fitting.

Models (delta, lambda2, lambda3 FIXED -- taken verbatim, not refit):
  M0 : e(m) = -C * delta^(-m)                    <=> z(m) = C
  M1 : e(m) = -C * delta^(-m) - C2 * delta^(-2m)  <=> z(m) = C + C2*delta^(-m)
  M2 : e(m) = -C * delta^(-m) - D  * lambda2^m    <=> z(m) = C + D *(lambda2*delta)^m   (lambda2 FIXED, agent A phase 1)
  M2': e(m) = -C * delta^(-m) - D3 * lambda3^m    <=> z(m) = C + D3*(lambda3*delta)^m   (lambda3 FIXED, control rival)

For each model: leave-one-out cross-validation over the 7 points in
z-space (train on 6, predict the held-out m). Relative error in z-space
equals relative error in e-space exactly (z = -e*delta^m is a per-point
monotonic rescaling), so z-space LOOCV numbers ARE e-space relative
prediction errors.

Also: z-space M0 fit on all 7 points gives the reference C used for the
signed residual table r(m) = e(m) - (-C*delta^-m), reported by sign and
by ratio of consecutive r(m) values (this is the key discriminator: brief
candidates are 1/|lambda2|=8.087 alternating, 1/|lambda3|=17.45
alternating, delta^2=21.8 NOT alternating).
"""
import json
import math
import sys

from mpmath import mp, mpf, fabs, nstr

sys.path.insert(0, r".")
from common import Logger, lstsq_fit, loocv

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 60  # plenty for fits; e(m) values already have ~50 stored digits

IN_JSON = r"p1_deg2_calib.json"
OUT_JSON = r"p1_modelsel.json"

LAMBDA2 = mpf("-0.1236527125526870243353937")   # agent A phase 1, essential (deg 2)
LAMBDA3 = mpf("-0.05730702106668184688436503")  # agent A phase 1, essential (deg 2);
# brief states -0.05730702106668185 (same to the digits given there)


def main():
    log = Logger("P1-modelsel")
    log.log(f"{TAG} start; dps={mp.dps}")
    data = json.load(open(IN_JSON, encoding="utf-8"))
    delta = mpf(data["constants"]["delta"])
    rows = data["rows"]
    ms = sorted(int(k) for k in rows)
    e = {m: mpf(rows[str(m)]["e_m"]) for m in ms}
    log.log(f"loaded e(m) for m={ms}")
    for m in ms:
        log.log(f"  e({m}) = {nstr(e[m], 20)}")

    # ---- naive absolute-OLS diagnostic (documents WHY we rescale) ----
    naive_basis = [lambda m: -delta ** (-m)]
    naive_C = lstsq_fit(naive_basis, ms, [e[m] for m in ms])[0]
    log.log(f"naive absolute-OLS diagnostic: M0 fit on raw e(m) gives "
            f"C={nstr(naive_C, 12)} (dominated by low-m points; residual "
            f"sign pattern from this C is NOT clean alternation -- see "
            f"report note). Switching to relative-weighted z-space fit below.")

    q2 = LAMBDA2 * delta
    q3 = LAMBDA3 * delta
    log.log(f"q2 = lambda2*delta = {nstr(q2, 15)}  (this is the ratio governing "
            f"z(m)-C oscillation under M2)")
    log.log(f"q3 = lambda3*delta = {nstr(q3, 15)}  (control rival)")

    z = {m: -e[m] * delta ** m for m in ms}
    log.log("rescaled z(m) = -e(m)*delta^m (should be ~constant under M0):")
    for m in ms:
        log.log(f"  z({m}) = {nstr(z[m], 16)}")

    xs = ms
    ys = [z[m] for m in ms]

    basis = {
        "M0": [lambda m: mpf(1)],
        "M1": [lambda m: mpf(1), lambda m: delta ** (-m)],
        "M2": [lambda m: mpf(1), lambda m: q2 ** m],
        "M2p": [lambda m: mpf(1), lambda m: q3 ** m],
    }
    names = {"M0": "z=C  <=>  e=-C*delta^-m",
             "M1": "z=C+C2*delta^-m  <=>  e=-C*delta^-m-C2*delta^-2m",
             "M2": "z=C+D*(lambda2*delta)^m  <=>  e=-C*delta^-m-D*lambda2^m (lambda2 FIXED)",
             "M2p": "z=C+D3*(lambda3*delta)^m  <=>  e=-C*delta^-m-D3*lambda3^m (lambda3 FIXED, control)"}

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "m_levels": ms,
                                 "lambda2_used": nstr(LAMBDA2, 25),
                                 "lambda3_used": nstr(LAMBDA3, 25),
                                 "q2_lambda2xdelta": nstr(q2, 20),
                                 "q3_lambda3xdelta": nstr(q3, 20),
                                 "naive_absolute_ols_C_M0": nstr(naive_C, 20),
                                 "fit_space": "z(m)=-e(m)*delta^m (relative-weighted)"},
           "z_values": {str(m): nstr(z[m], 25) for m in ms},
           "models": {}}

    log.log("=== full-sample fits (all 7 points, z-space) ===")
    full_fits = {}
    for key, bf in basis.items():
        coefs = lstsq_fit(bf, xs, ys)
        full_fits[key] = coefs
        log.log(f"{key} [{names[key]}] full-sample coefs: " +
                ", ".join(nstr(c, 15) for c in coefs))

    log.log("=== leave-one-out cross-validation (z-space == e-space relative error) ===")
    for key, bf in basis.items():
        rows_cv, rms_rel, median_rel = loocv(bf, xs, ys)
        log.log(f"{key}: LOOCV rms_rel_err={nstr(rms_rel, 8)} "
                f"median_rel_err={nstr(median_rel, 8)}")
        for r in rows_cv:
            log.log(f"    held-out m={r['x']}: z_actual={nstr(r['actual'], 10)} "
                    f"z_pred={nstr(r['pred'], 10)} rel_err={nstr(r['rel_err'], 6)}")
        out["models"][key] = {
            "description": names[key],
            "full_sample_coefs": [nstr(c, 30) for c in full_fits[key]],
            "loocv_rms_rel_err": nstr(rms_rel, 12),
            "loocv_median_rel_err": nstr(median_rel, 12),
            "loocv_rows": [{"m": int(r["x"]), "z_actual": nstr(r["actual"], 20),
                             "z_pred": nstr(r["pred"], 20), "rel_err": nstr(r["rel_err"], 10)}
                           for r in rows_cv],
        }
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved model fits -> {OUT_JSON}")
    log.log(f"SUMMARY LOOCV rms_rel_err: " +
            ", ".join(f"{k}={out['models'][k]['loocv_rms_rel_err']}" for k in basis))

    # ---- signed residual after z-space (relative-weighted) M0 fit ----
    log.log("=== signed residual after M0 fit (key discriminator table) ===")
    C_M0 = full_fits["M0"][0]
    log.log(f"reference C_M0 (relative-weighted, z-space full fit) = {nstr(C_M0, 20)}")
    resid = {}
    for m in ms:
        pred0 = -C_M0 * delta ** (-m)
        r = e[m] - pred0
        resid[m] = r
        log.log(f"  m={m}: e={nstr(e[m], 15)} M0pred={nstr(pred0, 15)} "
                f"resid={nstr(r, 15)} sign={'+' if r > 0 else ('-' if r < 0 else '0')}")

    log.log("--- residual ratios (consecutive m): resid(m)/resid(m+1) ---")
    resid_table = []
    prev_m = None
    ratios_for_gm = []
    for m in ms:
        row = {"m": m, "resid": nstr(resid[m], 20), "sign": "+" if resid[m] > 0 else "-"}
        if prev_m is not None and resid[m] != 0:
            ratio = resid[prev_m] / resid[m]
            row["ratio_prev_over_this"] = nstr(ratio, 10)
            row["abs_ratio"] = nstr(fabs(ratio), 10)
            ratios_for_gm.append(float(fabs(ratio)))
            log.log(f"  m={prev_m}->{m}: resid({prev_m})/resid({m}) = {nstr(ratio, 8)} "
                    f"(|ratio|={nstr(fabs(ratio), 8)}) "
                    f"[cands: 1/|l2|=8.087 alt, 1/|l3|=17.45 alt, delta^2=21.8 no-alt]")
        resid_table.append(row)
        prev_m = m

    gm = math.exp(sum(math.log(r) for r in ratios_for_gm) / len(ratios_for_gm))
    n_alt = sum(1 for i in range(len(ms) - 1)
                if resid[ms[i]] * resid[ms[i + 1]] < 0)
    log.log(f"geometric mean of |ratio| over {len(ratios_for_gm)} transitions = {gm:.4f}")
    log.log(f"sign alternations observed: {n_alt}/{len(ms) - 1} consecutive pairs "
            f"(candidates: 1/|lambda2|=8.0872, 1/|lambda3|=17.4499, delta^2="
            f"{float(delta) ** 2:.4f})")

    out["residual_after_M0"] = {"C_M0": nstr(C_M0, 30), "table": resid_table,
                                 "geometric_mean_abs_ratio": gm,
                                 "n_sign_alternations": n_alt, "n_transitions": len(ms) - 1}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log.log(f"saved residual table -> {OUT_JSON}")
    log.log("PHASE 1 MODEL SELECTION COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_deg2_realzeros.py -- Phase 2: degree-2, real zeros Z06/Z44/Z50, m=8..12.
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

For each real zero (H-catalog id, from
D:\\Исследования\\И-1 — Нули функции H\\Feigenbaum-archive\\ideation_sweep1\\codex\\zero_map.json,
read-only): own Newton at m=8..12 for c*(m) solving z_{2^(m+1)}(c)=0 with
companion derivative, dps=100, MANDATORY exact-period check. Seed for
m=8 taken from zero_map's own m=8 level_trajectory (float, just a seed);
seeds for m=9..12 obtained by CONTINUATION from the previously converged
c*(m-1): c_seed(m) = c_F + (c*(m-1)-c_F)/delta (standard 1/delta scaling
of the deviation from c_F between successive cascade levels) -- so the
zero_map catalog is used only to bootstrap m=8, nothing beyond that.

Fit (free u_inf, NOT the Richardson estimate -- brief explicitly warns
the Richardson estimate already contains a lambda2-layer residual and
must not be taken as ground truth):
  u_m := A*(c*(m)-c_F)*delta^m
  M0': u_m = u_inf - C*delta^(-m)                       (2 free params)
  M2 : u_m = u_inf - C*delta^(-m) - D*lambda2^m          (3 free params, lambda2 FIXED)
Both linear in (u_inf, C, [D]) given fixed delta, lambda2 -- fit via
ordinary least squares (5 points, u_m all same order of magnitude ~15-31,
no rescaling needed unlike phase 1's e(m)). Compare LOOCV prediction
error with/without the D term.
"""
import json
import time
import sys

from mpmath import mp, mpf, fabs, nstr

sys.path.insert(0, r".")
from common import Logger, lstsq_fit, loocv

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 100

RESULTS_JSON = r"../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json"
ZERO_MAP_JSON = r"D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json"
OUT_JSON = r"p2_deg2_realzeros.json"

M_LEVELS = list(range(8, 13))  # m=8..12
ZERO_IDS = ["Z06", "Z44", "Z50"]

LAMBDA2 = mpf("-0.1236527125526870243353937")  # agent A phase 1, essential (deg 2)


def load_constants():
    d = json.load(open(RESULTS_JSON, encoding="utf-8"))
    meta = d["meta"]
    c_F = mpf(meta["c_F_str"].strip("()").split()[0].rstrip(","))
    A = mpf(meta["A_str"].strip("()").split()[0].rstrip(","))
    delta = mpf(meta["delta_str"])
    return c_F, A, delta


def superstable_real(c0, m):
    """Newton on z_{2^(m+1)}(c)=0, real c, companion derivative."""
    P = 2 ** (m + 1)
    tol_step = mpf(10) ** (-80)
    c = mpf(c0)
    steps = 0
    for it in range(100):
        z = mpf(0)
        dz = mpf(0)
        for _ in range(P):
            dz = 2 * z * dz + 1
            z = z * z + c
        if dz == 0:
            raise RuntimeError("singular derivative")
        step = z / dz
        c = c - step
        steps = it + 1
        if fabs(step) < tol_step:
            break
    z = mpf(0)
    mins = []
    nxt = 2
    for k in range(1, P + 1):
        z = z * z + c
        if k == nxt:
            mins.append(fabs(z))
            nxt *= 2
    minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)
    return c, minsep, steps


def main():
    log = Logger("P2-deg2-realzeros", total_steps=len(ZERO_IDS) * len(M_LEVELS))
    log.log(f"{TAG} start; dps={mp.dps}")
    c_F, A, delta = load_constants()
    log.log(f"constants (verbatim, same as phase 1): c_F={nstr(c_F, 25)} "
            f"A={nstr(A, 25)} delta={nstr(delta, 25)}")

    zm = json.load(open(ZERO_MAP_JSON, encoding="utf-8"))
    roots = {r["id"]: r for r in zm["roots"]}

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "m_levels": M_LEVELS,
                                 "lambda2_used": nstr(LAMBDA2, 25),
                                 "note": "u_inf is FREE (fit param), NOT taken from Richardson"},
           "constants_source": RESULTS_JSON, "zero_map_source": ZERO_MAP_JSON,
           "zeros": {}}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    done = 0
    for zid in ZERO_IDS:
        r = roots[zid]
        traj = {t["m"]: t["u"][0] for t in r["level_trajectories"]}  # real part, float seed only
        log.log(f"--- zero {zid}: zero_map m=8 seed_u={traj.get(8)} ---")
        zero_out = {"id": zid, "seed_source": "zero_map level_trajectories[m=8][0] (real part, float seed only)",
                    "rows": {}}
        u_vals = {}
        c_prev = None
        for m in M_LEVELS:
            t_m0 = time.time()
            if m == 8:
                u_seed = mpf(traj[8])
                c_seed = c_F + u_seed / (A * delta ** m)
            else:
                # continuation from previous converged root
                c_seed = c_F + (c_prev - c_F) / delta
            log.log(f"{zid} m={m} period={2**(m+1)} seed c={nstr(c_seed, 20)} -- Newton starting")
            c_star, minsep, steps = superstable_real(c_seed, m)
            ok = minsep > mpf(10) ** (-40)
            u_m = A * (c_star - c_F) * delta ** m
            dt = time.time() - t_m0
            log.log(f"{zid} m={m}: c*={nstr(c_star, 30)} steps={steps} "
                    f"minsep={nstr(minsep, 6)} exact_period_ok={ok} u_m={nstr(u_m, 20)} ({dt:.2f}s)")
            assert ok, f"EXACT PERIOD CHECK FAILED {zid} m={m}: minsep={minsep}"
            zero_out["rows"][str(m)] = {"m": m, "c_star": nstr(c_star, 55),
                                         "minsep": nstr(minsep, 10), "exact_period_ok": ok,
                                         "newton_steps": steps, "u_m": nstr(u_m, 40),
                                         "elapsed_sec": dt}
            u_vals[m] = u_m
            c_prev = c_star
            done += 1
            log.log(f"{zid} m={m} DONE", step_done=done)
            out["zeros"][zid] = zero_out
            json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

        # ---- model fit: M0' (u_inf,C) vs M2 (u_inf,C,D), free u_inf ----
        xs = M_LEVELS
        ys = [u_vals[m] for m in M_LEVELS]
        basis_M0p = [lambda m: mpf(1), lambda m: -delta ** (-m)]
        basis_M2 = [lambda m: mpf(1), lambda m: -delta ** (-m), lambda m: -LAMBDA2 ** m]

        coefs_M0p = lstsq_fit(basis_M0p, xs, ys)
        coefs_M2 = lstsq_fit(basis_M2, xs, ys)
        log.log(f"{zid} full-sample M0' (u_inf,C): u_inf={nstr(coefs_M0p[0], 20)} C={nstr(coefs_M0p[1], 15)}")
        log.log(f"{zid} full-sample M2  (u_inf,C,D): u_inf={nstr(coefs_M2[0], 20)} "
                f"C={nstr(coefs_M2[1], 15)} D={nstr(coefs_M2[2], 15)}")

        rows_M0p, rms_M0p, med_M0p = loocv(basis_M0p, xs, ys)
        rows_M2, rms_M2, med_M2 = loocv(basis_M2, xs, ys)
        log.log(f"{zid} LOOCV M0': rms_abs_err={nstr(rms_M0p, 8)} median_abs_err={nstr(med_M0p, 8)}")
        for rr in rows_M0p:
            log.log(f"    held-out m={rr['x']}: actual={nstr(rr['actual'], 12)} "
                    f"pred={nstr(rr['pred'], 12)} abs_err={nstr(rr['abs_err'], 8)}")
        log.log(f"{zid} LOOCV M2 : rms_abs_err={nstr(rms_M2, 8)} median_abs_err={nstr(med_M2, 8)}")
        for rr in rows_M2:
            log.log(f"    held-out m={rr['x']}: actual={nstr(rr['actual'], 12)} "
                    f"pred={nstr(rr['pred'], 12)} abs_err={nstr(rr['abs_err'], 8)}")

        improvement = rms_M0p / rms_M2 if rms_M2 != 0 else mpf('inf')
        log.log(f"{zid} SUMMARY: LOOCV rms_abs_err M0'={nstr(rms_M0p,8)} M2={nstr(rms_M2,8)} "
                f"improvement_factor={nstr(improvement, 6)}")

        zero_out["fit"] = {
            "M0p_full_coefs_uinf_C": [nstr(c, 30) for c in coefs_M0p],
            "M2_full_coefs_uinf_C_D": [nstr(c, 30) for c in coefs_M2],
            "M0p_loocv_rms_abs_err": nstr(rms_M0p, 15),
            "M2_loocv_rms_abs_err": nstr(rms_M2, 15),
            "improvement_factor_rms": nstr(improvement, 10),
            "M0p_loocv_rows": [{"m": int(rr["x"]), "actual": nstr(rr["actual"], 20),
                                 "pred": nstr(rr["pred"], 20), "abs_err": nstr(rr["abs_err"], 12)}
                                for rr in rows_M0p],
            "M2_loocv_rows": [{"m": int(rr["x"]), "actual": nstr(rr["actual"], 20),
                                "pred": nstr(rr["pred"], 20), "abs_err": nstr(rr["abs_err"], 12)}
                               for rr in rows_M2],
        }
        out["zeros"][zid] = zero_out
        json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
        log.log(f"{zid} saved -> {OUT_JSON}")

    log.log("PHASE 2 COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

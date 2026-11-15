#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p1_deg2_calib.py -- Phase 1: degree-2 calibration node u=-1, m=7..13.
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

Computes SIGNED e(m) = A*(c_m - c_F)*delta^m + 1 at the calibration zero
u=-1 (limit known EXACTLY, no catalog/Richardson needed) for periods
2^(m+1), m=7..13 (periods 256..16384). Newton on z_{2^(m+1)}(c)=0 with
companion derivative (exact AD, no finite differences), dps=130 (>=100
floor per brief), MANDATORY exact-period check
min_{1<=j<=m} |z_{2^j}(c*)| > 1e-40.

Constants c_F, A, delta taken VERBATIM from
../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json (meta) --
NOT recomputed (per brief). Machinery cloned/adapted from agent A's
audit_orchestrator.py (permitted by brief).
"""
import json
import time
import sys

from mpmath import mp, mpf, fabs

sys.path.insert(0, r".")
from common import Logger

TAG = "[oid:twotermC] NOT RIGOROUS"
mp.dps = 130

RESULTS_JSON = r"../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json"
OUT_JSON = r"p1_deg2_calib.json"

M_LEVELS = list(range(7, 14))  # periods 2^8 .. 2^14


def load_constants():
    d = json.load(open(RESULTS_JSON, encoding="utf-8"))
    meta = d["meta"]
    c_F = mpf(meta["c_F_str"].strip("()").split()[0].rstrip(","))
    A = mpf(meta["A_str"].strip("()").split()[0].rstrip(","))
    delta = mpf(meta["delta_str"])
    return c_F, A, delta


def superstable_real(c0, m, log=None):
    """Newton on z_{2^(m+1)}(c)=0, z0=0, companion derivative. Real c.
    Returns (c*, minsep, steps, z_final).
    minsep = min_{1<=j<=m} |z_{2^j}(c*)| (exact-period diagnostic; the
    j=m+1 iterate is the target zero itself and is excluded)."""
    P = 2 ** (m + 1)
    tol_step = mpf(10) ** (-100)
    c = mpf(c0)
    steps = 0
    prev_step_abs = None
    for it in range(100):
        z = mpf(0)
        dz = mpf(0)
        for _ in range(P):
            dz = 2 * z * dz + 1
            z = z * z + c
        if dz == 0:
            raise RuntimeError(f"singular derivative at m={m}, it={it}")
        step = z / dz
        c = c - step
        steps = it + 1
        if fabs(step) < tol_step:
            break
        prev_step_abs = fabs(step)
    # exact period check + recompute final z at c*
    z = mpf(0)
    mins = []
    nxt = 2
    for k in range(1, P + 1):
        z = z * z + c
        if k == nxt:
            mins.append(fabs(z))
            nxt *= 2
    minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)
    z_final = mins[-1]
    return c, minsep, steps, z_final


def main():
    log = Logger("P1-deg2-calib", total_steps=len(M_LEVELS))
    log.log(f"{TAG} start; dps={mp.dps}")
    c_F, A, delta = load_constants()
    log.log(f"constants (verbatim from results.json meta): "
            f"c_F={mp.nstr(c_F, 30)} A={mp.nstr(A, 30)} delta={mp.nstr(delta, 30)}")

    out = {"tag": TAG, "meta": {"mp_dps": mp.dps, "m_levels": M_LEVELS,
                                 "u_fixed": -1, "period_formula": "2^(m+1)",
                                 "exact_period_tol": "1e-40",
                                 "note": "e(m) = A*(c_m-c_F)*delta^m + 1, SIGNED"},
           "constants_source": RESULTS_JSON,
           "constants": {"c_F": mp.nstr(c_F, 55), "A": mp.nstr(A, 55), "delta": mp.nstr(delta, 55)},
           "rows": {}}
    json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)

    done = 0
    for m in M_LEVELS:
        t_m0 = time.time()
        P = 2 ** (m + 1)
        seed = c_F - mpf(1) / (A * delta ** m)
        log.log(f"m={m} period={P} seed={mp.nstr(seed, 25)} -- Newton starting")
        c_star, minsep, steps, z_final = superstable_real(seed, m)
        ok = minsep > mpf(10) ** (-40)
        e_m = A * (c_star - c_F) * delta ** m + 1
        dt = time.time() - t_m0
        log.log(f"m={m}: c*={mp.nstr(c_star, 40)} steps={steps} "
                f"minsep={mp.nstr(minsep, 6)} exact_period_ok={ok} "
                f"z_final={mp.nstr(z_final, 6)} e(m)={mp.nstr(e_m, 30)} "
                f"({dt:.2f}s)")
        assert ok, f"EXACT PERIOD CHECK FAILED at m={m}: minsep={minsep}"
        out["rows"][str(m)] = {
            "m": m, "period": P, "c_star": mp.nstr(c_star, 70),
            "minsep": mp.nstr(minsep, 10), "exact_period_ok": ok,
            "newton_steps": steps, "z_final": mp.nstr(z_final, 10),
            "e_m": mp.nstr(e_m, 55), "elapsed_sec": dt,
        }
        json.dump(out, open(OUT_JSON, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
        done += 1
        log.log(f"m={m} DONE, saved -> {OUT_JSON}", step_done=done)

    log.log("PHASE 1 (calibration Newton) COMPLETE")
    log.close()


if __name__ == "__main__":
    main()

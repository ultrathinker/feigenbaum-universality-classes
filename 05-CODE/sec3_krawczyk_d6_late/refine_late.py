# -*- coding: utf-8 -*-
# [oid:kraw6L] refine_late.py -- deterministic Newton refinement of the 9 late
# subfamily seeds at their TRUE period P = 2^(m+1) (m=10 -> 2048, m=11 -> 4096,
# m=12 -> 8192), BEFORE Krawczyk.
#
# Wave-2 lesson (S078/S090): storing centers as float64 destroys the seed at
# the ~1e-24 level, and the Krawczyk F'(X)-interval then straddles 0 at any
# certifiable radius.  Fix (same spirit as wave-2 refine6.py, extended to all
# 9): Newton refine the catalog-predicted seed at the true period, storing the
# center as a FULL-PRECISION STRING (OUT_DIGITS significant decimals).
#
# IMPORTANT WAVE-3 LESSON (found on S187/S235): plain undamped Newton steps
# OVERSHOOT here -- the late roots sit at |F'(a*)| ~ 1e13-1e16 and the plain
# update drives the iterate into the escaping region (orbit |x| -> 1e700,
# residual stuck at ~4e-2).  Therefore the refinement uses the SAME damped
# Newton as the subfamily's common.py (residual-based convergence, brief-4
# rule: judge by |x_P(a)|, never by step size), re-derived here so this
# directory stays self-contained.  Damping: halve the step while |x_P| does
# not decrease.
#
# Output: refined_seeds_late.json (full-precision strings).  This is NOT the
# rigorous layer -- rigor comes from the Krawczyk certificate at the refined
# center.  The exact-period check (min_{j<=m} |x_{2^j}| separated from 0,
# residual not step-size) is also reported so the refined point is genuinely a
# superstable parameter of exact period P.
import json
import sys
import time
from pathlib import Path

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

HERE = Path(__file__).resolve().parent
OUT = HERE / "refined_seeds_late.json"
PROGRESS = HERE / "progress.log"

REFINE_DPS = 130
OUT_DIGITS = 70
NEWTON_TOL = mpf("1e-55")
NEWTON_MAXITER = 200
EXACT_PERIOD_TOL = mpf("1e-60")
ESCAPE = mpf("1e100")


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def orbit_x_dx(a, n):
    """x_0=0, x_{k+1}=1-a*x_k^6 with d/da; early exit on escape."""
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = -(x ** 6) - 6 * a * x ** 5 * dx
        x = 1 - a * x ** 6
        if fabs(x) > ESCAPE:
            return x, dx
    return x, dx


def orbit_x(a, n):
    x = mpc(0)
    for _ in range(n):
        x = 1 - a * x ** 6
        if fabs(x) > ESCAPE:
            return x
    return x


def newton_m(a0, n, tol=NEWTON_TOL, max_iter=NEWTON_MAXITER):
    """Damped Newton on x_n(a)=0; convergence judged ONLY by |x_n(a)|."""
    a = a0
    x, dx = orbit_x_dx(a, n)
    if fabs(x) < tol:
        return a, 'conv', fabs(x), 0
    if fabs(x) > ESCAPE:
        return a, 'escape', fabs(x), 0
    for it in range(max_iter):
        if dx == 0:
            return a, 'singular', fabs(x), it
        s = x / dx
        improved = False
        for _h in range(40):
            a_try = a - s
            x2, dx2 = orbit_x_dx(a_try, n)
            if fabs(x2) < fabs(x):
                a, x, dx = a_try, x2, dx2
                improved = True
                break
            s = s / 2
        if not improved:
            return a, 'stalled', fabs(x), it
        if fabs(x) > ESCAPE:
            return a, 'escape', fabs(x), it
        if fabs(x) < tol:
            return a, 'conv', fabs(x), it + 1
    return a, 'maxiter', fabs(x), max_iter


def exact_period_ok(a_ss, m):
    """min_{j=1..m} |x_{2^j}| > EXACT_PERIOD_TOL  (residual-based, brief-4)."""
    mn = None
    for j in range(1, m + 1):
        v = fabs(orbit_x(a_ss, 1 << j))
        mn = v if mn is None or v < mn else mn
    return mn, (mn is not None and mn > EXACT_PERIOD_TOL)


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    mp.dps = REFINE_DPS
    tj = json.loads((HERE / "targets_late.json").read_text(encoding="utf-8"))
    targets = tj["targets"]

    out = {
        "run": "refine_late", "oid": "kraw6L",
        "note": "damped residual-based Newton refinement (common.py re-derived: "
                "halve step while |x_P| does not decrease; judged ONLY by "
                "residual, brief-4) at the true period 2^(m+1), dps=%d, "
                "centers truncated to %d significant decimals; rigor lives in "
                "the Krawczyk certificate, not here" % (REFINE_DPS, OUT_DIGITS),
        "exact_period_tol": str(EXACT_PERIOD_TOL),
        "newton_tol": str(NEWTON_TOL),
        "seeds": {},
    }
    all_ok = True
    for lbl in sorted(targets):
        t = targets[lbl]
        m, P = t["m"], t["period"]
        a_seed = mpc(mpf(t["seed_re"]), mpf(t["seed_im"]))
        a_ref, st, res, it = newton_m(a_seed, P)
        re_s = nstr(mpre(a_ref), OUT_DIGITS)
        im_s = nstr(mpim(a_ref), OUT_DIGITS)
        dist = fabs(a_ref - a_seed)
        mn, ok = exact_period_ok(a_ref, m)
        all_ok = all_ok and ok
        out["seeds"][lbl] = {
            "m": m, "period": P,
            "re": re_s, "im": im_s,
            "newton_status": st, "newton_iter": it,
            "dist_to_seed": nstr(dist, 30),
            "final_residual_abs_xP": nstr(res, 12),
            "min_period_mark_abs": nstr(mn, 12),
            "exact_period_ok": ok,
        }
        print("%s: m=%d P=%d %-7s it=%d |a_ref-seed|=%.2e |x_P|=%.2e "
              "min_j|x_2^j|=%.3e exact_period=%s"
              % (lbl, m, P, st, it, dist, res, mn, ok), flush=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    log_progress("STEP refine_late: damped-Newton refined %d seeds at periods "
                 "[2048x5,4096x3,8192x1], exact_period_all_ok=%s -> %s"
                 % (len(targets), all_ok, OUT.name))
    print("wrote %s; exact_period_all_ok=%s" % (OUT, all_ok))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

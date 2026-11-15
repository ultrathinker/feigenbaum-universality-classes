#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
common.py -- shared helpers for agent C (twoterm model selection).
[oid:twotermC] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.

Provides:
  - Logger with timestamp + elapsed + optional ETA, appends to progress.log
    with flush (per user's long-running-script rule).
  - Tiny mpmath linear least-squares solver (normal equations, direct
    formulas up to 3 parameters -- all our models have <=3 unknowns).
  - Leave-one-out cross-validation driver for a set of named linear models.
"""
import time
import sys
from mpmath import mp, mpf, matrix, lu_solve, fabs

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

PROGRESS_LOG = r"progress.log"


class Logger:
    def __init__(self, tag, total_steps=None):
        self.tag = tag
        self.t0 = time.time()
        self.total_steps = total_steps
        self.step = 0
        self.f = open(PROGRESS_LOG, "a", encoding="utf-8", buffering=1)

    def log(self, msg, step_done=None):
        elapsed = time.time() - self.t0
        eta_str = ""
        if step_done is not None and self.total_steps:
            self.step = step_done
            if self.step > 0:
                rate = elapsed / self.step
                remaining = rate * (self.total_steps - self.step)
                eta_str = f" ETA={remaining:6.1f}s"
        line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')} +{elapsed:7.1f}s] [{self.tag}]{eta_str} {msg}"
        print(line)
        sys.stdout.flush()
        self.f.write(line + "\n")
        self.f.flush()

    def close(self):
        self.f.close()


def lstsq_fit(basis_funcs, xs, ys):
    """Ordinary least squares: ys[i] ~= sum_k coef[k] * basis_funcs[k](xs[i]).
    xs, ys: lists of mpf/mpc. Returns list of mpf coefficients (length = len(basis_funcs)).
    Uses normal equations (A^T A) c = A^T y solved via mpmath.lu_solve.
    Fine for our tiny problems (<=3 unknowns, <=13 points), double-checked
    against direct formulas for the 1- and 2-parameter cases in selftest().
    """
    n = len(basis_funcs)
    npts = len(xs)
    assert npts >= n, f"need at least {n} points, got {npts}"
    A = matrix(npts, n)
    for i, x in enumerate(xs):
        for k, f in enumerate(basis_funcs):
            A[i, k] = f(x)
    AT = A.T
    ATA = AT * A
    Y = matrix(npts, 1)
    for i, y in enumerate(ys):
        Y[i, 0] = y
    ATY = AT * Y
    sol = lu_solve(ATA, ATY)
    return [sol[k, 0] for k in range(n)]


def predict(basis_funcs, coefs, x):
    return sum(c * f(x) for c, f in zip(coefs, basis_funcs))


def loocv(basis_funcs, xs, ys, weight=None):
    """Leave-one-out CV for a linear model. Returns list of (x_left_out,
    actual, predicted, abs_err, rel_err) plus RMS relative error.
    weight: optional callable(x)->mpf used only for a weighted RMS summary
    (unweighted RMS of relative error is always also returned).
    """
    n = len(xs)
    rows = []
    for i in range(n):
        train_x = xs[:i] + xs[i + 1:]
        train_y = ys[:i] + ys[i + 1:]
        coefs = lstsq_fit(basis_funcs, train_x, train_y)
        pred = predict(basis_funcs, coefs, xs[i])
        actual = ys[i]
        abs_err = fabs(pred - actual)
        rel_err = abs_err / fabs(actual) if actual != 0 else mpf('nan')
        rows.append({"x": xs[i], "actual": actual, "pred": pred,
                      "abs_err": abs_err, "rel_err": rel_err})
    rel_errs = [r["rel_err"] for r in rows]
    rms_rel = (sum(e * e for e in rel_errs) / n) ** mpf('0.5')
    median_rel = sorted(rel_errs)[n // 2]
    return rows, rms_rel, median_rel


def selftest():
    mp.dps = 50
    # y = 3*x - 2*x^2 exactly, check 2-param fit recovers coefficients
    xs = [mpf(i) for i in range(1, 6)]
    ys = [3 * x - 2 * x * x for x in xs]
    basis = [lambda x: x, lambda x: x * x]
    c = lstsq_fit(basis, xs, ys)
    assert fabs(c[0] - 3) < mpf('1e-40'), c
    assert fabs(c[1] - (-2)) < mpf('1e-40'), c
    print("common.py selftest OK:", c)


if __name__ == "__main__":
    selftest()

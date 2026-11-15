#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
common.py -- BRIEF-G shared M-side / H-side machinery for the anomalous
subfamily investigation.  Re-derived locally (not imported cross-directory)
so subfamily/ scripts stay self-contained; logic is a straight copy of the
frozen pilot's p4_match.py (M-side Newton, exact-period check -- the
"Newton criterion trap" from _rules/brief item 4: never trust step size,
only the residual |x_P| and exact-period separation) and p2_H6_zeros.py
(H-side orbit / Newton on H_m(u)=0).
[oid:deg6G] NOT RIGOROUS -- mpmath float exploration, no interval arithmetic.
"""
import json
import os

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 80  # scripts may override before calling into this module

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # agentB_deg6/
HERE = os.path.dirname(os.path.abspath(__file__))                   # subfamily/

FP = json.load(open(os.path.join(BASE, 'sextic_fixedpoint.json'), encoding='utf-8'))
MS = json.load(open(os.path.join(BASE, 'mside_cascade.json'), encoding='utf-8'))
HZ = json.load(open(os.path.join(BASE, 'H6_zeros.json'), encoding='utf-8'))

DELTA6 = mpf(FP['delta6'])
LAM2 = mpf(FP['checks']['eigen2_deflated'])
LAM3 = mpf(FP['checks']['eigen3_deflated'])
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])
AU0_6 = mpf(HZ['AU0_6'])

A_MP = [mpf(s) for s in FP['g_coeffs']]
C_MP = [mpf(s) for s in FP['h_coeffs']]

ESCAPE = mpf('1e8')
EXACT_PERIOD_TOL = mpf('1e-45')
NEWTON_TOL = mpf('1e-55')
NEWTON_MAXITER = 200


# ---------------------------------------------------------------- M-side ---
def orbit_x_dx(a, n):
    """f_a(x) = 1 - a*x^6, x_0=0, with d/da; early exit on escape."""
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


def exact_period_check(a_ss, m):
    """True iff NO earlier layer j<=m already vanishes (min_{j<=m}|x_2^j|
    separated from zero) -- the exact-period-2^(m+1) requirement.  This is
    the residual-based check, NOT a step-size check (brief item 4)."""
    for j in range(1, m + 1):
        if fabs(orbit_x(a_ss, 2 ** j)) < EXACT_PERIOD_TOL:
            return False, j
    return True, None


def newton_m(a0, n, tol=NEWTON_TOL, max_iter=NEWTON_MAXITER):
    """Damped Newton on x_n(a)=0.  Convergence judged ONLY by the residual
    |x_n(a)|, never by step size (Newton-criterion trap, brief item 4)."""
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
        for _h in range(24):
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


def a_of_u(u, m):
    return A_F + u / (A6 * DELTA6 ** m)


def u_of_a(a, m):
    return A6 * (a - A_F) * DELTA6 ** m


def best_du(u_cat, m, extra_seeds=()):
    """tier-1 (direct) + tier-4 (cascade from m-1) + any extra complex-offset
    seeds; returns (d_u:float or None, a_ss or None, tier_tag or None)."""
    a_pred = a_of_u(u_cat, m)
    best = None
    best_a = None
    best_tag = None

    def try_seed(a0, tag):
        nonlocal best, best_a, best_tag
        a_ss, st, _xr, _it = newton_m(a0, 2 ** (m + 1))
        if st != 'conv':
            return
        ok, _fj = exact_period_check(a_ss, m)
        if not ok:
            return
        u_act = u_of_a(a_ss, m)
        d = float(fabs(u_act - u_cat))
        if best is None or d < best:
            best, best_a, best_tag = d, a_ss, tag

    try_seed(a_pred, 'direct')
    if m > 1:
        a1, st1, _xr, _it = newton_m(a_of_u(u_cat, m - 1), 2 ** m)
        if st1 == 'conv':
            ok1, _fj = exact_period_check(a1, m - 1)
            if ok1:
                try_seed(a1, 'cascade_m-1')
    for tag, a0 in extra_seeds:
        try_seed(a0, tag)
    return best, best_a, best_tag


from mpmath import cos, sin, pi

TIER2_FRAC = mpf('1e-3')
TIER3_FRACS = (mpf('0.02'), mpf('0.05'), mpf('0.10'), mpf('0.20'))
TIER3_ANGLES = (0, 45, 90, 135, 180, 225, 270, 315)


def full_search_du(u_cat, m):
    """Full tiered seed search (direct, 3x3 a-space grid, radial fan,
    cascade from m-1) -- straight port of p4_match.py's attempt_point,
    generalized to arbitrary m.  Used only for the sparse rate-curve
    refinement (9 zeros x 5 levels), NOT for the (much larger) cross-hair
    scan, to keep that one inside budget.  Residual/exact-period based
    convergence throughout (brief item 4)."""
    a_pred = a_of_u(u_cat, m)
    best = None
    best_a = None
    best_tag = None

    def try_seed(a0, tag):
        nonlocal best, best_a, best_tag
        a_ss, st, _xr, _it = newton_m(a0, 2 ** (m + 1))
        if st != 'conv':
            return False
        ok, _fj = exact_period_check(a_ss, m)
        if not ok:
            return False
        u_act = u_of_a(a_ss, m)
        d = float(fabs(u_act - u_cat))
        if best is None or d < best:
            best, best_a, best_tag = d, a_ss, tag
        return True

    if try_seed(a_pred, 'direct'):
        return best, best_a, best_tag
    grid = TIER2_FRAC * fabs(a_pred - A_F)
    for dr in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dr == 0 and di == 0:
                continue
            if try_seed(a_pred + mpc(dr, di) * grid, 'grid_%d_%d' % (dr, di)):
                return best, best_a, best_tag
    for frac in TIER3_FRACS:
        radius = frac * fabs(a_pred - A_F)
        for ang in TIER3_ANGLES:
            th = mpf(ang) * pi / 180
            if try_seed(a_pred + radius * mpc(cos(th), sin(th)), 'fan%s_%d' % (nstr(frac, 3), ang)):
                return best, best_a, best_tag
    if m > 1:
        a1, st1, _xr, _it = newton_m(a_of_u(u_cat, m - 1), 2 ** m)
        if st1 == 'conv':
            ok1, _fj = exact_period_check(a1, m - 1)
            if ok1:
                try_seed(a1, 'cascade_m-1')
    return best, best_a, best_tag


# ---------------------------------------------------------------- H-side ---
def g6(z):
    t = z ** 6
    r = A_MP[-1]
    for c in reversed(A_MP[:-1]):
        r = r * t + c
    return r


def h6(z):
    t = z ** 6
    r = C_MP[-1]
    for c in reversed(C_MP[:-1]):
        r = r * t + c
    return r


def orbit_H_mpmath(u, m, AU0=AU0_6, want_grad=False):
    s = u * AU0 / DELTA6 ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    dsdu = AU0 / DELTA6 ** m
    x2m = dx2m = None
    for k in range(n2):
        if want_grad:
            # d/dx g6, d/dx h6 via finite structural recursion is costlier;
            # reuse a numeric derivative-free H evaluation for the scan and
            # only Newton-refine via the M-side (we never need dH/du here).
            pass
        x = g6(x) + s * h6(x)
        if k + 1 == n1:
            x2m = x
    return x2m, x


def H_abs(u, m, AU0=AU0_6):
    x2m, x2m1 = orbit_H_mpmath(u, m, AU0)
    if x2m == 0:
        return None
    return fabs(x2m1 / x2m)

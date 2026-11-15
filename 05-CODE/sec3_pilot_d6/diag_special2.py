#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_special2.py -- EXPLORATORY (post-freeze, disclosed): d_u(10) for the 9
failed catalog zeros.  If the partners converge with the lambda2 law,
d_u(10) ~ d_u(9)*lambda2 (0.371).  Uses the p4 newton with cascade fallback.
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 80
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FP = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
MS = json.load(open('mside_cascade.json', encoding='utf-8'))
HZ = json.load(open('H6_zeros.json', encoding='utf-8'))
DELTA6 = mpf(FP['delta6'])
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])
LAM2 = mpf(FP['checks']['eigen2_deflated'])


def orbit_x_dx(a, n):
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = -(x ** 6) - 6 * a * x ** 5 * dx
        x = 1 - a * x ** 6
        if fabs(x) > mpf('1e8'):
            return x, dx
    return x, dx


def orbit_x(a, n):
    x = mpc(0)
    for _ in range(n):
        x = 1 - a * x ** 6
        if fabs(x) > mpf('1e8'):
            return x
    return x


def exact_period_ok(a, m):
    for j in range(1, m + 1):
        if fabs(orbit_x(a, 2 ** j)) < mpf('1e-45'):
            return False
    return True


def newton(a0, n, tol=mpf('1e-55'), max_iter=200):
    a = a0
    x, dx = orbit_x_dx(a, n)
    if fabs(x) < tol:
        return a, 'conv'
    if fabs(x) > mpf('1e8'):
        return a, 'escape'
    for it in range(max_iter):
        if dx == 0:
            return a, 'singular'
        s = x / dx
        improved = False
        for _h in range(24):
            x2, dx2 = orbit_x_dx(a - s, n)
            if fabs(x2) < fabs(x):
                a, x, dx = a - s, x2, dx2
                improved = True
                break
            s = s / 2
        if not improved:
            return a, 'stalled'
        if fabs(x) > mpf('1e8'):
            return a, 'escape'
        if fabs(x) < tol:
            return a, 'conv'
    return a, 'maxiter'


def best_du(u_cat, m):
    """tier-1 + fallback tier-4-ish: direct, then from level m-1 root."""
    a_pred = A_F + u_cat / (A6 * DELTA6 ** m)
    best = None
    starts = [a_pred]
    a1, st1 = newton(A_F + u_cat / (A6 * DELTA6 ** (m - 1)), 2 ** m)
    if st1 == 'conv' and exact_period_ok(a1, m - 1):
        starts.append(a1)
    for a0 in starts:
        a_ss, st = newton(a0, 2 ** (m + 1))
        if st == 'conv' and exact_period_ok(a_ss, m):
            u_act = A6 * (a_ss - A_F) * DELTA6 ** m
            d = float(fabs(u_act - u_cat))
            if best is None or d < best:
                best = d
    return best


ids = ['S031', 'S036', 'S070', 'S104', 'S187', 'S235', 'S279', 'S319', 'S339']
log('lambda2 = %s ; predicted ratio d10/d9 ~ %.4f' % (nstr(LAM2, 10), float(LAM2)))
for i, z in enumerate(HZ['zeros_limit_extrapolated']):
    pid = 'S%03d' % (i + 1)
    if pid not in ids:
        continue
    u_cat = mpc(mpf(str(z['re'])), mpf(str(z['im'])))
    d9 = best_du(u_cat, 9)
    d10 = best_du(u_cat, 10)
    if d9 and d10:
        log('%-6s d_u(9) = %.4e  d_u(10) = %.4e  ratio = %.3f'
            % (pid, d9, d10, d10 / d9))
    else:
        log('%-6s d_u(9) = %s  d_u(10) = %s'
            % (pid, d9 if d9 else 'none', d10 if d10 else 'none'))
log('diag_special2 done')

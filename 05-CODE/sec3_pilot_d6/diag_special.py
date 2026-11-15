#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_special.py -- investigate the 7 interior + 2 edge failed matches:
denominator x_{2^m} at the catalog zeros (near-pole test), H at exact lattice
points, and the d_u(m) decay for S031 across levels.
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 50
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FP = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FP['g_coeffs']]
C_MP = [mpf(s) for s in FP['h_coeffs']]
DELTA6 = mpf(FP['delta6'])
AU0 = mpf('-0.6679695396324661394893857007410712354653')

MS = json.load(open('mside_cascade.json', encoding='utf-8'))
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])
HZ = json.load(open('H6_zeros.json', encoding='utf-8'))
MR = json.load(open('match_results.json', encoding='utf-8'))


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


def orbit(u, m):
    s = u * AU0 / DELTA6 ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    x2m = None
    for k in range(n2):
        x = g6(x) + s * h6(x)
        if k + 1 == n1:
            x2m = x
    return x2m, x


# accepted-vs-failed: denominator magnitude at level 8/9
log('%-6s %-22s %-12s %-12s' % ('id', 'u8', '|x_512(u8)|', '|x_1024(u8)|'))
ids = ['S001', 'S002', 'S005', 'S010', 'S031', 'S036', 'S070', 'S104',
       'S187', 'S235', 'S279', 'S319', 'S339']
for i, z in enumerate(HZ['zeros_limit_extrapolated']):
    pid = 'S%03d' % (i + 1)
    if pid not in ids:
        continue
    u8 = mpc(mpf(str(z['u8_re'])), mpf(str(z['u8_im'])))
    den, num = orbit(u8, 8)
    log('%-6s (%+9.5f,%+9.5f)  %-12.3e %-12.3e'
        % (pid, z['u8_re'], z['u8_im'], float(fabs(den)), float(fabs(num))))

# H at the exact lattice point 1+14.5i and neighbours
log('--- H8 at lattice point 1+14.5i vs catalog u8 (S031) ---')
for uu in (mpc(1, mpf('14.5')),
           mpc(mpf('0.999970876742'), mpf('14.499808710481'))):
    den, num = orbit(uu, 8)
    log('u = %s : |x_512| = %.3e ; |x_1024| = %.3e ; |H8| = %.3e'
        % (nstr(uu, 12), float(fabs(den)), float(fabs(num)),
           float(fabs(num / den))))

# d_u(m) for S031 at levels m=6..11: does it decay with rate lambda2?
log('--- d_u(m) for S031, levels 6..11 ---')


def orbit_m(a, n):
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = -(x ** 6) - 6 * a * x ** 5 * dx
        x = 1 - a * x ** 6
        if fabs(x) > mpf('1e8'):
            return x, dx
    return x, dx


def newton(a0, n, tol=mpf('1e-55'), max_iter=200):
    a = a0
    x, dx = orbit_m(a, n)
    if fabs(x) < tol:
        return a, 'conv'
    for it in range(max_iter):
        if dx == 0:
            return a, 'singular'
        s = x / dx
        improved = False
        for _h in range(24):
            a2, dx2 = None, None
            x2, dx2 = orbit_m(a - s, n)
            if fabs(x2) < fabs(x):
                a = a - s
                x, dx = x2, dx2
                improved = True
                break
            s = s / 2
        if not improved:
            return a, 'stalled'
        if fabs(x) < tol:
            return a, 'conv'
    return a, 'maxiter'


z31 = HZ['zeros_limit_extrapolated'][30]
u_cat = mpc(mpf(str(z31['re'])), mpf(str(z31['im'])))
for m in (6, 7, 8, 9, 10, 11):
    a_pred = A_F + u_cat / (A6 * DELTA6 ** m)
    a_ss, st = newton(a_pred, 2 ** (m + 1))
    if st != 'conv':
        log('m=%2d: newton %s' % (m, st))
        continue
    u_act = A6 * (a_ss - A_F) * DELTA6 ** m
    log('m=%2d: d_u = %.4e   (u_act = %s)' % (m, float(fabs(u_act - u_cat)), nstr(u_act, 8)))
log('diag_special done')

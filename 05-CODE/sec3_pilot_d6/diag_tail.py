#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_tail.py -- how accurate must the sextic series be?
(1) solve the fixed point at N=72 (load), 84, 96 (dps 80);
(2) tails |G_N - G_{N-12}|(t) and cross-differences at t<=1.15;
(3) SENSITIVITY: real universal-cascade superstable s*_m (m=8,9,10) computed
    with each series (h fixed at N=72 -- s*h error negligible), gauged
    u* = s* delta^m / (-Ls_series): the shift of u* between series lengths is
    the error the H6 zero positions inherit from series truncation.
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr, pi, cos, matrix, lu_solve

mp.dps = 80
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FPJ = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A72 = [mpf(s) for s in FPJ['g_coeffs']]
C72 = [mpf(s) for s in FPJ['h_coeffs']]
AL = mpf(FPJ['alpha6'])
D6 = mpf(FPJ['delta6'])


def G(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def colloc_points(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def Gp(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def residuals_and_jacobian(a, alpha, pts):
    N = len(a) - 1
    npts = len(pts)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    a6 = alpha ** 6
    for i, t in enumerate(pts):
        ta = t / a6
        U = G(a, ta)
        Up = Gp(a, ta)
        W = U ** 6
        G_W = G(a, W)
        Gp_W = Gp(a, W)
        F[i] = G(a, t) + alpha * G_W
        base = alpha * Gp_W * 6 * U ** 5
        pow_t = mpf(1)
        pow_W = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            J[i, j - 1] = pow_t + alpha * pow_W + base * ta ** j
        J[i, N] = G_W - 36 * U ** 5 * t * Gp_W * Up / a6
    return F, J


def solve_fixedpoint(N, a_init, alpha_init):
    pts = colloc_points(N + 1)
    a = list(a_init) + [mpf(0)] * (N - (len(a_init) - 1))
    alpha = mpf(alpha_init)
    for it in range(100):
        F, J = residuals_and_jacobian(a, alpha, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < mpf('1e-60'):
            return a, alpha
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        for _h in range(30):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            F2, _ = residuals_and_jacobian(a_new, al_new, pts)
            nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            raise RuntimeError('stall')
    return a, alpha


log('solving N=84 and N=96 (dps 80)...')
A84, al84 = solve_fixedpoint(84, A72, AL)
log('N=84 done (alpha6=%s)' % nstr(al84, 30))
A96, al96 = solve_fixedpoint(96, A84, al84)
log('N=96 done (alpha6=%s)' % nstr(al96, 30))

log('--- tails ---')
for t in (mpf('0.9'), mpf('1.0'), mpf('1.05'), mpf('1.1'), mpf('1.15')):
    log('t=%s: |G72-G60|=%.2e |G84-G72|=%.2e |G96-G84|=%.2e |G96-G72|=%.2e'
        % (nstr(t, 5), float(fabs(G(A72, t) - G(A72[:61], t))),
           float(fabs(G(A84, t) - G(A72, t))),
           float(fabs(G(A96, t) - G(A84, t))),
           float(fabs(G(A96, t) - G(A72, t)))))

log('--- sensitivity: real universal cascade u*_m per series ---')


def h6(x, co):
    t = x ** 6
    r = co[-1]
    for c in reversed(co[:-1]):
        r = r * t + c
    return r


def s_star(a_ser, m):
    n = 2 ** (m + 1)

    def fres(s):
        x = mpf(0)
        ds = mpf(0)
        for _ in range(n):
            gp = 6 * x ** 5 * Gp(a_ser, x ** 6)
            hp = 6 * x ** 5 * Gp(C72, x ** 6)
            ds = (gp + s * hp) * ds + h6(x, C72)
            x = G(a_ser, x ** 6) + s * h6(x, C72)
        return x, ds
    s = mpf('-1e-8')
    for _ in range(60):
        x, ds = fres(s)
        if fabs(x) < mpf('1e-55'):
            break
        ss = x / ds
        improved = False
        for _h in range(30):
            x2, _ = fres(s - ss)
            if fabs(x2) < fabs(x):
                s = s - ss
                improved = True
                break
            ss = ss / 2
        if not improved:
            break
    return s


for m in (8, 9, 10):
    rows = {}
    for tag, ser in (('N72', A72), ('N84', A84), ('N96', A96)):
        sm1 = s_star(ser, m - 1)
        sm = s_star(ser, m)
        Ls_loc = sm * D6 ** m
        u_star = -sm1 * D6 ** (m - 1) / Ls_loc * 1  # gauge: u(s*_{m-1}) at level m-1
        rows[tag] = (sm1, Ls_loc, u_star)
    (s72, L72, u72) = rows['N72']
    (s84, L84, u84) = rows['N84']
    (s96, L96, u96) = rows['N96']
    log('m=%d: Ls N72=%s' % (m, nstr(L72, 20)))
    log('   u*(N72)=%.12f  u*(N84)=%.12f  u*(N96)=%.12f' % (u72, u84, u96))
    log('   shifts: |u72-u84|=%.2e |u84-u96|=%.2e |u72-u96|=%.2e'
        % (float(fabs(u72 - u84)), float(fabs(u84 - u96)), float(fabs(u72 - u96))))
log('diag_tail done')

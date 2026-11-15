#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag4_spectrum.py -- systematic mpmath sigma-scan of the exact DT collocation
matrix (N=48 and N=56) for the sextic class: locate ALL eigenvalues in
(0.01, 9.3), check N-stability, and certify the leading candidates by
function-space eigenpair residuals (dense grid, no matrix involved).
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, matrix, lu_solve, fabs, nstr, pi, cos

mp.dps = 50

T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FPJ = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FPJ['g_coeffs']]
ALPHA = mpf(FPJ['alpha6'])


def colloc_points(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def G_eval(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def g_of_x(a, x):
    return G_eval(a, x ** 6)


def gp_of_x(a, x):
    return 6 * x ** 5 * Gp_eval(a, x ** 6)


def DT_scalar(a, alpha, c, x):
    """(DT phi)(x) for phi(x) = sum_j c[j] x^{6j} (c[0] assumed 0)."""
    G1 = G_eval(a, mpf(1))

    def phi(z):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * (z ** 6) + cc
        return r

    phi1 = sum(c)
    gx = g_of_x(a, x)
    gpx = gp_of_x(a, x)
    u = g_of_x(a, x / alpha)
    gpu = gp_of_x(a, u)
    return (-(phi1 / G1) * (gx - x * gpx)
            - alpha * phi(u) - alpha * gpu * phi(x / alpha))


def DT_monomial(j, x):
    a, alpha = A_MP, ALPHA
    t = x ** 6
    G1 = G_eval(a, mpf(1))
    ta = t / alpha ** 6
    U = G_eval(a, ta)
    W = U ** 6
    gx = G_eval(a, t)
    gpx = 6 * x ** 5 * Gp_eval(a, t)
    gpu = 6 * U ** 5 * Gp_eval(a, W)
    return (-(1 / G1) * (gx - x * gpx)
            - alpha * (U ** (6 * j)) - alpha * gpu * ta ** j)


def build_A(N, a_series):
    pts = colloc_points(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        x = t ** (mpf(1) / 6)
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            F[i, j - 1] = DT_monomial(j, x)
    cols = []
    for jcol in range(N):
        b = matrix(N, 1)
        for i in range(N):
            b[i] = F[i, jcol]
        xj = lu_solve(V, b)
        cols.append([xj[i] for i in range(N)])
    A = matrix(N, N)
    for jcol in range(N):
        for i in range(N):
            A[i, jcol] = cols[jcol][i]
    return A


def matvec(M, v):
    n = len(v)
    w = matrix(n, 1)
    for i in range(n):
        s = mpc(0)
        for j in range(n):
            s += M[i, j] * v[j]
        w[i] = s
    return w


def inverse_iter(A, n, sigma, iters=50):
    sig = mpc(sigma)
    M = matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = A[i, j]
    for i in range(n):
        M[i, i] = A[i, i] - sig
    v = matrix(n, 1)
    for i in range(n):
        v[i] = mpc(1) / (i + 2)
    for it in range(iters):
        try:
            w = lu_solve(M, v)
        except Exception:  # noqa: BLE001
            return None, None, None
        nrm = max(abs(w[i]) for i in range(n))
        if nrm == 0 or not (nrm < mpf('1e40')):
            return None, None, None
        for i in range(n):
            w[i] = w[i] / nrm
        v = w
        if it in (8, 20):
            Av = matvec(A, v)
            sig = sum(v[i].conjugate() * Av[i] for i in range(n)) / sum(abs(v[i]) ** 2 for i in range(n))
            for i in range(n):
                M[i, i] = A[i, i] - sig
    Av = matvec(A, v)
    k = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[k] / v[k]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, v, resid


SIGMAS = [0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 4.5, 7.0]

log('building A48...')
A48 = build_A(48, A_MP)
log('building A56...')
A56 = build_A(56, A_MP)

found = {}
for N, A in ((48, A48), (56, A56)):
    eigs = []
    for sig in SIGMAS:
        lam, v, resid = inverse_iter(A, N, mpf(sig))
        if lam is None:
            log('N=%d sigma=%.2f: failed' % (N, sig))
            continue
        eigs.append((lam, v, resid, sig))
        log('N=%d sigma=%5.2f -> lam = %s  (resid %.1e)'
            % (N, sig, nstr(lam, 22), float(resid)))
    found[N] = eigs

# dedupe per N
def dedupe(eigs, tol_rel=0.01):
    out = []
    for lam, v, resid, sig in eigs:
        if all(abs(lam - l2) > tol_rel * max(abs(lam), abs(l2)) for l2, _v, _r, _s in out):
            out.append((lam, v, resid, sig))
    return out


log('--- deduped eigenvalues per N ---')
d48 = dedupe(found[48])
d56 = dedupe(found[56])
for tag, dd in (('N=48', d48), ('N=56', d56)):
    for lam, _v, _r, _s in sorted(dd, key=lambda z: -abs(z[0])):
        log('%s  %s' % (tag, nstr(lam, 22)))

log('--- N-stability pairing ---')
stable = []
for lam48, v48, r48, s48 in d48:
    best = min(d56, key=lambda z: abs(z[0] - lam48))
    diff = abs(best[0] - lam48)
    log('lam48=%s  nearest56=%s  |diff|=%.3e' % (nstr(lam48, 20), nstr(best[0], 20), float(diff)))
    stable.append((lam48, best[0], diff, v48))

log('--- function-space eigenpair residuals for |lam|>1e-3 pairs with |diff|/|lam|<1e-6 ---')
for lam48, lam56, diff, v48 in stable:
    if abs(lam48) < mpf('1e-3') or diff > abs(lam48) * mpf('1e-6'):
        continue
    # eigenvector -> t-series (monomial) -> function phi(x)=sum c_j x^{6j}
    c = [mpf(0)] + [v48[i] for i in range(48)]
    worst = mpf(0)
    for t in colloc_points(201):
        x = t ** (mpf(1) / 6)
        lhs = DT_scalar(A_MP, ALPHA, c, x)
        rhs = lam48 * sum(c[j] * x ** (6 * j) for j in range(1, 49))
        worst = max(worst, abs(lhs - rhs))
    log('lam=%s : dense-grid residual = %.3e' % (nstr(lam48, 20), float(worst)))

log('diag4 done')

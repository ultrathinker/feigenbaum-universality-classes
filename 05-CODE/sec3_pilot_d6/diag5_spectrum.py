#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag5_spectrum.py -- definitive eigenvalue hunt for the sextic DT operator:
(a) full mpmath eig of the collocation matrix for a ladder of small N;
(b) long targeted mpmath inverse iterations at N=48,56,64 for candidates;
(c) dense-grid FUNCTION-SPACE residual certification (no matrix involved) for
    every surviving eigenvector; N-stability across 48/56/64.
Also checks the rescaling-family numerology (alpha6^{-12} etc.).
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, matrix, lu_solve, fabs, nstr, pi, cos

mp.dps = 45

T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FPJ = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FPJ['g_coeffs']]
ALPHA = mpf(FPJ['alpha6'])


def colloc_points(npt):
    return [(1 - cos(pi * (i + mpf('0.5')) / npt)) / 2 for i in range(npt)]


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


def DT_scalar_x(a, alpha, c, x):
    G1 = G_eval(a, mpf(1))

    def phi(z):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * (z ** 6) + cc
        return r

    phi1 = sum(c)
    t = x ** 6
    gx = G_eval(a, t)
    gpx = 6 * x ** 5 * Gp_eval(a, t)
    u = G_eval(a, (x / alpha) ** 6)
    W = u ** 6
    gpu = 6 * u ** 5 * Gp_eval(a, W)
    return (-(phi1 / G1) * (gx - x * gpx)
            - alpha * phi(u) - alpha * gpu * phi(x / alpha))


def DT_monomial_col(j, x, a, alpha):
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


def build_A(N, a, alpha):
    pts = colloc_points(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        x = t ** (mpf(1) / 6)
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            F[i, j - 1] = DT_monomial_col(j, x, a, alpha)
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
        s = v[0] * 0
        for j in range(n):
            s += M[i, j] * v[j]
        w[i] = s
    return w


def inverse_iter_long(A, n, sigma, iters=250, refresh_at=(120, 180)):
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
        if nrm == 0 or not (nrm < mpf('1e50')):
            return None, None, None
        for i in range(n):
            w[i] = w[i] / nrm
        v = w
        if it in refresh_at:
            Av = matvec(A, v)
            sig = sum(v[i].conjugate() * Av[i] for i in range(n)) / sum(abs(v[i]) ** 2 for i in range(n))
            for i in range(n):
                M[i, i] = A[i, i] - sig
    Av = matvec(A, v)
    k = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[k] / v[k]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, v, resid


def dense_residual(c, lam, a, alpha, npts=161):
    worst = mpf(0)
    for t in colloc_points(npts):
        x = t ** (mpf(1) / 6)
        lhs = DT_scalar_x(a, alpha, c, x)
        rhs = lam * sum(c[j] * x ** (6 * j) for j in range(1, len(c)))
        worst = max(worst, abs(lhs - rhs))
    return worst


# ---------------- (a) full mpmath eig ladder ---------------------------------- #
log('--- (a) full eig ladder (mpmath) ---')
ladder_specs = {}
for N in (12, 16, 20, 24):
    A = build_A(N, A_MP, ALPHA)
    t_eig0 = time.time()
    try:
        E, _ = mp.eig(A)
        evs = sorted(E, key=lambda z: -abs(z))
        ladder_specs[N] = evs
        log('N=%d eig (%.0fs): %s'
            % (N, time.time() - t_eig0, '  '.join(nstr(e, 10) for e in evs[:8])))
    except Exception as exc:  # noqa: BLE001
        log('N=%d eig failed: %s' % (N, exc))

# ---------------- (b)+(c) targeted long inverse iterations --------------------- #
log('--- (b) long targeted inverse iterations at N=48,56,64 ---')
A48 = build_A(48, A_MP, ALPHA)
A56 = build_A(56, A_MP, ALPHA)
A64 = build_A(64, A_MP, ALPHA)
CANDS = [mpf('0.5820118020232834'), mpf('0.2078466910490824'),
         mpf('0.04654172867954866'), mpf('0.0215642988199406'),
         mpf('0.01000464627051392'), mpf('0.00215573818')]
results = {}
for cand in CANDS:
    row = {}
    for N, A in ((48, A48), (56, A56), (64, A64)):
        lam, v, resid = inverse_iter_long(A, N, cand)
        if lam is None:
            row[N] = None
            log('cand %.10f N=%d: FAILED' % (cand, N))
            continue
        c = [mpf(0)] + [v[i] for i in range(N)]
        dr = dense_residual(c, lam, A_MP, ALPHA)
        row[N] = (lam, resid, dr, v)
        log('cand %.10f N=%d: lam=%s matresid=%.1e dense_resid=%.1e'
            % (cand, N, nstr(lam, 20), float(resid), float(dr)))
    results[str(cand)] = row

log('--- summary: candidates with clean dense residuals at all N ---')
for cand, row in results.items():
    ok = all(row.get(N) is not None and float(row[N][2]) < 1e-15 for N in (48, 56, 64))
    if ok and all(row[N] is not None for N in (48, 56, 64)):
        d4856 = abs(row[48][0] - row[56][0])
        d5664 = abs(row[56][0] - row[64][0])
        log('cand %s: lam48=%s |d48-56|=%.2e |d56-64|=%.2e'
            % (cand, nstr(row[48][0], 20), float(d4856), float(d5664)))

log('--- numerology: alpha6^{-6k} family ---')
for k in (1, 2, 3, 4):
    log('alpha6^{-%d} = %s' % (6 * k, nstr(ALPHA ** (-6 * k), 20)))
log('alpha6^{-16} = %s' % nstr(ALPHA ** (-16), 20))
log('diag5 done')

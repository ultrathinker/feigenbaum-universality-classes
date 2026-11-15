#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag3_spectrum.py -- diagnostic round 2:
(1) DT collocation in the CHEBYSHEV basis T_j(2t-1) (well-conditioned V),
    float eig across N -- locate the physical N-stable eigenvalues;
(2) mpmath SHIFTED INVERSE iteration on the monomial-basis matrix (exact at
    dps 50) for the candidates -> polished eigenvalues + N-stability.
[oid:deg6B] NOT RIGOROUS -- float/mpmath diagnostic.
"""
import json
import time

import numpy as np
from mpmath import mp, mpf, matrix, lu_solve, fabs, nstr

mp.dps = 50

T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FPJ = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FPJ['g_coeffs']]
ALPHA = mpf(FPJ['alpha6'])
A_CO = np.array([float(s) for s in FPJ['g_coeffs']])
ALPHA_F = float(ALPHA)


def cheb_nodes(n):
    i = np.arange(n) + 0.5
    return (1.0 - np.cos(np.pi * i / n)) / 2.0


def cheb_nodes_mp(n):
    from mpmath import cos as mcos, pi as mpi
    return [(1 - mcos(mpi * (i + mpf('0.5')) / n)) / 2 for i in range(n)]


# ---------------- float Chebyshev-basis matrix ------------------------------- #
def G_eval(a, t):
    r = np.full_like(t, a[-1])
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    s = np.zeros_like(t)
    for j in range(1, len(a)):
        s = s * t + j * a[j]
    return s


def DT_of_phi(a, alpha, phi_vals_at, x):
    """DT(phi)(x), phi given as a callable; vectorised float."""
    t = x ** 6
    G1 = G_eval(a, np.array([1.0]))[0]
    ta = t / alpha ** 6
    U = G_eval(a, ta)
    W = U ** 6
    gx = G_eval(a, t)
    gpx = 6.0 * x ** 5 * Gp_eval(a, t)
    gpu = 6.0 * U ** 5 * Gp_eval(a, W)
    one = np.ones_like(x)
    return (-(phi_vals_at(one) / G1) * (gx - x * gpx)
            - alpha * phi_vals_at(U) - alpha * gpu * phi_vals_at((x / alpha)))


def build_B(N):
    """B = W^{-1} F, W_{ij} = T_j(2 t_i - 1) (Chebyshev basis), F_{ij} = DT(T_j)(t_i^{1/6})."""
    pts = cheb_nodes(N)
    x = pts ** (1.0 / 6.0)
    W = np.polynomial.chebyshev.chebvander(2 * pts - 1, N - 1)
    F = np.empty((N, N))
    for j in range(N):
        ph = np.polynomial.chebyshev.Chebyshev([0.0] * j + [1.0], domain=[-1, 1])
        F[:, j] = DT_of_phi(A_CO, ALPHA_F, lambda z, ph=ph: ph(2 * z - 1), x)
    return np.linalg.solve(W, F)


log('--- (1) Chebyshev-basis float spectra ---')
spectra = {}
for N in (24, 32, 40, 48, 56, 64):
    B = build_B(N)
    ev = np.linalg.eigvals(B)
    ev = ev[np.argsort(-np.abs(ev))]
    spectra[N] = ev
    log('N=%d top: %s' % (N, '  '.join('%.5f%+.5fi' % (e.real, e.imag) for e in ev[:7])))
    log('N=%d |e|<1.05: %s' % (N, '  '.join('%.5f%+.5fi' % (e.real, e.imag)
                                             for e in ev if abs(e) < 1.05)))

# ---------------- mpmath shifted inverse iteration ---------------------------- #
def G_eval_mp(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval_mp(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def DT_monomial_mp(j, x):
    """DT(phi)(x) for phi(x) = x^{6j}, mpmath scalar-exact, monomial basis."""
    a, alpha = A_MP, ALPHA
    t = x ** 6
    G1 = G_eval_mp(a, mpf(1))
    ta = t / alpha ** 6
    U = G_eval_mp(a, ta)
    W = U ** 6
    gx = G_eval_mp(a, t)
    gpx = 6 * x ** 5 * Gp_eval_mp(a, t)
    gpu = 6 * U ** 5 * Gp_eval_mp(a, W)
    return (-(1 / G1) * (gx - x * gpx)
            - alpha * (U ** (6 * j)) - alpha * gpu * ta ** j)


def build_A_mp(N):
    pts = cheb_nodes_mp(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        x = t ** (mpf(1) / 6)
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            F[i, j - 1] = DT_monomial_mp(j, x)
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
        s = mpf(0)
        for j in range(n):
            s += M[i, j] * v[j]
        w[i] = s
    return w


def inverse_iter(A, n, sigma, iters=60):
    """Shifted inverse power iteration with one Rayleigh refresh; complex-capable
    via mpc sigma. Returns (eigval, eigvec)."""
    from mpmath import mpc
    sig = sigma if not isinstance(sigma, (float, int)) else mpc(mpf(sigma))
    M = matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = A[i, j] - (sig if i == j else 0) * (1 if i == j else 0)
    for i in range(n):
        M[i, i] = A[i, i] - sig
    v = matrix(n, 1)
    for i in range(n):
        v[i] = mpc(1) / (i + 2)
    for it in range(iters):
        try:
            w = lu_solve(M, v)
        except Exception:  # noqa: BLE001
            return None, None
        nrm = max(abs(w[i]) for i in range(n))
        if nrm == 0:
            return None, None
        for i in range(n):
            w[i] = w[i] / nrm
        v = w
        if it in (10, 25):
            Av = matvec(A, v)
            # Rayleigh refresh
            sig = sum(v[i].conjugate() * Av[i] for i in range(n)) / sum(abs(v[i]) ** 2 for i in range(n))
            M = matrix(n, n)
            for i in range(n):
                for j in range(n):
                    M[i, j] = A[i, j]
            for i in range(n):
                M[i, i] = A[i, i] - sig
    Av = matvec(A, v)
    # eigval from componentwise ratios (median-ish): use max-|.| component
    k = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[k] / v[k]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, (v, resid)


log('--- (2) mpmath shifted inverse iteration, candidates from Chebyshev float N=48 ---')
ev48 = spectra[48]
cands = sorted([e for e in ev48 if 0.05 < abs(e) < 1.05], key=lambda e: -abs(e))
log('candidates (float N=48 Chebyshev): %s'
    % '  '.join('%.6f%+.6fi' % (e.real, e.imag) for e in cands[:8]))

A48 = build_A_mp(48)
log('built A48 (mpmath)')
seen = []
for e in cands[:8]:
    if any(abs(e - s) < 0.02 for s in seen):
        continue
    seen.append(e)
    from mpmath import mpc
    lam, vr = inverse_iter(A48, 48, mpc(mpf(float(e.real)), mpf(float(e.imag))))
    if lam is None:
        log('cand %.4f%+.4fi: inverse iteration failed' % (e.real, e.imag))
        continue
    log('cand %.4f%+.4fi -> lam = %s   resid=%s'
        % (e.real, e.imag, nstr(lam, 25), nstr(vr[1], 3)))
log('diag3 done')

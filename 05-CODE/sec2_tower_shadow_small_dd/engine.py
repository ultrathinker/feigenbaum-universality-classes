#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
engine_own.py -- [s3_F5_d8] independent numerics for BRIEF-S3-F5.
NOT RIGOROUS -- mpmath float, no interval arithmetic.

Own re-typed implementation of the SAME degree-d Cvitanovic-Feigenbaum
equation/linearisation used across the project (wave1_spectra_d/spectra_grid.py,
wave2_vector_homotopy/engine.py, wave2_check_d205/core.py, and MiniMax's
b1_4_F5_depth/engine.py which this task audits). Written independently here
(own function bodies, own optimisation: cached LU decomposition of the
collocation Vandermonde matrix reused across the N linear solves, instead of
calling lu_solve() N separate times -- turns an O(N^4) matrix build into
O(N^3), which is what makes N=56 at dps=40 cheap enough for a foreground run).

Equation (real degree d, t = |x|^d):
    g(x) = -alpha*g(g(x/alpha)),  g(0)=1,  G(t) = sum a_j t^j, a_0=1.
Linearised doubling operator on Phi(t), Phi(0)=0 (t-series c[0]=0):
    (DT Phi)(t) = -(Phi(1)/G(1))*K(t) - alpha*Phi(W)
                  - alpha*d*sign(U)*|U|^(d-1)*G'(W)*Phi(t/alpha^d)
    K(t) = G(t) - d*t*G'(t),  U = G(t/alpha^d),  W = |U|^d.
"""
import json
import os

from mpmath import mp, mpf, pi, cos, matrix
from mpmath import eig as mp_eig

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE1_JSON = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')


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


def colloc_points(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def build_DT_matrix(a, alpha, d, N):
    """A = V^{-1} F on span{t,...,t^N}. Single LU decomposition of V, reused
    for all N right-hand sides (F columns) via L_solve/U_solve -- O(N^3)
    instead of O(N^4) from calling lu_solve(V, .) N times."""
    pts = colloc_points(N)
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        W = aU ** d
        Kt = G_eval(a, t) - d * t * Gp_eval(a, t)
        gauge = -(1 / G1) * Kt
        c3 = -alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W)
        pow_t = mpf(1)
        pow_W = mpf(1)
        pow_ta = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            pow_ta *= ta
            V[i, j - 1] = pow_t
            F[i, j - 1] = gauge - alpha * pow_W + c3 * pow_ta
    LU, p = mp.LU_decomp(V)
    A = matrix(N, N)
    for jcol in range(N):
        b = matrix(N, 1)
        for i in range(N):
            b[i] = F[i, jcol]
        y = mp.L_solve(LU, b, p)
        x = mp.U_solve(LU, y)
        for i in range(N):
            A[i, jcol] = x[i]
    return A


def DT_scalar_t(a, alpha, d, c, t):
    """Scalar-exact (DT Phi)(t) for a t-series c (c[0]=0); ground truth for
    dense functional residual checks of any (lam, vec) candidate."""
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    ta = t / ad
    U = G_eval(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    Kt = G_eval(a, t) - d * t * Gp_eval(a, t)

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    Phi1 = sum(c)
    return (-(Phi1 / G1) * Kt - alpha * horner0(c, W)
            - alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W) * horner0(c, ta))


def eigenpair_dense_residual(a, alpha, d, lam, vec, npts=401):
    """Functional (non-truncated) residual of a candidate (lam, vec) pair,
    evaluated at npts collocation points independent of matrix size N used
    to produce (lam, vec) -- catches truncation artifacts."""
    c = [mpf(0)] + list(vec)
    worst = mpf(0)
    vnorm = mpf(0)
    for t in colloc_points(npts):
        lhs = DT_scalar_t(a, alpha, d, c, t)

        def horner0(cc, tt):
            r = cc[-1]
            for x in reversed(cc[:-1]):
                r = r * tt + x
            return r
        rhs = lam * horner0(c, t)
        worst = max(worst, abs(lhs - rhs))
        vnorm = max(vnorm, abs(horner0(c, t)))
    return worst / max(vnorm, mpf('1e-30'))


def load_wave1_points():
    with open(WAVE1_JSON, encoding='utf-8') as f:
        data = json.load(f)
    return {r['d']: r for r in data['points']}


def residuals_and_jacobian_fp(a, alpha, d, pts):
    N = len(a) - 1
    ad = alpha ** d
    F = [mpf(0)] * len(pts)
    J = []
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval(a, ta)
        Up = Gp_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        dWdU = d * sU * aU ** (d - 1)
        W = aU ** d
        G_W = G_eval(a, W)
        Gp_W = Gp_eval(a, W)
        F[i] = G_eval(a, t) + alpha * G_W
        base = alpha * Gp_W * dWdU
        row = [mpf(0)] * (N + 1)
        pow_t = mpf(1)
        pow_W = mpf(1)
        pow_ta = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            pow_ta *= ta
            row[j - 1] = pow_t + alpha * pow_W + base * pow_ta
        row[N] = G_W - d * t * dWdU * Gp_W * Up / ad
        J.append(row)
    return F, J


def newton_fp(d, N, a0, al0, tol=mpf('1e-32'), max_it=60):
    from mpmath import lu_solve as mp_lu_solve
    pts = colloc_points(N + 1)
    a = list(a0) + [mpf(0)] * (N + 1 - len(a0))
    alpha = mpf(al0)
    for it in range(max_it):
        F, J = residuals_and_jacobian_fp(a, alpha, d, pts)
        nrm = max(abs(x) for x in F)
        if nrm < tol:
            return a, alpha, it, nrm
        Jm = matrix(J)
        rhs = matrix([-x for x in F])
        dx = mp_lu_solve(Jm, rhs)
        lam = mpf(1)
        for _h in range(30):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = residuals_and_jacobian_fp(a_new, al_new, d, pts)
                nrm2 = max(abs(x) for x in F2)
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            if nrm < mpf('1e-22'):
                return a, alpha, it, nrm
            raise RuntimeError('newton stalled, res %s' % nrm)
    if nrm < mpf('1e-22'):
        return a, alpha, max_it, nrm
    raise RuntimeError('newton no convergence, res %s' % nrm)


def solve_fp_continuation(d, N, a_seed, al_seed, ladder=(10, 18, 28, 40)):
    a_cur, al_cur = a_seed, al_seed
    for Nl in list(ladder) + [N]:
        Nl_use = min(Nl, N)
        a_trunc = a_cur[:Nl_use + 1] + [mpf(0)] * max(0, Nl_use + 1 - len(a_cur))
        a_cur, al_cur, it, r = newton_fp(mpf(d), Nl_use, a_trunc, al_cur)
        if Nl_use == N:
            break
    if len(a_cur) < N + 1:
        a_cur, al_cur, it, r = newton_fp(mpf(d), N, a_cur, al_cur)
    return a_cur, al_cur, it, r


def fp_residual_dense(a, alpha, d, npts=801):
    ad = alpha ** d
    worst = mpf(0)
    for t in colloc_points(npts):
        U = G_eval(a, t / ad)
        r = G_eval(a, t) + alpha * G_eval(a, abs(U) ** d)
        worst = max(worst, abs(r))
    return worst


def full_eig_window(a, alpha, d, N, window=0.5, imag_tol=mpf('1e-10')):
    """Full mpmath eig (right AND left eigenvectors), filtered to
    |Re(lam)| <= window. Complex pairs with |Im| > |lam|*imag_tol are KEPT
    (flagged complex=True) rather than silently dropped, so genuine complex
    modes are not lost by a real-only filter.
    Returns list of dicts sorted by |lam| descending."""
    A = build_DT_matrix(a, alpha, d, N)
    E, EL, ER = mp_eig(A, left=True, right=True)
    out = []
    for k in range(N):
        lam = E[k]
        is_complex = abs(lam.imag) > abs(lam) * imag_tol + mpf('1e-28')
        lamval = lam if is_complex else lam.real
        if abs(lam.real) > window and not is_complex:
            continue
        if is_complex and abs(lam) > window:
            continue
        rvec = [ER[i, k] for i in range(N)]
        lvec = [EL[i, k] for i in range(N)]
        out.append({'lam': lamval, 'rvec': rvec, 'lvec': lvec, 'complex': is_complex})
    out.sort(key=lambda r: -abs(r['lam']))
    return out


def cosine_overlap(v1, v2):
    """Cosine between two REAL coefficient vectors in the shared {t,...,t^N}
    monomial basis (pad shorter with zeros)."""
    n = max(len(v1), len(v2))
    x1 = list(v1) + [mpf(0)] * (n - len(v1))
    x2 = list(v2) + [mpf(0)] * (n - len(v2))
    n1 = sum(x * x for x in x1) ** mpf('0.5')
    n2 = sum(x * x for x in x2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(x1[i] * x2[i] for i in range(n))
    return abs(dot / (n1 * n2))

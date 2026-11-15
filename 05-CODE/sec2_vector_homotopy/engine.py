#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
engine.py -- [oid:vecM] BRIEF-M shared numerics.
NOT RIGOROUS -- mpmath float, no interval arithmetic.

Reuses the SAME equation/linearisation as wave1_spectra_d/spectra_grid.py
(read-only reuse of the method; code re-typed here so this folder writes only
to itself). Fixed points (a48, alpha48) are loaded from the already-computed,
already-certified wave1_spectra_d/spectra_grid.json (no re-solving of Newton's
method -- saves the expensive part; only the LINEAR spectral step is redone
here, independently, with full eig (right+left) instead of targeted inverse
iteration, so that eigenVECTOR continuity can be checked explicitly).

Equation (degree-d Cvitanovic-Feigenbaum, real collocation, t=|x|^d):
    g(x) = -alpha*g(g(x/alpha)),  g(0)=1,  G(t)=sum a_j t^j, a_0=1.
Linearised doubling operator on Phi(0)=0 perturbations (t-series c[0]=0):
    (DT Phi)(t) = -(Phi(1)/G(1))*K(t) - alpha*Phi(W)
                  - alpha*d*sign(U)*|U|^(d-1)*G'(W)*Phi(t/alpha^d)
    K(t) = G(t) - d*t*G'(t),  U = G(t/alpha^d),  W = |U|^d.
"""
import json
import os

from mpmath import mp, mpf, pi, cos, matrix, lu_solve
from mpmath import eig as mp_eig

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE1_JSON = os.path.join(HERE, '..', 'wave1_spectra_d', 'spectra_grid.json')


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
    """A = V^{-1} F on span{t,...,t^N}; coefficient vectors of Phi are directly
    comparable across different d (same abstract monomial basis)."""
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
    A = matrix(N, N)
    for jcol in range(N):
        b = matrix(N, 1)
        for i in range(N):
            b[i] = F[i, jcol]
        x = lu_solve(V, b)
        for i in range(N):
            A[i, jcol] = x[i]
    return A


def DT_scalar_t(a, alpha, d, c, t):
    """Scalar-exact (DT Phi)(t) for a t-series c (c[0]=0); ground truth used
    for dense functional residuals of any (lam, vec) candidate."""
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


def dense_eigen_residual(a, alpha, d, hcoef, lam, npts=121):
    worst = mpf(0)
    phmax = mpf(0)
    for t in colloc_points(npts):

        def horner0(cc, tt):
            r = cc[-1]
            for x in reversed(cc[:-1]):
                r = r * tt + x
            return r

        lhs = DT_scalar_t(a, alpha, d, hcoef, t)
        rhs = lam * horner0(hcoef, t)
        worst = max(worst, abs(lhs - rhs))
        phmax = max(phmax, abs(horner0(hcoef, t)))
    return worst, phmax


def load_wave1_points():
    with open(WAVE1_JSON, encoding='utf-8') as f:
        data = json.load(f)
    return {r['d']: r for r in data['points']}


def residuals_and_jacobian_fp(a, alpha, d, pts):
    """Newton residual+Jacobian for the fixed-point equation (copied from
    wave1_spectra_d/spectra_grid.py's residuals_and_jacobian; needed here so
    this folder can solve NEW d values -- e.g. bridge refinement -- without
    writing into wave1's folder)."""
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


def newton_fp(d, N, a0, al0, tol=mpf('1e-28'), max_it=60):
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
        dx = lu_solve(Jm, rhs)
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
            if nrm < mpf('1e-20'):
                return a, alpha, it, nrm  # accept: at the dps noise floor already
            raise RuntimeError('newton stalled, res %s' % nrm)
    if nrm < mpf('1e-20'):
        return a, alpha, max_it, nrm
    raise RuntimeError('newton no convergence, res %s' % nrm)


def solve_fp_continuation(d, N, a_seed, al_seed):
    """Warm-started Newton continuation from a nearby (a_seed, al_seed) fixed
    point to the new target d, ladder N: 10 -> 18 -> 28 -> N (mirrors wave1's
    solve_fp warm-start ladder)."""
    a_cur, al_cur = a_seed, al_seed
    for Nl in (10, 18, 28, N):
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


def fixedpoint_from_wave1(rec):
    """Full-precision (49-coeff, N=48) fixed point a(t), alpha, as solved and
    dense-residual-certified by wave1 (fp_dense_resid_N48 typically <=1e-30).
    Used as-is for ALL matrix truncations N here (24/32/40/48): the nonlinear
    background solution is already converged; only the LINEAR operator matrix
    is truncated at various N -- isolates N-dependence of the spectrum from
    N-dependence of the fixed point (simplification vs wave1, noted in report)."""
    a48 = [mpf(s) for s in rec['_a48']]
    al48 = mpf(rec['alpha48'])
    return a48, al48


def full_eig_window(a, alpha, d, N, window=0.42, imag_tol=1e-12):
    """Full mpmath eig (right AND left eigenvectors), filtered to
    |Re(lam)| <= window and |Im(lam)| <= |lam|*imag_tol (near-real only).
    Returns list of dicts sorted by |lam| descending, each with
    lam (mpf), rvec (list[mpf]), lvec (list[mpf]) both length N, raw
    (unnormalised) coefficient vectors in the {t,...,t^N} monomial basis."""
    A = build_DT_matrix(a, alpha, d, N)
    E, EL, ER = mp_eig(A, left=True, right=True)
    out = []
    for k in range(N):
        lam = E[k]
        if abs(lam.imag) > abs(lam) * mpf(imag_tol) + mpf('1e-30'):
            continue
        lamr = lam.real
        if abs(lamr) > window:
            continue
        rvec = [ER[i, k].real for i in range(N)]
        lvec = [EL[i, k].real if hasattr(EL[i, k], 'real') else EL[i, k]
                for i in range(N)]
        out.append({'lam': lamr, 'rvec': rvec, 'lvec': lvec})
    out.sort(key=lambda r: -abs(r['lam']))
    return out


def cosine_overlap(v1, v2):
    n1 = sum(x * x for x in v1) ** mpf('0.5')
    n2 = sum(x * x for x in v2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(v1[i] * v2[i] for i in range(len(v1)))
    return abs(dot / (n1 * n2))


def biorth_coeff(lvec_i, vec_target, rvec_i):
    """Expansion coefficient of vec_target along right-eigvec rvec_i using the
    LEFT eigenvector lvec_i as dual (biorthogonal normalisation <l_i,r_i>=1
    applied here), appropriate for a NON-NORMAL operator: c_i = <l_i, target> /
    <l_i, r_i>. |c_i| close to 1 (with other targets' |c_j| small) signals
    'target IS (mostly) mode i' even when raw cosine is distorted by
    non-orthogonality."""
    denom = sum(lvec_i[k] * rvec_i[k] for k in range(len(rvec_i)))
    if denom == 0:
        return None
    numer = sum(lvec_i[k] * vec_target[k] for k in range(len(vec_target)))
    return numer / denom

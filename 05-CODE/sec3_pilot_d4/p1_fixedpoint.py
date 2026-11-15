#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p1_fixedpoint.py -- Phase 1 of the quartic universality-class float pilot.
[oid:quartK5]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

Solves the quartic Cvitanovic-Feigenbaum equation (BRIEF convention)

    g(x) = -alpha4 * g(g(x/alpha4)),   g(0) = 1,

where g is an even function of x^4:  g(x) = G(t), t = x^4,  G(t) = sum_j a_j t^j,
a_0 = 1.  NOTE: since g depends on x^4 only, g(-z)=g(z) and the sign of the inner
argument is irrelevant; the sign of alpha4 is fixed by g(1) = -1/alpha4.

Method: collocation at Chebyshev-type nodes on (0,1] + Newton with analytic
Jacobian (mpmath, dps=60); continuation in the number of modes N.

Then: linearisation of the doubling operator on normalised quartic-symmetric
perturbations  R(f)(x) = -alpha_f f(f(x/alpha_f)), alpha_f = -1/f(1):

    (DT phi)(x) = -(phi(1)/g(1)) [g(x) - x g'(x)]
                  - alpha4 * phi(g(x/alpha4)) - alpha4 * g'(g(x/alpha4)) * phi(x/alpha4)

which preserves phi(0)=0 and the "function of x^4" subspace.  Built as a matrix on
t-series coefficients; leading eigenvalue -> delta4, eigenfunction -> h4
(normalised h(1)=1).  Validated by dense-grid SCALAR residuals (no series-matrix
involved in the check).

In t = x^4 variable the operator reads (Phi = phi as t-series):
    K(t)  = G(t) - 4 t G'(t)
    U(t)  = G(t/alpha4^4)
    W(t)  = U(t)^4
    (DT Phi)(t) = -(Phi(1)/G(1)) K(t) - alpha4 * Phi(W(t))
                  - 4 alpha4 U(t)^3 G'(W(t)) * Phi(t/alpha4^4)

Self-consistency checks reported (BRIEF: "самосогласованностью своего счёта"):
  (1) fixed-point residual on dense grid,
  (2) N-convergence ladder of alpha4,
  (3) eigenpair scalar residual  max |DT(h)(x) - delta4 h(x)| on dense grid,
  (4) literature sanity anchors delta4 ~ 7.2847, alpha4 ~ 1.6903 (anchors only,
      the operative check is the M-side cascade in p3).

Output: quartic_fixedpoint.json
"""
import json
import time

import numpy as np
from mpmath import mp, mpf, mpmathify, nstr, pi, cos, matrix, lu_solve, fabs

mp.dps = 60

TAG = "[oid:quartK5] NOT RIGOROUS float pilot, phase 1"
T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


# --------------------------------------------------------------------------- #
#  series / evaluation helpers (t-series of length N+1, a[0]=1)
# --------------------------------------------------------------------------- #
def G_eval(a, t):
    """Horner evaluation of G(t) with coefficient list a[0..N]."""
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    """G'(t)."""
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def colloc_points(npts):
    """npts Chebyshev-type nodes inside (0,1), clustered at both ends."""
    pts = []
    for i in range(npts):
        pts.append((1 - cos(pi * (i + mpf('0.5')) / npts)) / 2)
    return pts


# --------------------------------------------------------------------------- #
#  fixed point: unknowns a[1..N] (+ a0=1 implicit) and alpha
# --------------------------------------------------------------------------- #
def residuals_and_jacobian(a, alpha, pts):
    """F_i = G(t_i) + alpha*G(W_i), W_i = U_i^4, U_i = G(t_i/alpha^4).
    Returns (F, J) with J (npts x npts+1) over unknowns (a_1..a_N, alpha).
    npts = N+1."""
    N = len(a) - 1
    npts = len(pts)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    a4 = alpha ** 4
    for i, t in enumerate(pts):
        ta = t / a4
        U = G_eval(a, ta)
        Up = Gp_eval(a, ta)
        W = U ** 4
        G_W = G_eval(a, W)
        Gp_W = Gp_eval(a, W)
        F[i] = G_eval(a, t) + alpha * G_W
        # d/da_j : t^j + alpha*( W^j + Gp_W * 4 U^3 * ta^j )
        base = alpha * Gp_W * 4 * U ** 3
        pow_t = mpf(1)
        pow_W = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            J[i, j - 1] = pow_t + alpha * pow_W + base * ta ** j
        # d/dalpha : G(W) - 16 U^3 t G'(W) G'(ta) / alpha^4
        J[i, N] = G_W - 16 * U ** 3 * t * Gp_W * Up / a4
    return F, J


def solve_fixedpoint(N, a_init=None, alpha_init=None):
    """Newton with damping; returns (a list, alpha) or raises."""
    pts = colloc_points(N + 1)
    if a_init is None:
        a = [mpf(1)] + [mpf(0)] * N
        # crude quartic-family-motivated start: g ~ 1 - 1.59 x^4 -> a_1 = -1.59
        if N >= 1:
            a[1] = mpf('-1.59')
        alpha = mpf('1.69')
    else:
        a = list(a_init) + [mpf(0)] * (N - (len(a_init) - 1))
        alpha = mpf(alpha_init)
    for it in range(200):
        F, J = residuals_and_jacobian(a, alpha, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < mpf('1e-48'):
            return a, alpha, it, nrm
        # solve J dx = -F (square system npts x npts with last col alpha)
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        try:
            dx = lu_solve(J, rhs)
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError('lu_solve failed: %s' % exc)
        # damped update
        lam = mpf(1)
        for _halve in range(30):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = residuals_and_jacobian(a_new, al_new, pts)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:  # noqa: BLE001
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            raise RuntimeError('Newton stalled, residual %s' % nstr(nrm, 5))
    raise RuntimeError('Newton did not converge, residual %s' % nstr(nrm, 5))


def fp_residual_dense(a, alpha, npts=4001):
    pts = colloc_points(npts)
    a4 = alpha ** 4
    worst = mpf(0)
    for t in pts:
        U = G_eval(a, t / a4)
        r = G_eval(a, t) + alpha * G_eval(a, U ** 4)
        worst = max(worst, abs(r))
    return worst


# --------------------------------------------------------------------------- #
#  scalar DT evaluation (ground truth for checks; x complex-capable)
# --------------------------------------------------------------------------- #
def g_of_x(a, x):
    """g(x) for x complex/plain scalar: g(x) = G(x^4)."""
    return G_eval(a, x ** 4)


def gp_of_x(a, x):
    """g'(x) = 4 x^3 G'(x^4)."""
    return 4 * x ** 3 * Gp_eval(a, x ** 4)


def DT_scalar(a, alpha, c, x):
    """(DT phi)(x) for phi(x) = sum_j c[j] x^{4j} (c[0] assumed 0).
    DT phi = -(phi(1)/g(1)) [g(x) - x g'(x)] - alpha phi(g(x/alpha))
             - alpha g'(g(x/alpha)) phi(x/alpha)."""
    G1 = G_eval(a, mpf(1))

    def phi(z):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * (z ** 4) + cc
        return r

    phi1 = sum(c)
    gx = g_of_x(a, x)
    gpx = gp_of_x(a, x)
    u = g_of_x(a, x / alpha)
    gpu = gp_of_x(a, u)
    return (-(phi1 / G1) * (gx - x * gpx)
            - alpha * phi(u) - alpha * gpu * phi(x / alpha))


def DT_unprojected_scalar(a, alpha, c, x):
    """Bare (no normalisation projection) linearisation term; used only for the
    coordinate-mode diagnostic."""
    def phi(z):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * (z ** 4) + cc
        return r

    u = g_of_x(a, x / alpha)
    gpu = gp_of_x(a, u)
    return -alpha * phi(u) - alpha * gpu * phi(x / alpha)


def R_of_f(f, alpha_f, x):
    """Doubling operator R(f)(x) = -alpha_f f(f(x/alpha_f)) for plain scalar f."""
    return -alpha_f * f(f(x / alpha_f))


# --------------------------------------------------------------------------- #
#  collocation matrix of DT on the span{t, t^2, ..., t^N} (t = x^4)
# --------------------------------------------------------------------------- #
def build_DT_collocation(a, alpha, N, xpower=mpf('0.9')):
    """A = V^{-1} F,  V_{ij} = t_i^j  (i=1..N nodes, j=1..N),  F_{ij} = DT(e_j)(t_i^{1/4}).
    Scalar-exact operator evaluation at the nodes (no series truncation of the
    operator image); V is inverted in mpmath (conditioning ~ (1+sqrt2)^N is fine
    at dps 60).  Eigenproblem of A -> spectrum of DT on normalised quartic modes."""
    pts_t = colloc_points(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts_t):
        x = t ** (mpf(1) / 4)
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            e_j = [mpf(0)] * (j + 1)
            e_j[j] = mpf(1)
            F[i, j - 1] = DT_scalar(a, alpha, e_j, x)
    # A = V^{-1} F, columnwise (mpmath lu_solve takes vector RHS)
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


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def main():
    log(TAG)
    log('mp.dps = %d' % mp.dps)

    # ---- continuation ladder in N ------------------------------------------------
    ladder = [4, 8, 14, 20, 28, 36, 44]
    sols = []
    a_cur, al_cur = None, None
    for N in ladder:
        try:
            a, al, iters, res = solve_fixedpoint(N, a_cur, al_cur)
        except RuntimeError as exc:
            # fall back to crude start
            log('N=%d continuation failed (%s); retrying crude start' % (N, exc))
            a, al, iters, res = solve_fixedpoint(N)
        sols.append((N, a, al))
        log('N=%2d: alpha4 = %s   (newton iters %d, colloc res %.2e)'
            % (N, mpmathify(al), iters, float(res)))
        a_cur, al_cur = a, al

    # ---- final solve at N=48 (and 56 for convergence estimate) -------------------
    N1 = 48
    a1, al1, it1, res1 = solve_fixedpoint(N1, a_cur, al_cur)
    N2 = 56
    a2, al2, it2, res2 = solve_fixedpoint(N2, a1, al1)
    log('N=%d alpha4=%s' % (N1, mpmathify(al1)))
    log('N=%d alpha4=%s' % (N2, mpmathify(al2)))
    d_alpha = abs(al1 - al2)
    d_coef = max(abs(a1[j] - a2[j]) for j in range(1, N1 + 1))
    log('N-convergence: |alpha(N48)-alpha(N56)| = %.3e , max|da| = %.3e'
        % (float(d_alpha), float(d_coef)))

    dense_res = fp_residual_dense(a1, al1, 4001)
    dense_res2 = fp_residual_dense(a2, al2, 4001)
    log('dense residual (N=48): %.3e   (N=56): %.3e'
        % (float(dense_res), float(dense_res2)))

    # radius-of-convergence estimate (t-plane)
    Rat = abs(a1[-1]) ** (-mpf(1) / N1)
    tail = [abs(a1[j]) for j in range(N1 - 8, N1 + 1)]
    Rat2 = (tail[0] / tail[-1]) ** (mpf(1) / 8)
    log('series radius estimate (t-plane): last-coeff %.6f, tail-ratio %.6f'
        % (float(Rat), float(Rat2)))

    g1 = G_eval(a1, mpf(1))
    log('g(1) = %s  (should equal -1/alpha4 = %s)' % (mpmathify(g1), mpmathify(-1 / al1)))

    # ---- validate DT formula by finite differences of R -------------------------
    eps = mpf(10) ** (-25)
    c_test = [mpf(0)] + [mpf(1) / (j + 1) for j in range(1, 12)]

    def f_eps(x):
        return g_of_x(a1, x) + eps * (sum(c_test[j] * x ** (4 * j) for j in range(1, 12)))

    alpha_eps = -1 / f_eps(mpf(1))
    fd_worst = mpf(0)
    for t in colloc_points(31):
        x = t ** (mpf(1) / 4)
        fd = (R_of_f(f_eps, alpha_eps, x) - g_of_x(a1, x)) / eps
        an = DT_scalar(a1, al1, c_test, x)
        fd_worst = max(fd_worst, abs(fd - an))
    log('DT formula vs finite difference of R (eps=1e-25): max diff = %.3e' % float(fd_worst))

    # ---- coordinate-mode diagnostic (unprojected DT on N = G - x g') -------------
    Kcoef = [a1[0]] + [a1[j] * (1 - 4 * j) for j in range(1, N1 + 1)]
    coord_ratios = []
    for t in colloc_points(11):
        x = t ** (mpf(1) / 4)
        num = DT_unprojected_scalar(a1, al1, Kcoef, x)
        den = sum(Kcoef[j] * x ** (4 * j) for j in range(0, N1 + 1))
        coord_ratios.append(num / den)
    coord_c = coord_ratios[len(coord_ratios) // 2]
    coord_spread = max(coord_ratios) - min(coord_ratios)
    log('coordinate mode DT_unproj(K)/K ~ %s (spread %.2e); for reference alpha4^2 = %.9f, alpha4^4 = %.9f'
        % (nstr(coord_c, 18), float(coord_spread), float(al1 ** 2), float(al1 ** 4)))

    # ---- linearisation spectrum (collocation matrix) -----------------------------
    # NOTE: float eig of A is ill-conditioned (Vandermonde basis); the authoritative
    # values are the mpmath power-iteration eigenpairs below, validated by the
    # dense-grid scalar residual.  The float spectrum is kept out of the JSON.
    A = build_DT_collocation(a1, al1, N1)

    # ---- leading eigenpairs by mpmath power iteration with deflation --------------
    def matvec(M, v):
        w = matrix(len(v), 1)
        for i in range(len(v)):
            s = mpf(0)
            for j in range(len(v)):
                s += M[i, j] * v[j]
            w[i] = s
        return w

    def power_iter(M, n, iters=600, tol=mpf('1e-45')):
        v = matrix(n, 1)
        for i in range(n):
            v[i] = mpf(1) / (i + 3)
        lam = mpf(1)
        for it in range(iters):
            w = matvec(M, v)
            lam_new = max(abs(w[i]) for i in range(n))
            if lam_new == 0:
                return mpf(0), v
            for i in range(n):
                v[i] = w[i] / lam_new
            if abs(lam_new - lam) < tol * abs(lam_new):
                lam = lam_new
                return lam, v, it
            lam = lam_new
        return lam, v, iters

    def deflate(M, lam1, v1, n):
        """M - lam1 * v1 w^T/(w^T v1) using left eigenvector w from power iter on M^T."""
        MT = matrix(n, n)
        for i in range(n):
            for j in range(n):
                MT[i, j] = M[j, i]
        _lw, w, _it = power_iter(MT, n)
        denom = sum(w[i] * v1[i] for i in range(n))
        D = matrix(n, n)
        for i in range(n):
            for j in range(n):
                D[i, j] = M[i, j] - lam1 * v1[i] * w[j] / denom
        return D

    lam1, v1, it1_ = power_iter(A, N1)
    log('power iteration: delta4 = %s  (%d iters)' % (nstr(lam1, 40), it1_))
    # normalise h: h(1) = 1
    s1 = sum(v1[i] for i in range(N1))
    hcoef = [mpf(0)] + [v1[i] / s1 for i in range(N1)]
    log('h4 normalised h(1)=1; h coeffs 1..5: %s'
        % [nstr(hcoef[j], 12) for j in range(1, 6)])

    D1 = deflate(A, lam1, v1, N1)
    lam2, _v2, _it2 = power_iter(D1, N1)
    log('second eigenvalue (deflated): %s' % nstr(lam2, 20))
    D2 = deflate(D1, lam2, _v2, N1)
    lam3, _v3, _it3 = power_iter(D2, N1)
    log('third eigenvalue (deflated): %s' % nstr(lam3, 20))

    # ---- eigenpair scalar residual check -----------------------------------------
    worst = mpf(0)
    for t in colloc_points(201):
        x = t ** (mpf(1) / 4)
        lhs = DT_scalar(a1, al1, hcoef, x)
        rhs = lam1 * sum(hcoef[j] * x ** (4 * j) for j in range(1, N1 + 1))
        worst = max(worst, abs(lhs - rhs))
    log('eigenpair scalar residual max on 201 pts: %.3e' % float(worst))

    # ---- literature anchors -------------------------------------------------------
    lit_delta = mpf('7.284686217')
    lit_alpha = mpf('1.690302971')
    log('anchors: |delta4-7.284686217| = %.3e ; |alpha4-1.690302971| = %.3e'
        % (float(abs(lam1 - lit_delta)), float(abs(al1 - lit_alpha))))

    # ---- output --------------------------------------------------------------------
    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'mp_dps': mp.dps, 'equation': 'g(x) = -alpha4*g(g(x/alpha4)), g(0)=1, g in x^4',
                 'elapsed_sec': time.time() - T0},
        'alpha4': nstr(al1, 50),
        'delta4': nstr(lam1, 50),
        'N': N1,
        'g_coeffs': [nstr(x, 50) for x in a1],
        'h_coeffs': [nstr(x, 50) for x in hcoef],  # h(x) = sum h_coeffs[j] * x^{4j}
        'g1_check': nstr(g1, 30),
        'checks': {
            'alpha4_N48': nstr(al1, 40), 'alpha4_N56': nstr(al2, 40),
            'alpha4_N_diff': float(d_alpha), 'g_coeffs_N_diff_max': float(d_coef),
            'residual_dense_N48': float(dense_res), 'residual_dense_N56': float(dense_res2),
            'DT_formula_vs_finite_difference': float(fd_worst),
            'coordinate_mode_ratio': nstr(coord_c, 20),
            'coordinate_mode_spread': float(coord_spread),
            'eigenpair_scalar_residual': float(worst),
            'eigen2_deflated': nstr(lam2, 20), 'eigen3_deflated': nstr(lam3, 20),
            'radius_t_lastcoeff': float(Rat), 'radius_t_tailratio': float(Rat2),
            'alpha_ladder': {str(N): nstr(al, 30) for N, _a, al in sols},
            'note_float_spectrum': 'float eig of collocation matrix ill-conditioned; only mpmath power-iteration eigenvalues reported',
            'anchor_delta4_diff': float(abs(lam1 - lit_delta)),
            'anchor_alpha4_diff': float(abs(al1 - lit_alpha)),
        },
    }
    with open('quartic_fixedpoint.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote quartic_fixedpoint.json')


if __name__ == '__main__':
    main()

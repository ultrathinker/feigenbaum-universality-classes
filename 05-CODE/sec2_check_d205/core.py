#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
core.py -- [K-D205-CHECK] NOT RIGOROUS float exploration (mpmath), no interval
arithmetic. Blind independent re-derivation (no access to wave1_spectra_d).

Cvitanovic-Feigenbaum equation for the family f_a(x) = 1 - a|x|^d, generalized
to real (not necessarily integer) degree d:

    g(x) = -alpha_d * g(g(x/alpha_d)),   g(0) = 1,   g depends on |x|^d.

Write t = |x|^d,  g(x) = G(t),  G(t) = sum_{j=0}^N a_j t^j,  a_0 = 1.

Method mirrors the (permitted, read) quartic_pilot/p1_fixedpoint.py convention
(degree-4 pilot), generalized from integer degree 4 to real degree d by
replacing x^4 -> |x|^d, 4x^3 -> d*sign(x)*|x|^(d-1) everywhere derivatives of
the absolute-power appear. Derivation redone from scratch and checked against
the d=2.0 literature anchors (see calibrate.py) -- NOT copied from any degree-2
result file.

Collocation + Newton (mpmath, dps=60), continuation in N and in d.

Linearisation of the doubling operator on normalised perturbations phi (phi(0)=0):
  R(f)(x) = -alpha_f f(f(x/alpha_f)),  alpha_f fixed by f(0)=1 self-consistency.
  (DT phi)(x) = -(phi(1)/g(1)) [g(x) - x g'(x)]
                - alpha*phi(g(x/alpha)) - alpha*g'(g(x/alpha))*phi(x/alpha)
In t-representation (Phi_t(w) = sum_j c_j w^j, c_0=0):
  K(t)  = G(t) - d*t*G'(t)
  ta    = t/alpha^d
  U     = G(ta)                      (real, can be negative)
  W     = |U|^d
  (DT Phi)(t) = -(Phi_t(1)/G(1))*K(t) - alpha*Phi_t(W)
                - d*alpha*sign(U)*|U|^(d-1)*G'(W)*Phi_t(ta)
"""
import time
from mpmath import mp, mpf, nstr, pi, cos, matrix, lu_solve, fabs, eig

mp.dps = 60


def log(fh, t0, msg):
    line = "[%8.1fs] %s" % (time.time() - t0, msg)
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def sgn(x):
    if x > 0:
        return mpf(1)
    if x < 0:
        return mpf(-1)
    return mpf(0)


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
    pts = []
    for i in range(npts):
        pts.append((1 - cos(pi * (i + mpf('0.5')) / npts)) / 2)
    return pts


def residuals_and_jacobian(a, alpha, pts, d):
    N = len(a) - 1
    npts = len(pts)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    ad = alpha ** d
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval(a, ta)
        Up = Gp_eval(a, ta)
        sU = sgn(U)
        absU = fabs(U)
        W = absU ** d
        G_W = G_eval(a, W)
        Gp_W = Gp_eval(a, W)
        F[i] = G_eval(a, t) + alpha * G_W
        base = alpha * Gp_W * d * sU * (absU ** (d - 1) if absU != 0 else mpf(0))
        pow_t = mpf(1)
        pow_W = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            J[i, j - 1] = pow_t + alpha * pow_W + base * ta ** j
        J[i, N] = G_W - d * d * sU * (absU ** (d - 1) if absU != 0 else mpf(0)) * Gp_W * Up * t / ad
    return F, J


def solve_fixedpoint(N, d, a_init=None, alpha_init=None, tol=mpf('1e-50'), maxit=200):
    pts = colloc_points(N + 1)
    if a_init is None:
        a = [mpf(1)] + [mpf(0)] * N
        if N >= 1:
            a[1] = mpf('-1.5')
        alpha = mpf('2.5')
    else:
        a = list(a_init) + [mpf(0)] * (N - (len(a_init) - 1))
        alpha = mpf(alpha_init)
    nrm = None
    for it in range(maxit):
        F, J = residuals_and_jacobian(a, alpha, pts, d)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < tol:
            return a, alpha, it, nrm
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        for _h in range(40):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = residuals_and_jacobian(a_new, al_new, pts, d)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            raise RuntimeError('Newton stalled N=%d d=%s res=%s' % (N, d, nstr(nrm, 5)))
    raise RuntimeError('Newton not converged N=%d d=%s res=%s' % (N, d, nstr(nrm, 5)))


def fp_residual_dense(a, alpha, d, npts=2001):
    pts = colloc_points(npts)
    ad = alpha ** d
    worst = mpf(0)
    for t in pts:
        ta = t / ad
        U = G_eval(a, ta)
        W = fabs(U) ** d
        r = G_eval(a, t) + alpha * G_eval(a, W)
        worst = max(worst, abs(r))
    return worst


# --------------------------------------------------------------------------- #
#  scalar DT (ground truth, for checks and coordinate-mode diagnostic)
# --------------------------------------------------------------------------- #
def DT_scalar_t(a, alpha, d, c, t):
    """(DT Phi)(t), Phi_t(w) = sum_j c[j] w^j, c[0]=0 assumed."""
    G1 = G_eval(a, mpf(1))

    def Phi_t(w):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * w + cc
        return r

    Phi1 = sum(c)
    Gt = G_eval(a, t)
    Gpt = Gp_eval(a, t)
    K = Gt - d * t * Gpt
    ad = alpha ** d
    ta = t / ad
    U = G_eval(a, ta)
    sU = sgn(U)
    absU = fabs(U)
    W = absU ** d
    GpW = Gp_eval(a, W)
    term1 = -(Phi1 / G1) * K
    term2 = -alpha * Phi_t(W)
    term3 = -d * alpha * sU * (absU ** (d - 1) if absU != 0 else mpf(0)) * GpW * Phi_t(ta)
    return term1 + term2 + term3


def DT_unprojected_t(a, alpha, d, c, t):
    """Bare (unprojected) linearisation term; used only for coordinate-mode test."""
    def Phi_t(w):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * w + cc
        return r

    ad = alpha ** d
    ta = t / ad
    U = G_eval(a, ta)
    sU = sgn(U)
    absU = fabs(U)
    W = absU ** d
    GpW = Gp_eval(a, W)
    term2 = -alpha * Phi_t(W)
    term3 = -d * alpha * sU * (absU ** (d - 1) if absU != 0 else mpf(0)) * GpW * Phi_t(ta)
    return term2 + term3


def build_DT_matrix(a, alpha, d, N):
    """A = V^{-1} F ; V_ij = t_i^j, F_ij = DT(e_j)(t_i); i,j = 1..N.
    Eigenproblem of A gives the (projected) spectrum on normalised modes."""
    pts_t = colloc_points(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts_t):
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            e_j = [mpf(0)] * (j + 1)
            e_j[j] = mpf(1)
            F[i, j - 1] = DT_scalar_t(a, alpha, d, e_j, t)
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


def power_iter(M, n, iters=800, tol=mpf('1e-45')):
    v = matrix(n, 1)
    for i in range(n):
        v[i] = mpf(1) / (i + 3)
    lam = mpf(1)
    for it in range(iters):
        w = matvec(M, v)
        lam_new = max(abs(w[i]) for i in range(n))
        if lam_new == 0:
            return mpf(0), v, it
        # preserve sign relative to dominant component
        idx = max(range(n), key=lambda i: abs(w[i]))
        signed = w[idx] / abs(w[idx])
        for i in range(n):
            v[i] = w[i] / lam_new
        # signed eigenvalue estimate via Rayleigh-like ratio on dominant comp
        lam_signed = w[idx] / (v[idx] * lam_new) * lam_new if v[idx] != 0 else lam_new
        if abs(lam_new - lam) < tol * abs(lam_new):
            lam = lam_new
            break
        lam = lam_new
    # recompute true signed eigenvalue via Rayleigh quotient v^T M v / v^T v
    Mv = matvec(M, v)
    num = sum(v[i] * Mv[i] for i in range(n))
    den = sum(v[i] * v[i] for i in range(n))
    lam_signed = num / den
    return lam_signed, v, it


def deflate(M, lam1, v1, n):
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


def eigenpair_residual(a, alpha, d, hcoef, lam, npts=201):
    pts = colloc_points(npts)
    worst = mpf(0)
    for t in pts:
        lhs = DT_scalar_t(a, alpha, d, hcoef, t)
        rhs = lam * sum(hcoef[j] * t ** j for j in range(1, len(hcoef)))
        worst = max(worst, abs(lhs - rhs))
    return worst

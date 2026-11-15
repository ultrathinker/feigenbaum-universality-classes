#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""diag_cond.py -- condition-number ladder of the EIGEN collocation Jacobian
(J2 in the d=2 template: (DT B_j)(t_i) - lam B_j(t_i) rows + norm row) in the
Chebyshev B-basis for d=4,6 vs N, to confirm the conditioning grows slowly
(polynomial) instead of exponentially like the monomial basis did in wave2."""
import time

from mpmath import mp, mpf, matrix, lu_solve, cos, sin, pi, nstr, eig

mp.dps = 90

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


def cheb_vals(npts):
    return [(cos(pi * (i + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


def clenshaw(b, s):
    d1 = mpf(0)
    d2 = mpf(0)
    for k in range(len(b) - 1, 0, -1):
        d1, d2 = 2 * s * d1 - d2 + b[k], d1
    return s * d1 - d2 + b[0]


def dcheb_coeffs(b):
    N = len(b) - 1
    c = [mpf(0)] * (N + 2)
    for k in range(N, 0, -1):
        c[k - 1] = c[k + 1] + 2 * k * b[k]
    c[0] = c[0] / 2
    return c[:N]


class ChebFun:
    def __init__(self, b):
        self.b = list(b)
        self._db = None

    @classmethod
    def from_free(cls, x):
        N = len(x)
        b = [mpf(0)] + [mpf(v) for v in x]
        b[0] = mpf(1) - sum(((-1) ** k) * b[k] for k in range(1, N + 1))
        return cls(b)

    def G(self, t):
        return clenshaw(self.b, 2 * t - 1)

    def Gp(self, t):
        if self._db is None:
            self._db = dcheb_coeffs(self.b)
        return 2 * clenshaw(self._db, 2 * t - 1)


def residual_jacobian(x, alpha, d, pts):
    cf = ChebFun.from_free(x)
    b = cf.b
    db = dcheb_coeffs(b)
    N = len(x)
    npts = len(pts)
    ad = alpha ** d
    sgn = [(-1) ** k for k in range(N + 1)]
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    for i, t in enumerate(pts):
        ta = t / ad
        U = clenshaw(b, 2 * ta - 1)
        W = U ** d
        sW = 2 * W - 1
        GW = clenshaw(b, sW)
        GpW = 2 * clenshaw(db, sW)
        Gp_ta = 2 * clenshaw(db, 2 * (t / ad) - 1)
        F[i] = clenshaw(b, 2 * t - 1) + alpha * GW
        s_t = 2 * t - 1
        s_ta = 2 * ta - 1
        Tprev_t, Tcur_t = mpf(1), s_t
        Tprev_W, Tcur_W = mpf(1), sW
        Tprev_ta, Tcur_ta = mpf(1), s_ta
        for j in range(1, N + 1):
            if j == 1:
                dG_t = Tcur_t - sgn[1]
                dG_W = Tcur_W - sgn[1]
                dG_ta = Tcur_ta - sgn[1]
            else:
                Tnext_t = 2 * s_t * Tcur_t - Tprev_t
                Tnext_W = 2 * sW * Tcur_W - Tprev_W
                Tnext_ta = 2 * s_ta * Tcur_ta - Tprev_ta
                dG_t = Tnext_t - sgn[j]
                dG_W = Tnext_W - sgn[j]
                dG_ta = Tnext_ta - sgn[j]
                Tprev_t, Tcur_t = Tcur_t, Tnext_t
                Tprev_W, Tcur_W = Tcur_W, Tnext_W
                Tprev_ta, Tcur_ta = Tcur_ta, Tnext_ta
            J[i, j - 1] = dG_t + alpha * (dG_W + d * U ** (d - 1) * GpW * dG_ta)
        J[i, N] = GW - d * d * t * U ** (d - 1) * GpW * Gp_ta / ad
    return F, J


def solve_fixedpoint(d, N, x_init=None, alpha_init=None, tol=mpf('1e-70'), maxit=300):
    pts = cheb_vals(N + 1)
    if x_init is None:
        x = [mpf(0)] * N
        x[0] = mpf('-0.5')
        alpha = mpf('7.5') if d == 4 else mpf('9.5')
    else:
        x = [mpf(v) for v in x_init] + [mpf(0)] * (N - len(x_init))
        alpha = mpf(alpha_init)
    for it in range(maxit):
        F, J = residual_jacobian(x, alpha, d, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < tol:
            return x, alpha, it, nrm
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        improved = False
        for _halve in range(40):
            x_new = [x[j] + lam * dx[j] for j in range(N)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = residual_jacobian(x_new, al_new, d, pts)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                x, alpha, nrm = x_new, al_new, nrm2
                improved = True
                break
            lam /= 2
        if not improved:
            raise RuntimeError('stalled res=%s' % nstr(nrm, 5))
    raise RuntimeError('no converge res=%s' % nstr(nrm, 5))


def eig_jacobian_mp(a, alpha, c, lam, d, pts):
    """Eigen collocation Jacobian J2 (N+1 x N+1) at point (c,lam), general d.
    Returns J and also the residual rows M[i,j] = (DT B_j)(t_i) - lam B_j(t_i)."""
    N = len(c)
    sgn = [(-1) ** k for k in range(N + 1)]
    ad = alpha ** d
    b = list(a)  # full Cheb coeffs incl b0
    db = dcheb_coeffs(b)
    J = matrix(N + 1, N + 1)

    def Bfun(j, t):
        s = 2 * t - 1
        Tprev, Tcur = mpf(1), s
        for _k in range(1, j):
            Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
        return Tcur - sgn[j]

    def DT_Bj(j, t):
        phi1 = Bfun(j, mpf(1))
        U = clenshaw(b, 2 * (t / ad) - 1)
        W = U ** d
        GpW = 2 * clenshaw(db, 2 * W - 1)
        Gp_ta = 2 * clenshaw(db, 2 * (t / ad) - 1)
        return (alpha * phi1 * clenshaw(b, 2 * t - 1)
                + alpha ** (2 - d) * d * d * t * phi1 * U ** (d - 1) * GpW * Gp_ta
                - alpha * Bfun(j, W)
                - alpha * d * U ** (d - 1) * GpW * Bfun(j, t / ad))

    for i, t in enumerate(pts):
        for j in range(1, N + 1):
            J[i, j - 1] = DT_Bj(j, t) - lam * Bfun(j, t)
        ph = sum(c[j - 1] * Bfun(j, t) for j in range(1, N + 1))
        J[i, N] = -ph
    for j in range(1, N + 1):
        J[N, j - 1] = 2 if j % 2 == 1 else 0
    J[N, N] = 0
    return J


def cond(M, n):
    """cond = ||M||_1 * ||M^-1||_1 (mpmath, not rigorous, for diagnostics)."""
    n1 = max(sum(abs(M[i, j]) for i in range(n)) for j in range(n))
    cols = []
    for jc in range(n):
        e = matrix(n, 1)
        e[jc, 0] = mpf(1)
        x = lu_solve(M, e)
        cols.append([x[i, 0] for i in range(n)])
    n1inv = max(sum(abs(cols[jc][i]) for jc in range(n)) for i in range(n))
    return n1 * n1inv


def main():
    ladder = [16, 24, 32, 40, 48]
    for d in [4, 6]:
        log('===== d=%d =====' % d)
        x_cur, al_cur = None, None
        for N in ladder:
            x, al, it, res = solve_fixedpoint(d, N, x_cur, al_cur)
            x_cur, al_cur = x, al
            cf = ChebFun.from_free(x)
            bfull = cf.b
            # solve eigenpair via Newton-like: first get lambda2 from DT eig matrix
            epts = cheb_vals(N)
            sgn = [(-1) ** k for k in range(N + 1)]

            def Bfun(j, t):
                s = 2 * t - 1
                Tprev, Tcur = mpf(1), s
                for _k in range(1, j):
                    Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
                return Tcur - sgn[j]

            def DT_apply(phival, t):
                phi1 = phival(mpf(1))
                U = cf.G(t / al ** d)
                W = U ** d
                return (al * phi1 * cf.G(t)
                        + al ** (2 - d) * d * d * t * phi1 * U ** (d - 1) * cf.Gp(W) * cf.Gp(t / al ** d)
                        - al * phival(W)
                        - al * d * U ** (d - 1) * cf.Gp(W) * phival(t / al ** d))

            V = matrix(N, N)
            for i, t in enumerate(epts):
                s = 2 * t - 1
                Tprev, Tcur = mpf(1), s
                for j in range(1, N + 1):
                    if j == 1:
                        V[i, 0] = Tcur - sgn[1]
                    else:
                        Tnext = 2 * s * Tcur - Tprev
                        V[i, j - 1] = Tnext - sgn[j]
                        Tprev, Tcur = Tcur, Tnext
            cols = []
            for j in range(1, N + 1):
                vals = matrix(N, 1)
                for i, t in enumerate(epts):
                    def phival(tt, j=j):
                        return Bfun(j, tt)
                    vals[i] = DT_apply(phival, t)
                coef = lu_solve(V, vals)
                cols.append([coef[i] for i in range(N)])
            A = matrix(N, N)
            for jc in range(N):
                for i in range(N):
                    A[i, jc] = cols[jc][i]
            E, _ER = eig(A)
            anch = mpf('0.29183840898354510453') if d == 4 else mpf('0.37107146768193013247')
            best = max((E[i].real for i in range(N)), key=lambda z: abs(z))  # placeholder
            # pick eigenvalue nearest anchor
            lam = min((E[i].real for i in range(N)), key=lambda z: abs(z - anch))
            # eigenvector: build from eigen decomposition right vectors
            ER = eig(A, left=False, right=True)[1]
            idx = min(range(N), key=lambda i: abs(E[i].real - lam))
            v = [ER[i, idx].real for i in range(N)]
            # normalize Phi(1)=1: sum_{j odd} 2 c_j = 1
            s = sum(2 * v[j - 1] for j in range(1, N + 1, 2))
            v = [c / s for c in v]
            J2 = eig_jacobian_mp(bfull, al, v, lam, d, epts)
            k = cond(J2, N + 1)
            err_anchor = abs(lam - anch)
            log('  N=%2d: lambda2=%.18f  |diff anchor|=%.3e  cond(J2)=%.3e'
                % (N, float(lam), float(err_anchor), float(k)))


if __name__ == '__main__':
    main()

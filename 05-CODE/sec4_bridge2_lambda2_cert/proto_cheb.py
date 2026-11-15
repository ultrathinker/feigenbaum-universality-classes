#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
proto_cheb.py -- validation prototype for MOST-2: general-d (d=4,6) lambda2
certificate in the TRIGONOMETRIC Chebyshev basis T_k(2t-1)=cos(k*theta),
t=x^d in [0,1].  Float layer (mpmath, NOT rigorous): solves the collocation
fixed point (b_1..b_N, alpha) and the chained eigenpair (c_1..c_N, lambda)
and checks lambda2 against the float anchors.

Conventions (from wave2 + d=2 template):
  g(x) = G(x^d), t = x^d, G(t) = sum_{k=0}^N b_k T_k(2t-1), b0 locked by G(0)=1.
  g(x) = -alpha*g(g(x/alpha))  ->  G(t) + alpha*G(W) = 0,  W = U^d, U = G(t/alpha^d).
  alpha = the "delta" constant (7.28468621707334 d=4; 9.29624683277137 d=6).
  Linearisation DT (includes variation of alpha, alpha(g) = -1/g(1)):
  DT Phi = a Phi(1) G(t) + a^{2-d} d^2 t Phi(1) U^{d-1} G'(W) G'(t/a^d)
           - a Phi(W) - a d U^{d-1} G'(W) Phi(t/a^d)
  Eigenpair: (DT Phi)(t) = lambda Phi(t) at N nodes, Phi(1)=1, Phi(0)=0.
"""
import time

from mpmath import mp, mpf, mpmathify, matrix, lu_solve, cos, sin, pi, nstr, sqrt

mp.dps = 80

ANCHORS = {
    4: mpf('0.29183840898354510453'),
    6: mpf('0.37107146768193013247'),
}


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


T0 = time.time()


# --------------------------- Chebyshev utilities ----------------------------

def cheb_vals(npts):
    """Gauss-Chebyshev roots s_i = cos(pi*(i+1/2)/npts), t=(s+1)/2 in (0,1)."""
    return [(cos(pi * (i + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


def clenshaw(b, s):
    """sum_k b[k] T_k(s) via Clenshaw."""
    d1 = mpf(0)
    d2 = mpf(0)
    for k in range(len(b) - 1, 0, -1):
        d1, d2 = 2 * s * d1 - d2 + b[k], d1
    return s * d1 - d2 + b[0]


def dcheb_coeffs(b):
    """coeffs of d/ds (sum b_k T_k) = sum c_k T_k."""
    N = len(b) - 1
    c = [mpf(0)] * (N + 2)
    for k in range(N, 0, -1):
        c[k - 1] = c[k + 1] + 2 * k * b[k]
    c[0] = c[0] / 2
    return c[:N]


class ChebFun:
    """G(t) = sum_k b_k T_k(2t-1); b[0] locked by G(0)=1."""

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
        """dG/dt = 2 dG/ds."""
        if self._db is None:
            self._db = dcheb_coeffs(self.b)
        return 2 * clenshaw(self._db, 2 * t - 1)


# --------------------------- fixed point ----------------------------

def residual_jacobian(x, alpha, d, pts):
    """x=(b_1..b_N); F (npts), J (npts x N+1)."""
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
        F[i] = clenshaw(b, 2 * t - 1) + alpha * GW
        # dG/db_j at point p = T_j(2p-1) - (-1)^j  (b0 locked)
        # recurrence for T_j at s_t, s_W, s_ta
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
        J[i, N] = GW - d * d * t * U ** (d - 1) * GpW * Gp_ta_cf(t, alpha, ad, db) / ad
    return F, J


def Gp_ta_cf(t, alpha, ad, db):
    """G'(t/alpha^d) computed as 2*clenshaw(db, 2*(t/ad)-1)."""
    return 2 * clenshaw(db, 2 * (t / ad) - 1)


def solve_fixedpoint(d, N, x_init=None, alpha_init=None, tol=mpf('1e-60'), maxit=300):
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
            raise RuntimeError('Newton stalled (fixed point) res=%s' % nstr(nrm, 5))
    raise RuntimeError('fixed point no converge res=%s' % nstr(nrm, 5))


# --------------------------- linearisation DT ----------------------------

def DT_apply(cf, alpha, phival, t, d):
    """(DT Phi)(t) in t-space; phival callable t -> Phi(t), Phi(0)=0, Phi(1) used."""
    a = alpha
    ad = a ** d
    U = cf.G(t / ad)
    W = U ** d
    GpW = cf.Gp(W)
    Gp_ta = cf.Gp(t / ad)
    phi1 = phival(mpf(1))
    return (a * phi1 * cf.G(t)
            + a ** (2 - d) * d * d * t * phi1 * U ** (d - 1) * GpW * Gp_ta
            - a * phival(W)
            - a * d * U ** (d - 1) * GpW * phival(t / ad))


def build_eig_matrix(cf, alpha, d, N):
    """DT matrix in B-basis: B_j(t)=T_j(2t-1)-(-1)^j, j=1..N."""
    epts = cheb_vals(N)
    sgn = [(-1) ** k for k in range(N + 1)]

    def Bfun(j, t):
        s = 2 * t - 1
        Tprev, Tcur = mpf(1), s
        for _k in range(1, j):
            Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
        return Tcur - sgn[j]

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
            vals[i] = DT_apply(cf, alpha, phival, t, d)
        coef = lu_solve(V, vals)
        cols.append([coef[i] for i in range(N)])
    A = matrix(N, N)
    for jc in range(N):
        for i in range(N):
            A[i, jc] = cols[jc][i]
    return A, epts, Bfun


def main():
    for d in [4, 6]:
        log('=== d=%d ===' % d)
        ladder = [8, 12, 16, 20, 24]
        x_cur, al_cur = None, None
        sols = {}
        for N in ladder:
            x, al, it, res = solve_fixedpoint(d, N, x_cur, al_cur)
            sols[N] = (x, al)
            log('  FP N=%2d: alpha=%s (it %d res %.2e)' % (N, nstr(al, 30), it, float(res)))
            x_cur, al_cur = x, al
        alpha = al_cur
        # expected alpha (delta anchors)
        exp_alpha = mpf('7.2846862170733433643089305679955530694780466189573') if d == 4 else mpf('9.2962468327713700828344765663674575503066887803593')
        log('  alpha vs delta-anchor: |diff|=%.3e' % float(abs(alpha - exp_alpha)))
        cf = ChebFun.from_free(x_cur)
        # eigen matrix at final N
        N = ladder[-1]
        A, epts, Bfun = build_eig_matrix(cf, alpha, d, N)
        from mpmath import eig as mp_eig
        E, ER = mp_eig(A, right=True)
        order = sorted(range(N), key=lambda i: -abs(E[i]))
        log('  spectrum (top 6 by |lam|):')
        for rank in range(min(6, N)):
            log('    #%d %s' % (rank, nstr(E[order[rank]], 26)))
        # find lambda2 nearest to anchor
        anch = ANCHORS[d]
        best = min((E[i].real for i in range(N)), key=lambda z: abs(z - anch))
        log('  lambda2 near anchor: %s ; |diff anchor|=%.3e' % (nstr(best, 30), float(abs(best - anch))))
        # delta should be rank 0
        log('  delta(eig, rank0)=%s' % nstr(E[order[0]].real, 24))


if __name__ == '__main__':
    main()

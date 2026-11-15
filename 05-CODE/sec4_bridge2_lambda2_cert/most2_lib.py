#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
most2_lib.py -- MOST-2 certificate library (general-d, d=4,6), Chebyshev basis.

Basis: G(t) = sum_{k=0}^N b_k T_k(2t-1),  T_k(2t-1) = cos(k*theta),
theta = arccos(2t-1),  t = x^d in [0,1],  b_0 locked by G(0) = 1:
    b_0 = 1 - sum_{k>=1} (-1)^k b_k.
Fixed point (K1):  G(t) + alpha*G(W) = 0 at N+1 Gauss-Chebyshev nodes,
    W = U^d,  U = G(t/alpha^d).   Unknowns b_1..b_N, alpha.
Linearisation (K2): derivative of R(g) = -alpha(g) g(g(x/alpha(g))),
    alpha(g) = -1/g(1), evaluated on functions of x^d:
    (DT Phi)(t) = a Phi(1) G(t) + a^{2-d} d^2 t Phi(1) U^{d-1} G'(W) G'(t/a^d)
                - a Phi(W) - a d U^{d-1} G'(W) Phi(t/a^d)
    Eigenpair: (DT Phi)(t) = lambda Phi(t),  Phi(1) = 1,  Phi(0) = 0
    (basis B_j(t) = T_j(2t-1) - (-1)^j,  j = 1..N).

Float layer: mpmath (dps configurable) continuation + Newton -> midpoints.
Rigorous layer: python-flint arb ball Krawczyk, two chained stages (K1 box
flows into K2), plus independent interval Gauss-Seidel verifier.
"""
import json
import time

from flint import arb, arb_mat, ctx
from mpmath import (mp, mpf, mpmathify, matrix, lu_solve, cos, sin, pi as mppi,
                    nstr, eig as mp_eig)

# --------------------------------------------------------------------------- #
#  precision
# --------------------------------------------------------------------------- #
mp.dps = 170
ctx.dps = 160

SGN = {}


def sgn_int(N):
    if N not in SGN:
        SGN[N] = [(-1) ** k for k in range(N + 1)]
    return SGN[N]


# --------------------------------------------------------------------------- #
#  mpmath Chebyshev utilities
# --------------------------------------------------------------------------- #
def cheb_vals_mp(npts):
    """Gauss-Chebyshev roots of T_npts: s_i = cos(pi*(i+1/2)/npts), t=(s+1)/2."""
    return [(cos(mppi * (mpf(i) + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


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


# --------------------------------------------------------------------------- #
#  mpmath: fixed point
# --------------------------------------------------------------------------- #
def fp_residual_jacobian_mp(x, alpha, d, pts):
    cf = ChebFun.from_free(x)
    b = cf.b
    db = dcheb_coeffs(b)
    N = len(x)
    npts = len(pts)
    ad = alpha ** d
    sgn = sgn_int(N)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    for i, t in enumerate(pts):
        ta = t / ad
        U = clenshaw(b, 2 * ta - 1)
        W = U ** d
        sW = 2 * W - 1
        GW = clenshaw(b, sW)
        GpW = 2 * clenshaw(db, sW)
        Gp_ta = 2 * clenshaw(db, 2 * ta - 1)
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


def solve_fixedpoint_mp(d, N, x_init=None, alpha_init=None, tol=mpf('1e-80'),
                        maxit=300, log=None):
    pts = cheb_vals_mp(N + 1)
    if x_init is None:
        x = [mpf(0)] * N
        x[0] = mpf('-0.5')
        alpha = mpf('7.5') if d == 4 else mpf('9.5')
    else:
        x = [mpf(v) for v in x_init] + [mpf(0)] * (N - len(x_init))
        alpha = mpf(alpha_init)
    for it in range(maxit):
        F, J = fp_residual_jacobian_mp(x, alpha, d, pts)
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
                F2, _ = fp_residual_jacobian_mp(x_new, al_new, d, pts)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                x, alpha, nrm = x_new, al_new, nrm2
                improved = True
                break
            lam /= 2
        if not improved:
            raise RuntimeError('FP Newton stalled res=%s' % nstr(nrm, 5))
    raise RuntimeError('FP no converge res=%s' % nstr(nrm, 5))


# --------------------------------------------------------------------------- #
#  mpmath: eigenpair on the collocation system (c_1..c_N, lambda)
# --------------------------------------------------------------------------- #
def bfun_at_s_list(s, N, sgn):
    """T_j(s) - sgn[j] for j=0..N, by recurrence. Returns list index j."""
    out = [mpf(1) - sgn[0]]
    if N >= 1:
        out.append(s - sgn[1])
    Tprev, Tcur = mpf(1), s
    for j in range(2, N + 1):
        Tnext = 2 * s * Tcur - Tprev
        out.append(Tnext - sgn[j])
        Tprev, Tcur = Tcur, Tnext
    return out


def eigen_F_J_mp(a, alpha, c, lam, d, pts):
    """a: full Cheb coeffs (b0..bN). c: length N. Returns F (N+1), J (N+1 x N+1)."""
    N = len(c)
    sgn = sgn_int(N)
    ad = alpha ** d
    b = list(a)
    db = dcheb_coeffs(b)
    F = matrix(N + 1, 1)
    J = matrix(N + 1, N + 1)
    for i, t in enumerate(pts):
        ta = t / ad
        U = clenshaw(b, 2 * ta - 1)
        W = U ** d
        Gt = clenshaw(b, 2 * t - 1)
        GpW = 2 * clenshaw(db, 2 * W - 1)
        Gp_ta = 2 * clenshaw(db, 2 * ta - 1)
        Bl_t = bfun_at_s_list(2 * t - 1, N, sgn)
        Bl_W = bfun_at_s_list(2 * W - 1, N, sgn)
        Bl_ta = bfun_at_s_list(2 * ta - 1, N, sgn)
        Ud1 = U ** (d - 1)
        coef = alpha ** (2 - d) * d * d * t * Ud1 * GpW * Gp_ta
        fac = alpha * d * Ud1 * GpW
        Phit = mpf(0)
        for j in range(1, N + 1):
            Phit += c[j - 1] * Bl_t[j]
        F[i] = -lam * Phit
        for j in range(1, N + 1):
            phi1 = mpf(2) if j % 2 == 1 else mpf(0)  # B_j(1) = 1-(-1)^j
            DT_Bj = (alpha * phi1 * Gt + coef * phi1
                     - alpha * Bl_W[j] - fac * Bl_ta[j])
            F[i] += c[j - 1] * DT_Bj
            J[i, j - 1] = DT_Bj - lam * Bl_t[j]
        J[i, N] = -Phit
    # normalization Phi(1) = 1 : sum_{j odd} 2 c_j = 1
    F[N] = mpf(-1)
    for j in range(1, N + 1, 2):
        F[N] += 2 * c[j - 1]
        J[N, j - 1] = mpf(2)
    for j in range(2, N + 1, 2):
        J[N, j - 1] = mpf(0)
    J[N, N] = mpf(0)
    return F, J


def build_eig_seed(d, a, alpha, N, anchor, dps_eig=90):
    """Build DT matrix in B-basis, eig, return (c, lam) nearest anchor, Phi(1)=1."""
    from mpmath import mp as _mp
    _mp.dps = dps_eig
    sgn = sgn_int(N)
    ad = alpha ** d
    b = list(a)
    db = dcheb_coeffs(b)
    epts = cheb_vals_mp(N)

    def Bfun(j, t):
        s = 2 * t - 1
        Tprev, Tcur = mpf(1), s
        for _k in range(1, j):
            Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
        return Tcur - sgn[j]

    # interpolation matrix V[i,j] = B_j(t_i)
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
            vals[i] = DT_apply_mp(b, db, alpha, d, phival, t)
        coef = lu_solve(V, vals)
        cols.append([coef[i] for i in range(N)])
    A = matrix(N, N)
    for jc in range(N):
        for i in range(N):
            A[i, jc] = cols[jc][i]
    E, ER = mp_eig(A, right=True)
    idx = min(range(N), key=lambda i: abs(E[i].real - anchor))
    lam = E[idx].real
    v = [ER[i, idx].real for i in range(N)]
    ss = sum(2 * v[j - 1] for j in range(1, N + 1, 2))
    v = [x / ss for x in v]
    _mp.dps = mp.dps
    return v, lam


def DT_apply_mp(b, db, alpha, d, phival, t):
    ad = alpha ** d
    U = clenshaw(b, 2 * (t / ad) - 1)
    W = U ** d
    phi1 = phival(mpf(1))
    return (alpha * phi1 * clenshaw(b, 2 * t - 1)
            + alpha ** (2 - d) * d * d * t * phi1 * U ** (d - 1) * 2 * clenshaw(db, 2 * W - 1) * 2 * clenshaw(db, 2 * (t / ad) - 1)
            - alpha * phival(W)
            - alpha * d * U ** (d - 1) * 2 * clenshaw(db, 2 * W - 1) * phival(t / ad))


def solve_eigen_mp(d, a, alpha, N, anchor, c_init=None, lam_init=None,
                   tol=mpf('1e-70'), maxit=200, log=None):
    pts = cheb_vals_mp(N)
    if c_init is None or lam_init is None:
        c, lam = build_eig_seed(d, a, alpha, N, anchor)
    else:
        c = list(c_init) + [mpf(0)] * (N - len(c_init))
        lam = mpf(lam_init)
    for it in range(maxit):
        F, J = eigen_F_J_mp(a, alpha, c, lam, d, pts)
        nrm = max(abs(F[i]) for i in range(N + 1))
        if nrm < tol:
            return c, lam, it, nrm
        rhs = matrix(N + 1, 1)
        for i in range(N + 1):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        step = mpf(1)
        improved = False
        for _halve in range(40):
            c_new = [c[k] + step * dx[k] for k in range(N)]
            lam_new = lam + step * dx[N]
            try:
                F2, _ = eigen_F_J_mp(a, alpha, c_new, lam_new, d, pts)
                nrm2 = max(abs(F2[i]) for i in range(N + 1))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                c, lam, nrm = c_new, lam_new, nrm2
                improved = True
                break
            step /= 2
        if not improved:
            raise RuntimeError('EIG Newton stalled res=%s' % nstr(nrm, 5))
    raise RuntimeError('EIG no converge res=%s' % nstr(nrm, 5))


# --------------------------------------------------------------------------- #
#  arb helpers
# --------------------------------------------------------------------------- #
def ball(rstr):
    return arb('-' + rstr).union(arb(rstr))


def arb_eye(n):
    M = arb_mat(n, n)
    for i in range(n):
        M[i, i] = arb(1)
    return M


def gauss_nodes_arb(npts):
    ts = []
    for i in range(npts):
        t = (cos(mppi * (mpf(i) + mpf('0.5')) / npts) + 1) / 2
        ts.append(arb(nstr(t, ctx.dps - 30)))
    return ts


def mp_to_arb(x, dps=None):
    return arb(nstr(x, dps or ctx.dps - 30))


def mp_vec_to_arb(v, dps=None):
    return [arb(nstr(x, dps or ctx.dps - 30)) for x in v]


def rad_ub(x):
    return float(x.rad().upper())


class TrigNode:
    """Chebyshev data at ball argument t in (0,1): T_k(2t-1) = cos(k*th)."""

    def __init__(self, t, b, N):
        s = 2 * t - 1
        self.t = t
        th = s.acos()
        self.th = th
        self.T = [arb(1)]
        for k in range(1, N + 1):
            self.T.append((arb(k) * th).cos())
        st = th.sin()
        self.sinth = st
        self.G = b[0]
        dsum = arb(0)
        for k in range(1, N + 1):
            self.G = self.G + b[k] * self.T[k]
            dsum = dsum + arb(k) * b[k] * (arb(k) * th).sin()
        self.Gp = 2 * dsum / st


def b0_of(bs, N):
    s = arb(1)
    for k in range(1, N + 1):
        s = s - arb(sgn_int(N)[k]) * bs[k - 1]
    return s


# --------------------------------------------------------------------------- #
#  arb: K1 fixed point F and J
# --------------------------------------------------------------------------- #
def fp_F_J_arb(d, N, bX, aX, pts, want_J=True):
    """bX: N+1 arb balls (incl b0), aX: arb alpha ball. Returns (F, rows or None)."""
    ad = aX ** d
    sgn = sgn_int(N)
    F = []
    Jrows = []
    for t in pts:
        ta = t / ad
        nta = TrigNode(ta, bX, N)
        U = nta.G
        W = U ** d
        nW = TrigNode(W, bX, N)
        nt = TrigNode(t, bX, N)
        F.append(nt.G + aX * nW.G)
        if want_J:
            row = []
            Ud1 = U ** (d - 1)
            GpW = nW.Gp
            Gp_ta = nta.Gp
            for j in range(1, N + 1):
                dG_t = nt.T[j] - arb(sgn[j])
                dG_W = nW.T[j] - arb(sgn[j])
                dG_ta = nta.T[j] - arb(sgn[j])
                row.append(dG_t + aX * (dG_W + arb(d) * Ud1 * GpW * dG_ta))
            row.append(nW.G - arb(d * d) * t * Ud1 * GpW * Gp_ta / ad)
            Jrows.append(row)
    return F, Jrows


# --------------------------------------------------------------------------- #
#  arb: K2 eigen F and J (chained on K1 box)
# --------------------------------------------------------------------------- #
def DTcol_ball(d, N, sgn, aX, node_t, node_ta, node_W, pts, j, ti):
    nt, nta, nW = node_t[ti], node_ta[ti], node_W[ti]
    t = pts[ti]
    U = nta.G
    Ud1 = U ** (d - 1)
    GpW = nW.Gp
    Gp_ta = nta.Gp
    phi1 = arb(1 - sgn[j])
    return (aX * phi1 * nt.G
            + aX ** (2 - d) * arb(d * d) * t * phi1 * Ud1 * GpW * Gp_ta
            - aX * (nW.T[j] - arb(sgn[j]))
            - aX * arb(d) * Ud1 * GpW * (nta.T[j] - arb(sgn[j])))


def eigen_F_J_arb(d, N, bX, aX, c, lam, pts):
    """c: length-N arb list, lam: arb. Returns (F list N+1, Jrows N+1 x N+1)."""
    sgn = sgn_int(N)
    ad = aX ** d
    node_t = [TrigNode(t, bX, N) for t in pts]
    node_ta = [TrigNode(t / ad, bX, N) for t in pts]
    node_W = []
    for t in pts:
        nta = TrigNode(t / ad, bX, N)
        node_W.append(TrigNode(nta.G ** d, bX, N))
    F = []
    Jrows = []
    for ti in range(N):
        nt = node_t[ti]
        dt = arb(0)
        ph = arb(0)
        for j in range(1, N + 1):
            dt = dt + c[j - 1] * DTcol_ball(d, N, sgn, aX, node_t, node_ta, node_W, pts, j, ti)
            ph = ph + c[j - 1] * (nt.T[j] - arb(sgn[j]))
        F.append(dt - lam * ph)
        row = []
        for j in range(1, N + 1):
            row.append(DTcol_ball(d, N, sgn, aX, node_t, node_ta, node_W, pts, j, ti)
                       - lam * (nt.T[j] - arb(sgn[j])))
        row.append(-ph)
        Jrows.append(row)
    Fn = arb(-1)
    for j in range(1, N + 1, 2):
        Fn = Fn + 2 * c[j - 1]
    F.append(Fn)
    normrow = [arb(2) if j % 2 == 1 else arb(0) for j in range(1, N + 1)] + [arb(0)]
    Jrows.append(normrow)
    return F, Jrows


# --------------------------------------------------------------------------- #
#  generic Krawczyk + interval Gauss-Seidel (same as bridge_lib)
# --------------------------------------------------------------------------- #
def mp_mat_inverse(M, n):
    cols = []
    for jcol in range(n):
        e = matrix(n, 1)
        e[jcol, 0] = mpf(1)
        x = lu_solve(M, e)
        cols.append([x[i, 0] for i in range(n)])
    Y = matrix(n, n)
    for jcol in range(n):
        for i in range(n):
            Y[i, jcol] = cols[jcol][i]
    return Y


def krawczyk_generic(xhat_mp_list, Jmid_mp, F_at_point, F_J_on_box, dps, box_r_str):
    n = len(xhat_mp_list)
    xhat_arb = mp_vec_to_arb(xhat_mp_list, dps)
    X_vec = [x + ball(box_r_str) for x in xhat_arb]
    Fh = F_at_point(xhat_arb)
    FX, JXrows = F_J_on_box(X_vec)
    Y_mp = mp_mat_inverse(Jmid_mp, n)
    Y = arb_mat([[mp_to_arb(Y_mp[i, j], dps) for j in range(n)] for i in range(n)])
    JX = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            JX[i, j] = JXrows[i][j]
    Fhm = arb_mat(n, 1)
    for i in range(n):
        Fhm[i, 0] = Fh[i]
    Xm = arb_mat(n, 1)
    for i in range(n):
        Xm[i, 0] = X_vec[i]
    xhm = arb_mat(n, 1)
    for i in range(n):
        xhm[i, 0] = xhat_arb[i]
    Kimg = xhm - Y * Fhm + (arb_eye(n) - Y * JX) * (Xm - xhm)
    contained = True
    fail_idx = []
    for i in range(n):
        if not X_vec[i].contains_interior(Kimg[i, 0]):
            contained = False
            fail_idx.append(i)
    halfwidths = [rad_ub(Kimg[i, 0]) for i in range(n)]
    return {'contained': contained, 'fail_idx': fail_idx, 'K': Kimg, 'X': X_vec,
            'xhat_arb': xhat_arb, 'halfwidths': halfwidths}


def gauss_seidel_verify(xhat_mp_list, Jmid_mp, F_at_point, F_J_on_box, dps, box_r_str,
                        sweeps=3):
    n = len(xhat_mp_list)
    xhat_arb = mp_vec_to_arb(xhat_mp_list, dps)
    X0 = [x + ball(box_r_str) for x in xhat_arb]
    Fh = F_at_point(xhat_arb)
    Y_mp = mp_mat_inverse(Jmid_mp, n)
    Y = [[mp_to_arb(Y_mp[i, j], dps) for j in range(n)] for i in range(n)]
    YFh = []
    for i in range(n):
        s = arb(0)
        for j in range(n):
            s = s + Y[i][j] * Fh[j]
        YFh.append(s)
    X = list(X0)
    contained_each = [False] * n
    for _sweep in range(sweeps):
        _FX, JXrows = F_J_on_box(X)
        B = [[arb(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                s = arb(0)
                for k in range(n):
                    s = s + Y[i][k] * JXrows[k][j]
                B[i][j] = s
        for i in range(n):
            rhs = -YFh[i]
            for j in range(n):
                if j == i:
                    continue
                rhs = rhs - B[i][j] * (X[j] - xhat_arb[j])
            bii = B[i][i]
            if bii.contains(0):
                continue
            zi = xhat_arb[i] + rhs / bii
            try:
                newXi = X[i].intersection(zi)
            except ValueError:
                contained_each[i] = False
                X[i] = X0[i]
                continue
            contained_each[i] = X0[i].contains_interior(zi)
            X[i] = newXi
    all_contained = all(contained_each)
    halfwidths = [rad_ub(X[i]) for i in range(n)]
    return {'contained': all_contained, 'contained_each': contained_each, 'X': X,
            'halfwidths': halfwidths}


# --------------------------------------------------------------------------- #
#  certification orchestrators
# --------------------------------------------------------------------------- #
def cond_of(Jmid_mp, n):
    """||M||_1 * ||M^-1||_1 (mpmath diagnostic, NOT rigorous)."""
    n1 = max(sum(abs(Jmid_mp[i, j]) for i in range(n)) for j in range(n))
    cols = []
    for jc in range(n):
        e = matrix(n, 1)
        e[jc, 0] = mpf(1)
        x = lu_solve(Jmid_mp, e)
        cols.append([x[i, 0] for i in range(n)])
    n1inv = max(sum(abs(cols[jc][i]) for jc in range(n)) for i in range(n))
    return n1 * n1inv


DEFAULT_RADII = ['1e-50', '1e-45', '1e-40', '1e-38', '1e-36', '1e-34', '1e-32',
                 '1e-30', '1e-28', '1e-26', '1e-24', '1e-22', '1e-20', '1e-18',
                 '1e-16', '1e-14', '1e-12', '1e-10', '1e-8', '1e-6']


def certify_fixedpoint(d, N, a_mid, alpha_mid, dps, radii=None, log=None):
    pts_mp = cheb_vals_mp(N + 1)
    pts_arb = gauss_nodes_arb(N + 1)
    xhat_mp = list(a_mid[1:]) + [alpha_mid]
    Fmp, Jmp = fp_residual_jacobian_mp(list(a_mid[1:]), alpha_mid, d, pts_mp)
    sgn = sgn_int(N)

    def F_at_point(v):
        b = [b0_of(v[:N], N)] + v[:N]
        F, _ = fp_F_J_arb(d, N, b, v[N], pts_arb, want_J=False)
        return F

    def F_J_on_box(v):
        b = [b0_of(v[:N], N)] + v[:N]
        return fp_F_J_arb(d, N, b, v[N], pts_arb, want_J=True)

    for r in (radii or DEFAULT_RADII):
        res = krawczyk_generic(xhat_mp, Jmp, F_at_point, F_J_on_box, dps, r)
        if res['contained']:
            res['box_r_used'] = r
            if log:
                log('  K1(d=%d,N=%d) OK box_r=%s alpha_hw=%.3e'
                    % (d, N, r, res['halfwidths'][N]))
            return res
    raise RuntimeError('K1 fixed-point Krawczyk failed for all tried box radii')


def certify_eigen(d, N, a_mid, alpha_mid, c_mid, lam_mid, dps, aX_box, alphaX_box,
                  radii=None, log=None):
    pts_mp = cheb_vals_mp(N)
    pts_arb = gauss_nodes_arb(N)
    xhat_mp = list(c_mid) + [lam_mid]
    Fmp, Jmp = eigen_F_J_mp([b0_of(a_mid[1:], N)] + list(a_mid[1:]), alpha_mid,
                            c_mid, lam_mid, d, pts_mp)

    def F_at_point(v):
        c = v[:N]
        lam = v[N]
        F, _ = eigen_F_J_arb(d, N, aX_box, alphaX_box, c, lam, pts_arb)
        return F

    def F_J_on_box(v):
        c = v[:N]
        lam = v[N]
        return eigen_F_J_arb(d, N, aX_box, alphaX_box, c, lam, pts_arb)

    for r in (radii or DEFAULT_RADII):
        res = krawczyk_generic(xhat_mp, Jmp, F_at_point, F_J_on_box, dps, r)
        if res['contained']:
            res['box_r_used'] = r
            if log:
                log('  K2(d=%d,N=%d) OK box_r=%s lambda_hw=%.3e'
                    % (d, N, r, res['halfwidths'][N]))
            return res
    raise RuntimeError('K2 eigenpair Krawczyk failed for all tried box radii')


def run_certification(d, ladder, anchor, dps=None, log=None, save_path=None):
    """Full ladder: fixed point continuation + K1 cert + K2 cert, per N.
    Returns dict N -> result (with all strings)."""
    if dps is None:
        dps = ctx.dps - 30
    out = {}
    x_cur, al_cur = None, None
    prev_lam = None
    for N in ladder:
        x_cur, al_cur, it, res = solve_fixedpoint_mp(d, N, x_cur, al_cur,
                                                     log=log)
        a_full = [mpf(1) - sum(((-1) ** k) * x_cur[k - 1] for k in range(1, N + 1))] + list(x_cur)
        if log:
            log('  FP N=%d converged it=%d res=%.2e alpha=%s'
                % (N, it, float(res), nstr(al_cur, 30)))
        # eigen midpoint
        if prev_lam is None:
            c, lam, ite, rese = solve_eigen_mp(d, a_full, al_cur, N, anchor, log=log)
        else:
            c, lam, ite, rese = solve_eigen_mp(d, a_full, al_cur, N, anchor,
                                               c_init=None, lam_init=prev_lam, log=log)
        if log:
            log('  EIG N=%d converged it=%d res=%.2e lambda2=%s'
                % (N, ite, float(rese), nstr(lam, 40)))
        prev_lam = lam
        # K1
        res1 = certify_fixedpoint(d, N, a_full, al_cur, dps, log=log)
        aX_box = [b0_of(res1['X'][:N], N)] + res1['X'][:N]
        alphaX_box = res1['X'][N]
        # K2
        res2 = certify_eigen(d, N, a_full, al_cur, c, lam, dps, aX_box, alphaX_box,
                             log=log)
        condJ2 = float(cond_of(eigen_F_J_mp(a_full, al_cur, c, lam, d, cheb_vals_mp(N))[1], N + 1))
        rec = {
            'N': N,
            'alpha_mid': nstr(al_cur, 60),
            'alpha_rad': '%.3e' % res1['halfwidths'][N],
            'k1_box_r': res1['box_r_used'],
            'k1_contained': res1['contained'],
            'lambda2_mid': nstr(lam, 80),
            'lambda2_rad': '%.3e' % res2['halfwidths'][N],
            'k2_box_r': res2['box_r_used'],
            'k2_contained': res2['contained'],
            'cond_J2': condJ2,
            'c_mid': [nstr(v, 60) for v in c],
            'b_mid': [nstr(v, 60) for v in a_full],
            'anchor': nstr(mpmathify(anchor), 80),
            'anchor_diff': '%.3e' % abs(lam - mpmathify(anchor)),
        }
        out[N] = rec
        if save_path:
            with open(save_path, 'w', encoding='utf-8') as f:
                json.dump(out, f, indent=2, ensure_ascii=False)
        if log:
            log('  == N=%d: lambda2=%s +- %s  (anchor diff %s)  cond_J2=%.2e'
                % (N, nstr(lam, 50), rec['lambda2_rad'], rec['anchor_diff'], condJ2))
    return out

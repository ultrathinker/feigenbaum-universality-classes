#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bridge_lib.py -- shared generic-degree (d=4,6,...) functions for the "bridge"
lambda2 certificate: collocation fixed point g_N + eigenpair (Phi,lambda) +
normalization, solved first at high precision in mpmath (real, to get a good
midpoint) then certified rigorously in python-flint arb ball arithmetic via
Krawczyk (two chained stages, same pattern as the d=2 template
kilo_verify/v3_krawczyk.py) plus an independent interval Gauss-Seidel check.

Equation (BRIEF/template convention):
    g(x) = -alpha * g(g(x/alpha)),  g(0) = 1,  g even function of x^d,
    g(x) = G(t), t = x^d, G(t) = sum_{j=0}^N a_j t^j, a_0 = 1.

Linearisation (DT) on the "function of x^d" subspace, in t-variable
(Phi(t) = sum_{j=1}^N c_j t^j, c_0 = 0 -- perturbations must vanish at 0):
    K(t)  = G(t) - d*t*G'(t)
    U(t)  = G(t/alpha^d)
    W(t)  = U(t)^d
    (DT Phi)(t) = -(Phi(1)/G(1)) K(t) - alpha*Phi(W(t))
                  - d*alpha*U(t)^(d-1)*G'(W(t))*Phi(t/alpha^d)

Derived by hand from the scalar x-space formula in p1_fixedpoint.py (verified
against its d=4 Jacobian formulas, which are the d=4 special case of the
general-d formulas used here). See progress.log for the derivation notes.
"""
import time

from mpmath import mp, mpf, matrix as mpmatrix, lu_solve, pi as mppi, cos as mpcos, nstr
from flint import arb, arb_mat, ctx


# --------------------------------------------------------------------------- #
#  collocation nodes (Chebyshev-type, clustered at both ends of (0,1))
# --------------------------------------------------------------------------- #
def colloc_points_mp(npts):
    return [(1 - mpcos(mppi * (mpf(i) + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def colloc_points_arb(npts, dps):
    # high-precision mpmath cos, converted to arb via decimal string (thin balls
    # -- rigorous because arb() parses the decimal string with correct rounding
    # to a ball that CONTAINS the true real value represented by the string;
    # since the string itself is only an approximation of the true node (an
    # irrational number in general), we inflate by 1 ulp at the string's
    # decimal precision to guarantee containment of the mathematically exact node.
    # NOTE: exact node values are irrelevant to soundness -- collocation just
    # needs *some* fixed real points; using an arb ball that CONTAINS a nearby
    # real number is fine, we do not need it to contain a specific irrational.
    pts = colloc_points_mp(npts)
    return [arb(nstr(t, dps)) for t in pts]


# --------------------------------------------------------------------------- #
#  mpmath (real, high precision) Horner helpers
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
#  fixed point: mpmath Newton (unknowns a_1..a_N, alpha)
# --------------------------------------------------------------------------- #
def fp_residual_jacobian_mp(a, alpha, pts, d):
    N = len(a) - 1
    npts = len(pts)
    F = mpmatrix(npts, 1)
    J = mpmatrix(npts, N + 1)
    ad = alpha ** d
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval_mp(a, ta)
        Up = Gp_eval_mp(a, ta)
        W = U ** d
        G_W = G_eval_mp(a, W)
        Gp_W = Gp_eval_mp(a, W)
        F[i] = G_eval_mp(a, t) + alpha * G_W
        base = alpha * Gp_W * d * U ** (d - 1)
        pow_t = mpf(1)
        pow_W = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            J[i, j - 1] = pow_t + alpha * pow_W + base * ta ** j
        J[i, N] = G_W - d * d * U ** (d - 1) * t * Gp_W * Up / ad
    return F, J


def solve_fixedpoint_mp(d, N, a_init=None, alpha_init=None, tol=mpf('1e-70'), maxit=200, log=None):
    pts = colloc_points_mp(N + 1)
    if a_init is None:
        a = [mpf(1)] + [mpf(0)] * N
        alpha = mpf('1.5')
    else:
        a = list(a_init) + [mpf(0)] * (N - (len(a_init) - 1))
        alpha = mpf(alpha_init)
    nrm = None
    for it in range(maxit):
        F, J = fp_residual_jacobian_mp(a, alpha, pts, d)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < tol:
            return a, alpha, it, nrm
        rhs = mpmatrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        for _halve in range(40):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = fp_residual_jacobian_mp(a_new, al_new, pts, d)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            raise RuntimeError('Newton stalled (fixed point), residual %s' % nstr(nrm, 5))
        if log:
            log('    FP Newton it=%d resid=%.3e' % (it, float(nrm)))
    raise RuntimeError('Newton did not converge (fixed point), residual %s' % nstr(nrm, 5))


def fp_residual_dense_mp(a, alpha, d, npts=2001):
    pts = colloc_points_mp(npts)
    ad = alpha ** d
    worst = mpf(0)
    for t in pts:
        U = G_eval_mp(a, t / ad)
        r = G_eval_mp(a, t) + alpha * G_eval_mp(a, U ** d)
        worst = max(worst, abs(r))
    return worst


# --------------------------------------------------------------------------- #
#  eigenpair: mpmath Newton (unknowns c_1..c_N, lambda) chained on (a, alpha)
# --------------------------------------------------------------------------- #
def eig_residual_jacobian_mp(a, alpha, c, lam, pts, d):
    """pts: N collocation t-nodes. c: length-N list, c[j-1] = coeff of t^j."""
    N = len(c)
    G1 = G_eval_mp(a, mpf(1))
    ad = alpha ** d
    F = mpmatrix(N + 1, 1)
    J = mpmatrix(N + 1, N + 1)
    Phi1 = sum(c)
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval_mp(a, ta)
        Up = Gp_eval_mp(a, ta)
        W = U ** d
        Gp_W = Gp_eval_mp(a, W)
        Kt = G_eval_mp(a, t) - d * t * Gp_eval_mp(a, t)
        Ud1 = U ** (d - 1)
        base = alpha * d * Ud1 * Gp_W
        PhiW = mpf(0)
        Phita = mpf(0)
        Phit = mpf(0)
        pow_W = mpf(1)
        pow_ta = mpf(1)
        pow_t = mpf(1)
        for k in range(1, N + 1):
            pow_W *= W
            pow_ta *= ta
            pow_t *= t
            PhiW += c[k - 1] * pow_W
            Phita += c[k - 1] * pow_ta
            Phit += c[k - 1] * pow_t
        DT_val = -(Phi1 / G1) * Kt - alpha * PhiW - base * Phita
        F[i] = DT_val - lam * Phit
        pow_W = mpf(1)
        pow_ta = mpf(1)
        pow_t = mpf(1)
        for k in range(1, N + 1):
            pow_W *= W
            pow_ta *= ta
            pow_t *= t
            dDT_dk = -(mpf(1) / G1) * Kt - alpha * pow_W - base * pow_ta
            J[i, k - 1] = dDT_dk - lam * pow_t
        J[i, N] = -Phit
    F[N] = Phi1 - 1
    for k in range(1, N + 1):
        J[N, k - 1] = mpf(1)
    J[N, N] = mpf(0)
    return F, J


def solve_eigen_mp(d, a, alpha, N, c_init=None, lam_init=None, tol=mpf('1e-70'), maxit=200, log=None):
    pts = colloc_points_mp(N)
    if c_init is None:
        c = [mpf(1) / N] * N
    else:
        c = list(c_init) + [mpf(0)] * (N - len(c_init))
    lam = mpf(lam_init)
    nrm = None
    for it in range(maxit):
        F, J = eig_residual_jacobian_mp(a, alpha, c, lam, pts, d)
        nrm = max(abs(F[i]) for i in range(N + 1))
        if nrm < tol:
            return c, lam, it, nrm
        rhs = mpmatrix(N + 1, 1)
        for i in range(N + 1):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lstep = mpf(1)
        for _halve in range(40):
            c_new = [c[k] + lstep * dx[k] for k in range(N)]
            lam_new = lam + lstep * dx[N]
            try:
                F2, _ = eig_residual_jacobian_mp(a, alpha, c_new, lam_new, pts, d)
                nrm2 = max(abs(F2[i]) for i in range(N + 1))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                c, lam, nrm = c_new, lam_new, nrm2
                break
            lstep /= 2
        else:
            raise RuntimeError('Newton stalled (eigenpair), residual %s' % nstr(nrm, 5))
        if log:
            log('    EIG Newton it=%d resid=%.3e' % (it, float(nrm)))
    raise RuntimeError('Newton did not converge (eigenpair), residual %s' % nstr(nrm, 5))


def eig_residual_dense_mp(a, alpha, c, lam, d, npts=2001):
    """Scalar (function-space, no series-Jacobian) residual check in x-space,
    independent sanity check of the eigen collocation solution."""
    pts_t = colloc_points_mp(npts)
    ad = alpha ** d
    G1 = G_eval_mp(a, mpf(1))
    Phi1 = sum(c)
    worst = mpf(0)
    for t in pts_t:
        ta = t / ad
        U = G_eval_mp(a, ta)
        W = U ** d
        Gp_W = Gp_eval_mp(a, W)
        Kt = G_eval_mp(a, t) - d * t * Gp_eval_mp(a, t)
        PhiW = G_eval_mp([mpf(0)] + c, W)
        Phita = G_eval_mp([mpf(0)] + c, ta)
        Phit = G_eval_mp([mpf(0)] + c, t)
        DT_val = -(Phi1 / G1) * Kt - alpha * PhiW - d * alpha * U ** (d - 1) * Gp_W * Phita
        worst = max(worst, abs(DT_val - lam * Phit))
    return worst


# --------------------------------------------------------------------------- #
#  arb (ball) helpers
# --------------------------------------------------------------------------- #
def ball(rstr):
    return arb('-' + rstr).union(arb(rstr))


def arb_eye(n):
    M = arb_mat(n, n)
    for i in range(n):
        M[i, i] = arb(1)
    return M


def rad_ub(x):
    """rigorous half-width upper bound of an arb ball, as a python float
    (rounded UP so it remains a valid upper bound)."""
    return float(x.rad().upper())


def mp_to_arb(x, dps):
    return arb(nstr(x, dps))


def mp_vec_to_arb(v, dps):
    return [arb(nstr(x, dps)) for x in v]


def mp_mat_inverse(M, n):
    """mpmath inverse via n LU solves (float/mpf, not rigorous -- used only to
    build the Krawczyk preconditioner Y ~ J(mid)^{-1}, whose only job is to make
    the Krawczyk map contractive; any invertible approximation works, exactness
    is not required for soundness)."""
    cols = []
    for jcol in range(n):
        e = mpmatrix(n, 1)
        e[jcol, 0] = mpf(1)
        x = lu_solve(M, e)
        cols.append([x[i, 0] for i in range(n)])
    Y = mpmatrix(n, n)
    for jcol in range(n):
        for i in range(n):
            Y[i, jcol] = cols[jcol][i]
    return Y


def G_eval_arb(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval_arb(a, t):
    s = arb(0)
    for j in range(1, len(a)):
        s = s + arb(j) * a[j] * t ** (j - 1)
    return s


def fp_F_arb(a, alpha, pts, d):
    """residual vector only (arb), for boxes or points."""
    ad = alpha ** d
    F = []
    for t in pts:
        ta = t / ad
        U = G_eval_arb(a, ta)
        W = U ** d
        G_W = G_eval_arb(a, W)
        F.append(G_eval_arb(a, t) + alpha * G_W)
    return F


def fp_F_J_arb(a, alpha, pts, d):
    N = len(a) - 1
    ad = alpha ** d
    F = []
    Jrows = []
    for t in pts:
        ta = t / ad
        U = G_eval_arb(a, ta)
        Up = Gp_eval_arb(a, ta)
        W = U ** d
        G_W = G_eval_arb(a, W)
        Gp_W = Gp_eval_arb(a, W)
        F.append(G_eval_arb(a, t) + alpha * G_W)
        base = alpha * Gp_W * arb(d) * U ** (d - 1)
        row = []
        pow_t = arb(1)
        pow_W = arb(1)
        for j in range(1, N + 1):
            pow_t = pow_t * t
            pow_W = pow_W * W
            row.append(pow_t + alpha * pow_W + base * ta ** j)
        row.append(G_W - arb(d * d) * U ** (d - 1) * t * Gp_W * Up / ad)
        Jrows.append(row)
    return F, Jrows


def eig_F_J_arb(a, alpha, c, lam, pts, d, G1):
    """a, alpha: arb (possibly ball-valued, chained from stage-1 box).
    c: length-N arb list, lam: arb. pts: N arb t-nodes."""
    N = len(c)
    ad = alpha ** d
    Phi1 = arb(0)
    for x in c:
        Phi1 = Phi1 + x
    F = []
    Jrows = []
    for t in pts:
        ta = t / ad
        U = G_eval_arb(a, ta)
        W = U ** d
        Gp_W = Gp_eval_arb(a, W)
        Kt = G_eval_arb(a, t) - arb(d) * t * Gp_eval_arb(a, t)
        Ud1 = U ** (d - 1)
        base = alpha * arb(d) * Ud1 * Gp_W
        PhiW = arb(0)
        Phita = arb(0)
        Phit = arb(0)
        pow_W = arb(1)
        pow_ta = arb(1)
        pow_t = arb(1)
        for k in range(1, N + 1):
            pow_W = pow_W * W
            pow_ta = pow_ta * ta
            pow_t = pow_t * t
            PhiW = PhiW + c[k - 1] * pow_W
            Phita = Phita + c[k - 1] * pow_ta
            Phit = Phit + c[k - 1] * pow_t
        DT_val = -(Phi1 / G1) * Kt - alpha * PhiW - base * Phita
        F.append(DT_val - lam * Phit)
        row = []
        pow_W = arb(1)
        pow_ta = arb(1)
        pow_t = arb(1)
        for k in range(1, N + 1):
            pow_W = pow_W * W
            pow_ta = pow_ta * ta
            pow_t = pow_t * t
            dDT_dk = -(arb(1) / G1) * Kt - alpha * pow_W - base * pow_ta
            row.append(dDT_dk - lam * pow_t)
        row.append(-Phit)
        Jrows.append(row)
    Frow_norm = arb(0)
    for x in c:
        Frow_norm = Frow_norm + x
    F.append(Frow_norm - arb(1))
    normrow = [arb(1)] * N + [arb(0)]
    Jrows.append(normrow)
    return F, Jrows


def timer(t0):
    return time.time() - t0


# --------------------------------------------------------------------------- #
#  continuation ladders (mpmath, real high precision) up to a target N
# --------------------------------------------------------------------------- #
DEFAULT_LADDER = [8, 14, 20, 28, 36, 44]


def continue_fixedpoint_mp(d, target_Ns, tol=mpf('1e-70'), log=None):
    """Solve the fixed point by continuation through DEFAULT_LADDER then through
    each N in target_Ns (ascending). Returns dict N -> (a, alpha, iters, resid)."""
    a_cur = [mpf(1), mpf('-1.6')]
    al_cur = mpf('1.5')
    out = {}
    ladder = sorted(set(DEFAULT_LADDER + list(target_Ns)))
    for N in ladder:
        a_cur, al_cur, it, res = solve_fixedpoint_mp(d, N, a_cur, al_cur, tol=tol, maxit=300, log=log)
        if log:
            log('  FP N=%d converged it=%d resid=%.3e alpha=%s'
                % (N, it, float(res), nstr(al_cur, 30)))
        if N in target_Ns:
            out[N] = (list(a_cur), al_cur, it, res)
    return out


def continue_eigen_mp(d, a_by_N, alpha_by_N, target_Ns, lam_seed, tol=mpf('1e-70'), log=None):
    """Solve the eigenpair by continuation through DEFAULT_LADDER then each N in
    target_Ns, using the fixed point (a,alpha) AT THE SAME N (a_by_N, alpha_by_N
    dicts must contain all ladder N's and target N's)."""
    c_cur = None
    lam_cur = mpf(lam_seed)
    out = {}
    ladder = sorted(set(DEFAULT_LADDER + list(target_Ns)))
    for N in ladder:
        a_N = a_by_N[N]
        al_N = alpha_by_N[N]
        c_cur, lam_cur, it, res = solve_eigen_mp(d, a_N, al_N, N, c_cur, lam_cur, tol=tol, maxit=300, log=log)
        if log:
            log('  EIG N=%d converged it=%d resid=%.3e lambda2=%s'
                % (N, it, float(res), nstr(lam_cur, 30)))
        if N in target_Ns:
            out[N] = (list(c_cur), lam_cur, it, res)
    return out


# --------------------------------------------------------------------------- #
#  Krawczyk (rigorous): K(X) = xhat - Y*F(xhat) + (I - Y*J(X))*(X - xhat)
#  K(X) subset int(X)  =>  unique zero of F in X (Krawczyk/Moore theorem).
# --------------------------------------------------------------------------- #
def krawczyk_generic(xhat_mp_list, Jmid_mp, F_at_point, F_J_on_box, dps, box_r_str):
    """xhat_mp_list: list of mpf (the unknowns' converged real values).
    Jmid_mp: mpmath (N x N) Jacobian at that point (real, used only to build
    the preconditioner -- need not be exact).
    F_at_point(xhat_arb_list) -> list of arb (residual at the point, arb balls
    possibly non-thin if chained parameters are ball-valued).
    F_J_on_box(X_arb_list) -> (F list, J rows) evaluated with X ball-valued.
    Returns dict with contained(bool), halfwidths, K image, X box."""
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
    return {
        'contained': contained, 'fail_idx': fail_idx, 'K': Kimg, 'X': X_vec,
        'xhat_arb': xhat_arb, 'halfwidths': halfwidths, 'Y': Y, 'JX': JX,
    }


def certify_fixedpoint(d, a_mid, alpha_mid, N, dps, box_r_str, log=None):
    """K1: fixed point. Unknowns = a_1..a_N, alpha (n = N+1)."""
    pts_mp = colloc_points_mp(N + 1)
    pts_arb = colloc_points_arb(N + 1, dps)
    xhat_mp = list(a_mid[1:]) + [alpha_mid]
    Fmp, Jmp = fp_residual_jacobian_mp(a_mid, alpha_mid, pts_mp, d)

    def F_at_point(v):
        a = [arb(1)] + v[:N]
        alpha = v[N]
        F, _ = fp_F_J_arb(a, alpha, pts_arb, d)
        return F

    def F_J_on_box(v):
        a = [arb(1)] + v[:N]
        alpha = v[N]
        return fp_F_J_arb(a, alpha, pts_arb, d)

    res = krawczyk_generic(xhat_mp, Jmp, F_at_point, F_J_on_box, dps, box_r_str)
    if log:
        log('  K1(d=%d,N=%d): contained=%s box_r=%s alpha_halfwidth=%.3e'
            % (d, N, res['contained'], box_r_str, res['halfwidths'][N] if res['contained'] else float('nan')))
    return res


DEFAULT_RADII = ['1e-45', '1e-40', '1e-35', '1e-30', '1e-25', '1e-20', '1e-15', '1e-10', '1e-7', '1e-4']


def certify_fixedpoint_auto(d, a_mid, alpha_mid, N, dps, log=None, radii=None):
    for r in (radii or DEFAULT_RADII):
        res = certify_fixedpoint(d, a_mid, alpha_mid, N, dps, r, log=log)
        if res['contained']:
            res['box_r_used'] = r
            return res
    raise RuntimeError('K1 fixed-point Krawczyk failed for all tried box radii')


def certify_eigen_auto(d, a_mid, alpha_mid, c_mid, lam_mid, N, dps, aX_box, alphaX_box,
                        log=None, radii=None):
    for r in (radii or DEFAULT_RADII):
        res = certify_eigen(d, a_mid, alpha_mid, c_mid, lam_mid, N, dps, r, aX_box, alphaX_box, log=log)
        if res['contained']:
            res['box_r_used'] = r
            return res
    raise RuntimeError('K2 eigenpair Krawczyk failed for all tried box radii')


def certify_eigen(d, a_mid, alpha_mid, c_mid, lam_mid, N, dps, box_r_str,
                   aX_box, alphaX_box, log=None):
    """K2: eigenpair, CHAINED on stage-1 box (aX_box: length N+1 arb list
    including a0=1 exactly; alphaX_box: arb, both ball-valued from K1's proven
    box). Unknowns = c_1..c_N, lambda (n = N+1)."""
    pts_mp = colloc_points_mp(N)
    pts_arb = colloc_points_arb(N, dps)
    xhat_mp = list(c_mid) + [lam_mid]
    Fmp, Jmp = eig_residual_jacobian_mp(a_mid, alpha_mid, c_mid, lam_mid, pts_mp, d)
    G1_box = G_eval_arb(aX_box, arb(1))

    def F_at_point(v):
        c = v[:N]
        lam = v[N]
        F, _ = eig_F_J_arb(aX_box, alphaX_box, c, lam, pts_arb, d, G1_box)
        return F

    def F_J_on_box(v):
        c = v[:N]
        lam = v[N]
        return eig_F_J_arb(aX_box, alphaX_box, c, lam, pts_arb, d, G1_box)

    res = krawczyk_generic(xhat_mp, Jmp, F_at_point, F_J_on_box, dps, box_r_str)
    if log:
        log('  K2(d=%d,N=%d): contained=%s box_r=%s lambda2_halfwidth=%.3e'
            % (d, N, res['contained'], box_r_str, res['halfwidths'][N] if res['contained'] else float('nan')))
    return res


# --------------------------------------------------------------------------- #
#  Independent verifier: preconditioned interval GAUSS-SEIDEL (Hansen-Sengupta).
#  Different containment criterion/code path than the Krawczyk MATRIX formula:
#  solves the preconditioned linear system Y*J(X)*(X-xhat) = -Y*F(xhat)
#  component-by-component with interval division, sequentially tightening each
#  component and re-substituting -- NOT the single matrix expression above.
# --------------------------------------------------------------------------- #
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
        # B = Y * JX  (n x n)
        B = [[arb(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                s = arb(0)
                for k in range(n):
                    s = s + Y[i][k] * JXrows[k][j]
                B[i][j] = s
        for i in range(n):
            # row i:  B[i,i]*(x_i - xhat_i) = -YFh[i] - sum_{j!=i} B[i,j]*(X_j - xhat_i_j)
            rhs = -YFh[i]
            for j in range(n):
                if j == i:
                    continue
                rhs = rhs - B[i][j] * (X[j] - xhat_arb[j])
            bii = B[i][i]
            if bii.contains(0):
                # cannot divide rigorously -- Gauss-Seidel step inconclusive for
                # this component this sweep; keep previous X[i]
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

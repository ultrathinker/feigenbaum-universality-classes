#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v4_tail.py -- TAIL STAGE: the N->infty gap for lambda2 (BRIEF-EXEC-4-KILO-TAIL).

Part A (float, mpmath dps 60): N=64 continuation solve, Chebyshev tail b_49..b_64,
  geometric decay rate, lambda2(N) ladder 56/64 via deflation chain + Newton refine
  (eig-precision insufficient for 1e-53-scale ladder diffs), sup-diffs ||G_N-G_M||_inf.
  -> hypothesis-grade numbers for the true fixed point tail (NOT a certificate).

Part B (balls, python-flint arb/acb, dps 130): rigorous pieces of the reduction theorem.
  B2: F_N(t) = G_N(t) + alpha*G_N(W_N(t)) as a FUNCTION on [0,1]:
      8192-point Gauss-Chebyshev; exact discrete-orthogonality identity
      (2D <= 2M-2) gives sum a_k^2 = (2/M)sum F_i^2, hence
      ||F||_inf <= sqrt(D+1) * ((2/M) sum F_i^2)^(1/2), D = 4608.
  B3: (i) mini-K2 Krawczyk with tight boxes (1e-30/1e-31, fallback 1e-18 certified
      config of v3) -> certified tight (v, lambda) box;
      (ii) r_phi = DT[phi~] - lambda~ phi~ as a FUNCTION on [0,1]:
      16384-point Gauss-Chebyshev, same Parseval bound, D = 9168.
  B4: tail sensitivity R(tau1,tau2): rigorous decomposed bound of
      DT_{G*}[phi~] - DT_{G_N}[phi~] under |G*-G_N| <= tau1, |G*'-G_N'| <= tau2
      (uses Markov bound sup|G_N''| <= 4 sum k^3 |b_k| and sup|phi~'| <= 2 sum j|v_j|).
  B5: discrete pencil resolvent S_disc(d) = ||(W - z V)^{-1} V||_inf on circles
      |z - lambda2~| = d -- rigorous in balls for the DISCRETE N=48 pencil
      (guide constant for hypothesis H2 of the reduction theorem).

Part C (float probe): Bernstein-ellipse containment W(E_rho), |G|,|G'| on E_rho --
  feasibility documentation for a future Lanford/EKW-style analytic tail validation.

Usage: python v4_tail.py [A|B|C|all]
"""
import json
import math
import statistics
import time

from mpmath import (mp, mpf, mpc, mpmathify, matrix as mpmatrix, lu_solve, cos as mcos,
                    sin as msin, acos as macos, pi as mppi, nstr, sqrt as msqrt)

mp.dps = 60

TAG = '[kid:verify] TAIL stage: N->inf gap for lambda2 (BRIEF-EXEC-4-KILO-TAIL)'
T0 = time.time()
LOG_PATH = 'v4_run.log'
ARG = 'all'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


# =====================================================================================
# shared float machinery (from v1_cheb.py, unchanged semantics)
# =====================================================================================

def cheb_vals(Nn, npts):
    return [(mcos(mppi * (i + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


def clenshaw(b, s):
    d1 = mpf(0) if not isinstance(b[0], mpc) else mpc(0)
    d2 = d1
    for k in range(len(b) - 1, 0, -1):
        d1, d2 = 2 * s * d1 - d2 + b[k], d1
    return s * d1 - d2 + b[0]


def dcheb_coeffs(b):
    Nn = len(b) - 1
    c = [mpf(0)] * (Nn + 2)
    for k in range(Nn, 0, -1):
        c[k - 1] = c[k + 1] + 2 * k * b[k]
    c[0] = c[0] / 2
    return c[:Nn]


class ChebFun:
    def __init__(self, b):
        self.b = list(b)
        self._db = None

    def G(self, t):
        return clenshaw(self.b, 2 * t - 1)

    def Gp(self, t):
        if self._db is None:
            self._db = dcheb_coeffs(self.b)
        return 2 * clenshaw(self._db, 2 * t - 1)


def set_from_free(x):
    Nn = len(x)
    b = [mpf(0)] + [mpf(v) for v in x]
    b[0] = mpf(1) - sum(((-1) ** k) * b[k] for k in range(1, Nn + 1))
    return b


def residual_jacobian(x, alpha, pts):
    b = set_from_free(x)
    db = dcheb_coeffs(b)
    Nn = len(x)
    npts = len(pts)
    F = mpmatrix(npts, 1)
    J = mpmatrix(npts, Nn + 1)
    a2 = alpha ** 2
    sgn = [(-1) ** k for k in range(Nn + 1)]
    for i, t in enumerate(pts):
        ta = t / a2
        U = clenshaw(b, 2 * ta - 1)
        W = U ** 2
        sW = 2 * W - 1
        GW = clenshaw(b, sW)
        GpW = 2 * clenshaw(db, sW)
        Gp_ta = 2 * clenshaw(db, 2 * ta - 1)
        F[i] = clenshaw(b, 2 * t - 1) + alpha * GW
        Tprev_t = mpf(1)
        Tcur_t = 2 * t - 1
        Tprev_W = mpf(1)
        Tcur_W = sW
        s_ta = 2 * ta - 1
        Tprev_ta = mpf(1)
        Tcur_ta = s_ta
        for j in range(1, Nn + 1):
            if j == 1:
                dG_t = Tcur_t - sgn[1]
                dG_W = Tcur_W - sgn[1]
                dG_ta = Tcur_ta - sgn[1]
            else:
                Tnext_t = 2 * (2 * t - 1) * Tcur_t - Tprev_t
                Tnext_W = 2 * sW * Tcur_W - Tprev_W
                Tnext_ta = 2 * s_ta * Tcur_ta - Tprev_ta
                dG_t = Tnext_t - sgn[j]
                dG_W = Tnext_W - sgn[j]
                dG_ta = Tnext_ta - sgn[j]
                Tprev_t, Tcur_t = Tcur_t, Tnext_t
                Tprev_W, Tcur_W = Tcur_W, Tnext_W
                Tprev_ta, Tcur_ta = Tcur_ta, Tnext_ta
            J[i, j - 1] = dG_t + alpha * (dG_W + 2 * U * GpW * dG_ta)
        J[i, Nn] = GW - 4 * t * U * GpW * Gp_ta / a2
    return F, J


def solve_fixedpoint(Nn, x_init=None, alpha_init=None, tol='1e-48'):
    pts = cheb_vals(Nn, Nn + 1)
    if x_init is None:
        x = [mpf(0)] * Nn
        x[0] = mpf('-0.75')
        alpha = mpf('2.5')
    else:
        x = [mpf(v) for v in x_init] + [mpf(0)] * (Nn - len(x_init))
        alpha = mpf(alpha_init)
    tol = mpf(tol)
    for it in range(300):
        F, J = residual_jacobian(x, alpha, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < tol:
            return x, alpha, it, nrm
        rhs = mpmatrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        improved = False
        for _halve in range(30):
            x_new = [x[j] + lam * dx[j] for j in range(Nn)]
            al_new = alpha + lam * dx[Nn]
            try:
                F2, _ = residual_jacobian(x_new, al_new, pts)
                nrm2 = max(abs(F2[i]) for i in range(len(pts)))
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                x, alpha, nrm = x_new, al_new, nrm2
                improved = True
                break
            lam /= 2
        if not improved:
            raise RuntimeError('Newton stalled at res=%s' % nstr(nrm, 5))
    raise RuntimeError('Newton no-converge res=%s' % nstr(nrm, 5))


def DT_apply(cf, alpha, phival, t):
    a = alpha
    a2 = a ** 2
    U = cf.G(t / a2)
    W = U ** 2
    GpW = cf.Gp(W)
    Gp_ta = cf.Gp(t / a2)
    return (a * phival(1) * cf.G(t)
            + 4 * t * phival(1) * U * GpW * Gp_ta
            - a * phival(W)
            - 2 * a * U * GpW * phival(t / a2))


def Bfun(j, t):
    s = 2 * t - 1
    Tprev, Tcur = mpf(1), s
    for k in range(1, j):
        Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
    return Tcur - ((-1) ** j)


def build_eig_matrix(Nn, x, alpha):
    cf = ChebFun(set_from_free(x))
    epts = cheb_vals(0, Nn)
    V = mpmatrix(Nn, Nn)
    for i, t in enumerate(epts):
        s = 2 * t - 1
        Tprev = mpf(1)
        Tcur = s
        for j in range(1, Nn + 1):
            if j == 1:
                V[i, 0] = Tcur + 1
            else:
                Tnext = 2 * s * Tcur - Tprev
                V[i, j - 1] = Tnext - ((-1) ** j)
                Tprev, Tcur = Tcur, Tnext
    cols = []
    for j in range(1, Nn + 1):
        vals = mpmatrix(Nn, 1)
        for i, t in enumerate(epts):
            def phival(tt, j=j):
                return Bfun(j, tt)
            vals[i] = DT_apply(cf, alpha, phival, t)
        coef = lu_solve(V, vals)
        cols.append([coef[i] for i in range(Nn)])
    A = mpmatrix(Nn, Nn)
    for jc in range(Nn):
        for i in range(Nn):
            A[i, jc] = cols[jc][i]
    return A


def matvec(M, v):
    n = len(v)
    w = mpmatrix(n, 1)
    for i in range(n):
        s = mpf(0)
        for j in range(n):
            s += M[i, j] * v[j]
        w[i] = s
    return w


def power_iter(M, n, iters=400, tol=mpf('1e-12')):
    v = mpmatrix(n, 1)
    for i in range(n):
        v[i] = mpf(1) / (i + 3)
    lam = mpf(0)
    for it in range(iters):
        w = matvec(M, v)
        lw = max(abs(w[i]) for i in range(n))
        if lw == 0:
            return mpf(0), v, it
        for i in range(n):
            v[i] = w[i] / lw
        if abs(lw - lam) < tol * lw:
            return lw, v, it
        lam = lw
    return lam, v, iters


def wielandt(M, lam1, v1, n):
    r = max(range(n), key=lambda i: abs(v1[i]))
    D = mpmatrix(n, n)
    vr = v1[r]
    for i in range(n):
        for j in range(n):
            D[i, j] = M[i, j] - (lam1 * v1[i] / vr) * (1 if j == r else 0)
    return D


def cached_DT_rows(Nn, b, alpha, epts):
    """rows[i][j-1] = DT[B_j](t_i); Bj[i][j-1] = B_j(t_i) -- one pass, node-cached.
    Formula = v3 DTcol_mp (fd-validated in v1); cross-checked vs DT_apply in part_a."""
    a2 = alpha ** 2
    sgn = [(-1) ** k for k in range(Nn + 1)]
    db = dcheb_coeffs(b)
    DTB = []
    BB = []
    for t in epts:
        ta = t / a2
        U = clenshaw(b, 2 * ta - 1)
        W = U * U
        sW = 2 * W - 1
        GpW = 2 * clenshaw(db, sW)
        Gp_ta = 2 * clenshaw(db, 2 * ta - 1)
        Gt = clenshaw(b, 2 * t - 1)

        def tlist(s):
            L = [mpf(1), s]
            Tp, Tc = mpf(1), s
            for j in range(2, Nn + 1):
                Tp, Tc = Tc, 2 * s * Tc - Tp
                L.append(Tc)
            return L

        Lt = tlist(2 * t - 1)
        LW = tlist(sW)
        Lta = tlist(2 * ta - 1)
        row = []
        rowb = []
        for j in range(1, Nn + 1):
            phi1 = mpf(1 - sgn[j])
            row.append(alpha * phi1 * Gt + 4 * t * phi1 * U * GpW * Gp_ta
                       - alpha * (LW[j] - sgn[j]) - 2 * alpha * U * GpW * (Lta[j] - sgn[j]))
            rowb.append(Lt[j] - sgn[j])
        DTB.append(row)
        BB.append(rowb)
    return DTB, BB


def eigen_newton_refine(Nn, cf, alpha, v0, lam0, DTB, BB):
    """Newton on F(v,lam) = [sum_j v_j DT[B_j](t_i) - lam*B-vals ; odd-norm-1].
    Uses cached DT rows. Returns (v, lam, res)."""
    v = [mpf(u) for u in v0]
    lam = mpf(lam0)
    n = Nn
    for it in range(30):
        F = mpmatrix(n + 1, 1)
        for i in range(n):
            dtv = sum(v[j - 1] * DTB[i][j - 1] for j in range(1, n + 1))
            phv = sum(v[j - 1] * BB[i][j - 1] for j in range(1, n + 1))
            F[i] = dtv - lam * phv
        F[n] = sum(2 * v[j - 1] for j in range(1, n + 1, 2)) - 1
        nrm = max(abs(F[i]) for i in range(n + 1))
        if nrm < mpf('1e-90'):
            return v, lam, nrm
        J = mpmatrix(n + 1, n + 1)
        for i in range(n):
            for j in range(1, n + 1):
                J[i, j - 1] = DTB[i][j - 1] - lam * BB[i][j - 1]
            J[i, n] = -sum(v[j - 1] * BB[i][j - 1] for j in range(1, n + 1))
        for j in range(1, n + 1):
            J[n, j - 1] = mpf(2) if j % 2 == 1 else mpf(0)
        J[n, n] = mpf(0)
        dx = lu_solve(J, -F)
        v = [v[j] + dx[j] for j in range(n)]
        lam = lam + dx[n]
    return v, lam, nrm


def part_a():
    log('=== PART A (float dps=120 deep): N=64 continuation, tail rate, lambda2 ladder ===')
    mp.dps = 120
    DEEP_TOL = '1e-100'
    ladder = [4, 8, 14, 20, 28, 36, 44, 48, 52, 56, 60, 64]
    sols = {}
    x_cur, al_cur = None, None
    for Nn in ladder:
        x, al, it, res = solve_fixedpoint(Nn, x_cur, al_cur, DEEP_TOL)
        sols[Nn] = (x, al)
        log('N=%2d: alpha=%s (it %d, res %.2e)' % (Nn, nstr(al, 60), it, float(res)))
        x_cur, al_cur = x, al
    x40, al40, _, _ = solve_fixedpoint(40, sols[36][0], sols[36][1], DEEP_TOL)
    sols[40] = (x40, al40)

    # ---- tail coefficients at N=64 -------------------------------------------------
    x64, al64 = sols[64]
    b64 = set_from_free(x64)
    log('--- Chebyshev tail |b_k|, N=64 deep solution ---')
    tail_tab = {}
    for k in range(30, 65):
        tail_tab[k] = float(abs(b64[k]))
    log('  |b_41|=%.3e |b_48|=%.3e |b_49|=%.3e |b_56|=%.3e |b_64|=%.3e'
        % (tail_tab[41], tail_tab[48], tail_tab[49], tail_tab[56], tail_tab[64]))
    nz = [k for k in range(33, 65) if abs(b64[k]) > 0]
    k0, k1 = min(nz), max(nz)
    lrs = [float(abs(b64[k]) / abs(b64[k + 1])) for k in range(k0, k1) if abs(b64[k + 1]) > 0]
    if not lrs:
        raise RuntimeError('no nonzero tail ratios -- deep solve failed to populate tail')
    rho_min = min(lrs)
    rho_ls = math.exp(statistics.fmean(math.log(r) for r in lrs))
    log('per-mode ratio |b_k/b_{k+1}| k=%d..%d: min=%.4f mean(geo)=%.4f' % (k0, k1, rho_min, rho_ls))
    rho_env = 0.9 * rho_min
    B0 = max(float(abs(b64[k])) * rho_env ** k for k in nz)
    log('strict geometric envelope on k=%d..%d: |b_k| <= %.3e * %.4f^{-k}' % (k0, k1, B0, rho_env))
    tail_sum = B0 * rho_env ** (-(k1 + 1)) / (1 - 1 / rho_env)
    log('envelope tail after k=%d: sum <= %.3e (sup-norm on [0,1])' % (k1, tail_sum))

    # ---- lambda2 ladder: deflation chain level 2 + Newton refine ---------------------
    lam2 = {}
    for Nn in [40, 48, 56, 64]:
        x, al = sols[Nn]
        cf = ChebFun(set_from_free(x))
        epts = cheb_vals(0, Nn)
        bN = set_from_free(x)
        DTB, BB = cached_DT_rows(Nn, bN, al, epts)
        # cross-check cached DT formula vs DT_apply at 2 sample nodes (guard vs typos)
        for i_chk in [Nn // 3, 2 * Nn // 3]:
            j_chk = Nn // 2
            ref = DT_apply(cf, al, (lambda tt, jc=j_chk: Bfun(jc, tt)), epts[i_chk])
            diff = abs(DTB[i_chk][j_chk - 1] - ref)
            if diff > mpf('1e-45'):
                raise RuntimeError('cached DT formula mismatch: %s' % nstr(diff, 5))
        # matrix A = V^{-1} DTB (V = BB)
        DTBm = mpmatrix(Nn, Nn)
        Vm_ = mpmatrix(Nn, Nn)
        for i in range(Nn):
            for j in range(Nn):
                DTBm[i, j] = DTB[i][j]
                Vm_[i, j] = BB[i][j]
        colsA = []
        for jc in range(Nn):
            colv = mpmatrix(Nn, 1)
            for i in range(Nn):
                colv[i] = DTBm[i, jc]
            cc = lu_solve(Vm_, colv)
            colsA.append([cc[i, 0] for i in range(Nn)])
        A = mpmatrix(Nn, Nn)
        for jc in range(Nn):
            for i in range(Nn):
                A[i, jc] = colsA[jc][i]
        Mcur = A
        chain = []
        for level in range(3):
            lw, v, _it = power_iter(Mcur, Nn)
            w = matvec(Mcur, v)
            sign = sum(v[i] * w[i] for i in range(Nn)) / sum(v[i] * v[i] for i in range(Nn))
            chain.append((lw * (1 if sign > 0 else -1), v))
            Mcur = wielandt(Mcur, sign, v, Nn)
        lam_chain, v_chain = chain[2]
        # seed vector in B-basis: A v_chain = lam v_chain means V^{-1}DTB v = lam v,
        # i.e. phi = sum (v)_j B_j interpolates DT phi = lam phi at nodes -> use eigvec of A
        vals = mpmatrix(Nn, 1)
        for i in range(Nn):
            vals[i] = v_chain[i]
        v0 = lu_solve(Vm_, vals)
        s_odd = sum(v0[j - 1] for j in range(1, Nn + 1, 2))
        v0 = [u / (2 * s_odd) for u in v0]
        v_ref, lam_ref, res_ref = eigen_newton_refine(Nn, cf, al, v0, lam_chain, DTB, BB)
        lam2[Nn] = lam_ref
        log('lambda2(N=%d) = %s (newton res %.2e, chain %.12f)'
            % (Nn, nstr(lam_ref, 50), float(res_ref), float(lam_chain)))
    d48_40 = float(abs(lam2[48] - lam2[40]))
    d56_48 = float(abs(lam2[56] - lam2[48]))
    d64_56 = float(abs(lam2[64] - lam2[56]))
    d64_48 = float(abs(lam2[64] - lam2[48]))
    log('lambda2 ladder diffs: |48-40|=%.3e |56-48|=%.3e |64-56|=%.3e |64-48|=%.3e'
        % (d48_40, d56_48, d64_56, d64_48))

    # ---- sup-diffs G, G' ---------------------------------------------------------------
    grid = cheb_vals(0, 4001)
    diffs = {}
    for (Na, Nb) in [(40, 48), (48, 56), (56, 64), (48, 64)]:
        ca = ChebFun(set_from_free(sols[Na][0]))
        cb = ChebFun(set_from_free(sols[Nb][0]))
        dG = max(abs(ca.G(t) - cb.G(t)) for t in grid)
        dGp = max(abs(ca.Gp(t) - cb.Gp(t)) for t in grid)
        diffs['%d_%d' % (Na, Nb)] = {'dG': float(dG), 'dGp': float(dGp)}
        log('sup ||G%d-G%d||=%.3e ; ||Gp%d-Gp%d||=%.3e' % (Na, Nb, dG, Na, Nb, dGp))

    out_a = {
        'alpha_ladder': {str(k): nstr(sols[k][1], 50) for k in sorted(sols)},
        'tail_b_abs': {str(k): tail_tab[k] for k in tail_tab},
        'ratio_min_k40_63': rho_min, 'ratio_geomean': rho_ls,
        'envelope': {'B0': B0, 'rho_env': rho_env, 'tail_after_64': tail_sum},
        'lambda2_ladder': {str(k): nstr(lam2[k], 50) for k in lam2},
        'lambda2_diffs': {'48_40': d48_40, '56_48': d56_48, '64_56': d64_56, '64_48': d64_48},
        'sup_diffs': diffs,
        'b64_1_to_64': [nstr(v, 50) for v in x64],
        'alpha64': nstr(al64, 50),
    }
    return out_a


# =====================================================================================
# Part B -- balls
# =====================================================================================

def part_b():
    from flint import arb, acb, arb_mat, acb_mat, ctx
    ctx.dps = 130
    mp.dps = 140

    _v1 = json.load(open('v1_cheb.json', encoding='utf-8'))
    N = 48
    ALPHA_STR = _v1['alpha']
    X48_STRS = _v1['x48_cheb_b1_to_b48']
    LAM2_STR = _v1['lambda2']
    V2_STRS = _v1['v_lambda2_B']

    sgn_int = [(-1) ** k for k in range(N + 1)]

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
            t = (mcos(mppi * (mpf(i) + mpf('0.5')) / npts) + 1) / 2
            ts.append(arb(nstr(t, 120)))
        return ts

    class TrigNode:
        def __init__(self, t, b):
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

    def b0_of(bs):
        s = arb(1)
        for k in range(1, len(bs) + 1):
            s = s - arb(sgn_int[k]) * bs[k - 1]
        return s

    # ---- K1 box (same construction as v3, radius 1e-44; v3 K1 certificate covers it) --
    bhat = [arb(s) for s in X48_STRS]
    ahat = arb(ALPHA_STR)
    RBOX = '1e-44'
    Xb = [b + ball(RBOX) for b in bhat]
    Xa = ahat + ball(RBOX)
    bX = [b0_of(Xb)] + Xb
    a2X = Xa * Xa

    # ---------------- B2: F_N residual as a function, M=8192 -------------------------
    log('=== B2: F_N dense residual as FUNCTION (M=8192 Gauss-Chebyshev, balls) ===')
    M1 = 8192
    pts1 = gauss_nodes_arb(M1)
    degF = 48 * 96  # deg G_N(W_N) = 48*96 ; deg G_N = 48
    q_total = arb(0)  # sum of POSITIVE upper bounds of F_i^2
    fmax = arb(0)
    for it_, t in enumerate(pts1):
        if it_ % 1024 == 0 and it_ > 0:
            log('  B2 progress %d/%d' % (it_, M1))
        nta = TrigNode(t / a2X, bX)
        nW = TrigNode(nta.G * nta.G, bX)
        nt = TrigNode(t, bX)
        F = nt.G + Xa * nW.G
        aF = abs(F).upper()
        q_total = q_total + aF * aF
        if bool(aF > fmax):
            fmax = aF
    s2_total = q_total * (arb(2) / M1)
    supF = (arb(degF + 1) * s2_total).sqrt()
    log('  max_i |F(t_i)| ub = %s' % fmax.str(5))
    log('  sqrt((2/M) sum F^2) = %s' % s2_total.sqrt().upper().str(5))
    log('  => sup_[0,1] |F_N| <= sqrt(D+1)*.. = %s (D=%d, 2D<=%d=2M-2: %s)'
        % (supF.upper().str(5), degF, 2 * M1 - 2, 2 * degF <= 2 * M1 - 2))
    B2 = {'M': M1, 'degF': degF, 'max_pt': fmax.str(10), 'sup_bound': supF.upper().str(10)}

    # ---------------- B3(i): mini-K2 with tight boxes --------------------------------
    log('=== B3(i): mini-K2 Krawczyk, tight boxes ===')
    vh = [mpmathify(s).real for s in V2_STRS]
    s_odd = sum(vh[j - 1] for j in range(1, N + 1, 2))
    vh = [v / (2 * s_odd) for v in vh]
    lam2_hat = mpf(LAM2_STR)
    pts48 = gauss_nodes_arb(N)

    node_t = [TrigNode(t, bX) for t in pts48]
    node_ta = []
    node_W = []
    for t in pts48:
        nta = TrigNode(t / a2X, bX)
        node_ta.append(nta)
        node_W.append(TrigNode(nta.G * nta.G, bX))

    DTcol_cache = {}

    def DTcol_ball(j, ti):
        key = (j, ti)
        if key in DTcol_cache:
            return DTcol_cache[key]
        nt, nta_, nW = node_t[ti], node_ta[ti], node_W[ti]
        t = pts48[ti]
        U = nta_.G
        GpW = nW.Gp
        Gp_ta = nta_.Gp
        phi1 = arb(1 - sgn_int[j])
        val = (Xa * phi1 * nt.G
               + 4 * t * phi1 * U * GpW * Gp_ta
               - Xa * (nW.T[j] - arb(sgn_int[j]))
               - 2 * Xa * U * GpW * (nta_.T[j] - arb(sgn_int[j])))
        DTcol_cache[key] = val
        return val

    def eigen_F_ball(vb, lam):
        out = []
        for ti in range(N):
            t = pts48[ti]
            nt = node_t[ti]
            dt = arb(0)
            ph = arb(0)
            for j in range(1, N + 1):
                c = vb[j - 1]
                dt = dt + c * DTcol_ball(j, ti)
                ph = ph + c * (nt.T[j] - arb(sgn_int[j]))
            out.append(dt - lam * ph)
        norm = arb(-1)
        for j in range(1, N + 1, 2):
            norm = norm + 2 * vb[j - 1]
        out.append(norm)
        return out

    b_mp = [mpf(1) - sum(((-1) ** k) * mpf(X48_STRS[k - 1]) for k in range(1, N + 1))] + [mpf(s) for s in X48_STRS]
    a_mp = mpf(ALPHA_STR)
    a2_mp = a_mp ** 2
    t48_mp = [(mcos(mppi * (mpf(i) + mpf('0.5')) / N) + 1) / 2 for i in range(N)]

    class TrigNodeMp:
        def __init__(self, t, b):
            s = 2 * t - 1
            th = macos(s)
            self.T = [mpf(1)]
            for k in range(1, N + 1):
                self.T.append(mcos(k * th))
            st = msin(th)
            self.G = b[0]
            dsum = mpf(0)
            for k in range(1, N + 1):
                self.G = self.G + b[k] * self.T[k]
                dsum = dsum + k * b[k] * msin(k * th)
            self.Gp = 2 * dsum / st

    node_mp = []
    for t in t48_mp:
        nt = TrigNodeMp(t, b_mp)
        nta_ = TrigNodeMp(t / a2_mp, b_mp)
        nW = TrigNodeMp(nta_.G * nta_.G, b_mp)
        node_mp.append((nt, nta_, nW))

    def DTcol_mp(j, i):
        nt, nta_, nW = node_mp[i]
        t = t48_mp[i]
        U = nta_.G
        phi1 = mpf(1 - sgn_int[j])
        return (a_mp * phi1 * nt.G
                + 4 * t * phi1 * U * nW.Gp * nta_.Gp
                - a_mp * (nW.T[j] - sgn_int[j])
                - 2 * a_mp * U * nW.Gp * (nta_.T[j] - sgn_int[j]))

    def eigen_J_mp(v, lam):
        Mx = mpmatrix(N + 1, N + 1)
        for i in range(N):
            nt = node_mp[i][0]
            for j in range(1, N + 1):
                Mx[i, j - 1] = DTcol_mp(j, i) - lam * (nt.T[j] - sgn_int[j])
            ph = sum(v[j - 1] * (nt.T[j] - sgn_int[j]) for j in range(1, N + 1))
            Mx[i, N] = -ph
        for j in range(1, N + 1):
            Mx[N, j - 1] = mpf(2) if j % 2 == 1 else mpf(0)
        Mx[N, N] = mpf(0)
        return Mx

    def mp_inverse(Mx, n):
        cols = []
        for jcol in range(n):
            e = mpmatrix(n, 1)
            e[jcol, 0] = mpf(1)
            x = lu_solve(Mx, e)
            cols.append([x[i, 0] for i in range(n)])
        Y = mpmatrix(n, n)
        for jcol in range(n):
            for i in range(n):
                Y[i, jcol] = cols[jcol][i]
        return Y

    J2 = eigen_J_mp(vh, lam2_hat)
    Y2m = mp_inverse(J2, N + 1)
    Y2 = arb_mat([[arb(nstr(Y2m[i, j], 120)) for j in range(N + 1)] for i in range(N + 1)])

    def krawczyk_k2(RVs, RLs):
        Xv_ = [v + ball(RVs) for v in vhat_arb]
        Xlam_ = lamhat_arb + ball(RLs)
        Fh = eigen_F_ball(vhat_arb, lamhat_arb)
        JX = arb_mat(N + 1, N + 1)
        for i in range(N):
            nt = node_t[i]
            for j in range(1, N + 1):
                JX[i, j - 1] = DTcol_ball(j, i) - Xlam_ * (nt.T[j] - arb(sgn_int[j]))
            ph = arb(0)
            for j in range(1, N + 1):
                ph = ph + Xv_[j - 1] * (nt.T[j] - arb(sgn_int[j]))
            JX[i, N] = -ph
        for j in range(1, N + 1):
            JX[N, j - 1] = arb(2) if j % 2 == 1 else arb(0)
        JX[N, N] = arb(0)
        Fhm = arb_mat(N + 1, 1)
        for i in range(N + 1):
            Fhm[i, 0] = Fh[i]
        Xm = arb_mat(N + 1, 1)
        for i in range(N):
            Xm[i, 0] = Xv_[i]
        Xm[N, 0] = Xlam_
        xhm = arb_mat(N + 1, 1)
        for i in range(N):
            xhm[i, 0] = vhat_arb[i]
        xhm[N, 0] = lamhat_arb
        Km = xhm - Y2 * Fhm + (arb_eye(N + 1) - Y2 * JX) * (Xm - xhm)
        nc = 0
        rad = []
        for i in range(N + 1):
            Xvec_i = (Xv_ + [Xlam_])[i]
            if Xvec_i.contains_interior(Km[i, 0]):
                nc += 1
            rad.append((Km[i, 0].rad() + abs(Km[i, 0].mid() - Xvec_i.mid())).upper())
        lam_rad = rad[N]
        return nc == N + 1, lam_rad, Xv_, Xlam_

    vhat_arb = [arb(nstr(v, 105)) for v in vh]
    lamhat_arb = arb(nstr(lam2_hat, 110))
    k2_ok, lam_rad_tight, Xv, Xlam_cert = False, None, None, None
    for (RVs, RLs) in [('1e-33', '1e-34'), ('1e-30', '1e-31'), ('1e-24', '1e-25'), ('1e-18', '1e-21')]:
        ok, lrad, Xv_, Xlam_ = krawczyk_k2(RVs, RLs)
        log('  mini-K2 (r_v=%s, r_lam=%s): contained %d/%d, lam-radius %s'
            % (RVs, RLs, (N + 1) if ok else 0, N + 1, lrad.str(5)))
        if ok:
            k2_ok, lam_rad_tight, Xv, Xlam_cert = True, lrad, Xv_, Xlam_
            log('  tight certified box: r_v=%s ; lambda2 certified to +-%s' % (RVs, lrad.str(5)))
            break
    if not k2_ok:
        raise RuntimeError('mini-K2 failed at all box sizes -- unexpected (v3 certified 1e-18)')
    # certified lambda box (use max of this run's radius and v3's 3.02e-31)
    LAM_CERT_RAD = '3.02e-31'
    lam_box = lamhat_arb + ball(LAM_CERT_RAD)

    # ---------------- B3(ii): r_phi as a function, M=16384 ---------------------------
    log('=== B3(ii): r_phi = DT[phi~] - lam phi~ as FUNCTION (M=16384, balls) ===')
    M2 = 16384
    pts2 = gauss_nodes_arb(M2)

    def phi_at(node, vb):
        s = arb(0)
        for j in range(1, N + 1):
            s = s + vb[j - 1] * (node.T[j] - arb(sgn_int[j]))
        return s

    phi1_ball = arb(0)
    for j in range(1, N + 1, 2):
        phi1_ball = phi1_ball + 2 * Xv[j - 1]

    degR = 48 + 47 * 96 + 48 * 96  # deg(U G'(W) phi(t/a^2)) = 48 + 47*96 + 48*96 = 9168
    q_r = arb(0)
    rmax = arb(0)
    for it_, t in enumerate(pts2):
        if it_ % 2048 == 0 and it_ > 0:
            log('  B3 progress %d/%d' % (it_, M2))
        nta = TrigNode(t / a2X, bX)
        U = nta.G
        Wb = U * U
        nW = TrigNode(Wb, bX)
        nt = TrigNode(t, bX)
        phW = phi_at(nW, Xv)
        phta = phi_at(nta, Xv)
        pht = phi_at(nt, Xv)
        DTv = (Xa * phi1_ball * nt.G
               + 4 * t * phi1_ball * U * nW.Gp * nta.Gp
               - Xa * phW
               - 2 * Xa * U * nW.Gp * phta)
        r = DTv - lam_box * pht
        ar_ = abs(r).upper()
        q_r = q_r + ar_ * ar_
        if bool(ar_ > rmax):
            rmax = ar_
    s2_r = q_r * (arb(2) / M2)
    supR = (arb(degR + 1) * s2_r).sqrt()
    log('  max_i |r(t_i)| ub = %s' % rmax.str(5))
    log('  deg(r) bound D=%d ; 2D <= 2M-2: %s' % (degR, 2 * degR <= 2 * M2 - 2))
    log('  => sup_[0,1] |r| <= sqrt(D+1)*sqrt((2/M)sum r^2) = %s' % supR.upper().str(5))
    B3 = {'M': M2, 'degR': degR, 'max_pt': rmax.str(10), 'sup_bound': supR.upper().str(10),
          'mini_k2_lam_radius': str(lam_rad_tight), 'phi1_ub': phi1_ball.abs_upper().str(10)}

    # ---------------- B4: tail sensitivity R(tau1, tau2) -------------------------------
    log('=== B4: tail sensitivity R(tau1,tau2), decomposed (2048 nodes) ===')
    M4 = 2048
    pts4 = gauss_nodes_arb(M4)
    C2 = arb(0)
    for k in range(1, N + 1):
        C2 = C2 + arb(k ** 3) * abs(bX[k])
    C2 = 4 * C2
    Lipphi = arb(0)
    for j in range(1, N + 1):
        Lipphi = Lipphi + arb(j) * abs(Xv[j - 1])
    Lipphi = 2 * Lipphi
    alpha_rad = Xa.rad().upper()
    log('  C2 = sup|G_N"|" ub (4*sum k^3|b_k|) = %s' % C2.upper().str(8))
    log('  Lipphi = sup|phi~\'| ub (2*sum j|v_j|) = %s' % Lipphi.upper().str(8))
    log('  |dalpha| = %s' % alpha_rad.str(5))
    tau_grid = [('1e-12', '1e-11'), ('1e-20', '1e-19'), ('1e-30', '1e-29'),
                ('1e-40', '1e-39'), ('1e-50', '1e-49'), ('1e-57', '1e-56'),
                ('1e-66', '1e-65'), ('1e-76', '1e-75')]
    B4 = {'C2': C2.upper().str(10), 'Lipphi': Lipphi.upper().str(10), 'dalpha': alpha_rad.str(10)}
    for (t1s, t2s) in tau_grid:
        t1 = arb(t1s)  # POSITIVE magnitude ball (thin)
        t2 = arb(t2s)
        rmax4 = arb(0)
        for it_, t in enumerate(pts4):
            nta = TrigNode(t / a2X, bX)
            U = nta.G
            Wb = U * U
            nW = TrigNode(Wb, bX)
            nt = TrigNode(t, bX)
            phW = phi_at(nW, Xv)
            phta = phi_at(nta, Xv)
            aU = abs(U).upper()
            aGpW = abs(nW.Gp).upper()
            aGpta = abs(nta.Gp).upper()
            aph1 = abs(phi1_ball).upper()
            aphW = abs(phW).upper()
            aphta = abs(phta).upper()
            aGt = abs(nt.G).upper()
            dU = t1
            dW = 2 * aU * t1 + t1 * t1
            dGpW = t2 + C2.upper() * dW
            Ust = aU + t1
            termA = abs(Xa).upper() * aph1 * t1
            termB = 4 * t * aph1 * (dU * aGpW * aGpta + Ust * dGpW * aGpta + Ust * (aGpW + dGpW) * t2)
            termC = abs(Xa).upper() * Lipphi.upper() * dW
            termD = 2 * abs(Xa).upper() * aphta * (dU * aGpW + Ust * dGpW)
            termE = alpha_rad * (aph1 * aGt + 4 * t * aph1 * aU * aGpW * aGpta + aphW + 2 * aU * aGpW * aphta)
            tot = (termA + termB + termC + termD + termE).upper()
            if float(tot) > float(rmax4):
                rmax4 = tot
        # safety factor for float-based max selection across 2048 positive balls
        rmax4 = rmax4 * arb('1.000000000001')
        B4['tau_%s' % t1s] = rmax4.upper().str(10)
        log('  R(tau1=%s, tau2=%s) <= %s' % (t1s, t2s, rmax4.upper().str(5)))

    # ---------------- B5: discrete pencil resolvent -------------------------------------
    log('=== B5: discrete pencil resolvent S_disc(d) (balls, N=48 pencil) ===')
    Wb_mat = arb_mat(N, N)
    Vb_mat = arb_mat(N, N)
    for i in range(N):
        nt = node_t[i]
        for j in range(1, N + 1):
            Wb_mat[i, j - 1] = DTcol_ball(j, i)
            Vb_mat[i, j - 1] = nt.T[j] - arb(sgn_int[j])

    def mat_norm_ub(Mx):
        n = Mx.nrows()
        mx = arb(0)
        for i in range(n):
            rs = arb(0)
            for j in range(n):
                rs = rs + abs(Mx[i, j]).upper()
            mx = max(mx, rs)
        return mx

    Wm = mpmatrix(N, N)
    Vm = mpmatrix(N, N)
    for i in range(N):
        for j in range(1, N + 1):
            Wm[i, j - 1] = DTcol_mp(j, i)
            Vm[i, j - 1] = node_mp[i][0].T[j] - sgn_int[j]

    Vnub = mat_norm_ub(Vb_mat)
    log('  ||V||_inf ub = %s' % Vnub.upper().str(8))
    B5 = {'V_norm': Vnub.upper().str(10)}
    for ds in ['0.02', '0.05', '0.065']:
        d = mpf(ds)
        Sm_d = arb(0)
        worst = None
        n_ok = 0
        for k in range(12):
            th = 2 * mppi * k / 12
            zr = lam2_hat + d * mcos(th)
            zi = d * msin(th)
            zc = mpc(zr, zi)
            Az_mp = mpmatrix(N, N)
            for i in range(N):
                for j in range(N):
                    Az_mp[i, j] = Wm[i, j] - zc * Vm[i, j]
            cols = []
            for jcol in range(N):
                e = mpmatrix(N, 1)
                e[jcol, 0] = mpf(1)
                x = lu_solve(Az_mp, e)
                cols.append([x[i, 0] for i in range(N)])
            # ||Y||_inf in mpmath, exact-ish, then as string -> arb
            Ynorm_mp = max(sum(abs(cols[jcol][i]) for jcol in range(N)) for i in range(N))
            Ynorm = arb(nstr(Ynorm_mp * (1 + mpf('1e-30')), 60))
            z_ball = acb(arb(nstr(zr, 120)), arb(nstr(zi, 120)))
            Az_ball = acb_mat(N, N)
            for i in range(N):
                for j in range(N):
                    Az_ball[i, j] = acb(Wb_mat[i, j]) - z_ball * acb(Vb_mat[i, j])
            Y_ball = acb_mat(N, N)
            for i in range(N):
                for jcol in range(N):
                    c = cols[jcol][i]
                    Y_ball[i, jcol] = acb(arb(nstr(c.real, 120)), arb(nstr(c.imag, 120)))
            E_ball = Y_ball * Az_ball
            En = arb_mat(N, N)
            for i in range(N):
                for j in range(N):
                    En[i, j] = abs(E_ball[i, j]).upper()
            for i in range(N):
                En[i, i] = En[i, i] - 1
            Enub = mat_norm_ub(En)
            if bool(Enub.upper() >= arb('0.5')):
                log('  d=%s k=%d: ||I - Y A||_ub = %s >= 0.5 -- point skipped (documented)'
                    % (ds, k, Enub.upper().str(5)))
                continue
            n_ok += 1
            inv_bound = Ynorm * (arb(1) / (arb(1) - Enub.upper()))
            S_pt = inv_bound * Vnub
            if bool(S_pt > Sm_d):
                Sm_d = S_pt
                worst = k
        B5['d_%s' % ds] = {'S_ub': Sm_d.upper().str(10), 'n_points': n_ok, 'worst_k': worst}
        log('  S_disc(d=%s) <= %s (points ok %d/12, worst k=%s)' % (ds, Sm_d.upper().str(8), n_ok, worst))

    return {'B2_FN_residual': B2, 'B3_rphi_residual': B3, 'B4_tail_sensitivity': B4,
            'B5_discrete_resolvent': B5,
            'meta': {'arb_dps': 130, 'mp_dps': 140, 'K1_box_radius': RBOX,
                     'lam2_cert_radius_used': LAM_CERT_RAD}}


# =====================================================================================
# Part C -- float probe: Bernstein ellipse containment
# =====================================================================================

def part_c(out_a):
    log('=== PART C (float probe): Bernstein-ellipse containment / Lanford feasibility ===')
    b64 = [mpf(s) for s in out_a['b64_1_to_64']]
    bfull = [mpf(1) - sum(((-1) ** k) * b64[k - 1] for k in range(1, 65))] + b64
    al = mpf(out_a['alpha64'])
    a2 = al ** 2
    cf = ChebFun(bfull)
    db = dcheb_coeffs(bfull)

    def Gp_complex(t):
        return 2 * clenshaw(db, 2 * t - 1)

    def mu_of(s):
        w = s + msqrt(s ** 2 - 1)
        if abs(w) < 1:
            w = s - msqrt(s ** 2 - 1)
        return abs(w)

    res = {}
    for rho_str in ['1.05', '1.1', '1.2', '1.3', '1.5', '2.0', '3.0']:
        rho = mpf(rho_str)
        muW_max = mpf(0)
        aG_max = mpf(0)
        aGp_max = mpf(0)
        mu_inner_max = mpf(0)
        for k in range(96):
            phi = 2 * mppi * k / 96
            w = rho * mpc(mcos(phi), msin(phi))
            s = (w + 1 / w) / 2
            t = (1 + s) / 2
            ta = t / a2
            U = cf.G(ta)
            Wc = U * U
            sW = 2 * Wc - 1
            muW_max = max(muW_max, mu_of(sW))
            aG_max = max(aG_max, abs(U))
            aGp_max = max(aGp_max, abs(Gp_complex(ta)))
            mu_inner_max = max(mu_inner_max, mu_of(2 * ta - 1))
        res[rho_str] = {'mu_W_max': float(muW_max), 'absG_max': float(aG_max),
                        'absGp_max': float(aGp_max), 'mu_inner_max': float(mu_inner_max)}
        log('  rho=%4s: mu(W(E_rho))<=%.4f ; |G|<=%.3f ; |Gp|<=%.3f ; mu(inner)<=%.4f'
            % (rho_str, float(muW_max), float(aG_max), float(aGp_max), float(mu_inner_max)))
    return res


def main():
    log(TAG)
    out = {'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S')}}
    if ARG in ('A', 'all'):
        out['partA'] = part_a()
    if ARG in ('B', 'all'):
        out['partB'] = part_b()
    if ARG in ('C', 'all'):
        if 'partA' not in out:
            prev = json.load(open('v4_tail.json', encoding='utf-8'))
            out['partA'] = prev['partA']
            out['meta']['partA_reused'] = True
        out['partC'] = part_c(out['partA'])
    out['meta']['elapsed_sec'] = time.time() - T0
    try:
        merged = json.load(open('v4_tail.json', encoding='utf-8'))
    except Exception:
        merged = {}
    merged.update(out)
    merged['meta'].update(out['meta'])
    with open('v4_tail.json', 'w', encoding='utf-8') as f:
        json.dump(merged, f, indent=2, ensure_ascii=False)
    log('wrote v4_tail.json (merged parts: %s)' % sorted(k for k in merged if k != 'meta'))


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        ARG = sys.argv[1]
    main()

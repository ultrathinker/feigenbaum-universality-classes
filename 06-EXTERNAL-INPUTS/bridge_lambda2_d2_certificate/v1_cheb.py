#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v1_cheb.py -- INDEPENDENT Phase-1 verification of agent A's lambda2-degree2 spectrum.
[kid:verify] NOT RIGOROUS float (mpmath dps=60) for the solve itself; interval (python-flint)
part is done separately in v3_krawczyk.py. Code path deliberately DIFFERENT from agent A:

  * basis: Chebyshev T_k(2t-1) on t in [0,1]  (agent A: raw powers t^j, Vandermonde-like)
  * linearisation formula: derived independently by chain rule, kept UNSIMPLIFIED
    (uses the product G'(t/a^2)*G'(W); agent A's simplified form uses K(t)=G-2tG'),
    validated by a finite difference of the operator R itself
  * eigenvalues: mpmath.eig AND a self-written one-sided Wielandt deflation + power
    iteration chain (agent A: two-sided deflation that broke after 2-3 levels)
  * coordinate-mode test: finite coordinate substitution x->x(1+eps) applied to the
    equation/operator numerically (NO K(t) formula anywhere in this file)

Equation: g(x) = -alpha*g(g(x/alpha)), g(0)=1, g even; G(t)=g(x), t=x^2.
"""
import json
import time

from mpmath import mp, mpf, mpc, mpmathify, matrix, lu_solve, cos, pi, nstr, sqrt

mp.dps = 60

TAG = '[kid:verify] independent Chebyshev-basis re-verification of agent A phase 1 (degree-2 CF spectrum)'
T0 = time.time()
LOG_PATH = 'v1_run.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


# --------------------------- Chebyshev utilities (t in [0,1], s=2t-1) ---------


def cheb_vals(N, npts):
    """Gauss-Chebyshev roots s_i = cos(pi*(i+1/2)/npts), t=(s+1)/2, i=0..npts-1."""
    return [(cos(pi * (i + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


def clenshaw(b, s):
    """sum_k b[k] T_k(s) via Clenshaw (b: list k=0..N)."""
    d1 = mpf(0) if not isinstance(b[0], mpc) else mpc(0)
    d2 = d1
    for k in range(len(b) - 1, 0, -1):
        d1, d2 = 2 * s * d1 - d2 + b[k], d1
    return s * d1 - d2 + b[0]


def dcheb_coeffs(b):
    """coeffs c of d/ds (sum b_k T_k) = sum c_k T_k (len N)."""
    N = len(b) - 1
    c = [mpf(0)] * (N + 2)
    for k in range(N, 0, -1):
        c[k - 1] = c[k + 1] + 2 * k * b[k]
    c[0] = c[0] / 2
    return c[:N]  # coefficients for T_0..T_{N-1}


class ChebFun:
    """G(t) = sum_k b_k T_k(2t-1); b[0] enforces G(0)=1 if locked."""

    def __init__(self, b, lock_b0=True):
        self.b = list(b)
        self.lock = lock_b0  # if True, b[0] = 1 - sum_{k>=1} (-1)^k b_k

    def set_from_free(self, x):
        """x = (b_1..b_N) free part; b0 locked."""
        N = len(x)
        b = [mpf(0)] + [mpf(v) for v in x]
        if self.lock:
            b[0] = mpf(1) - sum(((-1) ** k) * b[k] for k in range(1, N + 1))
        self.b = b
        self._db = None

    def G(self, t):
        return clenshaw(self.b, 2 * t - 1)

    def Gp(self, t):
        """dG/dt = 2 d/ds."""
        if getattr(self, '_db', None) is None:
            self._db = dcheb_coeffs(self.b)
        return 2 * clenshaw(self._db, 2 * t - 1)


# --------------------------- fixed point, Chebyshev basis ---------------------


def residual_jacobian(x, alpha, pts):
    """x=(b_1..b_N); returns F (npts) and J (npts x N+1)."""
    cf = ChebFun([0] * (len(x) + 1))
    cf.set_from_free(x)
    b = cf.b
    db = dcheb_coeffs(b)
    N = len(x)
    npts = len(pts)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    a2 = alpha ** 2
    sgn = [(-1) ** k for k in range(N + 1)]
    for i, t in enumerate(pts):
        ta = t / a2
        U = clenshaw(b, 2 * ta - 1)
        W = U ** 2
        sW = 2 * W - 1
        GW = clenshaw(b, sW)
        GpW = 2 * clenshaw(db, sW)
        Gp_ta = 2 * clenshaw(db, 2 * ta - 1)
        F[i] = clenshaw(b, 2 * t - 1) + alpha * GW
        # d/db_j G(any point) = T_j(s) - (-1)^j (b0 locked); T_j by recurrence.
        # Nodes: t (outer G), W (composition G(W)), ta=t/a^2 (inner U, W=U^2).
        Tprev_t = mpf(1)          # T_0
        Tcur_t = 2 * t - 1        # T_1
        Tprev_W = mpf(1)
        Tcur_W = sW
        s_ta = 2 * ta - 1
        Tprev_ta = mpf(1)
        Tcur_ta = s_ta
        for j in range(1, N + 1):
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
        # d/dalpha = G(W) - 4 t U G'(W) G'(t/a^2) / alpha^2
        J[i, N] = GW - 4 * t * U * GpW * Gp_ta / a2
    return F, J


def solve_fixedpoint(N, x_init=None, alpha_init=None):
    pts = cheb_vals(N, N + 1)
    if x_init is None:
        x = [mpf(0)] * N
        x[0] = mpf('-0.75')  # g ~ 1 - 1.5 t = 0.25 - 0.75 T_1
        alpha = mpf('2.5')
    else:
        x = [mpf(v) for v in x_init] + [mpf(0)] * (N - len(x_init))
        alpha = mpf(alpha_init)
    for it in range(300):
        F, J = residual_jacobian(x, alpha, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < mpf('1e-48'):
            return x, alpha, it, nrm
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        dx = lu_solve(J, rhs)
        lam = mpf(1)
        improved = False
        for _halve in range(30):
            x_new = [x[j] + lam * dx[j] for j in range(N)]
            al_new = alpha + lam * dx[N]
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


def dense_residual(x, alpha, npts=4001):
    cf = ChebFun([0] * (len(x) + 1))
    cf.set_from_free(x)
    a2 = alpha ** 2
    worst = mpf(0)
    for t in cheb_vals(0, npts):
        U = cf.G(t / a2)
        r = cf.G(t) + alpha * cf.G(U ** 2)
        worst = max(worst, abs(r))
    return worst


# --------------------------- linearised DT (MY formula) -----------------------
# DT phi(t) = a*phi(1)*G(t) + 4 t phi(1) U(t) G'(W) G'(t/a^2)
#             - a phi(W) - 2 a U(t) G'(W) phi(t/a^2)
# (chain-rule form; NOT simplified through g'(x) = -g'(u) g'(x/a))


def DT_apply(cf, alpha, phival, t):
    """phival: callable t -> phi(t) (even perturbation, phi(0)=0 not required here)."""
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


def R_projected(fval, x):
    """R(f)(x) = -alpha_f f(f(x/alpha_f)), alpha_f = -1/f(1); fval: callable x->f(x)."""
    f1 = fval(1)
    af = -1 / f1
    return -af * fval(fval(x / af))


# --------------------------- main ---------------------------------------------

def main():
    log(TAG)
    log('mp.dps=%d, basis=T_k(2t-1) Chebyshev, lin.formula=chain-rule unsimplified' % mp.dps)

    ladder = [4, 8, 14, 20, 28, 36, 44]
    sols = {}
    x_cur, al_cur = None, None
    for N in ladder:
        try:
            x, al, it, res = solve_fixedpoint(N, x_cur, al_cur)
        except RuntimeError as exc:
            log('N=%d continuation failed (%s); crude restart' % (N, exc))
            x, al, it, res = solve_fixedpoint(N)
        sols[N] = (x, al)
        log('N=%2d: alpha=%s (it %d, res %.2e)' % (N, nstr(al, 50), it, float(res)))
        x_cur, al_cur = x, al

    Nfinal = 48
    x48, al48, it48, res48 = solve_fixedpoint(Nfinal, x_cur, al_cur)
    sols[Nfinal] = (x48, al48)
    x56, al56, it56, res56 = solve_fixedpoint(56, x48, al48)
    sols[56] = (x56, al56)
    x40, al40, it40, res40 = solve_fixedpoint(40, sols[36][0], sols[36][1])
    sols[40] = (x40, al40)
    log('N=40: alpha=%s' % nstr(al40, 50))
    log('N=48: alpha=%s' % nstr(al48, 50))
    log('N=56: alpha=%s' % nstr(al56, 50))
    log('N-convergence |a48-a56|=%.3e |a40-a48|=%.3e' % (float(abs(al48 - al56)), float(abs(al40 - al48))))

    cf48 = ChebFun([0] * (Nfinal + 1))
    cf48.set_from_free(x48)
    dr48 = dense_residual(x48, al48, 4001)
    dr56 = dense_residual(x56, al56, 4001)
    log('dense residual N=48: %.3e  N=56: %.3e' % (float(dr48), float(dr56)))

    g1 = cf48.G(mpf(1))
    log('G(1)=%s vs -1/alpha=%s  diff %.3e' % (nstr(g1, 40), nstr(-1 / al48, 40), float(abs(g1 + 1 / al48))))

    # anchors
    lit_alpha = mpf('2.502907875095892822283902873218215786381271376727149977336192056779235463179590')
    lit_delta = mpf('4.669201609102990671853203820466201617258185577475768632745651343004133555868548')
    agentA_alpha = mpf('2.50290787509589282228390287321821578638127137672714997756593')
    log('anchor |alpha-lit|=%.3e ; |alpha-agentA|=%.3e' % (float(abs(al48 - lit_alpha)), float(abs(al48 - agentA_alpha))))

    # ---- validate MY DT formula vs finite difference of R -----------------------
    eps = mpf(10) ** (-25)

    def phi_test(t):
        return sum(mpf(1) / (j + 3) * t ** j for j in range(1, 12))

    def gfun(x):
        return cf48.G(x ** 2)

    def f_eps(x):
        return cf48.G(x ** 2) + eps * phi_test(x ** 2)

    worst_fd = mpf(0)
    for t in cheb_vals(0, 31):
        x = sqrt(t)
        fd = (R_projected(f_eps, x) - R_projected(gfun, x)) / eps
        an = DT_apply(cf48, al48, phi_test, t)
        worst_fd = max(worst_fd, abs(fd - an))
    log('MY DT formula vs fd of R (eps=1e-25): max diff = %.3e' % float(worst_fd))

    # ---- coordinate mode: finite substitution x->x(1+eps) ON THE OPERATOR -------
    eps_c = mpf(10) ** (-20)
    lam_c = mpf(1) + eps_c

    def h_lam(x):
        return cf48.G((lam_c * x) ** 2) / lam_c

    psi_nodes = []
    for t in cheb_vals(0, 61):
        x = sqrt(t)
        psi_nodes.append((R_projected(h_lam, x) - gfun(x)) / eps_c)
    # psi_num(t) ~ tangent of the coordinate family (should be ~ xg'(x)-g(x) = -K)

    # ---- eigen collocation matrix in B_k(t) = T_k(2t-1) - (-1)^k -----------------
    Nmat = Nfinal
    epts = cheb_vals(0, Nmat)  # Nmat Gauss roots

    def Bfun(j, t):
        s = 2 * t - 1
        if j == 0:
            return mpf(0)
        Tprev, Tcur = mpf(1), s
        for k in range(1, j):
            Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
        return Tcur - ((-1) ** j)

    # interpolation matrix V[i,j] = B_j(t_i)
    V = matrix(Nmat, Nmat)
    for i, t in enumerate(epts):
        s = 2 * t - 1
        Tprev = mpf(1)
        Tcur = s
        for j in range(1, Nmat + 1):
            if j == 1:
                V[i, 0] = Tcur - (-1)
            else:
                Tnext = 2 * s * Tcur - Tprev
                V[i, j - 1] = Tnext - ((-1) ** j)
                Tprev, Tcur = Tcur, Tnext

    cols = []
    for j in range(1, Nmat + 1):
        vals = matrix(Nmat, 1)
        for i, t in enumerate(epts):
            def phival(tt, j=j):
                return Bfun(j, tt)
            vals[i] = DT_apply(cf48, al48, phival, t)
        coef = lu_solve(V, vals)
        cols.append([coef[i] for i in range(Nmat)])
    A = matrix(Nmat, Nmat)
    for jc in range(Nmat):
        for i in range(Nmat):
            A[i, jc] = cols[jc][i]
    log('built %dx%d DT collocation matrix in Chebyshev-B basis' % (Nmat, Nmat))

    # (a) full eig
    from mpmath import eig as mp_eig
    E, ER = mp_eig(A, right=True)
    order = sorted(range(Nmat), key=lambda i: -abs(E[i]))
    log('full eig(), top 12 by |lambda|:')
    for rank, idx in enumerate(order[:12]):
        log('  #%2d  %s' % (rank, nstr(E[idx], 28)))

    delta_eig = E[order[0]].real
    log('delta(eig) = %s ; |delta-lit|=%.3e' % (nstr(delta_eig, 50), float(abs(delta_eig - lit_delta))))

    # (b) power iteration + one-sided Wielandt deflation chain
    def matvec(M, v):
        n = len(v)
        w = matrix(n, 1)
        for i in range(n):
            s = mpf(0)
            for j in range(n):
                s += M[i, j] * v[j]
            w[i] = s
        return w

    def power_iter(M, n, iters=4000, tol=mpf('1e-40')):
        v = matrix(n, 1)
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
        D = matrix(n, n)
        vr = v1[r]
        for i in range(n):
            for j in range(n):
                D[i, j] = M[i, j] - (lam1 * v1[i] / vr) * (1 if j == r else 0)
        return D

    chain = []
    Mcur = A
    for level in range(6):
        lam, v, it = power_iter(Mcur, Nmat)
        # sign: use Rayleigh-like quotient with last iterate
        w = matvec(Mcur, v)
        nv = sum(v[i] * v[i] for i in range(Nmat))
        sign = (sum(v[i] * w[i] for i in range(Nmat)) / nv)
        chain.append((lam * (1 if sign > 0 else -1), v, it))
        log('power-iter level %d: |lambda|=%s signed=%s (%d iters)' % (
            level, nstr(lam, 24), nstr(lam * (1 if sign > 0 else -1), 24), it))
        Mcur = wielandt(Mcur, sign, v, Nmat)

    # compare chain vs eig
    log('deflation chain vs eig() cross-check:')
    for lv, (lamc, _v, _i) in enumerate(chain):
        best = min((E[i] for i in range(Nmat)), key=lambda z: abs(z - lamc))
        log('  level %d: chain %s vs eig %s diff %.3e' % (lv, nstr(lamc, 20), nstr(best, 20), float(abs(lamc - best))))

    def dot(u, w):
        return sum(u[i] * w[i] for i in range(Nmat))

    def norm(u):
        return sqrt(sum(abs(u[i]) ** 2 for i in range(Nmat)))

    # ---- coordinate-mode identification ------------------------------------------
    # THEORY (derived independently, see WORKLOG): the exact scaling family
    # h_lam(x)=g(lam x)/lam consists of fixed points of R only OFF the
    # normalisation sheet (h_lam(0)=1/lam != 1); its tangent psi=-K has psi(0)!=0,
    # i.e. lives off-sheet, where it has eigenvalue -alpha*g'(1) (computed below),
    # NOT 1/alpha^2. The on-sheet (phi(0)=0) tower 1/alpha^(2k) is a different
    # mechanism. On-sheet tests: (a) sheet projection of the finite-eps tangent,
    # (b) analytic tangents of the renormalised families g(lam x) [chi_a=2tG'(t)]
    # and g(x(1+eps x^2)) [chi_b=2t^2 G'(t)]; each -> Rayleigh quotient under M.
    # OFF-SHEET eigenvalue of the naive scaling generator:
    gp1 = cf48.Gp(mpf(1))          # G'(1)
    gpx1 = 2 * gp1                  # g'(1)
    offsheet_eig = -al48 * gpx1
    log('off-sheet scaling-generator eigenvalue -alpha*g\'(1) = %s (not an on-sheet eigenvalue)'
        % nstr(offsheet_eig, 30))

    def tocol(p):
        m = matrix(len(p), 1)
        for i, val in enumerate(p):
            m[i] = val
        return m

    def rayleigh(pvec_in, label):
        p = [v for v in pvec_in]
        w = matvec(A, tocol(p))
        np_ = norm(p)
        mu = sum(p[i] * w[i] for i in range(len(p))) / (np_ ** 2)
        r = max(abs(w[i] - mu * p[i]) for i in range(len(p))) / np_
        log('  Rayleigh(%s): mu*=%s  ||Mp-mu*p||/||p||=%.3e' % (label, nstr(mu, 18), float(r)))
        return mu, r, p

    # psi_num at 101 nodes incl extrapolation to t=0 for projection:
    psi_val = {}
    for t in cheb_vals(0, 61):
        x = sqrt(t)
        psi_val[t] = (R_projected(h_lam, x) - gfun(x)) / eps_c
    psi0 = psi_val[mpf(0)] if mpf(0) in psi_val else None
    # t=0 is not a Gauss node; evaluate psi_num at t=0 directly:
    psi0 = (R_projected(h_lam, mpf(0)) - gfun(mpf(0))) / eps_c
    psi_epts = matrix(Nmat, 1)
    for i, t in enumerate(epts):
        x = sqrt(t)
        psi_epts[i] = (R_projected(h_lam, x) - gfun(x)) / eps_c - psi0  # sheet projection
    pcoef = lu_solve(V, psi_epts)
    pvec = [pcoef[i] for i in range(Nmat)]
    log('finite-eps tangent (sheet-projected) tests:')
    # PREDICTION from my symmetry algebra: R(h_lam)=h_{phi(lam)}, phi'(1)=-alpha*g'(1)=alpha^2,
    # hence psi_num = [R(h)-g]/eps ~ phi'(1)*(-K) + O(eps). Testing this prediction against
    # the analytic K = G - 2tG' validates the algebra (formula used only as prediction target).
    log('  psi_num vs predicted phi\'(1)*(-K), phi\'(1)=alpha^2=%s:' % nstr(offsheet_eig, 16))
    for t in [mpf('0.05'), mpf('0.2'), mpf('0.5'), mpf('0.8')]:
        x = sqrt(t)
        pv = (R_projected(h_lam, x) - gfun(x)) / eps_c
        Kt = cf48.G(t) - 2 * t * cf48.Gp(t)
        log('    t=%.2f: psi_num/(phi1*-K) = %s' % (float(t), nstr(pv / (offsheet_eig * (-Kt)), 16)))
    mu_psi, r_psi, pvec = rayleigh(pvec, 'psi_num_tilde')

    # analytic tangents
    def chi_a(t):
        return 2 * t * cf48.Gp(t)

    def chi_b(t):
        return 2 * t * t * cf48.Gp(t)

    for nm, fn in [('chi_a=2tG\'', chi_a), ('chi_b=2t^2G\'', chi_b)]:
        vals = matrix(Nmat, 1)
        for i, t in enumerate(epts):
            vals[i] = fn(t)
        c = lu_solve(V, vals)
        rayleigh([c[i] for i in range(Nmat)], nm)

    # eigenvector factor checks (value-level independence of lambda2 from tower)
    v_rank1 = [ER[i, order[1]] for i in range(Nmat)]
    v_rank2 = [ER[i, order[2]] for i in range(Nmat)]
    v_rank4 = [ER[i, order[4]] for i in range(Nmat)]
    v_delta = [ER[i, order[0]] for i in range(Nmat)]
    for nm, vv in [('v_rank1', v_rank1), ('v_rank4', v_rank4), ('v_delta', v_delta)]:
        cs = dot(v_rank2, vv) / (norm(v_rank2) * norm(vv))
        log('cos(v_lambda2, %s) = %s' % (nm, nstr(cs, 12)))
    cos_r1 = dot(pvec, v_rank1) / (norm(pvec) * norm(v_rank1))
    cos_r2 = dot(pvec, v_rank2) / (norm(pvec) * norm(v_rank2))
    log('cos(psi_num_tilde, v_rank1[1/alpha^2]) = %s' % nstr(cos_r1, 12))
    log('cos(psi_num_tilde, v_rank2[lambda2])   = %s' % nstr(cos_r2, 12))

    # ---- candidate classification (my own list from MY alpha, delta) ------------
    l2 = E[order[2]].real
    l3 = E[order[3]].real
    log('lambda2 (rank2) = %s' % nstr(l2, 30))
    log('lambda3 (rank3) = %s' % nstr(l3, 30))
    cand = {
        '1/alpha': 1 / al48, '1/alpha^2': 1 / al48 ** 2, '1/alpha^3': 1 / al48 ** 3,
        '1/alpha^4': 1 / al48 ** 4, '1/alpha^6': 1 / al48 ** 6, '1/alpha^8': 1 / al48 ** 8,
        'alpha': al48, 'alpha^2': al48 ** 2, 'delta': delta_eig,
        '1/delta': 1 / delta_eig, 'sqrt(delta)': sqrt(delta_eig),
        'delta/alpha^2': delta_eig / al48 ** 2, 'delta/alpha^4': delta_eig / al48 ** 4,
    }
    log('--- candidate classification of ranks 1..8 ---')
    essentials = []
    for rank in range(1, 9):
        lv = E[order[rank]].real
        bn, bd = None, None
        for nm, cv in cand.items():
            dd = abs(abs(lv) - cv)
            if bd is None or dd < bd:
                bn, bd = nm, dd
        is_coord = bd < mpf('1e-30')
        log('  rank %d: %s  nearest=%s diff=%.3e -> %s' % (
            rank, nstr(lv, 24), bn, float(bd), 'COORDINATE' if is_coord else 'essential'))
        if not is_coord:
            essentials.append((rank, lv))

    # ---- continuous-operator eigenpair residuals (MY DT formula) ----------------
    def make_phi(coefvec):
        def phi(t):
            s = 2 * t - 1
            Tprev, Tcur = mpf(1), s
            acc = coefvec[0] * (Tcur - (-1))
            for k in range(2, Nmat + 1):
                Tprev, Tcur = Tcur, 2 * s * Tcur - Tprev
                acc = acc + coefvec[k - 1] * (Tcur - ((-1) ** k))
            return acc
        return phi

    for label, lamv, idx in [('delta', delta_eig, order[0]), ('lambda2', l2, order[2]), ('lambda3', l3, order[3])]:
        v = [ER[i, idx] for i in range(Nmat)]
        ssum = sum(v)
        v = [vv / ssum for vv in v]  # normalise phi(1)=1 -- B_k(1)=1-(-1)^k
        phi = make_phi(v)
        worst = mpf(0)
        for t in cheb_vals(0, 101):
            worst = max(worst, abs(DT_apply(cf48, al48, phi, t) - lamv * phi(t)))
        log('eigenpair %s continuous residual (101 pts): %.3e' % (label, float(worst)))

    # ---- N=40 cross-check of lambda2 --------------------------------------------
    cf40 = ChebFun([0] * 41)
    cf40.set_from_free(x40)
    epts40 = cheb_vals(0, 40)
    V40 = matrix(40, 40)
    for i, t in enumerate(epts40):
        s = 2 * t - 1
        Tprev, Tcur = mpf(1), s
        for j in range(1, 41):
            if j == 1:
                V40[i, 0] = Tcur + 1
            else:
                Tnext = 2 * s * Tcur - Tprev
                V40[i, j - 1] = Tnext - ((-1) ** j)
                Tprev, Tcur = Tcur, Tnext
    cols40 = []
    for j in range(1, 41):
        vals = matrix(40, 1)
        for i, t in enumerate(epts40):
            def phival(tt, j=j):
                return Bfun(j, tt)
            vals[i] = DT_apply(cf40, al40, phival, t)
        coef = lu_solve(V40, vals)
        cols40.append([coef[i] for i in range(40)])
    A40 = matrix(40, 40)
    for jc in range(40):
        for i in range(40):
            A40[i, jc] = cols40[jc][i]
    E40 = mp_eig(A40, right=False)
    o40 = sorted(E40, key=lambda z: -abs(z))
    l2_40 = min(o40[:8], key=lambda z: abs(z - l2))
    log('N=40 cross-check lambda2: %s vs N=48 %s (diff %.3e)' % (nstr(l2_40, 24), nstr(l2, 24), float(abs(l2_40 - l2))))

    # ---- comparison table with agent A -----------------------------------------
    agentA = {
        'alpha': agentA_alpha,
        'delta': mpf('4.669201609102990671853203820466201617258185577475768632745651343'),
        'lambda2': mpf('-0.123652712552687024335393652586'),
        'lambda3': mpf('-0.0573070210666818468843650322586'),
    }
    mine = {'alpha': al48, 'delta': delta_eig, 'lambda2': l2, 'lambda3': l3}
    log('--- comparison with agent A (theirs from p1_deg2.json) ---')
    comp = {}
    for k in agentA:
        d = abs(mine[k] - agentA[k])
        comp[k] = {'mine': nstr(mine[k], 45), 'agentA': nstr(agentA[k], 45), 'absdiff': float(d)}
        log('  %s: mine=%s agentA=%s |diff|=%.3e' % (k, nstr(mine[k], 40), nstr(agentA[k], 40), float(d)))

    # also power-iteration delta (independent of eig)
    lam0 = chain[0][0]
    log('power-iteration delta vs eig delta diff: %.3e' % float(abs(lam0 - delta_eig)))

    # spectrum table for report
    spec_top12 = [nstr(E[order[i]], 28) for i in range(12)]
    chain_vals = [nstr(c[0], 24) for c in chain]

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'mp_dps': mp.dps,
                 'basis': 'Chebyshev T_k(2t-1), b0 locked by G(0)=1',
                 'lin_formula': 'chain-rule unsimplified (4t phi(1) U Gp(W) Gp(t/a^2) form), fd-validated',
                 'not_rigorous': True, 'elapsed_sec': time.time() - T0},
        'alpha': nstr(al48, 55), 'delta': nstr(delta_eig, 55),
        'lambda2': nstr(l2, 40), 'lambda3': nstr(l3, 40),
        'alpha_ladder': {str(N): nstr(sols[N][1], 45) for N in sols},
        'checks': {
            'dense_residual_N48': float(dr48), 'dense_residual_N56': float(dr56),
            'alpha_N48_vs_N56': float(abs(al48 - al56)), 'alpha_N40_vs_N48': float(abs(al40 - al48)),
            'alpha_vs_lit_anchor': float(abs(al48 - lit_alpha)),
            'delta_vs_lit_anchor': float(abs(delta_eig - lit_delta)),
            'DT_formula_vs_fd': float(worst_fd),
            'G1_plus_invalpha': float(abs(g1 + 1 / al48)),
            'spectrum_top12': spec_top12,
            'deflation_chain': chain_vals,
            'offsheet_scaling_eigenvalue': nstr(offsheet_eig, 30),
            'rayleigh_psi_tilde': {'mu': nstr(mu_psi, 20), 'resid': float(r_psi)},
            'cos_psi_v_rank1': nstr(cos_r1, 12), 'cos_psi_v_rank2': nstr(cos_r2, 12),
            'cos_l2_v_rank1': nstr(dot(v_rank2, v_rank1) / (norm(v_rank2) * norm(v_rank1)), 12),
            'cos_l2_v_rank4': nstr(dot(v_rank2, v_rank4) / (norm(v_rank2) * norm(v_rank4)), 12),
            'lambda2_N40_crosscheck': nstr(l2_40, 28), 'lambda2_N40_diff': float(abs(l2_40 - l2)),
            'comparison_vs_agentA': comp,
            'poweriter_vs_eig_delta': float(abs(lam0 - delta_eig)),
        },
        # raw solution data for the interval (flint) stage:
        'x48_cheb_b1_to_b48': [nstr(v, 50) for v in x48],
        'v_delta_B': [nstr(ER[i, order[0]], 45) for i in range(Nmat)],
        'v_lambda2_B': [nstr(ER[i, order[2]], 45) for i in range(Nmat)],
    }
    with open('v1_cheb.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote v1_cheb.json')


if __name__ == '__main__':
    main()

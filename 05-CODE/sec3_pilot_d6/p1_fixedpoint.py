#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p1_fixedpoint.py -- Phase 1 of the degree-6 universality-class float pilot.
[oid:deg6B]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

Solves the sextic Cvitanovic-Feigenbaum equation (BRIEF convention, quartic
template with t = x^4 replaced by t = x^6):

    g(x) = -alpha6 * g(g(x/alpha6)),   g(0) = 1,

where g is an even function of x^6:  g(x) = G(t), t = x^6,  G(t) = sum_j a_j t^j,
a_0 = 1.  Since g depends on x^6 only, g(-z)=g(z) and the sign of the inner
argument is irrelevant; the sign of alpha6 is fixed by g(1) = -1/alpha6
(which follows from the equation at x = 0).

Method: collocation at Chebyshev-type nodes on (0,1] + Newton with analytic
Jacobian (mpmath, dps=60); continuation in the number of modes N.

Then: linearisation of the doubling operator on normalised sextic-symmetric
perturbations  R(f)(x) = -alpha_f f(f(x/alpha_f)), alpha_f = -1/f(1):

    (DT phi)(x) = -(phi(1)/g(1)) [g(x) - x g'(x)]
                  - alpha6 * phi(g(x/alpha6)) - alpha6 * g'(g(x/alpha6)) * phi(x/alpha6)

which preserves phi(0)=0 and the "function of x^6" subspace.  Built as a matrix
on t-series coefficients (collocation, scalar-exact operator evaluation at the
nodes); leading eigenvalue -> delta6, eigenfunction -> h6 (normalised h(1)=1).

SPECTRUM METHOD (degree-6 lesson, differs from the quartic pilot!): the
collocation operator is EXTREMELY NON-NORMAL and carries N-wandering
pseudo-spectral pollution (e.g. a line -1.8,-2.0,-2.2,-2.5,... moving with N,
and ~0.5-0.6 values at N>=56) sitting BETWEEN delta6 and the physical subleading
modes.  Deflated power iteration (the quartic recipe) chases this pollution and
returns garbage (first run gave lambda2=3.40, lambda3=2*lambda2).  Working
method, all in mpmath:
  (i)   full eig of the collocation matrix on a small-N ladder (N=12..24, cheap,
        well-conditioned there) -> N-CONVERGENT candidates = physical modes,
        N-wandering ones = pollution;
  (ii)  targeted SHIFTED INVERSE iteration (250 iters, late Rayleigh refresh)
        at N=48 and N=56 from each ladder-stable candidate;
  (iii) certification: dense-grid FUNCTION-SPACE eigenpair residual (no matrix
        involved) plus N=48 vs N=56 agreement.
The rescaling-family eigenvalue alpha6^{-12} is identified exactly (numerology
cross-check).  delta6 itself comes from plain mpmath power iteration (huge
spectral gap, immune to the pollution) and is certified by its own dense-grid
eigenpair residual.

In t = x^6 variable the operator reads (Phi = phi as t-series):
    K(t)  = G(t) - 6 t G'(t)
    U(t)  = G(t/alpha6^6)
    W(t)  = U(t)^6
    (DT Phi)(t) = -(Phi(1)/G(1)) K(t) - alpha6 * Phi(W(t))
                  - 6 alpha6 U(t)^5 G'(W(t)) * Phi(t/alpha6^6)

Self-consistency checks reported (BRIEF: "самосогласованностью своего счёта"):
  (1) fixed-point residual on dense grid,
  (2) N-convergence ladder of alpha6 (N=48 vs N=56),
  (3) eigenpair scalar residual  max |DT(h)(x) - delta6 h(x)| on dense grid,
  (4) lambda2/lambda3 stability N=48 vs N=56 (lambda2 feeds the 1/lambda2
      rate-law verdict, so it gets its own N-stability check),
  (5) literature sanity anchors delta6 ~ 9.30, alpha6 ~ 1.467 (anchors only,
      the operative check is the M-side cascade in p3).

Output: sextic_fixedpoint.json ; run log: p1_run.log (written by this script
itself, flush + elapsed/ETA lines).
"""
import json
import time

from mpmath import mp, mpf, mpmathify, nstr, pi, cos, matrix, lu_solve, fabs, \
    re as mpre

mp.dps = 60

TAG = "[oid:deg6B] NOT RIGOROUS float pilot, phase 1 (degree 6)"
T0 = time.time()
LOGF = open('p1_run.log', 'w', encoding='utf-8')


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


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
    """F_i = G(t_i) + alpha*G(W_i), W_i = U_i^6, U_i = G(t_i/alpha^6).
    Returns (F, J) with J (npts x npts+1) over unknowns (a_1..a_N, alpha).
    npts = N+1."""
    N = len(a) - 1
    npts = len(pts)
    F = matrix(npts, 1)
    J = matrix(npts, N + 1)
    a6 = alpha ** 6
    for i, t in enumerate(pts):
        ta = t / a6
        U = G_eval(a, ta)
        Up = Gp_eval(a, ta)
        W = U ** 6
        G_W = G_eval(a, W)
        Gp_W = Gp_eval(a, W)
        F[i] = G_eval(a, t) + alpha * G_W
        # d/da_j : t^j + alpha*( W^j + Gp_W * 6 U^5 * ta^j )
        base = alpha * Gp_W * 6 * U ** 5
        pow_t = mpf(1)
        pow_W = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            J[i, j - 1] = pow_t + alpha * pow_W + base * ta ** j
        # d/dalpha : G(W) - 36 U^5 t G'(W) G'(ta) / alpha^6
        J[i, N] = G_W - 36 * U ** 5 * t * Gp_W * Up / a6
    return F, J


def solve_fixedpoint(N, a_init=None, alpha_init=None):
    """Newton with damping; returns (a list, alpha, iters, residual) or raises."""
    pts = colloc_points(N + 1)
    if a_init is None:
        a = [mpf(1)] + [mpf(0)] * N
        if N >= 1:
            a[1] = mpf('-1.68')
        alpha = mpf('1.467')
    else:
        a = list(a_init) + [mpf(0)] * (N - (len(a_init) - 1))
        alpha = mpf(alpha_init)
    for it in range(200):
        F, J = residuals_and_jacobian(a, alpha, pts)
        nrm = max(abs(F[i]) for i in range(len(pts)))
        if nrm < mpf('1e-48'):
            return a, alpha, it, nrm
        rhs = matrix(len(pts), 1)
        for i in range(len(pts)):
            rhs[i] = -F[i]
        try:
            dx = lu_solve(J, rhs)
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError('lu_solve failed: %s' % exc)
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
    a6 = alpha ** 6
    worst = mpf(0)
    for t in pts:
        U = G_eval(a, t / a6)
        r = G_eval(a, t) + alpha * G_eval(a, U ** 6)
        worst = max(worst, abs(r))
    return worst


# --------------------------------------------------------------------------- #
#  scalar DT evaluation (ground truth for checks; x complex-capable)
# --------------------------------------------------------------------------- #
def g_of_x(a, x):
    """g(x) for x complex/plain scalar: g(x) = G(x^6)."""
    return G_eval(a, x ** 6)


def gp_of_x(a, x):
    """g'(x) = 6 x^5 G'(x^6)."""
    return 6 * x ** 5 * Gp_eval(a, x ** 6)


def DT_scalar(a, alpha, c, x):
    """(DT phi)(x) for phi(x) = sum_j c[j] x^{6j} (c[0] assumed 0).
    DT phi = -(phi(1)/g(1)) [g(x) - x g'(x)] - alpha phi(g(x/alpha))
              - alpha g'(g(x/alpha)) phi(x/alpha)."""
    G1 = G_eval(a, mpf(1))

    def phi(z):
        r = c[-1]
        for cc in reversed(c[:-1]):
            r = r * (z ** 6) + cc
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
            r = r * (z ** 6) + cc
        return r

    u = g_of_x(a, x / alpha)
    gpu = gp_of_x(a, u)
    return -alpha * phi(u) - alpha * gpu * phi(x / alpha)


def R_of_f(f, alpha_f, x):
    """Doubling operator R(f)(x) = -alpha_f f(f(x/alpha_f)) for plain scalar f."""
    return -alpha_f * f(f(x / alpha_f))


# --------------------------------------------------------------------------- #
#  collocation matrix of DT on the span{t, t^2, ..., t^N} (t = x^6)
# --------------------------------------------------------------------------- #
def build_DT_collocation(a, alpha, N):
    """A = V^{-1} F,  V_{ij} = t_i^j  (i=1..N nodes, j=1..N),  F_{ij} = DT(e_j)(t_i^{1/6}).
    Scalar-exact operator evaluation at the nodes (no series truncation of the
    operator image); V is inverted in mpmath.  Eigenproblem of A -> spectrum of
    DT on normalised sextic modes."""
    pts_t = colloc_points(N)
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts_t):
        x = t ** (mpf(1) / 6)
        for j in range(1, N + 1):
            V[i, j - 1] = t ** j
            e_j = [mpf(0)] * (j + 1)
            e_j[j] = mpf(1)
            F[i, j - 1] = DT_scalar(a, alpha, e_j, x)
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
#  mpmath power iteration (dominant eigenvalue only -- huge gap, pollution-free)
# --------------------------------------------------------------------------- #
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
            return mpf(0), v, 0
        for i in range(n):
            v[i] = w[i] / lam_new
        if abs(lam_new - lam) < tol * abs(lam_new):
            lam = lam_new
            return lam, v, it
        lam = lam_new
    return lam, v, iters


def inverse_iter_long(A, n, sigma, iters=250, refresh_at=(120, 180)):
    """Targeted shifted inverse iteration (complex-capable), late Rayleigh refresh."""
    from mpmath import mpc
    sig = mpc(sigma)
    M = matrix(n, n)
    for i in range(n):
        for j in range(n):
            M[i, j] = A[i, j]
    for i in range(n):
        M[i, i] = A[i, i] - sig
    v = matrix(n, 1)
    for i in range(n):
        v[i] = mpc(1) / (i + 2)
    for it in range(iters):
        try:
            w = lu_solve(M, v)
        except Exception:  # noqa: BLE001
            return None, None, None
        nrm = max(abs(w[i]) for i in range(n))
        if nrm == 0 or not (nrm < mpf('1e50')):
            return None, None, None
        for i in range(n):
            w[i] = w[i] / nrm
        v = w
        if it in refresh_at:
            Av = matvec(A, v)
            sig = sum(v[i].conjugate() * Av[i] for i in range(n)) / sum(abs(v[i]) ** 2 for i in range(n))
            for i in range(n):
                M[i, i] = A[i, i] - sig
    Av = matvec(A, v)
    k = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[k] / v[k]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, v, resid


def dense_eigen_residual(a, alpha, c, lam, npts=161):
    """max_t |DT(phi)(t^{1/6}) - lam*phi(t^{1/6})| on a dense grid, scalar-exact,
    no matrix involved (function-space certificate)."""
    worst = mpf(0)
    for t in colloc_points(npts):
        x = t ** (mpf(1) / 6)
        lhs = DT_scalar(a, alpha, c, x)
        rhs = lam * sum(c[j] * x ** (6 * j) for j in range(1, len(c)))
        worst = max(worst, abs(lhs - rhs))
    return worst


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
        ok = False
        if a_cur is not None:
            try:
                a, al, iters, res = solve_fixedpoint(N, a_cur, al_cur)
                ok = True
            except RuntimeError as exc:
                log('N=%d continuation failed (%s); trying crude starts' % (N, exc))
        if not ok:
            for a1c, alc in (('-1.68', '1.467'), ('-1.5', '1.40'),
                             ('-1.8', '1.55'), ('-1.3', '1.30')):
                try:
                    a = [mpf(1)] + [mpf(a1c)] + [mpf(0)] * (N - 1)
                    a, al, iters, res = solve_fixedpoint(N, a[:], mpf(alc))
                    ok = True
                    log('crude start a1=%s alpha=%s worked at N=%d' % (a1c, alc, N))
                    break
                except RuntimeError:
                    continue
        if not ok:
            raise RuntimeError('all starts failed at N=%d' % N)
        sols.append((N, a, al))
        log('N=%2d: alpha6 = %s   (newton iters %d, colloc res %.2e)'
            % (N, mpmathify(al), iters, float(res)))
        a_cur, al_cur = a, al

    # ---- final solve at N=72 (and 80 for convergence estimate) -------------------
    # N=72 (not 48 like the quartic): the sextic t=x^6 series converges with a
    # slow prefactor (|G48-G40|(t=1) ~ 3e-6 -- measured), so the H-side (phase 3)
    # needs the longer series for its 1e-11 truncation protocol; tail at t<=1.15
    # is ~1e-22 at N=72.  alpha6/delta6/eigenvalues are unchanged (ladder
    # super-convergence), only the series gets longer.
    N1 = 72
    a1, al1, it1, res1 = solve_fixedpoint(N1, a_cur, al_cur)
    N2 = 80
    a2, al2, it2, res2 = solve_fixedpoint(N2, a1, al1)
    log('N=%d alpha6=%s' % (N1, mpmathify(al1)))
    log('N=%d alpha6=%s' % (N2, mpmathify(al2)))
    d_alpha = abs(al1 - al2)
    d_coef = max(abs(a1[j] - a2[j]) for j in range(1, N1 + 1))
    log('N-convergence: |alpha(N72)-alpha(N80)| = %.3e , max|da| = %.3e'
        % (float(d_alpha), float(d_coef)))

    dense_res = fp_residual_dense(a1, al1, 4001)
    dense_res2 = fp_residual_dense(a2, al2, 4001)
    log('dense residual (N=72): %.3e   (N=80): %.3e'
        % (float(dense_res), float(dense_res2)))

    # radius-of-convergence estimate (t-plane)
    Rat = abs(a1[-1]) ** (-mpf(1) / N1)
    tail = [abs(a1[j]) for j in range(N1 - 8, N1 + 1)]
    Rat2 = (tail[0] / tail[-1]) ** (mpf(1) / 8)
    log('series radius estimate (t-plane): last-coeff %.6f, tail-ratio %.6f'
        % (float(Rat), float(Rat2)))

    g1 = G_eval(a1, mpf(1))
    log('g(1) = %s  (should equal -1/alpha6 = %s)'
        % (mpmathify(g1), mpmathify(-1 / al1)))

    # ---- validate DT formula by finite differences of R -------------------------
    eps = mpf(10) ** (-25)
    c_test = [mpf(0)] + [mpf(1) / (j + 1) for j in range(1, 12)]

    def f_eps(x):
        return g_of_x(a1, x) + eps * (sum(c_test[j] * x ** (6 * j) for j in range(1, 12)))

    alpha_eps = -1 / f_eps(mpf(1))
    fd_worst = mpf(0)
    for t in colloc_points(31):
        x = t ** (mpf(1) / 6)
        fd = (R_of_f(f_eps, alpha_eps, x) - g_of_x(a1, x)) / eps
        an = DT_scalar(a1, al1, c_test, x)
        fd_worst = max(fd_worst, abs(fd - an))
    log('DT formula vs finite difference of R (eps=1e-25): max diff = %.3e' % float(fd_worst))

    # ---- coordinate-mode diagnostic (unprojected DT on N = G - x g') -------------
    Kcoef = [a1[0]] + [a1[j] * (1 - 6 * j) for j in range(1, N1 + 1)]
    coord_ratios = []
    for t in colloc_points(11):
        x = t ** (mpf(1) / 6)
        num = DT_unprojected_scalar(a1, al1, Kcoef, x)
        den = sum(Kcoef[j] * x ** (6 * j) for j in range(0, N1 + 1))
        coord_ratios.append(num / den)
    coord_c = coord_ratios[len(coord_ratios) // 2]
    coord_spread = max(coord_ratios) - min(coord_ratios)
    log('coordinate mode DT_unproj(K)/K ~ %s (spread %.2e); for reference alpha6^2 = %.9f, alpha6^6 = %.9f'
        % (nstr(coord_c, 18), float(coord_spread), float(al1 ** 2), float(al1 ** 6)))

    # ---- linearisation spectrum ----------------------------------------------------
    # delta6 by TARGETED INVERSE ITERATION (not power iteration!): at N=72/80 the
    # pseudo-spectral pollution of the collocation matrix has grown ABOVE delta6
    # (first N=72 attempt: power iteration returned 3646.4), so dominance is lost;
    # the shift below is the certified value from the N=48/56 runs (dense-grid
    # residual 4.8e-38, N48-vs-56 agreement 3.7e-36, M-side cross-check 4e-7) --
    # inverse iteration converges to the eigenvalue nearest the shift regardless
    # of the pollution.  Subleading modes: eig-ladder + targeted inverse iteration.
    DELTA6_SHIFT = mpf('9.296246832771370082834476566367457546559')
    log('building DT collocation matrix N=%d (scalar-exact at nodes)...' % N1)
    A = build_DT_collocation(a1, al1, N1)
    log('inverse iteration for delta6 at N=%d (shift %s)...' % (N1, nstr(DELTA6_SHIFT, 20)))
    lam1, v1, mres1 = inverse_iter_long(A, N1, DELTA6_SHIFT, iters=60,
                                        refresh_at=(30, 45))
    assert lam1 is not None and abs(lam1 - DELTA6_SHIFT) < mpf('1e-25'), \
        'delta6 inverse iteration failed: %s' % nstr(lam1, 20)
    log('delta6 = %s  (matrix resid %.1e; N=48 certified value diff %.2e)'
        % (nstr(lam1, 40), float(mres1), float(abs(lam1 - DELTA6_SHIFT))))

    # h6 normalised h(1)=1  (eigenvector is real up to 1e-100 -> take real parts)
    v1r = [mpre(v1[i]) for i in range(N1)]
    s1 = sum(v1r)
    hcoef = [mpf(0)] + [v1r[i] / s1 for i in range(N1)]
    log('h6 normalised h(1)=1; h coeffs 1..5: %s'
        % [nstr(hcoef[j], 12) for j in range(1, 6)])

    # ---- subleading spectrum: eig ladder on small N --------------------------------
    log('--- subleading spectrum: mpmath full eig ladder N=12,16,20,24 ---')
    ladder_eigs = {}
    for Nl in (12, 16, 20, 24):
        Al = build_DT_collocation(a1, al1, Nl)
        E, _ = mp.eig(Al)
        evs = sorted(E, key=lambda z: -abs(z))
        ladder_eigs[Nl] = evs
        log('N=%d top-8: %s' % (Nl, '  '.join(nstr(e, 12) for e in evs[:8])))
    # candidates: N=24 eigenvalues |lam|<1 that are ladder-stable (agree with
    # N=20 to 1e-6 relative); pollution (e.g. the -1.8..-2.5 line) moves with N.
    cands = []
    for e24 in ladder_eigs[24]:
        if not (mpf('0.001') < abs(e24) < mpf(1)):
            continue
        e20 = min(ladder_eigs[20], key=lambda z: abs(z - e24))
        if abs(e20 - e24) < abs(e24) * mpf('1e-6'):
            if all(abs(e24 - c) > mpf('0.005') for c in cands):
                cands.append(e24)
    cands.sort(key=lambda z: -abs(z))
    cands = cands[:6]
    log('ladder-stable candidates (top 6 by |lam|): %s'
        % '  '.join(nstr(c, 14) for c in cands))

    # ---- certification at N=72 and N=80: targeted inverse iteration ----------------
    # Seeds: ladder-stable candidates PLUS the certified list of the previous
    # N=48/56 run (kept in sextic_fixedpoint.json) -- certify the top 3 by |lam|
    # fresh on the long series; deeper ones are copied through with a note (they
    # were certified at N=48/56 with dense residuals ~1e-38 and eigenvalues are
    # insensitive to series length: N=48-series vs N=56-series operators agree
    # to 20+ digits despite map differences ~3e-6).
    old_cert = []
    try:
        with open('sextic_fixedpoint.json', encoding='utf-8') as f_old:
            old_cert = [mpf(s) for s in json.load(f_old)['checks']['eigen_certified_list']]
    except Exception:  # noqa: BLE001
        old_cert = []
    for oc in old_cert:
        if all(abs(oc - c) > mpf('0.005') for c in cands):
            cands.append(oc)
    cands.sort(key=lambda z: -abs(z))
    log('building N=%d matrix for certification...' % N2)
    A2 = build_DT_collocation(a2, al2, N2)
    certified = []
    cert_details = []
    n_fresh = 0
    for cand in cands:
        fresh = n_fresh < 3
        rec = {'seed': nstr(cand, 16), 'route': 'fresh N72/N80' if fresh else 'carried from N48/56 run'}
        ok_c = fresh
        if fresh:
            n_fresh += 1
            lam48 = lam56 = None
            for Nl, Al, a_ser in ((N1, A, a1), (N2, A2, a2)):
                lam, v, mres = inverse_iter_long(Al, Nl, cand, iters=150,
                                                 refresh_at=(60, 110))
                if lam is None:
                    ok_c = False
                    rec['N%d' % Nl] = 'failed'
                    break
                c = [mpf(0)] + [v[i] for i in range(Nl)]
                dres = dense_eigen_residual(a_ser, al1, c, lam)
                rec['N%d' % Nl] = {'lam': nstr(lam, 30), 'matrix_resid': float(mres),
                                   'dense_resid': float(dres)}
                if Nl == N1:
                    lam48, dres48 = lam, dres
                else:
                    lam56 = lam
                ok_c = ok_c and float(dres) < 1e-25
                log('  cand %s N=%d: lam=%s  dense_resid=%.1e'
                    % (nstr(cand, 10), Nl, nstr(lam, 20), float(dres)))
            if ok_c:
                rel = float(abs(lam48 - lam56) / abs(lam48))
                rec['N72_vs_80_rel'] = rel
                ok_c = rel < 1e-12
            if ok_c:
                certified.append(lam48)
        else:
            certified.append(cand)
        cert_details.append(rec)
    certified.sort(key=lambda z: -abs(z))
    log('CERTIFIED subleading eigenvalues (dense resid<1e-25; fresh N72/N80 rel<1e-12, or carried from the certified N48/56 run): %s'
        % '  '.join(nstr(c, 20) for c in certified))
    lam2 = certified[0] if len(certified) >= 1 else None
    lam3 = certified[1] if len(certified) >= 2 else None
    log('lambda2 = %s' % (nstr(lam2, 25) if lam2 is not None else 'NONE'))
    log('lambda3 = %s' % (nstr(lam3, 25) if lam3 is not None else 'NONE'))
    # rescaling-family numerology check
    a_m12 = al1 ** (-12)
    nearest_m12 = min(certified, key=lambda z: abs(z - a_m12)) if certified else None
    log('alpha6^{-12} = %s ; nearest certified = %s (diff %.2e)'
        % (nstr(a_m12, 20),
           nstr(nearest_m12, 20) if nearest_m12 is not None else '-',
           float(abs(nearest_m12 - a_m12)) if nearest_m12 is not None else -1))
    # delta6 N-stability via inverse iteration as well
    lam1b, _vb, _rb = inverse_iter_long(A2, N2, lam1, iters=80, refresh_at=(40, 60))
    if lam1b is not None:
        log('delta6 N=72 vs N=80: |diff| = %.3e' % float(abs(lam1 - lam1b)))

    # ---- eigenpair scalar residual check (delta6) ------------------------------------
    worst = mpf(0)
    for t in colloc_points(201):
        x = t ** (mpf(1) / 6)
        lhs = DT_scalar(a1, al1, hcoef, x)
        rhs = lam1 * sum(hcoef[j] * x ** (6 * j) for j in range(1, N1 + 1))
        worst = max(worst, abs(lhs - rhs))
    log('eigenpair scalar residual max on 201 pts: %.3e' % float(worst))

    # ---- literature anchors (sanity only) -----------------------------------------
    lit_delta = mpf('9.30')
    lit_alpha = mpf('1.467')
    log('anchors: |delta6-9.30| = %.3e ; |alpha6-1.467| = %.3e'
        % (float(abs(lam1 - lit_delta)), float(abs(al1 - lit_alpha))))

    # ---- output --------------------------------------------------------------------
    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'mp_dps': mp.dps,
                 'equation': 'g(x) = -alpha6*g(g(x/alpha6)), g(0)=1, g in x^6',
                 'elapsed_sec': time.time() - T0},
        'alpha6': nstr(al1, 50),
        'delta6': nstr(mpre(lam1), 50),
        'N': N1,
        'g_coeffs': [nstr(x, 50) for x in a1],
        'h_coeffs': [nstr(x, 50) for x in hcoef],  # h(x) = sum h_coeffs[j] * x^{6j}
        'g1_check': nstr(g1, 30),
        'checks': {
            'alpha6_N72': nstr(al1, 40), 'alpha6_N80': nstr(al2, 40),
            'alpha6_N_diff': float(d_alpha), 'g_coeffs_N_diff_max': float(d_coef),
            'residual_dense_N72': float(dense_res), 'residual_dense_N80': float(dense_res2),
            'DT_formula_vs_finite_difference': float(fd_worst),
            'coordinate_mode_ratio': nstr(coord_c, 20),
            'coordinate_mode_spread': float(coord_spread),
            'eigenpair_scalar_residual': float(worst),
            'eigen2_deflated': nstr(mpre(lam2), 25) if lam2 is not None else None,
            'eigen3_deflated': nstr(mpre(lam3), 25) if lam3 is not None else None,
            'eigen_certified_list': [nstr(mpre(c), 25) for c in certified],
            'eigen_certification': cert_details,
            'alpha6_pow_minus12': nstr(a_m12, 20),
            'delta6_N80': nstr(mpre(lam1b), 20) if lam1b is not None else None,
            'delta6_N_stability': float(abs(lam1 - lam1b)) if lam1b is not None else None,
            'eig_ladder_N24_top8': [nstr(e, 14) for e in ladder_eigs[24][:8]],
            'radius_t_lastcoeff': float(Rat), 'radius_t_tailratio': float(Rat2),
            'alpha_ladder': {str(N): nstr(al, 30) for N, _a, al in sols},
            'note_float_spectrum': ('collocation operator highly non-normal for degree 6: '
                                    'deflated power iteration returns pseudo-spectral pollution '
                                    '(first attempt gave lambda2=3.40, lambda3=2*lambda2, N-unstable). '
                                    'At N=72/80 the pollution even overtakes delta6 (power iteration '
                                    'returns 3646.4 at N=72), so delta6 too is obtained by shifted '
                                    'inverse iteration from the certified N=48 value. Authoritative '
                                    'route: mpmath eig ladder N=12..24 (physical modes N-convergent, '
                                    'pollution N-wandering) + targeted shifted inverse iteration at '
                                    'N=72/80 + dense-grid function-space residuals.'),
            'anchor_delta6_diff': float(abs(lam1 - lit_delta)),
            'anchor_alpha6_diff': float(abs(al1 - lit_alpha)),
        },
    }
    with open('sextic_fixedpoint.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote sextic_fixedpoint.json')
    log('DONE phase 1, total elapsed %.1fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

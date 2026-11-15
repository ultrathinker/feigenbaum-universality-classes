#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectra_grid.py -- BRIEF-F: operator delta(d), lambda2(d), lambda3(d) on the grid
d = 2.0..3.6 step 0.1 (+ 4.0, 5.0, 6.0, 7.0 with bridges), then brackets for the
lambda2 zero and the crossover d* (delta(d*) = 1/|lambda2(d*)|) by bisection.
[oid:spectraF]  NOT RIGOROUS -- mpmath float exploration, no interval certificates.
The functional (dense-grid) residual reported next to every eigenvalue is a
quality measure, not a certificate.

Equation (degree-d Cvitanovic-Feigenbaum, REAL collocation in t = |x|^d):
    g(x) = -alpha*g(g(x/alpha)),  g(0)=1,  g(x)=G(t), t=|x|^d, G(0)=1,
    g(1) = -1/alpha.   With U = G(t/alpha^d), W = |U|^d:
    F(t) = G(t) + alpha*G(W) = 0.
For non-integer d everything stays real (no complexification), per BRIEF-F.

Linearisation of the normalised doubling operator R(f) = -alpha_f f(f(x/alpha_f)),
alpha_f = -1/f(1), on Phi(0)=0 perturbations (t-series c[0]=0):
    (DT Phi)(t) = -(Phi(1)/G(1))*K(t) - alpha*Phi(W)
                  - alpha*d*sU*|U|^(d-1)*G'(W)*Phi(t/alpha^d)
    K(t) = G(t) - d*t*G'(t),  sU = sign(U).

Spectrum route (degree-6 lessons of agent B, non-normal collocation operator):
  (i)   fixed point: Newton with damping + analytic Jacobian, continuation in N
        and in d (seed of the next d = solution of the previous one);
  (ii)  full mpmath eig ladder at small N (16 and 20) at selected d -> physical
        modes are ladder-stable, pseudo-spectral pollution is N-wandering;
  (iii) gauge tower test: eigenvalues equal to alpha^(-d*k) (k=1,2,...) are the
        rescaling-family tower (at d=2 this is 1/alpha^(2k) as in agent A; at
        d=6 agent B certified alpha6^(-12) exactly) -> rejected from the
        essential list; 1/alpha^(2k) is tested as well and the matched law is
        recorded;
  (iv)  certification of delta, lambda2, lambda3 (+ tower1) by targeted shifted
        inverse iteration at N=40 AND N=48 (own LU, factorised once per shift;
        Rayleigh refresh) + dense-grid function-space residual per eigenpair;
        N-drift = |lambda(N40) - lambda(N48)| goes into the tables;
  (v)   lambda2/lambda3 are tracked by CONTINUATION in d (shift = previous
        value), not by modulus rank: the order swaps near the lambda2 zero and
        when |lambda2| passes tower1.  cos(v2,v3) is recorded as a
        collision diagnostic; anchors d=2/4/6 are compared explicitly.

Outputs (this folder): spectra_grid.json (rewritten atomically after EVERY
point), progress.log (flush after every step).  CPU: one process, one thread.
"""
import json
import os
import time

from mpmath import mp, mpf, mpc, pi, cos, nstr
from mpmath import eig as mp_eig

mp.dps = 50

TAG = '[oid:spectraF] NOT RIGOROUS -- mpmath float exploration, no interval certificates'
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'spectra_grid.json')
PROG = os.path.join(HERE, 'progress.log')

N1 = 40          # first truncation order (series length N+1, matrix N)
N2 = 48          # second truncation order
LADDER_EIG_N = (16, 20)
NEWTON_TOL = mpf(10) ** (-42)
ITERS_INV = 60
REFRESH_AT = (25, 45)
ITERS_INV_LONG = 150
REFRESH_LONG = (60, 110, 140)
CERT_ABS = 1e-18     # dense residual certification threshold (absolute)


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


# --------------------------------------------------------------------------- #
#  series helpers (t-series, a[0] = 1) and collocation nodes
# --------------------------------------------------------------------------- #
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


# --------------------------------------------------------------------------- #
#  dense-list LU with reuse (partial pivoting, Doolittle)
# --------------------------------------------------------------------------- #
def lu_decomp(A):
    n = len(A)
    LU = [row[:] for row in A]
    piv = list(range(n))
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(LU[i][k]))
        if p != k:
            LU[k], LU[p] = LU[p], LU[k]
            piv[k], piv[p] = piv[p], piv[k]
        pivval = LU[k][k]
        if pivval == 0:
            raise ZeroDivisionError('singular LU pivot')
        rk = LU[k]
        for i in range(k + 1, n):
            m = LU[i][k] / pivval
            LU[i][k] = m
            if m != 0:
                ri = LU[i]
                for j in range(k + 1, n):
                    ri[j] -= m * rk[j]
    return LU, piv


def lu_solve_f(LU, piv, b):
    n = len(b)
    x = [b[piv[i]] for i in range(n)]
    for i in range(1, n):
        s = x[i]
        ri = LU[i]
        for j in range(i):
            s -= ri[j] * x[j]
        x[i] = s
    for i in range(n - 1, -1, -1):
        s = x[i]
        ri = LU[i]
        for j in range(i + 1, n):
            s -= ri[j] * x[j]
        x[i] = s / ri[i]
    return x


# --------------------------------------------------------------------------- #
#  fixed point for general real degree d (real collocation, t = |x|^d)
# --------------------------------------------------------------------------- #
def residuals_and_jacobian(a, alpha, d, pts):
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


def newton_fp(d, N, a0, al0, max_it=60):
    pts = colloc_points(N + 1)
    a = list(a0)
    a += [mpf(0)] * (N + 1 - len(a))
    alpha = mpf(al0)
    for it in range(max_it):
        F, J = residuals_and_jacobian(a, alpha, d, pts)
        nrm = max(abs(x) for x in F)
        if nrm < NEWTON_TOL:
            return a, alpha, it, nrm
        rhs = [-x for x in F]
        try:
            LU, piv = lu_decomp(J)
            dx = lu_solve_f(LU, piv, rhs)
        except ZeroDivisionError as exc:
            raise RuntimeError('lu failed: %s' % exc)
        lam = mpf(1)
        for _halve in range(30):
            a_new = [a[0]] + [a[j] + lam * dx[j - 1] for j in range(1, N + 1)]
            al_new = alpha + lam * dx[N]
            try:
                F2, _ = residuals_and_jacobian(a_new, al_new, d, pts)
                nrm2 = max(abs(x) for x in F2)
            except Exception:
                nrm2 = mpf('inf')
            if nrm2 < nrm:
                a, alpha, nrm = a_new, al_new, nrm2
                break
            lam /= 2
        else:
            raise RuntimeError('newton stalled, res %s' % nstr(nrm, 5))
    raise RuntimeError('newton no convergence, res %s' % nstr(nrm, 5))


def _fp_ok(a, al, d, thresh=mpf('1e-25'), npts=201):
    """Reject spurious collocation solutions: dense-grid residual check."""
    try:
        return fp_residual_dense(a, al, d, npts) < thresh
    except Exception:
        return False


def solve_fp(d, N, a_init=None, al_init=None):
    """Continuation solver: direct warm start, else N-ladder from warm/crude.
    Every candidate is dense-grid validated (spurious solutions of the projected
    equation at larger N have dense residuals ~1e-17, true ones <=1e-35)."""
    if a_init is not None and al_init is not None:
        try:
            a, al, it, r = newton_fp(d, N, a_init, al_init)
            if _fp_ok(a, al, d):
                return a, al, it, r
        except RuntimeError:
            pass
    if a_init is not None:
        a_cur, al_cur = a_init, al_init
        for Nl in (10, 18, 28, N):
            try:
                a_cur, al_cur, _it, _r = newton_fp(d, Nl,
                                                   a_cur[:Nl + 1] + [mpf(0)] * max(0, Nl + 1 - len(a_cur))
                                                   if len(a_cur) < Nl + 1 else a_cur[:Nl + 1],
                                                   al_cur)
            except RuntimeError:
                a_cur, al_cur = None, None
                break
        if a_cur is not None and _fp_ok(a_cur, al_cur, d):
            try:
                a, al, it, r = newton_fp(d, N, a_cur, al_cur)
                if _fp_ok(a, al, d):
                    return a, al, it, r
            except RuntimeError:
                pass
    for a1c, alc in (('−1.5', '2.5'), ('-1.6', '2.05'), ('-1.55', '1.83'),
                     ('-1.63', '1.56'), ('-1.68', '1.47'), ('-1.7', '1.42'),
                     ('-1.45', '2.25'), ('-1.75', '1.65')):
        a1c = mpf(a1c.replace('−', '-'))
        a_cur = [mpf(1), a1c] + [mpf(0)] * 10
        al_cur = mpf(alc)
        try:
            for Nl in (4, 6, 10, 16, 24, 32, N):
                a_cur, al_cur, _it, _r = newton_fp(d, Nl, a_cur, al_cur)
            if _fp_ok(a_cur, al_cur, d):
                return a_cur, al_cur, _it, _r
        except RuntimeError:
            continue
    raise RuntimeError('fixed point failed at d=%s N=%d' % (nstr(d, 6), N))


def fp_residual_dense(a, alpha, d, npts=801):
    ad = alpha ** d
    worst = mpf(0)
    for t in colloc_points(npts):
        U = G_eval(a, t / ad)
        r = G_eval(a, t) + alpha * G_eval(a, abs(U) ** d)
        worst = max(worst, abs(r))
    return worst


# --------------------------------------------------------------------------- #
#  linearisation DT in t-space (scalar-exact ground truth + matrix builder)
# --------------------------------------------------------------------------- #
def DT_scalar_t(a, alpha, d, c, t):
    """(DT Phi)(t), Phi t-series with c[0]=0; includes the gauge term."""
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    ta = t / ad
    U = G_eval(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    Kt = G_eval(a, t) - d * t * Gp_eval(a, t)
    PhiW = _horner0(c, W)
    Phita = _horner0(c, ta)
    Phi1 = sum(c)
    return (-(Phi1 / G1) * Kt - alpha * PhiW
            - alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W) * Phita)


def _horner0(c, t):
    r = c[-1]
    for cc in reversed(c[:-1]):
        r = r * t + cc
    return r


def DT_unproj_t(a, alpha, d, c, t):
    """Bare linearisation (no normalisation projection) on t-series c (c[0] may
    be nonzero) -- coordinate-mode diagnostic."""
    ad = alpha ** d
    ta = t / ad
    U = G_eval(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    return (-alpha * _horner0(c, W)
            - alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W) * _horner0(c, ta))


def build_DT_matrix(a, alpha, d, N):
    """A = V^{-1} F on span{t,...,t^N}; V factored once, all columns reused."""
    pts = colloc_points(N)
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    V = []
    F = []
    for t in pts:
        ta = t / ad
        U = G_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        W = aU ** d
        Kt = G_eval(a, t) - d * t * Gp_eval(a, t)
        gauge = -(1 / G1) * Kt
        c3 = -alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W)
        vrow = [mpf(0)] * N
        frow = [mpf(0)] * N
        pow_t = mpf(1)
        pow_W = mpf(1)
        pow_ta = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t
            pow_W *= W
            pow_ta *= ta
            vrow[j - 1] = pow_t
            frow[j - 1] = gauge - alpha * pow_W + c3 * pow_ta
        V.append(vrow)
        F.append(frow)
    LU, piv = lu_decomp(V)
    Acols = []
    for jcol in range(N):
        b = [F[i][jcol] for i in range(N)]
        Acols.append(lu_solve_f(LU, piv, b))
    A = [[Acols[j][i] for j in range(N)] for i in range(N)]
    return A


def matvec(A, v):
    n = len(v)
    w = []
    for i in range(n):
        s = v[0] * 0
        Ai = A[i]
        for j in range(n):
            s += Ai[j] * v[j]
        w.append(s)
    return w


def inverse_iter(A, n, sigma, iters=ITERS_INV, refresh_at=REFRESH_AT):
    """Shifted inverse iteration, real arithmetic; LU per shift, reused."""
    sig = mpf(sigma)
    M = [[A[i][j] - (sig if i == j else 0) for j in range(n)] for i in range(n)]
    v = [mpf(1) / (i + 2) for i in range(n)]
    LU = piv = None
    try:
        LU, piv = lu_decomp(M)
    except ZeroDivisionError:
        sig = sig * (1 + mpf('1e-20'))
        M = [[A[i][j] - (sig if i == j else 0) for j in range(n)] for i in range(n)]
        LU, piv = lu_decomp(M)
    for it in range(iters):
        try:
            w = lu_solve_f(LU, piv, v)
        except ZeroDivisionError:
            return None, None, None
        nrm = max(abs(x) for x in w)
        if nrm == 0 or not (nrm < mpf(10) ** 52):
            return None, None, None
        v = [x / nrm for x in w]
        if it in refresh_at:
            Av = matvec(A, v)
            sig = sum(v[i] * Av[i] for i in range(n)) / sum(v[i] * v[i] for i in range(n))
            M = [[A[i][j] - (sig if i == j else 0) for j in range(n)] for i in range(n)]
            try:
                LU, piv = lu_decomp(M)
            except ZeroDivisionError:
                pass
    Av = matvec(A, v)
    k = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[k] / v[k]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, v, resid


def dense_eigen_residual(a, alpha, d, c, lam, npts=161):
    worst = mpf(0)
    phmax = mpf(0)
    for t in colloc_points(npts):
        lhs = DT_scalar_t(a, alpha, d, c, t)
        rhs = lam * _horner0(c, t)
        worst = max(worst, abs(lhs - rhs))
        phmax = max(phmax, abs(_horner0(c, t)))
    return worst, phmax


# --------------------------------------------------------------------------- #
#  per-point computation
# --------------------------------------------------------------------------- #
ANCHORS = {
    '2.0': {'alpha': '2.5029078750958928222839028732182157863812713767271',
            'delta': '4.669201609102990671853203820466201617258185577464',
            'lam2': '-0.123652712552687024335393652586',
            'lam3': '-0.0573070210666818468843650322586'},
    '4.0': {'alpha': '1.6903029714052448533437801503241613482282780596328',
            'delta': '7.2846862170733433643089305679955530694780466189573',
            'lam2': '0.291838408984'},
    '6.0': {'alpha': '1.4677424503199009444538343151089737463687971293968',
            'delta': '9.2962468327713700828344765663674575503',
            'lam2': '0.3710714676819301324709834',
            'lam3': '-0.3419683225478021050998591'},
}


def tower_members(alpha, d, kmax=8):
    return [alpha ** (-d * k) for k in range(1, kmax + 1)]


def eig_ladder_enumerate(a, alpha, d):
    """Full mp.eig at LADDER_EIG_N; ladder-stable candidates -> dict."""
    from mpmath import matrix
    sets = {}
    for Nl in LADDER_EIG_N:
        Al = build_DT_matrix(a, alpha, d, Nl)
        E = mp_eig(matrix(Al), right=False)
        sets[Nl] = sorted(E, key=lambda z: -abs(z))
    n0, n1 = LADDER_EIG_N
    stable = []
    for e1 in sets[n1]:
        if abs(e1) < mpf('0.005'):
            continue
        e0 = min(sets[n0], key=lambda z: abs(z - e1))
        if abs(e0 - e1) < abs(e1) * mpf('1e-6'):
            if all(abs(e1 - s) > mpf('0.003') for s in stable):
                stable.append(e1)
    stable.sort(key=lambda z: -abs(z))
    return sets, stable


def certify_target(A40, A48, a40, al40, a48, al48, d, sigma):
    """Inverse iteration at both N; returns dict or None."""
    out = {'seed': nstr(sigma, 14)}
    lam40 = lam48 = None
    for tag, A, a, al, N in (('N40', A40, a40, al40, N1), ('N48', A48, a48, al48, N2)):
        lam, v, mres = inverse_iter(A, N, sigma)
        if lam is None or mres is None or float(mres) > 1e-20:
            lam, v, mres = inverse_iter(A, N, sigma, iters=ITERS_INV_LONG,
                                        refresh_at=REFRESH_LONG)
        if lam is None:
            out[tag] = 'FAILED'
            continue
        c = [mpf(0)] + list(v)
        ssum = sum(v)
        hcoef = [mpf(0)] + [x / ssum for x in v]
        dres, phmax = dense_eigen_residual(a, al, d, hcoef, lam)
        out[tag] = {'lam': nstr(lam, 30), 'lam_f': float(lam),
                    'mat_resid': float(mres), 'dense_resid': float(dres),
                    'dense_resid_rel': float(dres / max(phmax, mpf('1e-30')))}
        out[tag + '_vec'] = [nstr(x, 16) for x in v]
        if tag == 'N40':
            lam40 = lam
        else:
            lam48 = lam
    if lam40 is not None and lam48 is not None:
        out['drift'] = float(abs(lam40 - lam48))
        out['drift_rel'] = float(abs(lam40 - lam48) / max(abs(lam48), mpf('1e-12')))
        out['lam_f'] = float(lam48)
        out['lam'] = nstr(lam48, 30)
    return out


def compute_point(ds, seeds, full_eig=False):
    d = mpf(ds)
    t_start = time.time()
    rec = {'d': ds, 'role': seeds.get('role', 'grid')}

    # ---- fixed points at both truncations -------------------------------- #
    if seeds.get('a48') is not None:
        a40, al40, it40, _ = solve_fp(d, N1, seeds['a48'], seeds['al48'])
        a48, al48, it48, _ = solve_fp(d, N2, seeds['a48'], seeds['al48'])
    else:
        a48, al48, it48, _ = solve_fp(d, N2)
        a40, al40, it40, _ = solve_fp(d, N1, a48, al48)
    rec['alpha40'] = nstr(al40, 30)
    rec['alpha48'] = nstr(al48, 30)
    rec['alpha_drift'] = float(abs(al40 - al48))
    rec['newton_iters'] = [it40, it48]
    rec['fp_dense_resid_N48'] = float(fp_residual_dense(a48, al48, d))
    rec['fp_dense_resid_N40'] = float(fp_residual_dense(a40, al40, d))
    rec['g1_check'] = nstr(G_eval(a48, mpf(1)), 25)
    plog('d=%s: alpha48=%s (drift40-48 %.1e, dense res %.1e)'
         % (ds, nstr(al48, 18), rec['alpha_drift'], rec['fp_dense_resid_N48']))

    # ---- DT matrices ------------------------------------------------------ #
    A40 = build_DT_matrix(a40, al40, d, N1)
    A48 = build_DT_matrix(a48, al48, d, N2)

    # ---- full eig ladder (selected d) -------------------------------------- #
    if full_eig:
        sets, stable = eig_ladder_enumerate(a48, al48, d)
        tow = tower_members(al48, d)
        tow2 = [al48 ** (-2 * k) for k in range(1, 13)]
        enum = []
        for e in stable[:14]:
            if abs(e.imag) > abs(e) * mpf('1e-12'):
                enum.append({'lam': nstr(e, 20), 'kind': 'complex'})
                continue
            ev = e.real
            kind = 'essential'
            m_tow = min(abs(ev - t) / t for t in tow if t > mpf('0.001'))
            m_tow2 = min(abs(ev - t) / t for t in tow2 if t > mpf('0.001'))
            if m_tow < mpf('1e-6'):
                kind = 'tower_alpha^-dk'
            elif m_tow2 < mpf('1e-6'):
                kind = 'tower_alpha^-2k'
            enum.append({'lam': nstr(ev, 20), 'lam_f': float(ev),
                         'kind': kind,
                         'tower_match_rel': float(min(m_tow, m_tow2))})
        rec['eig_enum_N20'] = enum
        rec['eig_ladder_N16_top6'] = [nstr(e, 16) for e in sets[16][:6]]
        plog('d=%s: eig ladder enumeration (top 6): %s'
             % (ds, ' | '.join('%s:%s' % (x['lam'][:14], x['kind'][:14])
                               for x in enum[:6])))

    # ---- certify delta, lam2, lam3, tower1 --------------------------------- #
    targets = {}
    an0 = ANCHORS.get(ds, {})
    targets['delta'] = seeds.get('delta', mpf(an0.get('delta', '5.0')))
    targets['lam2'] = seeds.get('lam2', mpf(an0.get('lam2', '-0.1')))
    targets['lam3'] = seeds.get('lam3', mpf(an0.get('lam3', '-0.2')))
    tower1 = al48 ** (-d)
    targets['tower1'] = tower1
    rec['tower1_alpha_pow'] = nstr(tower1, 25)

    cert = {}
    for name in ('delta', 'lam2', 'lam3', 'tower1'):
        cert[name] = certify_target(A40, A48, a40, al40, a48, al48, d, targets[name])

    # tower classification on certified value
    for name in ('delta', 'lam2', 'lam3'):
        c = cert[name]
        if isinstance(c.get('N48'), dict):
            lam = mpf(c['N48']['lam'])
            m = min([abs(lam - t) / t for t in tower_members(al48, d, 10)
                     if t > mpf('0.0005')] or [mpf(1)])
            c['tower_match_rel'] = float(m)
            c['is_tower'] = bool(m < 1e-9)
    if isinstance(cert['tower1'].get('N48'), dict):
        cert['tower1']['match_rel'] = float(abs(mpf(cert['tower1']['N48']['lam']) - tower1) / tower1)
    rec['delta'] = cert['delta']
    rec['lam2'] = cert['lam2']
    rec['lam3'] = cert['lam3']
    rec['tower1'] = cert['tower1']

    # ---- coordinate-mode diagnostic ---------------------------------------- #
    Kcoef = [a48[0]] + [a48[j] * (1 - d * j) for j in range(1, N2 + 1)]
    ratios = []
    for t in colloc_points(11):
        num = DT_unproj_t(a48, al48, d, Kcoef, t)
        den = _horner0(Kcoef, t)
        ratios.append(num / den)
    rec['coord_mode_ratio'] = nstr(sorted(ratios)[5], 18)
    rec['coord_mode_spread'] = float(max(ratios) - min(ratios))

    # ---- eigenvector continuity + collision diagnostic ---------------------- #
    def vec_of(crt):
        if isinstance(crt.get('N40_vec'), list):
            return [mpf(s) for s in crt['N40_vec']]
        return None

    v2, v3 = vec_of(cert['lam2']), vec_of(cert['lam3'])
    if v2 is not None and v3 is not None:
        dot = sum(v2[i] * v3[i] for i in range(N1))
        n2 = sum(x * x for x in v2) ** mpf('0.5')
        n3 = sum(x * x for x in v3) ** mpf('0.5')
        rec['cos_v2_v3'] = float(abs(dot / (n2 * n3)))
    pv2 = seeds.get('v2_40')
    if v2 is not None and pv2 is not None:
        dot = sum(v2[i] * pv2[i] for i in range(N1))
        n2 = sum(x * x for x in v2) ** mpf('0.5')
        npv = sum(x * x for x in pv2) ** mpf('0.5')
        rec['cos_v2_prev'] = float(abs(dot / (n2 * npv)))

    # ---- anchor check ------------------------------------------------------- #
    if ds in ANCHORS:
        an = ANCHORS[ds]
        chk = {}
        for key, certkey in (('alpha', None), ('delta', 'delta'),
                             ('lam2', 'lam2'), ('lam3', 'lam3')):
            if key not in an:
                continue
            ref = mpf(an[key])
            if key == 'alpha':
                mine = al48
            else:
                cd = cert[certkey]
                mine = mpf(cd['N48']['lam']) if isinstance(cd.get('N48'), dict) else None
            if mine is not None:
                chk[key + '_diff'] = float(abs(mine - ref))
        rec['anchor_check'] = chk

    rec['elapsed_s'] = round(time.time() - t_start, 1)

    # ---- seeds out ----------------------------------------------------------- #
    def lam48_of(name):
        cd = cert[name]
        return mpf(cd['N48']['lam']) if isinstance(cd.get('N48'), dict) else targets[name]

    seeds_out = {'a48': a48, 'al48': al48, 'a40': a40, 'al40': al40,
                 'delta': lam48_of('delta'), 'lam2': lam48_of('lam2'),
                 'lam3': lam48_of('lam3'), 'v2_40': v2}
    return rec, seeds_out


# --------------------------------------------------------------------------- #
#  grid driver + bisections
# --------------------------------------------------------------------------- #
GRID = [('%0.1f' % (2.0 + 0.1 * k)) for k in range(17)] + \
       ['3.7', '3.8', '3.9', '4.0', '4.5', '5.0', '5.5', '6.0', '6.5', '7.0']
FULL_EIG_AT = {'2.0', '2.5', '3.0', '3.6', '4.0', '5.0', '6.0', '7.0'}
ROLE_BRIDGE = {'3.7', '3.8', '3.9', '4.5', '5.5', '6.5'}


def load_points():
    if os.path.exists(OUTJ):
        with open(OUTJ, encoding='utf-8') as f:
            data = json.load(f)
        return data
    return {'meta': {'tag': TAG, 'created': time.strftime('%Y-%m-%d %H:%M:%S'),
                     'mp_dps': mp.dps, 'N_pair': [N1, N2],
                     'grid': GRID, 'equation': 'g(x)=-alpha*g(g(x/alpha)), t=|x|^d, real collocation'},
            'points': [], 'brackets': {}}


def save_points(data):
    tmp = OUTJ + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, OUTJ)


def l2_of(rec):
    return rec['lam2']['lam_f']


def dl_of(rec):
    return rec['delta']['lam_f']


def crossover_f(rec):
    a = abs(l2_of(rec))
    return dl_of(rec) * a - 1.0


def nearest_seeds(data, dval):
    best = None
    for rec in data['points']:
        rd = abs(float(rec['d']) - dval)
        if best is None or rd < best[0]:
            best = (rd, rec)
    rec = best[1]
    a48 = [mpf(s) for s in rec['_a48']] if '_a48' in rec else None
    seeds = {'a48': a48, 'al48': mpf(rec['alpha48']),
             'delta': mpf(rec['delta']['lam']), 'lam2': mpf(rec['lam2']['lam']),
             'lam3': mpf(rec['lam3']['lam'])}
    return seeds


def run_grid(data):
    done = {rec['d'] for rec in data['points']}
    seeds = None
    plog('GRID START: %d points, N=%d/%d, dps=%d' % (len(GRID), N1, N2, mp.dps))
    times = []
    for i, ds in enumerate(GRID):
        if ds in done:
            plog('d=%s: already done, skip' % ds)
            continue
        role = 'bridge' if ds in ROLE_BRIDGE else 'grid'
        sd = dict(seeds or {})
        sd['role'] = role
        rec, seeds = compute_point(ds, sd, full_eig=(ds in FULL_EIG_AT))
        rec['_a48'] = [nstr(x, 40) for x in seeds['a48']]
        data['points'].append(rec)
        data['points'].sort(key=lambda r: float(r['d']))
        save_points(data)
        times.append(rec['elapsed_s'])
        eta = sum(times[-5:]) / len(times[-5:]) * (len(GRID) - i - 1)
        plog('d=%s DONE (%.0fs): delta=%s (drift %.1e, res %.1e)  '
             'lam2=%s (drift %.1e, res %.1e)  lam3=%s (drift %.1e, res %.1e)  '
             'tower1 match %.1e  cos(v2,v3)=%s  [ETA grid %.0fs]'
             % (ds, rec['elapsed_s'],
                rec['delta']['lam'][:16] if 'lam' in rec['delta'] else 'FAIL',
                rec['delta'].get('drift', -1), rec['delta']['N48'].get('dense_resid', -1) if isinstance(rec['delta'].get('N48'), dict) else -1,
                rec['lam2']['lam'][:16] if 'lam' in rec['lam2'] else 'FAIL',
                rec['lam2'].get('drift', -1), rec['lam2']['N48'].get('dense_resid', -1) if isinstance(rec['lam2'].get('N48'), dict) else -1,
                rec['lam3']['lam'][:16] if 'lam' in rec['lam3'] else 'FAIL',
                rec['lam3'].get('drift', -1), rec['lam3']['N48'].get('dense_resid', -1) if isinstance(rec['lam3'].get('N48'), dict) else -1,
                rec['tower1'].get('match_rel', -1),
                ('%.3f' % rec['cos_v2_v3']) if 'cos_v2_v3' in rec else '-',
                eta))


def bisect(data, f_rec, name, target_w=2e-3, max_steps=9):
    """Generic bisection over d using rec-level scalar functions."""
    pts = sorted(data['points'], key=lambda r: float(r['d']))
    bracket = None
    for r0, r1 in zip(pts, pts[1:]):
        f0, f1 = f_rec(r0), f_rec(r1)
        if f0 == 0:
            bracket = (r0, r0)
            break
        if f0 * f1 < 0:
            bracket = (r0, r1)
            break
    if bracket is None:
        plog('BISECT %s: no sign change found on grid' % name)
        return None
    lo, hi = bracket
    plog('BISECT %s: bracket [%.4f, %.4f] f=[%.6f, %.6f]'
         % (name, float(lo['d']), float(hi['d']), f_rec(lo), f_rec(hi)))
    steps = 0
    while float(hi['d']) - float(lo['d']) > target_w and steps < max_steps:
        dm = (float(lo['d']) + float(hi['d'])) / 2
        ds = ('%.6f' % dm).rstrip('0')
        if ds in {r['d'] for r in data['points']}:
            mid = next(r for r in data['points'] if r['d'] == ds)
        else:
            sd = nearest_seeds(data, dm)
            sd['role'] = 'bisect_' + name
            mid, seeds_out = compute_point(ds, sd, full_eig=False)
            mid['_a48'] = [nstr(x, 40) for x in seeds_out['a48']]
            data['points'].append(mid)
            data['points'].sort(key=lambda r: float(r['d']))
            save_points(data)
            plog('BISECT %s: d=%s -> f=%.3e (lam2=%s delta=%s, %.0fs)'
                 % (name, ds, f_rec(mid), mid['lam2']['lam'][:14],
                    mid['delta']['lam'][:14], mid['elapsed_s']))
        if f_rec(mid) == 0:
            lo = hi = mid
            break
        if f_rec(lo) * f_rec(mid) < 0:
            hi = mid
        else:
            lo = mid
        steps += 1
    dlo, dhi = float(lo['d']), float(hi['d'])
    out = {'bracket': [dlo, dhi], 'width': dhi - dlo,
           'f_lo': f_rec(lo), 'f_hi': f_rec(hi), 'steps': steps,
           'lam2_lo': l2_of(lo), 'lam2_hi': l2_of(hi)}
    if name == 'lam2_zero' and l2_of(hi) != l2_of(lo):
        out['linear_root'] = dlo - l2_of(lo) * (dhi - dlo) / (l2_of(hi) - l2_of(lo))
        out['slope_est'] = (l2_of(hi) - l2_of(lo)) / (dhi - dlo)
    if name == 'crossover':
        f_lo, f_hi = f_rec(lo), f_rec(hi)
        if f_hi != f_lo:
            out['linear_root'] = dlo - f_lo * (dhi - dlo) / (f_hi - f_lo)
    data['brackets'][name] = out
    save_points(data)
    plog('BISECT %s DONE: bracket [%.6f, %.6f] width %.1e  %s'
         % (name, dlo, dhi, dhi - dlo,
            ('linear root ~ %.6f' % out['linear_root']) if 'linear_root' in out else ''))
    return out


def main():
    plog('SPECTRA GRID START  ' + TAG)
    plog('N pair %d/%d, dps %d, eig ladder %s' % (N1, N2, mp.dps, LADDER_EIG_N))
    data = load_points()
    run_grid(data)
    plog('GRID COMPLETE (%d points). Brackets: lambda2 zero, crossover.' % len(data['points']))
    bisect(data, l2_of, 'lam2_zero')
    bisect(data, crossover_f, 'crossover')
    # repeat lam2-zero bisection refinement is implicit (target width 2e-3)
    plog('ALL DONE in %.0fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

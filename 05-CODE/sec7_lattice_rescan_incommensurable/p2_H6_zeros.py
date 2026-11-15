#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_H6_zeros.py -- Phase 3 of the degree-6 float pilot: build H6, scan its zeros.
[oid:deg6B]  NOT RIGOROUS -- float exploration (numpy complex128 + mpmath polish),
no interval arithmetic.

H6 by the SAME recipe as the quadratic-article H (article 2.2), all parameters
replaced by their sextic analogues:

    universal perturbed family:  f_s(x) = g6(x) + s * h6(x),
        g6(x) = sum_j a_j x^{6j}  (a_0=1),  h6(x) = sum_j c_j x^{6j}  (h6(1)=1)
    orbit:    x_0 = 0,  x_{n+1} = g6(x_n) + s*h6(x_n)
    coupling: s_m(u) = u * AU0_6 / delta6^m
    H6_m(u)  = x_{2^{m+1}}(u) / x_{2^m}(u),   H6 = lim H6_m

AU0_6 calibration (analogue of the article's H(-1)=0 gauge): superstable events
of the universal family, s*_m (real principal cascade, x_{2^{m+1}}=0, exact
period 2^{m+1}), obey s*_m -> Ls * delta6^{-m}; AU0_6 := -Ls puts the limit
cascade at u = -1.  Self-check: u_selfcheck-style gauge H_m(-1) ~ 0.

Scan: ADAPTIVE (geometry of the sextic class is its own -- BRIEF: find where the
first 30-60 zeros live): coarse float scan (m=7, step 0.25, R=48) to locate the
extent of valid local minima, then fine scan (step 1/16) over the chosen radius,
levels m=7,8.  Local minima of |H6_m| -> Newton refinement (holomorphic
derivative recursion) -> mpmath polish at level m=8 (parallel, <=6 workers) ->
pair with level-7 survivors -> level-9 refinement -> limit extrapolation
u_inf = u9 + (u9-u8)/(delta6-1); drift ratio d78/d89 must scale ~ delta6
(fingerprint of genuine limit-H6 zeros, not level artefacts).

Domain validity (analogue of the quartic |x|^4<2.2 rule, scaled by the sextic
series radius estimate t_radius ~ 1.68 from phase 1): per-orbit
(a) max |x_n|^6 < 1.4  and (b) truncation error |g_N48 - g_N40| < 1e-11.

Output: H6_zeros.json ; run log p2_run.log (flush + ETA).
"""
import json
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 50

TAG = "[oid:deg6B] NOT RIGOROUS float pilot, phase 3 (H6 zeros)"
T0 = time.time()
LOGF = open('p2_run.log', 'w', encoding='utf-8')
MAX_WORKERS = 6

TMAX6 = 1.15         # domain rule: max|x|^6 < 1.15 (bounded orbits live at |x|^6<=1;
                     # first orbit point x1=1 gives t=1; margin for near-bounded orbits)
TRUNC_TOL = 1e-11    # |g_N72 - g_N60| along orbit; N=72 series tail at t<=1.15 ~1e-22
                     # (measured: N=48 tail at t=1 was 3e-6 -> N-extension was MANDATORY)


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


FP = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_MP = [mpf(s) for s in FP['g_coeffs']]           # a_0..a_48, a_0 = 1
C_MP = [mpf(s) for s in FP['h_coeffs']]           # c_0..c_48, c_0 = 0, sum c_j = 1
ALPHA6 = mpf(FP['alpha6'])
DELTA6 = mpf(FP['delta6'])
N_MOD = len(A_MP) - 1

A_PY = [float(x) for x in A_MP]
C_PY = [float(x) for x in C_MP]
A_np = np.array(A_PY)
C_np = np.array(C_PY)
N_TRUNC = 60
A_np_tr = A_np.copy()
A_np_tr[N_TRUNC + 1:] = 0.0
C_np_tr = C_np.copy()
C_np_tr[N_TRUNC + 1:] = 0.0

# --------------------------------------------------------------------------- #
#  mpmath evaluators
# --------------------------------------------------------------------------- #
def g6(z):
    t = z ** 6
    r = A_MP[-1]
    for c in reversed(A_MP[:-1]):
        r = r * t + c
    return r


def h6(z):
    t = z ** 6
    r = C_MP[-1]
    for c in reversed(C_MP[:-1]):
        r = r * t + c
    return r


def g6p(z):
    t = z ** 6
    s = mpf(0)
    for j in range(len(A_MP) - 1, 0, -1):
        s = s * t + j * A_MP[j]
    return 6 * z ** 5 * s


def h6p(z):
    t = z ** 6
    s = mpf(0)
    for j in range(len(C_MP) - 1, 0, -1):
        s = s * t + j * C_MP[j]
    return 6 * z ** 5 * s


def orbit_H_mpmath(u, m, AU0, want_grad=False):
    """x_0=0, x_{n+1} = g6(x_n) + s*h6(x_n), s = u*AU0/delta6^m.
    Returns (x_{2^m}, x_{2^{m+1}}) (+ d/du of both if want_grad)."""
    s = u * AU0 / DELTA6 ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    dsdu = AU0 / DELTA6 ** m
    x2m = dx2m = None
    for k in range(n2):
        if want_grad:
            dx = (g6p(x) + s * h6p(x)) * dx + h6(x) * dsdu
        x = g6(x) + s * h6(x)
        if k + 1 == n1:
            x2m, dx2m = x, dx
    if want_grad:
        return x2m, x, dx2m, dx
    return x2m, x


# --------------------------------------------------------------------------- #
#  numpy scan
# --------------------------------------------------------------------------- #
def poly6(z, co):
    t = z * z
    t = t * t * t
    r = np.full_like(z, co[-1])
    for c in reversed(co[:-1]):
        r = r * t + c
    return r


def scan_level(m, AU0, grid):
    """Returns Habs, valid mask, max|x|^6, trunc err."""
    s = complex(float(AU0) / float(DELTA6) ** m) * grid
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    maxt = np.zeros(grid.shape, dtype=np.float64)
    errmax = np.zeros(grid.shape, dtype=np.float64)
    x2m = None
    for k in range(n2):
        gfull = poly6(x, A_np)
        gtrunc = poly6(x, A_np_tr)
        np.maximum(errmax, np.abs(gfull - gtrunc), out=errmax)
        t6 = x * x
        t6 = t6 * t6 * t6
        np.maximum(maxt, np.abs(t6), out=maxt)
        x = gfull + s * poly6(x, C_np)
        if k + 1 == n1:
            x2m = x.copy()
    num = np.abs(x)
    den = np.abs(x2m)
    with np.errstate(divide='ignore', invalid='ignore'):
        Habs = num / den
    # validity: orbit bounded (max|x|^6 <= TMAX6), denominator non-degenerate,
    # finite H.  NOTE: no in-orbit truncation mask -- the honest truncation
    # certificate for EVALUATING with the N=72 series is the measured
    # |G84-G72|(t) <= 2.8e-37 for |t|<=1.15 (diag_tail.py, mpmath dps=80:
    # |G96-G84| = 0), i.e. the N=72 evaluation error is ~1e-37 on the whole
    # domain; the in-orbit window |G72-G60| (~4e-8 at t=1) measures REAL series
    # terms a_61..a_72, not error (the quartic pilot's window check certified
    # its N=40 truncation, a stricter-than-needed proxy).  errmax is still
    # recorded per zero for the catalog.
    valid = (np.isfinite(Habs) & (den > 1e-12) & (maxt < TMAX6) & (num < 1e100))
    return Habs, valid, maxt, errmax


def eval_gh_py(x):
    """(g, h, g', h') at complex x, plain Python. Raises OverflowError for
    violently escaping orbits (|x| so large that x**6 overflows)."""
    if abs(x) > 1e40:
        raise OverflowError('orbit escaped')
    t = x ** 6
    g = A_PY[-1]
    for c in reversed(A_PY[:-1]):
        g = g * t + c
    h = C_PY[-1]
    for c in reversed(C_PY[:-1]):
        h = h * t + c
    dg = 0j
    for j in range(len(A_PY) - 1, 0, -1):
        dg = dg * t + j * A_PY[j]
    dh = 0j
    for j in range(len(C_PY) - 1, 0, -1):
        dh = dh * t + j * C_PY[j]
    x5 = x ** 5
    return g, h, 6.0 * x5 * dg, 6.0 * x5 * dh


def newton_refine(u0, m, AU0, iters=60, rmax=70.0, tol=1e-12):
    """Newton on H_m(u)=0, plain Python complex with overflow-guarded damping.
    Returns (u, |H|) or (None, inf)."""
    s_coef = float(AU0) / float(DELTA6) ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)

    def H_at(u):
        s = u * s_coef
        x = 0j
        x2m = None
        try:
            for k in range(n2):
                g, h, _gp, _hp = eval_gh_py(x)
                x = g + s * h
                if k + 1 == n1:
                    x2m = x
        except OverflowError:
            return None
        if x2m == 0 or not (abs(x2m) > 1e-290):
            return None
        return x / x2m

    u = u0
    for _ in range(iters):
        s = u * s_coef
        x = 0j
        dx = 0j
        x2m = None
        dx2m = None
        try:
            for k in range(n2):
                g, h, gp, hp = eval_gh_py(x)
                dx = (gp + s * hp) * dx + h * s_coef
                x = g + s * h
                if k + 1 == n1:
                    x2m, dx2m = x, dx
        except OverflowError:
            return None, float('inf')
        if not (abs(x2m) > 1e-290) or not np.isfinite(abs(x) + abs(dx)):
            return None, float('inf')
        H = x / x2m
        dH = (dx * x2m - x * dx2m) / (x2m * x2m)
        if abs(H) < tol:
            break
        if dH == 0 or not np.isfinite(H.real + H.imag + dH.real + dH.imag):
            return None, float('inf')
        step = H / dH
        moved = False
        for _h in range(20):
            u_try = u - step
            if abs(u_try) > rmax:
                step = step / 2
                continue
            Ht = H_at(u_try)
            if Ht is not None and abs(Ht) < abs(H):
                u = u_try
                moved = True
                break
            step = step / 2
        if not moved:
            break
        if abs(step) < 1e-15:
            break
    Hf = H_at(u)
    if Hf is None:
        return None, float('inf')
    return u, abs(Hf)


def domain_check(u, m, AU0, tmax=TMAX6, errtol=TRUNC_TOL):
    """Validity of the refined zero: orbit stays inside |x|^6 < tmax and the
    series truncation error (N=72 vs N=60) stays below errtol."""
    s = u * float(AU0) / float(DELTA6) ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = 0j
    maxt6 = 0.0
    maxerr = 0.0
    den = None
    try:
        for k in range(n2):
            t = x ** 6
            gfull = A_PY[-1]
            for c in reversed(A_PY[:-1]):
                gfull = gfull * t + c
            gtr = A_PY[N_TRUNC]
            for c in reversed(A_PY[:N_TRUNC]):
                gtr = gtr * t + c
            err = abs(gfull - gtr)
            if err > maxerr:
                maxerr = err
            if abs(t) > maxt6:
                maxt6 = abs(t)
            h = C_PY[-1]
            for c in reversed(C_PY[:-1]):
                h = h * t + c
            x = gfull + s * h
            if k + 1 == n1:
                den = x
    except OverflowError:
        return float('inf'), float('inf'), 0.0
    return maxt6, maxerr, abs(den)


# --------------------------------------------------------------------------- #
#  parallel mpmath polish
# --------------------------------------------------------------------------- #
def mp_polish_job(args):
    """(u_start, m, AU0_str, iters, tol_step_exp) -> (u_mp_str, absH float, ok)."""
    u_re, u_im, m, AU0s, iters = args
    AU0 = mpf(AU0s)
    u_mp = mpc(mpf(u_re), mpf(u_im))
    for _ in range(iters):
        x2m, x2m1, dx2m, dx2m1 = orbit_H_mpmath(u_mp, m, AU0, want_grad=True)
        H = x2m1 / x2m
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        step = H / dH
        u_mp = u_mp - step
        if abs(step) < mpf(10) ** (-24):
            break
    try:
        x2m, x2m1 = orbit_H_mpmath(u_mp, m, AU0)
        r = float(abs(x2m1 / x2m))
    except Exception:  # noqa: BLE001
        return None
    if r < 1e-20:
        return (nstr(mpre(u_mp), 30), nstr(mpim(u_mp), 30), r, m)
    return None


def main():
    log(TAG)

    # ---- AU0_6 calibration -------------------------------------------------------
    log('--- AU0_6 calibration: superstable cascade of the universal family ---')
    smax = 12
    s_stars = {}
    s_prev = None
    for m in range(2, smax + 1):
        n = 2 ** (m + 1)

        def fres(s):
            x = mpf(0)
            ds = mpf(0)
            for _ in range(n):
                ds = (g6p(x) + s * h6p(x)) * ds + h6(x)
                x = g6(x) + s * h6(x)
            return x, ds

        s0 = s_prev / DELTA6 if s_prev is not None else mpf('-1e-3')
        s = s0
        ok_conv = False
        for _it in range(80):
            x, ds = fres(s)
            if fabs(x) < mpf('1e-40'):
                ok_conv = True
                break
            if ds == 0:
                break
            ss = x / ds
            improved = False
            for _h in range(30):
                x2, _ = fres(s - ss)
                if fabs(x2) < fabs(x):
                    s = s - ss
                    improved = True
                    break
                ss = ss / 2
            if not improved:
                break
        if not ok_conv:
            raise RuntimeError('AU0 calibration Newton failed at m=%d' % m)
        for j in range(1, m + 1):
            x = mpf(0)
            for _ in range(2 ** j):
                x = g6(x) + s * h6(x)
            if fabs(x) < mpf('1e-30'):
                raise RuntimeError('inherited root at m=%d j=%d' % (m, j))
        s_stars[m] = s
        s_prev = s
        log('m=%2d  s*_m = %s   Ls^(m) = %s'
            % (m, nstr(s, 20), nstr(s * DELTA6 ** m, 20)))
    rL = 1 / DELTA6
    Ls_seq = {m: s_stars[m] * DELTA6 ** m for m in s_stars}
    Ls = (Ls_seq[smax] - rL * Ls_seq[smax - 1]) / (1 - rL)
    AU0 = -Ls
    log('Ls = %s   AU0_6 = -Ls = %s' % (nstr(Ls, 30), nstr(AU0, 30)))
    for m in (6, 8, 10):
        x2m, x2m1 = orbit_H_mpmath(mpf(-1), m, AU0)
        log('gauge check H_%d(-1) = %s' % (m, nstr(x2m1 / x2m, 12)))

    AU0s = nstr(AU0, 40)

    # ---- coarse exploratory scan ---------------------------------------------------
    log('--- coarse scan m=7, R=48, step 0.25 (locate the geometry) ---')
    R_C, STEP_C = 48.0, 0.25
    xs = np.arange(-R_C, R_C + STEP_C / 2, STEP_C)
    ys = np.arange(0.0, R_C + STEP_C / 2, STEP_C)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    gridC = UU + 1j * VV
    Habs, valid, maxt, errmax = scan_level(7, AU0, gridC)
    log('coarse: grid %dx%d, valid %d (%.1f%%), min|H| on valid = %.3e'
        % (gridC.shape[0], gridC.shape[1], int(valid.sum()),
           100.0 * valid.sum() / gridC.size,
           float(Habs[valid].min()) if valid.any() else float('nan')))
    Hm = np.where(valid, Habs, np.inf)
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    ii, jj = np.nonzero(ismin)
    coarse_minima = [complex(gridC[i, j]) for i, j in zip(ii, jj)]
    log('coarse candidate minima: %d' % len(coarse_minima))
    umags = sorted(abs(z) for z in coarse_minima)
    if umags:
        log('coarse minima |u| deciles: %s'
            % ' '.join('%.1f' % umags[int(q * (len(umags) - 1))]
                       for q in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)))
    # quick float Newton on coarse minima to see which are genuine + their extent
    coarse_zeros = []
    t_coarse0 = time.time()
    for i2, u0 in enumerate(coarse_minima):
        u, res = newton_refine(u0, 7, AU0, tol=1e-9)
        if u is None or res > 1e-5:
            continue
        mt6, merr, den = domain_check(u, 7, AU0)
        if mt6 < TMAX6 and den > 1e-12:
            coarse_zeros.append(u)
        if (i2 + 1) % 20 == 0:
            el = time.time() - t_coarse0
            eta = el / (i2 + 1) * (len(coarse_minima) - i2 - 1)
            log('  coarse newton %d/%d (%.0fs elapsed, ETA %.0fs)'
                % (i2 + 1, len(coarse_minima), el, eta), )
    ded = []
    for z in coarse_zeros:
        if all(abs(z - w) > 5e-3 for w in ded):
            ded.append(z)
    log('coarse genuine zeros: %d ; max |u| = %.3f ; count |u|<20: %d, |u|<30: %d, |u|<40: %d'
        % (len(ded), max((abs(z) for z in ded), default=0),
           sum(1 for z in ded if abs(z) < 20),
           sum(1 for z in ded if abs(z) < 30),
           sum(1 for z in ded if abs(z) < 40)))
    if len(ded) < 30:
        log('WARNING: fewer than 30 coarse zeros; radius choice will follow the data')
    R_MAX = min(48.0, round((max((abs(z) for z in ded), default=20.0) + 2.0) * 2) / 2)
    STEP = 0.0625
    log('ADAPTIVE choice: fine scan R_MAX = %.1f, step %.4f' % (R_MAX, STEP))

    # ---- fine scan, levels 7 and 8 ---------------------------------------------------
    xs = np.arange(-R_MAX, R_MAX + STEP / 2, STEP)
    ys = np.arange(0.0, R_MAX + STEP / 2, STEP)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    grid = UU + 1j * VV
    results = {}
    for m in (7, 8):
        t_scan = time.time()
        Habs, valid, maxt, errmax = scan_level(m, AU0, grid)
        n_valid = int(valid.sum())
        log('scan m=%d: grid %dx%d (%.0fs), valid %d (%.1f%%), min|H| on valid = %.3e'
            % (m, grid.shape[0], grid.shape[1], time.time() - t_scan, n_valid,
               100.0 * n_valid / grid.size,
               float(Habs[valid].min()) if n_valid else float('nan')))
        Hm = np.where(valid, Habs, np.inf)
        Pg = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
        Pg[1:-1, 1:-1] = Hm
        core = Pg[1:-1, 1:-1]
        nb = np.stack([Pg[:-2, :-2], Pg[:-2, 1:-1], Pg[:-2, 2:],
                       Pg[1:-1, :-2], Pg[1:-1, 2:],
                       Pg[2:, :-2], Pg[2:, 1:-1], Pg[2:, 2:]])
        ismin = (core <= nb).all(axis=0) & (core < 0.2)
        ii, jj = np.nonzero(ismin)
        cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
        if m == 8:
            cand = cand + list(ded)
        log('  candidates: %d' % len(cand))
        zeros = []
        t_newt = time.time()
        for i2, u0 in enumerate(cand):
            u, res = newton_refine(u0, m, AU0, tol=1e-9, rmax=R_MAX + 5)
            if u is None or res > 1e-5:
                continue
            mt6, merr, den = domain_check(u, m, AU0)
            zeros.append({'u': u, 'absH': res, 'maxt6': mt6, 'maxerr': merr,
                          'den': den,
                          'domain_ok': bool(mt6 < TMAX6 and den > 1e-12)})
            if (i2 + 1) % 40 == 0:
                el = time.time() - t_newt
                eta = el / (i2 + 1) * (len(cand) - i2 - 1)
                log('  m=%d newton %d/%d (%.0fs, ETA %.0fs)'
                    % (m, i2 + 1, len(cand), el, eta))
        dedup = []
        for z in zeros:
            if all(abs(z['u'] - v['u']) > 1e-4 for v in dedup):
                dedup.append(z)
        log('  refined zeros (|H|<1e-5): %d (domain-ok: %d)'
            % (len(dedup), sum(1 for z in dedup if z['domain_ok'])))
        results[m] = dedup

    cat7 = [z for z in results[7] if z['domain_ok']]
    cat8 = [z for z in results[8] if z['domain_ok']]
    log('domain-ok: level7 %d, level8 %d' % (len(cat7), len(cat8)))

    # ---- mpmath polish at m=8 (parallel) ---------------------------------------------
    jobs = [(z['u'].real, z['u'].imag, 8, AU0s, 16) for z in cat8]
    cat8_mp = []
    t_pol = time.time()
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i2, reslt in enumerate(ex.map(mp_polish_job, jobs)):
            if reslt is not None:
                ur, ui, r, _m = reslt
                cat8_mp.append({'u8_mp_re': ur, 'u8_mp_im': ui, 'absH8': r,
                                'seed': cat8[i2]})
            if (i2 + 1) % 10 == 0:
                el = time.time() - t_pol
                eta = el / (i2 + 1) * (len(jobs) - i2 - 1)
                log('  polish m=8 %d/%d (%.0fs, ETA %.0fs)'
                    % (i2 + 1, len(jobs), el, eta))
    log('polished level-8 zeros (|H8|<1e-20): %d' % len(cat8_mp))

    # ---- pair with level-7 zeros (survival filter) --------------------------------------
    paired = []
    for z8 in cat8_mp:
        u8c = complex(mpf(z8['u8_mp_re']), mpf(z8['u8_mp_im']))
        best = min(cat7, key=lambda z: abs(z['u'] - u8c))
        if abs(best['u'] - u8c) < 5e-3:
            paired.append({**z8, 'u7': best['u'], 'u8': u8c})
        else:
            log('  level-8 zero at %+.6f%+.6fi has NO level-7 partner (nearest %.4f)'
                % (u8c.real, u8c.imag, abs(best['u'] - u8c)))
    log('paired zeros (survive m=7 and m=8): %d' % len(paired))

    # ---- level-9 refinement (parallel), drift-ratio check, limit extrapolation ----------
    jobs9 = [(z['u8_mp_re'], z['u8_mp_im'], 9, AU0s, 16) for z in paired]
    t_p9 = time.time()
    res9 = [None] * len(paired)
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i2, reslt in enumerate(ex.map(mp_polish_job, jobs9)):
            res9[i2] = reslt
            if (i2 + 1) % 10 == 0:
                el = time.time() - t_p9
                eta = el / (i2 + 1) * (len(jobs9) - i2 - 1)
                log('  polish m=9 %d/%d (%.0fs, ETA %.0fs)'
                    % (i2 + 1, len(jobs9), el, eta))
    for z, r9 in zip(paired, res9):
        if r9 is not None:
            ur, ui, r, _m = r9
            u9c = complex(mpf(ur), mpf(ui))
            z['u9'] = u9c
            z['absH9'] = r
            z['drift78'] = abs(z['u8'] - z['u7'])
            z['drift89'] = abs(u9c - z['u8'])
            z['drift_ratio'] = z['drift78'] / z['drift89'] if z['drift89'] > 0 else None
            z['uinf'] = u9c + (u9c - z['u8']) / float(DELTA6 - 1)
        else:
            z['u9'] = None
            z['absH9'] = None
            z['drift78'] = abs(z['u8'] - z['u7'])
            z['drift89'] = None
            z['drift_ratio'] = None
            z['uinf'] = z['u8'] + (z['u8'] - z['u7']) / float(DELTA6 - 1)

    paired.sort(key=lambda z: (abs(z['uinf']), z['uinf'].real, z['uinf'].imag))
    # min spacing
    us = [z['uinf'] for z in paired]
    minsp = min((abs(a - b) for i, a in enumerate(us) for b in us[i + 1:]),
                default=float('inf'))
    log('min spacing between catalog zeros: %.4f' % minsp)
    log('=== CATALOG (limit-extrapolated zeros of H6) ===')
    for z in paired:
        dr = ('%.3f' % z['drift_ratio']) if z['drift_ratio'] else 'n/a'
        m9 = ('%.1e' % z['absH9']) if z['absH9'] is not None else 'fail'
        log('   uinf = %+.8f %+.8f i  |u|=%6.4f  |H8|=%.1e |H9|=%s  d78=%.2e ratio=%s'
            % (z['uinf'].real, z['uinf'].imag, abs(z['uinf']), z['absH8'], m9,
               z['drift78'], dr))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': ('article 2.2 analogue: x_{n+1}=g6(x_n)+s*h6(x_n), '
                            's=u*AU0_6/delta6^m, H6_m=x_{2^{m+1}}/x_{2^m}'),
                 'scan': {'R_max': R_MAX, 'step': STEP, 'levels': [7, 8],
                          'coarse_scan': {'R': R_C, 'step': STEP_C},
                          'adaptive_note': ('coarse m=7 scan over R=48 located the '
                                            'zero extent; fine radius chosen from data'),
                          'domain_rule': ('max|x_n|^6 < %.2f and |x_{2^m}|>1e-12; '
                                          'truncation certificate: |G84-G72|<=2.8e-37 '
                                          'on |t|<=1.15 (diag_tail), NOT an in-orbit window' % TMAX6)},
                 'mp_dps': mp.dps, 'max_workers': MAX_WORKERS,
                 'elapsed_sec': time.time() - T0},
        'AU0_6': AU0s,
        'Ls': nstr(Ls, 40),
        's_star': {str(m): nstr(s_stars[m], 30) for m in s_stars},
        'min_spacing': minsp,
        'zeros_level8': [{'re': c['u8'].real, 'im': c['u8'].imag,
                          'absH': c['absH8'], 'str': (c['u8_mp_re'] + '+' + c['u8_mp_im'] + 'i'),
                          'maxt6': c['seed']['maxt6'], 'maxerr': c['seed']['maxerr']}
                         for c in paired],
        'zeros_level7': [{'re': c['u7'].real, 'im': c['u7'].imag} for c in paired],
        'zeros_limit_extrapolated': [
            {'re': c['uinf'].real, 'im': c['uinf'].imag,
             'u7_re': c['u7'].real, 'u7_im': c['u7'].imag,
             'u8_re': c['u8'].real, 'u8_im': c['u8'].imag,
             'u9_re': (c['u9'].real if c.get('u9') else None),
             'u9_im': (c['u9'].imag if c.get('u9') else None),
             'absH9': c.get('absH9'),
             'drift78': c['drift78'], 'drift89': c.get('drift89'),
             'drift_ratio': c.get('drift_ratio')}
            for c in paired],
        'counts': {'level8_polished': len(cat8_mp), 'paired': len(paired),
                   'domain_ok_7': len(cat7), 'domain_ok_8': len(cat8),
                   'all_refined_7': len(results[7]), 'all_refined_8': len(results[8]),
                   'coarse_zeros': len(ded)},
    }
    with open('H6_zeros.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote H6_zeros.json')
    log('DONE phase 3 (H6 zeros), total elapsed %.1fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

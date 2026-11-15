#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Brief S -- rescan a 4x4 patch of u-plane at level m=8 with FOUR griddings:
   A) step 1/16 (control: reproduce catalog), B) step 1/24, C) step 1/10,
   D) 2000 random starts.  Each seed -> float Newton -> domain rule -> mpmath
   polish (dps=70, tol 1e-30).  Cluster by 1e-8.

NOT RIGOROUS.  Reuses sextic_fixedpoint.json (g/h coeffs, alpha6, delta6) and
H6_zeros.json (AU0_6) from agentB_deg6/.  Reimplements scan_level / newton_refine
/ domain_check / orbit_H_mpmath as in p2_H6_zeros.py.
"""
import json
import os
import sys
import time
import math

import numpy as np
from mpmath import mp, mpf, mpc, nstr, re as mpre, im as mpim, fabs

mp.dps = 70

ROOT = 'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
HERE = os.path.dirname(os.path.abspath(__file__))
DEG6 = os.path.join(os.path.dirname(HERE), 'agentB_deg6')

T0 = time.time()
LOGF = open(os.path.join(HERE, 'progress.log'), 'a', encoding='utf-8')


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


# ---------------- load constants ---------------- #
log('loading sextic_fixedpoint.json + H6_zeros.json (meta, AU0_6)')
fp = json.load(open(os.path.join(DEG6, 'sextic_fixedpoint.json'), encoding='utf-8'))
hz = json.load(open(os.path.join(DEG6, 'H6_zeros.json'), encoding='utf-8'))

A_MP = [mpf(s) for s in fp['g_coeffs']]
C_MP = [mpf(s) for s in fp['h_coeffs']]
DELTA6_MP = mpf(fp['delta6'])
AU0_MP = mpf(hz['AU0_6'])
AU0 = float(AU0_MP)
DELTA6 = float(DELTA6_MP)
M_LEVEL = 8
N_TRUNC = 60
TMAX6 = 1.15

A_PY = [float(x) for x in A_MP]
C_PY = [float(x) for x in C_MP]
A_np = np.array(A_PY, dtype=np.float64)
C_np = np.array(C_PY, dtype=np.float64)
A_np_tr = A_np.copy()
A_np_tr[N_TRUNC + 1:] = 0.0
C_np_tr = C_np.copy()
C_np_tr[N_TRUNC + 1:] = 0.0

log('AU0_6 = %s  (from H6_zeros.json)' % hz['AU0_6'])
log('DELTA6 = %.20f   N coeffs = %d' % (float(DELTA6_MP), len(A_MP)))


# ---------------- polynomial evaluators ---------------- #
def poly6(z, co):
    t = z * z
    t = t * t * t
    r = np.full(z.shape if hasattr(z, 'shape') else (), co[-1])
    for c in reversed(co[:-1]):
        r = r * t + c
    return r


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


def orbit_H_mpmath(u, m, AU0v, want_grad=False):
    s = u * AU0v / DELTA6_MP ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    dsdu = AU0v / DELTA6_MP ** m
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


# ---------------- numpy scan ---------------- #
def scan_level(grid, m=M_LEVEL):
    s = (AU0 / DELTA6 ** m) * grid
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
    valid = (np.isfinite(Habs) & (den > 1e-12) & (maxt < TMAX6) & (num < 1e100))
    return Habs, valid, maxt, errmax


def detect_local_minima(Habs, valid, xs, ys):
    Hm = np.where(valid, Habs, np.inf)
    if Hm.size < 9:
        return []
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    ii, jj = np.nonzero(ismin)
    return [(float(xs[i]), float(ys[j]), float(core[i, j]))
            for i, j in zip(ii, jj)]


# ---------------- float Newton refinement (matches p2_H6_zeros.newton_refine) ---------------- #
def eval_gh_py(x):
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


def newton_refine(u0, m=M_LEVEL, iters=80, rmax=70.0, tol=1e-9):
    s_coef = AU0 / DELTA6 ** m
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


def domain_check(u, m=M_LEVEL):
    s = u * AU0 / DELTA6 ** m
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


# ---------------- mpmath polish ---------------- #
def mp_polish_full(u_start_re, u_start_im, max_iter=200, tol=mpf('1e-40')):
    u_mp = mpc(mpf(u_start_re), mpf(u_start_im))
    last_H = None
    for _it in range(max_iter):
        x2m, x2m1, dx2m, dx2m1 = orbit_H_mpmath(u_mp, M_LEVEL, AU0_MP, want_grad=True)
        H = x2m1 / x2m
        last_H = abs(H)
        if last_H < tol:
            break
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        if dH == 0:
            break
        step = H / dH
        moved = False
        for _h in range(15):
            ut = u_mp - step
            x2mt, x2m1t = orbit_H_mpmath(ut, M_LEVEL, AU0_MP)
            if fabs(x2m1t / x2mt) < last_H:
                u_mp = ut
                moved = True
                break
            step = step / 2
        if not moved:
            break
        if fabs(step) < mpf('1e-45'):
            break
    x2m, x2m1 = orbit_H_mpmath(u_mp, M_LEVEL, AU0_MP)
    return u_mp, float(fabs(x2m1 / x2m))


def finalize_seed(re_seed, im_seed):
    """float Newton -> domain check -> mpmath polish to <=1e-30.  Returns
    tuple (re, im, absH, status) where re, im are mpmath nstr (35 digits)
    strings preserving full precision (otherwise json float64 would lose the
    digits needed to recover absH<=1e-40 at m=8 on re-eval)."""
    u, res = newton_refine(complex(re_seed, im_seed), M_LEVEL, tol=1e-9, rmax=60)
    if u is None or res > 1e-5 or not np.isfinite(res):
        return (str(re_seed), str(im_seed), float('inf'), 'refine_fail')
    mt6, merr, den = domain_check(u, M_LEVEL)
    if not (mt6 < TMAX6 and den > 1e-12):
        return (str(u.real), str(u.imag), res, 'domain_fail')
    try:
        u_mp, absH = mp_polish_full(str(re_seed), str(im_seed))
    except Exception:
        return (str(u.real), str(u.imag), res, 'mpmath_fail')
    re_s = nstr(mpre(u_mp), 35)
    im_s = nstr(mpim(u_mp), 35)
    if absH < 1e-30 and np.isfinite(absH):
        return (re_s, im_s, absH, 'polished')
    return (re_s, im_s, absH, 'float_only')


# ---------------- cluster ---------------- #
def cluster_zeros(zeros, eps=1e-8):
    """zeros: list of (re_str, im_str, absH, status).  Returns list of dict
    entries: {re, im, absH, status, count} with re/im as strings (full precision).
    Aggregate nearest duplicates within eps (using float comparison)."""
    out = []
    for re_str, im_str, absH, status in zeros:
        re_f = float(re_str)
        im_f = float(im_str)
        placed = False
        for o in out:
            if abs(re_f - float(o['re'])) < eps and abs(im_f - float(o['im'])) < eps:
                if absH < o['absH'] or o['status'] == 'float_only':
                    o['re'], o['im'], o['absH'], o['status'] = re_str, im_str, absH, status
                o['count'] += 1
                placed = True
                break
        if not placed:
            out.append({'re': re_str, 'im': im_str, 'absH': absH,
                        'status': status, 'count': 1})
    return out


def match_sets(setA, setB, eps=1e-6):
    """For each zero in setA find nearest in setB; pair if dist <= eps.
    Returns (pairs, a_only, b_only)."""
    pairs = []
    used = set()
    for i, a in enumerate(setA):
        a_re = float(a['re'])
        a_im = float(a['im'])
        best = (-1, 1e9)
        for j, b in enumerate(setB):
            if j in used:
                continue
            d = math.hypot(a_re - float(b['re']), a_im - float(b['im']))
            if d < best[1]:
                best = (j, d)
        if best[0] >= 0 and best[1] <= eps:
            pairs.append((i, best[0], best[1]))
            used.add(best[0])
    a_only = [i for i in range(len(setA))
              if not any(p[0] == i for p in pairs)]
    b_only = [j for j in range(len(setB)) if j not in used]
    return pairs, a_only, b_only


# ---------------- patch geometry ---------------- #
PATCH = dict(re_lo=18.0, re_hi=22.0, im_lo=17.0, im_hi=21.0)


def make_grid(step):
    xs = np.arange(PATCH['re_lo'], PATCH['re_hi'] + step / 2, step)
    ys = np.arange(PATCH['im_lo'], PATCH['im_hi'] + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return xs, ys, UU + 1j * VV


def scan_grid(step):
    xs, ys, grid = make_grid(step)
    log('grid step=1/%d: %d x %d = %d nodes'
        % (int(round(1 / step)), xs.size, ys.size, grid.size))
    Habs, valid, maxt, err = scan_level(grid, M_LEVEL)
    cands = detect_local_minima(Habs, valid, xs, ys)
    log('  candidates (local minima |H|<0.2): %d' % len(cands))
    return xs, ys, cands


def finalize_candidates(cands, label):
    zeros = []
    t0 = time.time()
    n = len(cands)
    for i2, (re, im, hmin) in enumerate(cands):
        z = finalize_seed(re, im)
        zeros.append(z)
        if (i2 + 1) % 25 == 0 or (i2 + 1) == n:
            el = time.time() - t0
            log('  %s finalize %d/%d  absH=%.1e  (%.0fs, ETA %.0fs)'
                % (label, i2 + 1, n, z[2], el,
                   el / (i2 + 1) * (n - i2 - 1)))
    return zeros


def random_scan(n_starts, seed=20260815):
    rng = np.random.default_rng(seed)
    xs = rng.uniform(PATCH['re_lo'], PATCH['re_hi'], n_starts)
    ys = rng.uniform(PATCH['im_lo'], PATCH['im_hi'], n_starts)
    return [(float(x), float(y)) for x, y in zip(xs, ys)]


# ---------------- main flow ---------------- #
def run_scan(name, step, cands):
    log('=== scan %s (step=%.6f) ===' % (name, step))
    zeros = finalize_candidates(cands, name)
    pol = [z for z in zeros if z[3] == 'polished']
    fo = [z for z in zeros if z[3] == 'float_only']
    df = [z for z in zeros if z[3] == 'domain_fail']
    rf = [z for z in zeros if z[3] == 'refine_fail']
    mp = [z for z in zeros if z[3] == 'mpmath_fail']
    log('   refined: %d  polished<=1e-30: %d  float_only: %d  domain_fail: %d  refine_fail: %d  mp_fail: %d'
        % (len(zeros), len(pol), len(fo), len(df), len(rf), len(mp)))
    clustered = cluster_zeros(pol + fo, eps=1e-8)
    log('   after cluster(1e-8): %d unique zeros' % len(clustered))
    return zeros, clustered


def summary_stats(clustered):
    qs = {'q16_re': [], 'q16_im': [], 'q24_re': [], 'q24_im': [],
          'q10_re': [], 'q10_im': []}
    for z in clustered:
        re, im = float(z['re']), float(z['im'])
        qs['q16_re'].append((16 * re) - math.floor(16 * re))
        qs['q16_im'].append((16 * im) - math.floor(16 * im))
        qs['q24_re'].append((24 * re) - math.floor(24 * re))
        qs['q24_im'].append((24 * im) - math.floor(24 * im))
        qs['q10_re'].append((10 * re) - math.floor(10 * re))
        qs['q10_im'].append((10 * im) - math.floor(10 * im))

    def dist_to_lattice(fractions, q):
        return [min(f, 1 - f) if f > 0 else 0.0 for f in fractions]

    out = {}
    for name, fracs in qs.items():
        q = int(name.split('_')[0][1:])
        ds = dist_to_lattice(fracs, q)
        if not ds:
            out[name] = dict(n=0, frac_near=0, median_d=0, min_d=0)
            continue
        out[name] = dict(
            n=len(ds),
            frac_near=sum(1 for d in ds if d < 0.004) / len(ds),
            median_d=float(np.median(ds)),
            min_d=float(min(ds)),
            max_d=float(max(ds)),
        )
    return out


def nearest_lattice_dist(re, im, q):
    return math.hypot(
        min((q * re) - math.floor(q * re), 1 - ((q * re) - math.floor(q * re))),
        min((q * im) - math.floor(q * im), 1 - ((q * im) - math.floor(q * im))),
    )


# ---------------- main ---------------- #
def main():
    log('--- rescan START (patch %s) ---' % str(PATCH))
    # catalog zeros in patch for cross-check
    cat_patch = [z for z in hz['zeros_level8']
                 if PATCH['re_lo'] <= z['re'] <= PATCH['re_hi']
                 and PATCH['im_lo'] <= z['im'] <= PATCH['im_hi']]
    log('CATALOG zeros inside patch: %d' % len(cat_patch))
    for z in cat_patch:
        log('  cat: re=%.10f  im=%.10f  absH=%.2e  str=%s'
            % (z['re'], z['im'], z['absH'], z['str']))

    # ---------------- scan A ---------------- #
    xsA, ysA, candsA = scan_grid(1.0 / 16)
    zA, cA = run_scan('A', 1.0 / 16, candsA)

    # ---------------- scan B ---------------- #
    xsB, ysB, candsB = scan_grid(1.0 / 24)
    zB, cB = run_scan('B', 1.0 / 24, candsB)

    # ---------------- scan C ---------------- #
    xsC, ysC, candsC = scan_grid(1.0 / 10)
    zC, cC = run_scan('C', 1.0 / 10, candsC)

    # ---------------- scan D (random) ---------------- #
    log('=== scan D (1000 random starts in patch) ===')
    startsD = random_scan(1000)
    zD = []
    t0 = time.time()
    n = len(startsD)
    for i2, (re, im) in enumerate(startsD):
        z = finalize_seed(re, im)
        zD.append(z)
        if (i2 + 1) % 200 == 0 or (i2 + 1) == n:
            el = time.time() - t0
            log('  D finalize %d/%d  absH=%.1e  (%.0fs, ETA %.0fs)'
                % (i2 + 1, n, z[2], el,
                   el / (i2 + 1) * (n - i2 - 1)))
    polD = [z for z in zD if z[3] == 'polished']
    foD = [z for z in zD if z[3] == 'float_only']
    dfD = [z for z in zD if z[3] == 'domain_fail']
    rfD = [z for z in zD if z[3] == 'refine_fail']
    mpD = [z for z in zD if z[3] == 'mpmath_fail']
    log('  D refined: %d  polished: %d  float_only: %d  domain_fail: %d  refine_fail: %d  mp_fail: %d'
        % (len(zD), len(polD), len(foD), len(dfD), len(rfD), len(mpD)))
    cD = cluster_zeros(polD + foD, eps=1e-8)
    log('  D cluster(1e-8): %d unique zeros' % len(cD))

    # ---------------- analyze ---------------- #
    log('--- analysis ---')
    raw_sets = dict(A=zA, B=zB, C=zC, D=zD)
    cluster_sets = dict(A=cA, B=cB, C=cC, D=cD)
    stats = {n: summary_stats(c) for n, c in cluster_sets.items()}

    for nm in 'ABCD':
        s = stats[nm]
        log('  %s: n=%d  q16 frac<0.004: re=%.2f im=%.2f  q24 re=%.2f im=%.2f  q10 re=%.2f im=%.2f  q16 med dist=%.3e (min=%.3e)'
            % (nm, s['q16_re']['n'],
               s['q16_re']['frac_near'], s['q16_im']['frac_near'],
               s['q24_re']['frac_near'], s['q24_im']['frac_near'],
               s['q10_re']['frac_near'], s['q10_im']['frac_near'],
               s['q16_re']['median_d'], s['q16_re']['min_d']))

    # pairwise intersections
    inter = {}
    keys = list(cluster_sets.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            n1, n2 = keys[i], keys[j]
            p, ao, bo = match_sets(cluster_sets[n1], cluster_sets[n2], 1e-6)
            inter[(n1, n2)] = dict(common=len(p), only_n1=len(ao), only_n2=len(bo),
                                   pairs=p)
            log('  intersect %s x %s (1e-6): common=%d  only_%s=%d  only_%s=%d'
                % (n1, n2, len(p), n1, len(ao), n2, len(bo)))

    # union of all four clustered sets
    union = []
    seen = []
    for nm in 'ABCD':
        for z in cluster_sets[nm]:
            placed = False
            for o in union:
                if abs(z['re'] - o['re']) < 1e-8 and abs(z['im'] - o['im']) < 1e-8:
                    o['sources'].add(nm)
                    placed = True
                    break
            if not placed:
                union.append({'re': z['re'], 'im': z['im'], 'absH': z['absH'],
                              'status': z['status'], 'sources': {nm}})
    log('UNION (clustered 1e-8): %d unique zeros total' % len(union))

    # per-source presence map: how many union zeros does each scan see?
    presence = {}
    for nm in 'ABCD':
        zs = {(round(float(z['re']), 12), round(float(z['im']), 12))
              for z in cluster_sets[nm]}
        presence[nm] = sum(1 for u in union
                            if (round(float(u['re']), 12),
                                round(float(u['im']), 12)) in zs)
    log('union zeros seen by each scan: %s' % presence)

# S104 (20.125013, 19.250015) check
    S104 = (20.125013, 19.250015)
    log('S104 reference: %s' % str(S104))
    found_s104 = {}
    for nm in 'ABCD':
        best = (None, 1e9)
        for z in cluster_sets[nm]:
            d = math.hypot(float(z['re']) - S104[0], float(z['im']) - S104[1])
            if d < best[1]:
                best = (z, d)
        if best[0] is not None:
            found_s104[nm] = (best[1], best[0])
            log('  scan %s nearest to S104: dist=%.3e re=%.14s im=%.14s absH=%.1e status=%s'
                % (nm, best[1], str(best[0]['re'])[:14], str(best[0]['im'])[:14],
                   best[0]['absH'], best[0]['status']))

    # S104 neighbors within 0.06
    s104_neighbors = [u for u in union
                       if math.hypot(float(u['re']) - S104[0], float(u['im']) - S104[1]) < 0.06]
    log('S104 neighbors (within 0.06): %d' % len(s104_neighbors))
    for u in sorted(s104_neighbors, key=lambda x: math.hypot(float(x['re']) - S104[0], float(x['im']) - S104[1])):
        d_to = math.hypot(float(u['re']) - S104[0], float(u['im']) - S104[1])
        log('  d=%.4f re=%s im=%s absH=%.1e sources=%s status=%s'
            % (d_to, str(u['re'])[:14], str(u['im'])[:14],
               u['absH'], sorted(u['sources']), u['status']))

    # catalog cross-check: which catalog zeros in patch appear in scan A?
    log('--- catalog vs scan A ---')
    cat_set_in_Ascan = 0
    for cz in cat_patch:
        best = (None, 1e9)
        for z in cluster_sets['A']:
            d = math.hypot(float(z['re']) - cz['re'], float(z['im']) - cz['im'])
            if d < best[1]:
                best = (z, d)
        if best[0] is not None and best[1] < 1e-9:
            cat_set_in_Ascan += 1
            log('  cat re=%.10f im=%.10f -> scanA re=%s im=%s d=%.3e'
                % (cz['re'], cz['im'], str(best[0]['re'])[:14], str(best[0]['im'])[:14], best[1]))
        else:
            log('  cat re=%.10f im=%.10f NOT recovered by scan A (best d=%.3e)'
                % (cz['re'], cz['im'],
                   best[1] if best[0] is not None else float('inf')))
    log('catalog zeros in patch recovered by scan A: %d / %d'
        % (cat_set_in_Ascan, len(cat_patch)))

    # catalog completeness: |A| / |A union B union C union D|
    seen_in_A = set()
    for z in cluster_sets['A']:
        seen_in_A.add((round(float(z['re']), 8), round(float(z['im']), 8)))
    total_union = set()
    for nm in 'ABCD':
        for z in cluster_sets[nm]:
            total_union.add((round(float(z['re']), 8), round(float(z['im']), 8)))
    if total_union:
        completeness = len(seen_in_A) / len(total_union)
        log('catalog completeness (A) / (A U B U C U D) in patch = %d / %d = %.4f'
            % (len(seen_in_A), len(total_union), completeness))
    else:
        completeness = None

    # ---------------- save outputs ---------------- #
    out = dict(
        meta=dict(
            patch=PATCH,
            level=M_LEVEL,
            mp_dps=mp.dps,
            polish_tol='1e-30',
            cluster_eps='1e-6 for cross, 1e-8 within scan',
            delta6=float(DELTA6_MP),
            AU0_6=hz['AU0_6'],
            elapsed_sec=time.time() - T0,
        ),
        catalog_patch=[dict(re=z['re'], im=z['im'], absH=z['absH'], str=z.get('str'))
                       for z in cat_patch],
        scans={
            'A': dict(step='1/16', raw=zA, clustered=[
                dict(re=z['re'], im=z['im'], absH=z['absH'], status=z['status'],
                     count=z['count']) for z in cA]),
            'B': dict(step='1/24', raw=zB, clustered=[
                dict(re=z['re'], im=z['im'], absH=z['absH'], status=z['status'],
                     count=z['count']) for z in cB]),
            'C': dict(step='1/10', raw=zC, clustered=[
                dict(re=z['re'], im=z['im'], absH=z['absH'], status=z['status'],
                     count=z['count']) for z in cC]),
            'D': dict(step='random2000', clustered=[
                dict(re=z['re'], im=z['im'], absH=z['absH'], status=z['status'],
                     count=z['count']) for z in cD]),
        },
        stats=stats,
        intersections={'%s_%s' % k: dict(common=v['common'], only_n1=v['only_n1'],
                                          only_n2=v['only_n2'])
                       for k, v in inter.items()},
        union=[dict(re=u['re'], im=u['im'], absH=u['absH'],
                    status=u['status'], sources=sorted(u['sources']))
               for u in union],
        presence_in_union=presence,
        s104_neighbors=[dict(re=u['re'], im=u['im'], absH=u['absH'],
                             sources=sorted(u['sources']),
                             d=math.hypot(u['re'] - S104[0], u['im'] - S104[1]))
                        for u in s104_neighbors],
        completeness_catalog=completeness,
        found_s104_in_each_scan={nm: dict(dist=v[0], re=v[1]['re'], im=v[1]['im'],
                                          absH=v[1]['absH'])
                                 for nm, v in found_s104.items()},
    )
    with open(os.path.join(HERE, 'rescan_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote rescan_results.json (elapsed %.1fs)' % (time.time() - T0))


if __name__ == '__main__':
    main()

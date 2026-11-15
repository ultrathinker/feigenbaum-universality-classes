#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1_d4_rescan.py -- R1 lattice-rescan driver for d=4.

Reuses the float p2_H4_zeros.py pipeline (H, newton_refine, domain_check,
mp-polish, scan_level) and only swaps the GRID generator and the OUTPUT.

Modes:
  sanity    : original 1/16 grid, no shift, R_max=20 (reproducibility check)
  step24    : 1/24 grid + random shift (seed=20260816), R_max=20
  random    : N random starts uniform in [-20,20] x [0,20] (seed=20260816)

Output: H4_<mode>.json with `zeros_level8` (only level 8, no level-7 pairing).
Dedup |deltau|<1e-8 per BRIEF (vs 1e-4 in original).  All other criteria as in
p2_H4_zeros.py.

Deviations from the original pipeline (called out explicitly):
  (1) grid generator + (2) dedup tolerance 1e-8 (was 1e-4) + (3) no level-7 pairing.
"""
import json
import math
import os
import sys
import time

import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 50

# Import the original pipeline as-is (without running its main()).
import p2_H4_zeros as p2

SEED = 20260816
R_MAX = 20.0
STEP16 = 0.0625
STEP24 = 1.0 / 24.0
TAG = "[oid:r1H4] R1 lattice-rescan driver (d=4)"

T0 = time.time()
LOGF = open('r1_d4_rescan.log', 'w', encoding='utf-8')


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    # force ASCII-safe output on Windows cp1252 console
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass
    try:
        print(line, flush=True)
    except UnicodeEncodeError:
        print(line.encode('ascii', 'replace').decode('ascii'), flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


# --------------------------------------------------------------------------- #
#  AU0 calibration -- copy of p2_H4_zeros.main() AU0 block (NOT modified)
# --------------------------------------------------------------------------- #
AU0_CACHE = 'H4_AU0_cache.json'


def calibrate_au0():
    """AU0 calibration -- same algorithm as p2_H4_zeros.main."""
    if os.path.exists(AU0_CACHE):
        cache = json.load(open(AU0_CACHE, encoding='utf-8'))
        AU0 = mpf(cache['AU0_4'])
        Ls = mpf(cache['Ls'])
        s_stars = {int(k): mpf(v) for k, v in cache['s_star'].items()}
        log('AU0_4 cache hit: AU0_4 = %s (gauge check skipped)' % nstr(AU0, 30))
        return AU0, Ls, s_stars
    log('--- AU0_4 calibration (copy of p2_H4_zeros.main) ---')
    smax = 12
    s_stars = {}
    s_prev = None
    for m in range(2, smax + 1):
        n = 2 ** (m + 1)

        def fres(s, m=m, n=n):
            x = mpf(0)
            ds = mpf(0)
            for _ in range(n):
                ds = (p2.g4p(x) + s * p2.h4p(x)) * ds + p2.h4(x)
                x = p2.g4(x) + s * p2.h4(x)
            return x, ds

        s0 = s_prev / p2.DELTA4 if s_prev is not None else mpf('-1e-3')
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
                x = p2.g4(x) + s * p2.h4(x)
            if fabs(x) < mpf('1e-30'):
                raise RuntimeError('inherited root at m=%d j=%d' % (j, m))
        s_stars[m] = s
        s_prev = s
        log('m=%2d  s*_m = %s   Ls^(m) = %s'
            % (m, nstr(s, 20), nstr(s * p2.DELTA4 ** m, 20)))
    rL = 1 / p2.DELTA4
    Ls_seq = {m: s_stars[m] * p2.DELTA4 ** m for m in s_stars}
    Ls = (Ls_seq[smax] - rL * Ls_seq[smax - 1]) / (1 - rL)
    AU0 = -Ls
    log('Ls = %s   AU0_4 = -Ls = %s' % (nstr(Ls, 30), nstr(AU0, 30)))
    for m in (6, 8, 10):
        x2m, x2m1 = p2.orbit_H_mpmath(mpf(-1), m, AU0)
        log('gauge check H_%d(-1) = %s' % (m, nstr(x2m1 / x2m, 12)))
    # save to cache
    cache = {'AU0_4': nstr(AU0, 40), 'Ls': nstr(Ls, 40),
             's_star': {str(m): nstr(s_stars[m], 30) for m in s_stars}}
    with open(AU0_CACHE, 'w', encoding='utf-8') as f:
        json.dump(cache, f, indent=2)
    log('cached AU0_4 -> %s' % AU0_CACHE)
    return AU0, Ls, s_stars


# --------------------------------------------------------------------------- #
#  mp-polish -- copy of p2_H4_zeros.mp_polish (was nested inside main())
# --------------------------------------------------------------------------- #
def mp_polish(u_start, m, AU0, iters=14, tol_step=mpf('1e-24')):
    u_mp = mpc(mpf(u_start.real), mpf(u_start.imag))
    for _ in range(iters):
        x2m, x2m1, dx2m, dx2m1 = p2.orbit_H_mpmath(u_mp, m, AU0, want_grad=True)
        H = x2m1 / x2m
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        step = H / dH
        u_mp = u_mp - step
        if abs(step) < tol_step:
            break
    x2m, x2m1 = p2.orbit_H_mpmath(u_mp, m, AU0)
    return u_mp, abs(x2m1 / x2m)


# --------------------------------------------------------------------------- #
#  grid generators
# --------------------------------------------------------------------------- #
def grid_orig(R, step, shift_re=0.0, shift_im=0.0):
    """The original 1/16 grid (no shift, R_max=20)."""
    xs = np.arange(-R, R + step / 2, step)
    ys = np.arange(0.0, R + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return UU + 1j * VV, (0.0, 0.0)


def grid_step24(R, step, seed=SEED):
    """1/24 grid + random shift in [0, step) per coordinate."""
    rng = np.random.default_rng(seed)
    sh_re = float(rng.uniform(0.0, step))
    sh_im = float(rng.uniform(0.0, step))
    xs = np.arange(-R + sh_re, R + step / 2, step)
    ys = np.arange(0.0 + sh_im, R + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return UU + 1j * VV, (sh_re, sh_im)


def random_starts(n, R, seed=SEED):
    """n random points uniform in [-R, R] x [0, R]."""
    rng = np.random.default_rng(seed)
    pts = np.empty(n, dtype=np.complex128)
    pts.real = rng.uniform(-R, R, n)
    pts.imag = rng.uniform(0.0, R, n)
    return pts, (None, None)


# --------------------------------------------------------------------------- #
#  scan + local minima on a (possibly shifted) grid
# --------------------------------------------------------------------------- #
def find_local_minima(grid):
    """Return complex values of cells that are local minima (≤ all 8 neighbours)."""
    if grid.ndim == 1:
        # arbitrary 1D points -- no local-min notion; treat all as candidates
        return [complex(z) for z in grid]
    Hm = grid  # callers pass |H| already
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    ii, jj = np.nonzero(ismin)
    return ii, jj


def run(mode, R=R_MAX, n_random=None):
    log('=== R1 d=4 rescan: mode=%s ===' % mode)
    AU0, Ls, s_stars = calibrate_au0()
    AU0s = nstr(AU0, 40)

    if mode == 'sanity':
        grid, shift = grid_orig(R, STEP16)
        step_used = STEP16
        seed_used = None
    elif mode == 'step24':
        grid, shift = grid_step24(R, STEP24)
        step_used = STEP24
        seed_used = SEED
    elif mode == 'random':
        # N = number of 1/16 nodes in same rectangle (per BRIEF)
        n_nodes = int((2 * R / STEP16)) * int((R / STEP16))
        N = n_random if n_random is not None else n_nodes
        grid, shift = random_starts(N, R)
        step_used = None
        seed_used = SEED
    else:
        raise ValueError('unknown mode %s' % mode)

    log('mode=%s: R=%.1f step=%s shift=(%.6f,%.6f) seed=%s'
        % (mode, R, step_used,
           (shift[0] if shift[0] is not None else -1),
           (shift[1] if shift[1] is not None else -1), seed_used))

    # ---- float scan ----------------------------------------------------------
    t_scan = time.time()
    if grid.ndim == 2:
        # grid-point scan: compute |H|, take local minima
        Habs, valid, maxt, errmax = p2.scan_level(8, AU0, grid)
        n_valid = int(valid.sum())
        log('scan m=8: grid %dx%d (%.0fs), valid %d (%.1f%%), min|H| on valid = %.3e'
            % (grid.shape[0], grid.shape[1], time.time() - t_scan, n_valid,
               100.0 * n_valid / grid.size,
               float(Habs[valid].min()) if n_valid else float('nan')))
        Hm = np.where(valid, Habs, np.inf)
        ii, jj = find_local_minima(Hm)
        cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
        log('  candidates: %d' % len(cand))
    else:
        # random starts: each start is a candidate (no local-min filtering)
        cand = [complex(z) for z in grid]
        log('  random starts (no local-min filter): %d' % len(cand))

    # ---- Newton refinement + domain check (level 8 only) -------------------
    zeros = []
    t_newt = time.time()
    for i2, u0 in enumerate(cand):
        try:
            u, res = p2.newton_refine(u0, 8, AU0, tol=1e-9, rmax=R + 5)
        except (OverflowError, ZeroDivisionError, ValueError):
            u, res = None, float('inf')
        if u is None or res > 1e-5:
            continue
        mt4, merr, den = p2.domain_check(u, 8, AU0)
        if mt4 < 2.2 and merr < 1e-11 and den > 1e-12:
            zeros.append({'u': u, 'absH': res, 'maxt4': mt4, 'maxerr': merr,
                          'den': den})
        if (i2 + 1) % 50 == 0:
            el = time.time() - t_newt
            eta = el / (i2 + 1) * (len(cand) - i2 - 1)
            log('  newton %d/%d (%.0fs, ETA %.0fs)' % (i2 + 1, len(cand), el, eta))
    log('  refined zeros (|H|<1e-5, domain-ok): %d' % len(zeros))

    # ---- mp polish (level 8) ------------------------------------------------
    cat8_mp = []
    t_pol = time.time()
    for i2, z in enumerate(zeros):
        try:
            u_mp, r = mp_polish(z['u'], 8, AU0)
        except Exception as exc:  # noqa: BLE001
            log('  polish failed near %.4f%+.4fi: %s' % (z['u'].real, z['u'].imag, exc))
            continue
        if float(r) < 1e-20:
            cat8_mp.append({'u8_mp': u_mp, 'absH8': float(r),
                            'u8_seed': z['u'], 'maxt4': z['maxt4'],
                            'maxerr': z['maxerr']})
        if (i2 + 1) % 50 == 0:
            el = time.time() - t_pol
            eta = el / (i2 + 1) * (len(zeros) - i2 - 1)
            log('  polish %d/%d (%.0fs, ETA %.0fs)' % (i2 + 1, len(zeros), el, eta))
    log('polished level-8 zeros (|H8|<1e-20): %d (%.0fs)' % (len(cat8_mp), time.time() - t_pol))

    # ---- dedup |deltau|<1e-8 ----------------------------------------------------
    dedup = []
    for z in cat8_mp:
        u8c = complex(mpre(z['u8_mp']), mpim(z['u8_mp']))
        if all(abs(u8c - v['u8c']) > 1e-8 for v in dedup):
            dedup.append({**z, 'u8c': u8c})
    log('dedup (|deltau|<1e-8): %d' % len(dedup))

    # ---- write JSON ---------------------------------------------------------
    out = {
        'meta': {'tag': TAG, 'mode': mode, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': ('article 2.2 analogue (d=4); uses p2_H4_zeros newton_refine, '
                            'domain_check, mp_polish as-is; differs in grid + dedup 1e-8'),
                 'scan': {'R_max': R, 'step': step_used,
                          'shift_re': shift[0], 'shift_im': shift[1],
                          'seed': seed_used, 'level': 8},
                 'pipeline_deviations': ['grid (this mode)', 'dedup 1e-8 (was 1e-4)',
                                         'no level-7 pairing (per BRIEF)'],
                 'mp_dps': mp.dps, 'elapsed_sec': time.time() - T0},
        'AU0_4': AU0s,
        'Ls': nstr(Ls, 40),
        'zeros_level8': [{'re': z['u8c'].real, 'im': z['u8c'].imag,
                          'absH': z['absH8'], 'str': nstr(z['u8_mp'], 25),
                          'maxt4': z['maxt4'], 'maxerr': z['maxerr']} for z in dedup],
        'counts': {'candidates': len(cand),
                   'refined': len(zeros),
                   'polished': len(cat8_mp),
                   'dedup': len(dedup)},
    }
    out_path = 'H4_%s.json' % mode
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote %s (n=%d)' % (out_path, len(dedup)))
    return out


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sanity'
    R = float(sys.argv[2]) if len(sys.argv) > 2 else R_MAX
    n_random = int(sys.argv[3]) if len(sys.argv) > 3 else None
    run(mode, R=R, n_random=n_random)
    log('DONE total %.1fs' % (time.time() - T0))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1_d6_rescan.py -- R1 lattice-rescan driver for d=6.

Reuses p2_H6_zeros.py pipeline (H, scan_level, newton_refine, domain_check,
mp_polish_job) and only swaps the GRID generator and the OUTPUT.

Modes:
  random    : 1/3 * num 1/16 nodes in R_max=48 (Im>=0), seed=20260816
  step24    : step 1/24 + random shift in [0,1/24), seed=20260816, R_max=48

Level 8 only (no level-7 pairing).  Dedup |deltau|<1e-8 per BRIEF.

Background log: progress_d6_<mode>.log with line every ~2% + ETA.
Checkpoints: H6_<mode>_ckpt.json every 5 minutes.

Deviations from original pipeline:
  (1) grid generator (no coarse R=48 step-0.25 pre-scan)
  (2) dedup 1e-8 (was 1e-4)
  (3) no level-7 pairing
  (4) no level-9 drift/refinement (limit extrapolation NOT done)
  (5) no late-9 paired survival filter
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 50

import p2_H6_zeros as p2

SEED = 20260816
R_MAX = 48.0
STEP16 = 0.0625
STEP24 = 1.0 / 24.0
MAX_WORKERS = 8  # per BRIEF cap
TAG = "[oid:r1H6] R1 lattice-rescan driver (d=6)"

T0 = time.time()
LOGF = None  # opened in main()


def log(msg):
    global LOGF
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass
    try:
        print(line, flush=True)
    except UnicodeEncodeError:
        print(line.encode('ascii', 'replace').decode('ascii'), flush=True)
    if LOGF is not None:
        LOGF.write(line + "\n")
        LOGF.flush()


# --------------------------------------------------------------------------- #
#  AU0 calibration (copy of p2_H6_zeros.main) -- not modified
# --------------------------------------------------------------------------- #
AU0_CACHE = 'H6_AU0_cache.json'


def calibrate_au0():
    if os.path.exists(AU0_CACHE):
        cache = json.load(open(AU0_CACHE, encoding='utf-8'))
        AU0 = mpf(cache['AU0_6'])
        Ls = mpf(cache['Ls'])
        s_stars = {int(k): mpf(v) for k, v in cache['s_star'].items()}
        log('AU0_6 cache hit: AU0_6 = %s (gauge check skipped)' % nstr(AU0, 30))
        return AU0, Ls, s_stars
    log('--- AU0_6 calibration (copy of p2_H6_zeros.main) ---')
    smax = 12
    s_stars = {}
    s_prev = None
    for m in range(2, smax + 1):
        n = 2 ** (m + 1)

        def fres(s, m=m, n=n):
            x = mpf(0)
            ds = mpf(0)
            for _ in range(n):
                ds = (p2.g6p(x) + s * p2.h6p(x)) * ds + p2.h6(x)
                x = p2.g6(x) + s * p2.h6(x)
            return x, ds

        s0 = s_prev / p2.DELTA6 if s_prev is not None else mpf('-1e-3')
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
                x = p2.g6(x) + s * p2.h6(x)
            if fabs(x) < mpf('1e-30'):
                raise RuntimeError('inherited root at m=%d j=%d' % (j, m))
        s_stars[m] = s
        s_prev = s
        log('m=%2d  s*_m = %s   Ls^(m) = %s'
            % (m, nstr(s, 20), nstr(s * p2.DELTA6 ** m, 20)))
    rL = 1 / p2.DELTA6
    Ls_seq = {m: s_stars[m] * p2.DELTA6 ** m for m in s_stars}
    Ls = (Ls_seq[smax] - rL * Ls_seq[smax - 1]) / (1 - rL)
    AU0 = -Ls
    log('Ls = %s   AU0_6 = -Ls = %s' % (nstr(Ls, 30), nstr(AU0, 30)))
    for m in (6, 8, 10):
        x2m, x2m1 = p2.orbit_H_mpmath(mpf(-1), m, AU0)
        log('gauge check H_%d(-1) = %s' % (m, nstr(x2m1 / x2m, 12)))
    cache = {'AU0_6': nstr(AU0, 40), 'Ls': nstr(Ls, 40),
             's_star': {str(m): nstr(s_stars[m], 30) for m in s_stars}}
    with open(AU0_CACHE, 'w', encoding='utf-8') as f:
        json.dump(cache, f, indent=2)
    log('cached AU0_6 -> %s' % AU0_CACHE)
    return AU0, Ls, s_stars


# --------------------------------------------------------------------------- #
#  grid generators
# --------------------------------------------------------------------------- #
def grid_step24(R, step, seed=SEED):
    rng = np.random.default_rng(seed)
    sh_re = float(rng.uniform(0.0, step))
    sh_im = float(rng.uniform(0.0, step))
    xs = np.arange(-R + sh_re, R + step / 2, step)
    ys = np.arange(0.0 + sh_im, R + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return UU + 1j * VV, (sh_re, sh_im)


def random_starts(n, R, seed=SEED):
    rng = np.random.default_rng(seed)
    pts = np.empty(n, dtype=np.complex128)
    pts.real = rng.uniform(-R, R, n)
    pts.imag = rng.uniform(0.0, R, n)
    return pts, (None, None)


def find_local_minima(grid_arr):
    Hm = grid_arr
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    ii, jj = np.nonzero(ismin)
    return ii, jj


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def run(mode, R=R_MAX, n_random=None, ckpt_seconds=300):
    global LOGF
    log_path = 'progress_d6_%s.log' % mode
    LOGF = open(log_path, 'w', encoding='utf-8')
    log('=== R1 d=6 rescan: mode=%s ===' % mode)
    AU0, Ls, s_stars = calibrate_au0()
    AU0s = nstr(AU0, 40)

    if mode == 'step24':
        grid, shift = grid_step24(R, STEP24)
        step_used = STEP24
        seed_used = SEED
    elif mode == 'random':
        # N = 1/3 * num 1/16 nodes in R_max=48, Im>=0 (per BRIEF)
        n_nodes = int((2 * R / STEP16)) * int((R / STEP16))
        N = n_random if n_random is not None else max(1, n_nodes // 3)
        grid, shift = random_starts(N, R)
        step_used = None
        seed_used = SEED
    else:
        raise ValueError('unknown mode %s' % mode)

    log('mode=%s: R=%.1f step=%s shift=(%.6f,%.6f) seed=%s'
        % (mode, R, step_used, (shift[0] if shift[0] is not None else -1),
           (shift[1] if shift[1] is not None else -1), seed_used))

    # ---- float scan ----------------------------------------------------------
    t_scan = time.time()
    if grid.ndim == 2:
        Habs, valid, maxt, errmax = p2.scan_level(8, AU0, grid)
        n_valid = int(valid.sum())
        log('scan m=8: grid %dx%d (%.0fs), valid %d (%.1f%%), min|H| on valid = %.3e'
            % (grid.shape[0], grid.shape[1], time.time() - t_scan, n_valid,
               100.0 * n_valid / grid.size,
               float(Habs[valid].min()) if n_valid else float('nan')))
        Hm = np.where(valid, Habs, np.inf)
        ii, jj = find_local_minima(Hm)
        cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
        log('  candidates (local minima): %d' % len(cand))
    else:
        cand = [complex(z) for z in grid]
        log('  random starts (no local-min filter): %d' % len(cand))

    # ---- Newton refinement + domain check (level 8 only) -------------------
    zeros = []
    t_newt = time.time()
    last_prog = -1
    last_ckpt = time.time()
    n_total = len(cand)
    for i2, u0 in enumerate(cand):
        try:
            u, res = p2.newton_refine(u0, 8, AU0, tol=1e-9, rmax=R + 5)
        except (OverflowError, ZeroDivisionError, ValueError):
            u, res = None, float('inf')
        if u is not None and res <= 1e-5:
            mt6, merr, den = p2.domain_check(u, 8, AU0)
            if mt6 < p2.TMAX6 and den > 1e-12:
                zeros.append({'u': u, 'absH': res, 'maxt6': mt6, 'maxerr': merr,
                              'den': den})
        # progress log every ~2%
        pct = int(100 * (i2 + 1) / n_total)
        if pct >= last_prog + 2:
            el = time.time() - t_newt
            eta = el / (i2 + 1) * (n_total - i2 - 1)
            log('  newton %d/%d (%d%%, %.0fs elapsed, ETA %.0fs), refined+ok: %d'
                % (i2 + 1, n_total, pct, el, eta, len(zeros)))
            last_prog = pct
        # checkpoint every ckpt_seconds
        if time.time() - last_ckpt > ckpt_seconds and zeros:
            write_ckpt(mode, zeros, i2 + 1, n_total, AU0s, Ls, s_stars,
                       step_used, shift, seed_used, R, tag='newton')
            last_ckpt = time.time()
    log('  refined zeros (|H|<1e-5, domain-ok): %d' % len(zeros))

    # ---- mp polish (parallel) ------------------------------------------------
    jobs = [(z['u'].real, z['u'].imag, 8, AU0s, 16) for z in zeros]
    cat8_mp = []
    seeds = list(zeros)
    t_pol = time.time()
    last_prog = -1
    last_ckpt = time.time()
    n_jobs = len(jobs)
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i2, reslt in enumerate(ex.map(p2.mp_polish_job, jobs)):
            if reslt is not None:
                ur, ui, r, _m = reslt
                cat8_mp.append({'u8_mp_re': ur, 'u8_mp_im': ui, 'absH8': r,
                                'seed': seeds[i2]})
            pct = int(100 * (i2 + 1) / n_jobs) if n_jobs else 100
            if pct >= last_prog + 2:
                el = time.time() - t_pol
                eta = el / (i2 + 1) * (n_jobs - i2 - 1) if (i2 + 1) else 0
                log('  polish %d/%d (%d%%, %.0fs elapsed, ETA %.0fs), ok: %d'
                    % (i2 + 1, n_jobs, pct, el, eta, len(cat8_mp)))
                last_prog = pct
            if time.time() - last_ckpt > ckpt_seconds and cat8_mp:
                write_ckpt(mode, zeros, i2 + 1, n_jobs, AU0s, Ls, s_stars,
                           step_used, shift, seed_used, R, tag='polish',
                           partial_mp=cat8_mp)
                last_ckpt = time.time()
    log('polished level-8 zeros (|H8|<1e-20): %d (%.0fs)' % (len(cat8_mp), time.time() - t_pol))

    # ---- dedup |deltau|<1e-8 --------------------------------------------------
    dedup = []
    for z in cat8_mp:
        u8c = complex(mpf(z['u8_mp_re']), mpf(z['u8_mp_im']))
        if all(abs(u8c - v['u8c']) > 1e-8 for v in dedup):
            dedup.append({**z, 'u8c': u8c})
    log('dedup (|deltau|<1e-8): %d' % len(dedup))

    # ---- write JSON ----------------------------------------------------------
    out = {
        'meta': {'tag': TAG, 'mode': mode, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': ('article 2.2 analogue (d=6); uses p2_H6_zeros newton_refine, '
                            'domain_check, mp_polish_job as-is; differs in grid + dedup 1e-8'),
                 'scan': {'R_max': R, 'step': step_used,
                          'shift_re': shift[0], 'shift_im': shift[1],
                          'seed': seed_used, 'level': 8},
                 'pipeline_deviations': ['grid (this mode)', 'dedup 1e-8 (was 1e-4)',
                                         'no level-7 pairing', 'no level-9 polish',
                                         'no late-9 survival filter'],
                 'mp_dps': mp.dps, 'max_workers': MAX_WORKERS,
                 'elapsed_sec': time.time() - T0},
        'AU0_6': AU0s,
        'Ls': nstr(Ls, 40),
        'zeros_level8': [{'re': z['u8c'].real, 'im': z['u8c'].imag,
                          'absH': z['absH8'], 'str': (z['u8_mp_re'] + '+' + z['u8_mp_im'] + 'i'),
                          'maxt6': z['seed']['maxt6'], 'maxerr': z['seed']['maxerr']}
                         for z in dedup],
        'counts': {'candidates': len(cand),
                   'refined': len(zeros),
                   'polished': len(cat8_mp),
                   'dedup': len(dedup)},
    }
    out_path = 'H6_%s.json' % mode
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote %s (n=%d)' % (out_path, len(dedup)))
    LOGF.close()
    return out


def write_ckpt(mode, zeros, n_done, n_total, AU0s, Ls, s_stars,
               step, shift, seed, R, tag='newton', partial_mp=None):
    out = {
        'meta': {'tag': TAG, 'mode': mode, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': 'ckpt mid-run', 'tag_stage': tag,
                 'scan': {'R_max': R, 'step': step,
                          'shift_re': shift[0], 'shift_im': shift[1],
                          'seed': seed, 'level': 8}},
        'AU0_6': AU0s,
        'Ls': nstr(Ls, 40),
        'progress': {'n_done': n_done, 'n_total': n_total,
                     'n_zeros_so_far': len(zeros)},
        'zeros_so_far': [{'re': z['u'].real, 'im': z['u'].imag,
                          'absH': z['absH'], 'maxt6': z['maxt6'],
                          'maxerr': z['maxerr']} for z in zeros],
    }
    if partial_mp is not None:
        out['mp_so_far'] = [{'re': complex(mpf(z['u8_mp_re'])).real,
                             'im': complex(mpf(z['u8_mp_im'])).imag,
                             'absH': z['absH8']} for z in partial_mp]
    path = 'H6_%s_ckpt.json' % mode
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('  ckpt @ %d/%d -> %s' % (n_done, n_total, path))


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'random'
    R = float(sys.argv[2]) if len(sys.argv) > 2 else R_MAX
    n_random = int(sys.argv[3]) if len(sys.argv) > 3 else None
    run(mode, R=R, n_random=n_random)
    log('DONE total %.1fs' % (time.time() - T0))

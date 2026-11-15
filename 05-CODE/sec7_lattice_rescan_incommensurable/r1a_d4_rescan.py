#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_d4_rescan.py -- Task 5 control: 'shift16' mode.  Copy of
r1_d4_rescan.py's driver (grid + Newton + domain_check + mp-polish + dedup),
extended with ONE new mode:

  shift16 : 1/16 grid + random shift in [0, 1/16) per coordinate
            (seed=20260817, DIFFERENT seed from the step24 run's 20260816),
            same R_max=20 region, level 8 only, dedup 1e-8 -- otherwise
            identical machinery to r1_d4_rescan.py's step24 mode.

Rationale (BRIEF item 5): if zeros found by a grid scan cluster near the
grid's OWN (possibly shifted) nodes regardless of where those nodes sit, that
demonstrates the clustering is an artifact of "start near a node -> Newton
converges to the nearest attracting root -> looks like a zero at the node",
not a real lattice property of H_4.  This driver reuses newton_refine,
domain_check, mp_polish, scan_level from p2_H4_zeros.py UNCHANGED, exactly
like r1_d4_rescan.py -- only the grid generator (shift16) is new.

Note on rmax: uses r1a_common.newton_from, whose default rmax=22.0 matches
p2_H4_zeros.py's OWN default (r1_d4_rescan.py instead used rmax=R+5=25.0 for
step24/sanity, an undeclared deviation -- see r1a_common.py docstring and
REPORT-R1A-D4-FORENSICS.md). We deliberately use p2's true default here so
shift16 is comparable to the ORIGINAL pipeline's reach, not to step24's.

Output: r1a_H4_shift16.json (same schema as H4_step24.json).
"""
import sys
import time

import numpy as np
from mpmath import mp, nstr, re as mpre, im as mpim

mp.dps = 50

import r1a_common as C

SEED16 = 20260817
STEP16 = C.STEP16
R_MAX = C.R_MAX

log = C.mk_logger('r1a_d4_rescan_shift16.log')


def grid_shift16(R, step, seed=SEED16):
    rng = np.random.default_rng(seed)
    sh_re = float(rng.uniform(0.0, step))
    sh_im = float(rng.uniform(0.0, step))
    xs = np.arange(-R + sh_re, R + step / 2, step)
    ys = np.arange(0.0 + sh_im, R + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return UU + 1j * VV, (sh_re, sh_im)


def find_local_minima(Hm):
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    return np.nonzero(ismin)


def run():
    log('=== r1a d=4 rescan: mode=shift16 (seed=%d) ===' % SEED16)
    AU0, Ls = C.load_au0()
    log('AU0_4 loaded from cache')

    grid, shift = grid_shift16(R_MAX, STEP16, seed=SEED16)
    log('grid=%s shift=(%.6f,%.6f)' % (grid.shape, shift[0], shift[1]))

    t_scan = time.time()
    Habs, valid, maxt, errmax = C.p2.scan_level(8, AU0, grid)
    n_valid = int(valid.sum())
    log('scan m=8: %.0fs, valid %d/%d (%.1f%%)'
        % (time.time() - t_scan, n_valid, grid.size, 100.0 * n_valid / grid.size))
    Hm = np.where(valid, Habs, np.inf)
    ii, jj = find_local_minima(Hm)
    cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
    log('candidates: %d' % len(cand))

    zeros = []
    t_newt = time.time()
    for i2, u0 in enumerate(cand):
        u, res = C.newton_from(u0, 8, AU0)  # p2's own default rmax=22.0
        if u is None or res > 1e-5:
            continue
        mt4, merr, den = C.p2.domain_check(u, 8, AU0)
        if mt4 < 2.2 and merr < 1e-11 and den > 1e-12:
            zeros.append({'u': u, 'absH': res, 'maxt4': mt4, 'maxerr': merr, 'den': den})
        if (i2 + 1) % 50 == 0:
            el = time.time() - t_newt
            eta = el / (i2 + 1) * (len(cand) - i2 - 1)
            log('  newton %d/%d (%.0fs, ETA %.0fs)' % (i2 + 1, len(cand), el, eta))
    log('refined zeros (domain-ok): %d' % len(zeros))

    cat8_mp = []
    t_pol = time.time()
    for i2, z in enumerate(zeros):
        try:
            u_mp, r = C.mp_polish(z['u'], 8, AU0)
        except Exception as exc:  # noqa: BLE001
            log('  polish failed near %.4f%+.4fi: %s' % (z['u'].real, z['u'].imag, exc))
            continue
        if float(r) < 1e-20:
            cat8_mp.append({'u8_mp': u_mp, 'absH8': float(r), 'u8_seed': z['u'],
                            'maxt4': z['maxt4'], 'maxerr': z['maxerr']})
        if (i2 + 1) % 50 == 0:
            log('  polish %d/%d (%.0fs)' % (i2 + 1, len(zeros), time.time() - t_pol))
    log('polished (|H8|<1e-20): %d' % len(cat8_mp))

    dedup = []
    for z in cat8_mp:
        u8c = complex(mpre(z['u8_mp']), mpim(z['u8_mp']))
        if all(abs(u8c - v['u8c']) > 1e-8 for v in dedup):
            dedup.append({**z, 'u8c': u8c})
    log('dedup |du|<1e-8: %d' % len(dedup))

    out = {
        'meta': {'tag': '[oid:r1aH4shift16] r1a shift16 control (d=4)', 'mode': 'shift16',
                 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'recipe': ('article 2.2 analogue (d=4); p2_H4_zeros newton_refine/domain_check/'
                            'mp_polish/scan_level used unmodified; rmax=22.0 (p2 default, NOT '
                            'r1_d4_rescan.py\'s 25.0); differs from p2 only in grid (shifted 1/16) '
                            'and dedup 1e-8 (was 1e-4) and no level-7 pairing'),
                 'scan': {'R_max': R_MAX, 'step': STEP16, 'shift_re': shift[0], 'shift_im': shift[1],
                          'seed': SEED16, 'level': 8, 'newton_rmax': C.P2_DEFAULT_RMAX},
                 'pipeline_deviations': ['grid (shifted 1/16, seed 20260817)',
                                         'dedup 1e-8 (was 1e-4)', 'no level-7 pairing'],
                 'mp_dps': mp.dps},
        'AU0_4': nstr(AU0, 40),
        'zeros_level8': [{'re': z['u8c'].real, 'im': z['u8c'].imag, 'absH': z['absH8'],
                          'str': nstr(z['u8_mp'], 25), 'maxt4': z['maxt4'], 'maxerr': z['maxerr']}
                         for z in dedup],
        'counts': {'candidates': len(cand), 'refined': len(zeros),
                   'polished': len(cat8_mp), 'dedup': len(dedup)},
    }
    C.dump_json(out, 'r1a_H4_shift16.json')
    log('wrote r1a_H4_shift16.json (n=%d)' % len(dedup))


if __name__ == '__main__':
    run()
    log('DONE')

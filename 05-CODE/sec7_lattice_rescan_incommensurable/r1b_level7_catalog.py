#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_level7_catalog.py -- reproduce the level-7 (m=7) zero catalog for d=6
EXACTLY as p2_H6_zeros.main()'s fine-scan stage does (grid 1/16, R_max=48,
unshifted): grid scan -> local minima -> newton_refine(tol=1e-9, p2's own
default rmax=70.0) -> domain_check(maxt6<1.15, den>1e-12; NO merr threshold,
matching p2_H6/r1_d6's own domain_ok rule) -> dedup(|du|<1e-4).

Needed because r1_d6_rescan.py (STEP24 driver) skips level-7 pairing, so its
641 zeros have no recorded m=7 partner.  Sanity check at the end: match
against the u7 already stored in H6_zeros_ORIG.json ('zeros_level7', paired
1:1 with 'zeros_level8') -- should reproduce all 342 to high precision.

Output: r1b_H6_level7.json
"""
import json
import time

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_level7_catalog.log')


def main():
    log('=== r1b level-7 catalog (grid 1/16, R=48, unshifted) ===')
    AU0, Ls = C.load_au0()
    log('AU0_6 loaded from cache')

    grid = C.grid_orig(C.R_MAX, C.STEP16)
    log('grid built: %s (%d points)' % (grid.shape, grid.size))

    t0 = time.time()
    Habs, valid, maxt, errmax = C.p2.scan_level(7, AU0, grid)
    n_valid = int(valid.sum())
    log('scan_level(7) done in %.1fs; valid %d/%d (%.1f%%)'
        % (time.time() - t0, n_valid, grid.size, 100.0 * n_valid / grid.size))

    Hm = np.where(valid, Habs, np.inf)
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < 0.2)
    ii, jj = np.nonzero(ismin)
    cand = [complex(grid[i, j]) for i, j in zip(ii, jj)]
    log('candidates (local minima): %d' % len(cand))

    zeros = []
    t1 = time.time()
    for i2, u0 in enumerate(cand):
        u, res = C.newton_from(u0, 7, AU0)  # p2's own default rmax=70.0
        if u is None or res > 1e-5:
            continue
        ok, dom = C.domain_ok(u, 7, AU0)
        zeros.append({'u': u, 'absH': res, 'domain_ok': ok, **dom})
        if (i2 + 1) % 100 == 0:
            el = time.time() - t1
            eta = el / (i2 + 1) * (len(cand) - i2 - 1)
            log('  newton %d/%d (%.0fs, ETA %.0fs)' % (i2 + 1, len(cand), el, eta))
    log('refined (|H|<1e-5): %d' % len(zeros))

    dedup = []
    for z in zeros:
        if all(abs(z['u'] - v['u']) > 1e-4 for v in dedup):
            dedup.append(z)
    log('dedup |du|<1e-4: %d' % len(dedup))

    cat7 = [z for z in dedup if z['domain_ok']]
    log('domain_ok subset (== p2 cat7): %d' % len(cat7))

    out = {
        'meta': {'tag': 'r1b level-7 catalog, reproduces p2_H6_zeros m=7 fine-scan stage',
                 'grid': {'R_max': C.R_MAX, 'step': C.STEP16, 'shift': [0.0, 0.0]},
                 'newton_rmax': C.P2_DEFAULT_RMAX,
                 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'elapsed_sec': time.time() - t0},
        'counts': {'candidates': len(cand), 'refined': len(zeros),
                   'dedup': len(dedup), 'domain_ok': len(cat7)},
        'zeros_level7_domain_ok': [{'re': z['u'].real, 'im': z['u'].imag,
                                    'absH': z['absH']} for z in cat7],
    }
    C.dump_json(out, 'r1b_H6_level7.json')
    log('wrote r1b_H6_level7.json (domain_ok=%d)' % len(cat7))

    orig = json.load(open(C.ORIG_PATH, encoding='utf-8'))
    u7_orig = [(z['re'], z['im']) for z in orig['zeros_level7']]
    my7 = np.array([(z['u'].real, z['u'].imag) for z in cat7]) if cat7 else np.empty((0, 2))
    matched = 0
    maxd = 0.0
    for (re, im) in u7_orig:
        if len(my7) == 0:
            continue
        d = np.abs(my7 - np.array([re, im])).max(axis=1)
        m = float(d.min())
        if m < 1e-6:
            matched += 1
        maxd = max(maxd, m)
    log('sanity: %d/%d ORIG-stored u7 values reproduced within 1e-6 (worst-among-matched-scan=%.2e)'
        % (matched, len(u7_orig), maxd))


if __name__ == '__main__':
    main()

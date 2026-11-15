#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_level7_catalog.py -- reproduce the level-7 (m=7) zero catalog EXACTLY as
p2_H4_zeros.main() does (original 1/16 grid, R_max=20, no shift): grid scan ->
local minima -> newton_refine(tol=1e-9) -> domain_check -> dedup(|du|<1e-4).

This is needed because r1_d4_rescan.py (STEP24 driver) deliberately skips
level-7 pairing (per BRIEF), so STEP24's level-8 zeros have no recorded m=7
partner.  To test "does this zero have a genuine m=7 partner" for BOTH the
62 ORIG and 146 STEP24 zeros on equal footing, we need our own independently
produced m=7 catalog -- built with the untouched p2 functions and the
UNSHIFTED 1/16 grid (this is what the audited pipeline itself uses to
discover level-7 zeros; we are testing the pipeline's own criterion, not
inventing a new one).

Sanity check performed at the end: match this catalog against the u7 values
already stored inside quartic_pilot/H4_zeros.json ('zeros_level7', paired
1:1 with 'zeros_level8') -- should reproduce all 62 to high precision.

Output: r1a_H4_level7.json
"""
import json
import time

import numpy as np

import r1a_common as C

log = C.mk_logger('r1a_level7_catalog.log')


def main():
    log('=== r1a level-7 catalog (grid 1/16, R=20, unshifted) ===')
    AU0, Ls = C.load_au0()
    log('AU0_4 loaded from cache')

    grid = C.grid_orig(C.R_MAX, C.STEP16)
    log('grid built: %s' % (grid.shape,))

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
        u, res = C.newton_from(u0, 7, AU0)
        if u is None or res > 1e-5:
            continue
        mt4, merr, den = C.p2.domain_check(u, 7, AU0)
        domain_ok = bool(mt4 < 2.2 and merr < 1e-11 and den > 1e-12)
        zeros.append({'u': u, 'absH': res, 'maxt4': mt4, 'maxerr': merr,
                      'den': den, 'domain_ok': domain_ok})
        if (i2 + 1) % 100 == 0:
            log('  newton %d/%d (%.0fs)' % (i2 + 1, len(cand), time.time() - t1))
    log('refined (|H|<1e-5): %d' % len(zeros))

    # dedup |du|<1e-4 exactly as p2_H4_zeros.main does for level 7 (applies to
    # ALL refined zeros, domain_ok or not; the domain filter is applied AFTER)
    dedup = []
    for z in zeros:
        if all(abs(z['u'] - v['u']) > 1e-4 for v in dedup):
            dedup.append(z)
    log('dedup |du|<1e-4: %d' % len(dedup))

    cat7 = [z for z in dedup if z['domain_ok']]
    log('domain_ok subset (== p2 cat7): %d' % len(cat7))

    out = {
        'meta': {'tag': 'r1a level-7 catalog, reproduces p2_H4_zeros m=7 stage exactly',
                 'grid': {'R_max': C.R_MAX, 'step': C.STEP16, 'shift': [0.0, 0.0]},
                 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'elapsed_sec': time.time() - t0},
        'counts': {'candidates': len(cand), 'refined': len(zeros),
                   'dedup': len(dedup), 'domain_ok': len(cat7)},
        'zeros_level7_domain_ok': [{'re': z['u'].real, 'im': z['u'].imag,
                                    'absH': z['absH']} for z in cat7],
        'zeros_level7_all_dedup': [{'re': z['u'].real, 'im': z['u'].imag,
                                    'absH': z['absH'], 'domain_ok': z['domain_ok']}
                                   for z in dedup],
    }
    C.dump_json(out, 'r1a_H4_level7.json')
    log('wrote r1a_H4_level7.json (domain_ok=%d)' % len(cat7))

    # ---- sanity check vs ORIG's stored zeros_level7 -------------------------
    orig = json.load(open(C.ORIG_PATH, encoding='utf-8'))
    u7_orig = [(z['re'], z['im']) for z in orig['zeros_level7']]
    my7 = np.array([(z.real, z.imag) for z in [zz['u'] for zz in cat7]])
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
    log('sanity: %d/%d ORIG-stored u7 values reproduced within 1e-6 (worst nearest-dist among matched-set scan = %.2e)'
        % (matched, len(u7_orig), maxd))


if __name__ == '__main__':
    main()

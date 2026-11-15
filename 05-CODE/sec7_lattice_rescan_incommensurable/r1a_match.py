#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_match.py -- ORIG(62) vs STEP24(146) matching at |du|<1e-6: BOTH / ORIG-ONLY
/ NEW-ONLY, with full index bookkeeping (so other r1a_ scripts can join back
to absH, high-precision string, etc).  Reimplements the same greedy
nearest-match logic as r1_stats.match_subsets (read-only reused elsewhere in
this folder), but keeps indices instead of just (re,im) tuples.

Output: r1a_match.json
"""
import time

import numpy as np

import r1a_common as C

log = C.mk_logger('r1a_match.log')

TOL = 1e-6


def main():
    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    step24_zs, step24_zl, _ = C.load_catalog('H4_step24.json')
    log('ORIG n=%d, STEP24 n=%d' % (len(orig_zs), len(step24_zs)))

    both_o, both_n, only_new, only_orig = C.match_indices(orig_zs, step24_zs, tol=TOL)
    log('BOTH=%d, ORIG-ONLY=%d, NEW-ONLY=%d (check: %d+%d=%d==62?)'
        % (len(both_o), len(only_orig), len(only_new), len(both_o), len(only_orig), len(both_o) + len(only_orig)))

    # nearest-distance diagnostics for ORIG-ONLY / NEW-ONLY (useful for task 2/3)
    orig_arr = np.array(orig_zs)
    step24_arr = np.array(step24_zs)

    only_orig_detail = []
    for j in only_orig:
        u = orig_zl[j]
        if len(step24_arr):
            d = np.abs(step24_arr - np.array(orig_zs[j])).max(axis=1)
            nearest_i = int(d.argmin())
            nearest_d = float(d.min())
        else:
            nearest_i, nearest_d = None, None
        only_orig_detail.append({'orig_idx': j, 're': u['re'], 'im': u['im'],
                                 'absH_stored': u.get('absH'),
                                 'nearest_step24_idx': nearest_i,
                                 'nearest_step24_dist': nearest_d})

    only_new_detail = []
    for i in only_new:
        u = step24_zl[i]
        if len(orig_arr):
            d = np.abs(orig_arr - np.array(step24_zs[i])).max(axis=1)
            nearest_j = int(d.argmin())
            nearest_d = float(d.min())
        else:
            nearest_j, nearest_d = None, None
        only_new_detail.append({'step24_idx': i, 're': u['re'], 'im': u['im'],
                                'absH_stored': u.get('absH'),
                                'nearest_orig_idx': nearest_j,
                                'nearest_orig_dist': nearest_d})

    both_detail = []
    for j, i in zip(both_o, both_n):
        both_detail.append({'orig_idx': j, 'step24_idx': i,
                            'orig_re': orig_zl[j]['re'], 'orig_im': orig_zl[j]['im'],
                            'step24_re': step24_zl[i]['re'], 'step24_im': step24_zl[i]['im'],
                            'dist': float(max(abs(orig_zl[j]['re'] - step24_zl[i]['re']),
                                              abs(orig_zl[j]['im'] - step24_zl[i]['im'])))})

    out = {
        'meta': {'tol': TOL, 'n_orig': len(orig_zs), 'n_step24': len(step24_zs),
                 'date': time.strftime('%Y-%m-%d %H:%M:%S')},
        'counts': {'both': len(both_o), 'only_orig': len(only_orig), 'only_new': len(only_new)},
        'both': both_detail,
        'only_orig': only_orig_detail,
        'only_new': only_new_detail,
    }
    C.dump_json(out, 'r1a_match.json')
    log('wrote r1a_match.json')


if __name__ == '__main__':
    main()

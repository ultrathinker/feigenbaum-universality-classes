#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""r1b_match.py -- ORIG(342) vs STEP24(641) matching at |du|<1e-6, d=6
analogue of r1a_match.py.  Output: r1b_match.json"""
import time

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_match.log')

TOL = 1e-6


def main():
    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    step24_zs, step24_zl, _ = C.load_catalog('H6_step24.json')
    log('ORIG n=%d, STEP24 n=%d' % (len(orig_zs), len(step24_zs)))

    both_o, both_n, only_new, only_orig = C.match_indices(orig_zs, step24_zs, tol=TOL)
    log('BOTH=%d, ORIG-ONLY=%d, NEW-ONLY=%d' % (len(both_o), len(only_orig), len(only_new)))

    orig_arr = np.array(orig_zs)
    step24_arr = np.array(step24_zs)

    only_orig_detail = []
    for j in only_orig:
        u = orig_zl[j]
        d = np.abs(step24_arr - np.array(orig_zs[j])).max(axis=1) if len(step24_arr) else None
        only_orig_detail.append({'orig_idx': j, 're': u['re'], 'im': u['im'],
                                 'absH_stored': u.get('absH'),
                                 'nearest_step24_dist': float(d.min()) if d is not None else None})

    only_new_detail = []
    for i in only_new:
        u = step24_zl[i]
        d = np.abs(orig_arr - np.array(step24_zs[i])).max(axis=1) if len(orig_arr) else None
        only_new_detail.append({'step24_idx': i, 're': u['re'], 'im': u['im'],
                                'absH_stored': u.get('absH'),
                                'nearest_orig_dist': float(d.min()) if d is not None else None})

    both_detail = []
    for j, i in zip(both_o, both_n):
        both_detail.append({'orig_idx': j, 'step24_idx': i,
                            'orig_re': orig_zl[j]['re'], 'orig_im': orig_zl[j]['im'],
                            'step24_re': step24_zl[i]['re'], 'step24_im': step24_zl[i]['im']})

    out = {'meta': {'tol': TOL, 'n_orig': len(orig_zs), 'n_step24': len(step24_zs),
                    'date': time.strftime('%Y-%m-%d %H:%M:%S')},
           'counts': {'both': len(both_o), 'only_orig': len(only_orig), 'only_new': len(only_new)},
           'both': both_detail, 'only_orig': only_orig_detail, 'only_new': only_new_detail}
    C.dump_json(out, 'r1b_match.json')
    log('wrote r1b_match.json')


if __name__ == '__main__':
    main()

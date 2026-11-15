#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_stats_lattice.py -- Task 4 (d=6): lattice statistics (q=16,24,32,48;
re/im/joint; binomial p) for ORIG, STEP24, BOTH, ORIG-ONLY, NEW-ONLY, plus
STEP24 vs its own shifted 1/24 grid vs the fixed unshifted grid.  Also
consolidates task1/2/3/late9 classification summaries -- this IS the file
the coordinator wants "all numbers" in.
Output: r1b_stats_d6.json
"""
import json
import os
import time

import r1b_common as C

log = C.mk_logger('r1b_stats_lattice.log')

Q_LIST = (16, 24, 32, 48)


def cat_stats(zs):
    re = [z[0] for z in zs]
    im = [z[1] for z in zs]
    return C.stats_block(re, im, q_list=Q_LIST)


def main():
    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'q_list': list(Q_LIST)}}

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    step24_zs, step24_zl, step24_j = C.load_catalog('H6_step24.json')
    match = json.load(open('r1b_match.json', encoding='utf-8'))

    both = [(d['orig_re'], d['orig_im']) for d in match['both']]
    only_orig = [(d['re'], d['im']) for d in match['only_orig']]
    only_new = [(d['re'], d['im']) for d in match['only_new']]

    log('n: ORIG=%d STEP24=%d BOTH=%d ORIG-ONLY=%d NEW-ONLY=%d'
        % (len(orig_zs), len(step24_zs), len(both), len(only_orig), len(only_new)))

    out['catalogs'] = {}
    out['catalogs']['ORIG_342'] = cat_stats(orig_zs)
    out['catalogs']['STEP24_641'] = cat_stats(step24_zs)
    out['catalogs']['BOTH_155'] = cat_stats(both)
    out['catalogs']['ORIG_ONLY_187'] = cat_stats(only_orig)
    out['catalogs']['NEW_ONLY_486'] = cat_stats(only_new)

    sc = step24_j['meta']['scan']
    sh_re, sh_im = sc['shift_re'], sc['shift_im']
    log('STEP24 shift = (%.6f, %.6f), step=%.6f' % (sh_re, sh_im, sc['step']))
    re24 = [z[0] for z in step24_zs]
    im24 = [z[1] for z in step24_zs]
    out['step24_vs_own_shifted_grid_q24'] = C.near_stats_shifted(re24, im24, 24, sh_re, sh_im)
    out['step24_vs_unshifted_grid_q24'] = C.near_stats(re24, im24, q=24)
    out['step24_vs_unshifted_grid_q16'] = C.near_stats(re24, im24, q=16)

    only_new_re = [z[0] for z in only_new]
    only_new_im = [z[1] for z in only_new]
    out['new_only_vs_own_shifted_grid_q24'] = C.near_stats_shifted(only_new_re, only_new_im, 24, sh_re, sh_im)
    only_orig_re = [z[0] for z in only_orig]
    only_orig_im = [z[1] for z in only_orig]
    out['orig_only_vs_unshifted_grid_q16'] = C.near_stats(only_orig_re, only_orig_im, q=16)

    out['task5_note'] = ('SKIPPED for d=6: a full fresh 1/16+shift rescan (shift16-equivalent) '
                         'would take ~the same order as STEP24, which took 5111s (~85 min) even '
                         'at 8 workers for the mp-polish stage alone -- far beyond the ~15 min / '
                         '<=6 workers budget given for this optional step. Not run.')

    if os.path.exists('r1b_authenticity.json'):
        auth = json.load(open('r1b_authenticity.json', encoding='utf-8'))
        out['task1_authenticity_summary'] = {
            k: auth.get(k) for k in ('ORIG_ONLY_sample_summary', 'NEW_ONLY_sample_summary', 'BOTH_full_summary')}
        out['task1_authenticity_summary']['sample_note'] = auth['meta']['sample_note']
    if os.path.exists('r1b_origonly_forensics.json'):
        oo = json.load(open('r1b_origonly_forensics.json', encoding='utf-8'))
        out['task2_origonly_classification'] = oo['classification_counts']
    if os.path.exists('r1b_newonly_forensics.json'):
        no = json.load(open('r1b_newonly_forensics.json', encoding='utf-8'))
        out['task3_newonly_classification_SAMPLE'] = no['classification_counts']
        out['task3_newonly_sample_meta'] = no['meta']
    if os.path.exists('r1b_late9.json'):
        l9 = json.load(open('r1b_late9.json', encoding='utf-8'))
        out['late9'] = {'n_found_by_step24': l9['n_found_by_step24'], 'rows': l9['rows']}
    if os.path.exists('r1b_H6_level7.json'):
        lv7 = json.load(open('r1b_H6_level7.json', encoding='utf-8'))
        out['level7_catalog_sanity'] = lv7['counts']

    C.dump_json(out, 'r1b_stats_d6.json')
    log('wrote r1b_stats_d6.json')

    print()
    print('%-24s %5s %7s %7s %8s %7s %7s %7s %7s' % ('catalog', 'n', 're16', 'im16', 'joint16', 're24', 'im24', 're32', 're48'))
    for nm, rec in out['catalogs'].items():
        n = rec['n']
        def v(q, k):
            r = rec['near'][str(q)]
            x = r.get(k)
            return '%.1f%%' % (100 * x) if x is not None else '-'
        print('%-24s %5d %7s %7s %8s %7s %7s %7s %7s'
              % (nm, n, v(16, 're_near'), v(16, 'im_near'), v(16, 'joint'),
                 v(24, 're_near'), v(24, 'im_near'), v(32, 're_near'), v(48, 're_near')))


if __name__ == '__main__':
    main()

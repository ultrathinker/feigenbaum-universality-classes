#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_stats_lattice.py -- Task 4: lattice statistics (q=16,24,32,48; re/im/joint;
binomial p) for ORIG, STEP24, BOTH, ORIG-ONLY, NEW-ONLY, and the 'genuine'
subsets from r1a_authenticity.json; plus STEP24 measured against its OWN
shifted 1/24 grid nodes; plus (once available) r1a_H4_shift16.json measured
against BOTH the unshifted-1/16 nodes and its own shifted-1/16 nodes.

Requires: r1a_match.json, r1a_authenticity.json, H4_step24.json, H4_sanity.json,
          quartic_pilot/H4_zeros.json, and (if present) r1a_H4_shift16.json.
Output: r1a_stats_d4.json  (this IS the consolidated numbers file the BRIEF
        for this sub-task asks for: "Все числа -- в r1a_stats_d4.json").
"""
import json
import os
import time

import numpy as np

import r1a_common as C

log = C.mk_logger('r1a_stats_lattice.log')

Q_LIST = (16, 24, 32, 48)


def cat_stats(zs):
    re = [z[0] for z in zs]
    im = [z[1] for z in zs]
    return C.stats_block(re, im, q_list=Q_LIST)


def main():
    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'q_list': list(Q_LIST),
                    'near_tol': '1/16 of node spacing (frac(q*x)<1/16)',
                    'p_method': 'binomial one-sided (scipy binomtest, alternative=greater)'}}

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    step24_zs, step24_zl, step24_j = C.load_catalog('H4_step24.json')
    match = json.load(open('r1a_match.json', encoding='utf-8'))
    auth = json.load(open('r1a_authenticity.json', encoding='utf-8'))

    both = [(d['orig_re'], d['orig_im']) for d in match['both']]
    only_orig = [(d['re'], d['im']) for d in match['only_orig']]
    only_new = [(d['re'], d['im']) for d in match['only_new']]

    log('n: ORIG=%d STEP24=%d BOTH=%d ORIG-ONLY=%d NEW-ONLY=%d'
        % (len(orig_zs), len(step24_zs), len(both), len(only_orig), len(only_new)))

    out['catalogs'] = {}
    out['catalogs']['ORIG_62'] = cat_stats(orig_zs)
    out['catalogs']['STEP24_146'] = cat_stats(step24_zs)
    out['catalogs']['BOTH_31'] = cat_stats(both)
    out['catalogs']['ORIG_ONLY_31'] = cat_stats(only_orig)
    out['catalogs']['NEW_ONLY_115'] = cat_stats(only_new)

    # ---- 'genuine' subsets from r1a_authenticity.json ----
    orig_genuine = [(r['re'], r['im']) for r in auth['ORIG'] if r['genuine']]
    orig_survive = [(r['re'], r['im']) for r in auth['ORIG'] if r['would_survive_orig_pipeline']]
    step24_genuine = [(r['re'], r['im']) for r in auth['STEP24'] if r['genuine']]
    step24_survive = [(r['re'], r['im']) for r in auth['STEP24'] if r['would_survive_orig_pipeline']]
    log('ORIG genuine=%d survive=%d ; STEP24 genuine=%d survive=%d'
        % (len(orig_genuine), len(orig_survive), len(step24_genuine), len(step24_survive)))
    out['catalogs']['ORIG_genuine'] = cat_stats(orig_genuine)
    out['catalogs']['ORIG_would_survive_pipeline'] = cat_stats(orig_survive)
    out['catalogs']['STEP24_genuine'] = cat_stats(step24_genuine)
    out['catalogs']['STEP24_would_survive_pipeline'] = cat_stats(step24_survive)

    # ---- STEP24 vs its OWN shifted grid ----
    sc = step24_j['meta']['scan']
    sh_re, sh_im = sc['shift_re'], sc['shift_im']
    log('STEP24 shift = (%.6f, %.6f), step=%.6f' % (sh_re, sh_im, sc['step']))
    re24 = [z[0] for z in step24_zs]
    im24 = [z[1] for z in step24_zs]
    out['step24_vs_own_shifted_grid_q24'] = C.near_stats_shifted(re24, im24, 24, sh_re, sh_im)
    out['step24_vs_unshifted_grid_q24'] = C.near_stats(re24, im24, q=24)
    out['step24_vs_unshifted_grid_q16'] = C.near_stats(re24, im24, q=16)

    # ---- shift16 control (task 5), if already run ----
    if os.path.exists('r1a_H4_shift16.json'):
        s16_zs, s16_zl, s16_j = C.load_catalog('r1a_H4_shift16.json')
        sc16 = s16_j['meta']['scan']
        sh16_re, sh16_im = sc16['shift_re'], sc16['shift_im']
        log('shift16: n=%d shift=(%.6f,%.6f)' % (len(s16_zs), sh16_re, sh16_im))
        out['catalogs']['SHIFT16'] = cat_stats(s16_zs)
        re16 = [z[0] for z in s16_zs]
        im16 = [z[1] for z in s16_zs]
        out['shift16_vs_unshifted_grid_q16'] = C.near_stats(re16, im16, q=16)
        out['shift16_vs_own_shifted_grid_q16'] = C.near_stats_shifted(re16, im16, 16, sh16_re, sh16_im)
        # match shift16 vs ORIG for reference
        both_o16, both_n16, only_new16, only_orig16 = C.match_indices(orig_zs, s16_zs, tol=1e-6)
        out['shift16_vs_orig_overlap'] = {'n_both': len(both_o16), 'n_only_shift16': len(only_new16),
                                          'n_only_orig': len(only_orig16)}
        log('shift16 overlap with ORIG: both=%d only_shift16=%d only_orig=%d'
            % (len(both_o16), len(only_new16), len(only_orig16)))
    else:
        out['shift16_note'] = 'r1a_H4_shift16.json not yet produced'
        log('r1a_H4_shift16.json not found -- run r1a_d4_rescan.py shift16 first')

    # ---- consolidate task 1/2/3 classification summaries (this IS the file
    # the BRIEF wants "all numbers" in) ----
    out['task1_authenticity_summary'] = {
        'ORIG': auth.get('ORIG_summary'), 'STEP24': auth.get('STEP24_summary'),
        'note': ('STEP24 zeros WITH a level-7 partner are EXACTLY the 31 BOTH '
                 'zeros (set-equality verified); all 115 NEW-ONLY zeros lack a '
                 'level-7 partner.')}
    if os.path.exists('r1a_origonly_forensics.json'):
        oo = json.load(open('r1a_origonly_forensics.json', encoding='utf-8'))
        out['task2_origonly_classification'] = oo['classification_counts']
    if os.path.exists('r1a_newonly_forensics.json'):
        no = json.load(open('r1a_newonly_forensics.json', encoding='utf-8'))
        out['task3_newonly_classification'] = no['classification_counts']
        beyond = sum(1 for z in no['zeros'] if z.get('beyond_orig_rmax'))
        out['task3_newonly_classification']['_beyond_orig_rmax22_of_which_not_genuine'] = beyond

    C.dump_json(out, 'r1a_stats_d4.json')
    log('wrote r1a_stats_d4.json')

    # ---- print compact summary table ----
    print()
    print('%-28s %5s %7s %7s %8s %7s %7s %7s %7s' % ('catalog', 'n', 're16', 'im16', 'joint16', 're24', 'im24', 're32', 're48'))
    for nm, rec in out['catalogs'].items():
        n = rec['n']
        def v(q, k):
            r = rec['near'][str(q)]
            x = r.get(k)
            return '%.1f%%' % (100 * x) if x is not None else '-'
        print('%-28s %5d %7s %7s %8s %7s %7s %7s %7s'
              % (nm, n, v(16, 're_near'), v(16, 'im_near'), v(16, 'joint'),
                 v(24, 're_near'), v(24, 'im_near'), v(32, 're_near'), v(48, 're_near')))


if __name__ == '__main__':
    main()

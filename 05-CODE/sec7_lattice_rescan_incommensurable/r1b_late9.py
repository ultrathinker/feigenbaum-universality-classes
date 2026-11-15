#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_late9.py -- d=6 "late 9" check (S031,S036,S070,S104,S187,S235,S279,S319,
S339 per 03-APPENDIX-B-DATA v0.4.1 SS B.2.3/B.3.2/B.6.3): these were REJECTED
by the standard float-pilot m=7/8 pairing pipeline (accepted=False in
agentB_deg6/match_results.json) and only certified later via a separate
RIGOROUS interval-Krawczyk re-open at m=10-12 -- i.e. they are NOT
"discovered the normal way" to begin with.  Publication table B.6.3 reports
100% (9/9) joint-near-q16 for this set, the most extreme lattice-adherence
number in the whole publication.

This script: (a) confirms they ARE present in H6_zeros_ORIG.json's 342
(matched within ~1e-6 of the Krawczyk-box u_cat), (b) checks whether STEP24
(1/24+shift rescan) recovers them at all, and if so at what distance and
whether that recovered point sits near the SHIFTED 1/24 grid or stays pinned
to the original 1/16 node, (c) reports their frac-to-1/16-node distance
precisely.

Output: r1b_late9.json
"""
import json
import time

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_late9.log')

LATE_IDS = ['S031', 'S036', 'S070', 'S104', 'S187', 'S235', 'S279', 'S319', 'S339']


def node_dist16(re, im):
    fr = C.frac([re], 16)[0] / 16.0
    fi = C.frac([im], 16)[0] / 16.0
    return fr, fi


def main():
    late = json.load(open('../sec3_pilot_d6/match_results.json', encoding='utf-8'))
    late9 = {r['id']: r for r in late['catalog_results'] if r['id'] in LATE_IDS}
    log('loaded %d/9 late records from match_results.json' % len(late9))

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    orig_arr = np.array(orig_zs)
    step24_zs, step24_zl, step24_j = C.load_catalog('H6_step24.json')
    step24_arr = np.array(step24_zs)
    sc = step24_j['meta']['scan']
    sh_re, sh_im = sc['shift_re'], sc['shift_im']
    step = sc['step']

    out_rows = []
    for lid in LATE_IDS:
        r = late9[lid]
        u_cat = complex(r['u_cat_re'], r['u_cat_im'])
        d_orig = np.abs(orig_arr - np.array([u_cat.real, u_cat.imag])).max(axis=1)
        j_orig = int(d_orig.argmin())
        u_orig = orig_zl[j_orig]

        d_step24 = np.abs(step24_arr - np.array([u_cat.real, u_cat.imag])).max(axis=1)
        j_step24 = int(d_step24.argmin())
        found_step24 = bool(d_step24.min() < 1e-4)  # looser tol -- these are hard cases

        fr16, fi16 = node_dist16(u_orig['re'], u_orig['im'])
        row = {'id': lid, 'u_cat_re': r['u_cat_re'], 'u_cat_im': r['u_cat_im'],
               'accepted_in_match_results': r.get('accepted'),
               'orig_match_dist': float(d_orig.min()), 'orig_idx': j_orig,
               'orig_re': u_orig['re'], 'orig_im': u_orig['im'],
               'dist_to_1_16_node': {'re_frac_of_step': float(fr16 * 16), 'im_frac_of_step': float(fi16 * 16),
                                     'abs_dist': float(np.hypot(fr16, fi16))},
               'step24_nearest_dist': float(d_step24.min()),
               'found_by_step24_at_1e-4': found_step24}
        if found_step24:
            u_s24 = step24_zl[j_step24]
            row['step24_re'] = u_s24['re']
            row['step24_im'] = u_s24['im']
            dre_shift = C.frac([u_s24['re'] - sh_re], 24)[0]
            dim_shift = C.frac([u_s24['im'] - sh_im], 24)[0]
            row['step24_dist_to_own_shifted_1_24_node'] = float(np.hypot(dre_shift, dim_shift) * 24 / 16)
            row['step24_dist_to_unshifted_1_16_node'] = float(np.hypot(*node_dist16(u_s24['re'], u_s24['im'])) * 16)
        out_rows.append(row)
        log('%s: |u_cat-ORIG|=%.2e dist_to_1/16_node=%.2e  STEP24: nearest=%.2e found@1e-4=%s'
            % (lid, row['orig_match_dist'], row['dist_to_1_16_node']['abs_dist'],
               row['step24_nearest_dist'], found_step24))

    n_found = sum(1 for r in out_rows if r['found_by_step24_at_1e-4'])
    log('SUMMARY: %d/9 late zeros recovered by STEP24 (tol 1e-4)' % n_found)

    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'ids': LATE_IDS,
                    'note': ('all 9 have accepted=False in match_results.json (rejected by '
                             'the standard m=7/8 float pipeline; certified only later via a '
                             'SEPARATE rigorous interval-Krawczyk re-open at m=10-12, per '
                             '03-APPENDIX-B-DATA SS B.2.3) -- they are present in '
                             'H6_zeros_ORIG.json (342) but were NOT discovered by the plain '
                             'grid-scan+Newton route this whole R1/R1B forensic exercise audits.')},
           'n_found_by_step24': n_found, 'rows': out_rows}
    C.dump_json(out, 'r1b_late9.json')
    log('wrote r1b_late9.json')


if __name__ == '__main__':
    main()

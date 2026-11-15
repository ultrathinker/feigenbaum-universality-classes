#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_origonly_forensics.py -- Task 2: for each of the 31 ORIG-ONLY zeros
(present in the original 1/16-grid catalog, NOT recovered by the 1/24+shift
rescan at |du|<1e-6), diagnose why the 1/24 scan missed it:
  - evaluate |H_8| on a 5x5 patch of the ACTUAL shifted 1/24 grid (shift read
    from H4_step24.json's own meta -- exact reproduction, not re-derived from
    the seed) centred on the node nearest the ORIG zero;
  - was that node a local minimum under the exact rule used by
    r1_d4_rescan.find_local_minima (<=8 neighbours AND <0.2)?
  - if yes: run the FULL acceptance pipeline (newton_refine tol=1e-9 -> res<=1e-5
    -> domain_check -> mp_polish -> |H|<1e-20) from that node and see where it
    converges (back to the ORIG zero, to a different zero, or nowhere/rejected);
  - if no: is the node merely *not a local minimum* (H too flat / neighbour
    lower) or is it >=0.2 (not even a scan candidate) or invalid (domain
    violation already at grid level)?

Classification per zero:
  'recovered_but_different_node'  -- pipeline from nearest node converges back
                                      to this exact zero (contradiction; would
                                      mean it should be in STEP24 -- investigate)
  'missed_not_local_min'          -- node existed, |H| valid, but not a local
                                      minimum on the shifted grid (a neighbour
                                      is lower) -- the SHIFT moved the discrete
                                      minimum onto/near a different candidate
  'missed_newton_diverged_or_rejected' -- was a local min, but Newton from it
                                      does not return to this zero / fails
                                      domain_check / fails mp_polish
  'missed_converged_elsewhere'    -- was a local min, Newton converges to a
                                      DIFFERENT zero (basin capture by a
                                      neighbouring root)
  'missed_invalid_node'           -- node itself invalid (domain violated in
                                      scan_level already)

Output: r1a_origonly_forensics.json
"""
import json
import time

import numpy as np

import r1a_common as C

log = C.mk_logger('r1a_origonly_forensics.log')


def main():
    AU0, Ls = C.load_au0()
    match = json.load(open('r1a_match.json', encoding='utf-8'))
    only_orig = match['only_orig']
    log('ORIG-ONLY zeros to analyze: %d' % len(only_orig))

    step24_meta = json.load(open('H4_step24.json', encoding='utf-8'))['meta']['scan']
    shift_re = step24_meta['shift_re']
    shift_im = step24_meta['shift_im']
    step = step24_meta['step']
    log('STEP24 grid: step=%.6f shift=(%.6f,%.6f) [read from H4_step24.json meta, exact]'
        % (step, shift_re, shift_im))

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    orig_arr = np.array(orig_zs)

    results = []
    t0 = time.time()
    for k, rec in enumerate(only_orig):
        u0 = complex(rec['re'], rec['im'])
        patch = C.local_min_patch(u0, AU0, step, shift_re, shift_im, m=8, half=2)
        entry = {'orig_idx': rec['orig_idx'], 're': rec['re'], 'im': rec['im'],
                 'abs_u': abs(u0),
                 'absH_stored': rec['absH_stored'],
                 'nearest_step24_dist': rec['nearest_step24_dist'],
                 'nearest_node': {'re': patch['node_re'], 'im': patch['node_im']},
                 'node_dist_from_zero': abs(complex(patch['node_re'], patch['node_im']) - u0),
                 'node_valid': patch['center_valid'],
                 'node_Habs': patch['center_Habs'],
                 'is_local_min_on_step24_grid': patch['is_local_min'],
                 'neighbour_Habs': patch['neighbour_Habs']}

        if not patch['center_valid']:
            entry['classification'] = 'missed_invalid_node'
        elif not patch['is_local_min']:
            entry['classification'] = 'missed_not_local_min'
        else:
            # it WAS a local-min candidate on the shifted grid -- run the pipeline
            node_u0 = complex(patch['node_re'], patch['node_im'])
            pipe = C.full_candidate_pipeline(node_u0, AU0, m=8)
            entry['pipeline_from_node'] = {
                'newton_u': [pipe['newton_u'].real, pipe['newton_u'].imag] if pipe.get('newton_u') else None,
                'newton_absH': pipe.get('newton_absH'),
                'domain_ok': pipe.get('domain', {}).get('ok'),
                'polish_absH': pipe.get('polish_absH'),
                'polish_u': [pipe['polish_u'].real, pipe['polish_u'].imag] if pipe.get('polish_u') else None,
                'stage_failed': pipe.get('stage_failed'),
            }
            if pipe.get('stage_failed') is not None:
                entry['classification'] = 'missed_newton_diverged_or_rejected'
            else:
                pu = pipe['polish_u']
                dist_to_this = abs(pu - u0)
                if dist_to_this < 1e-4:
                    entry['classification'] = 'recovered_but_different_node(contradiction-check)'
                else:
                    # which zero did it converge to?
                    d = np.abs(orig_arr - np.array([pu.real, pu.imag])).max(axis=1)
                    nearest_j = int(d.argmin())
                    entry['converged_to_orig_idx'] = nearest_j
                    entry['converged_to_dist'] = float(d.min())
                    entry['classification'] = 'missed_converged_elsewhere'
        results.append(entry)
        if (k + 1) % 10 == 0:
            log('  %d/%d (%.0fs)' % (k + 1, len(only_orig), time.time() - t0))

    from collections import Counter
    cls_counts = Counter(r['classification'] for r in results)
    log('classification counts: %s' % dict(cls_counts))

    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'step24_shift': [shift_re, shift_im], 'step24_step': step},
           'n': len(results), 'classification_counts': dict(cls_counts),
           'zeros': results}
    C.dump_json(out, 'r1a_origonly_forensics.json')
    log('wrote r1a_origonly_forensics.json')


if __name__ == '__main__':
    main()

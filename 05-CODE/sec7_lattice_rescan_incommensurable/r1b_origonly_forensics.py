#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_origonly_forensics.py -- Task 2 (d=6): for each of the 187 ORIG-ONLY
zeros, diagnose why the 1/24+shift rescan missed it (candidate-selection vs
domain_check vs Newton), same protocol as r1a_origonly_forensics.py for d=4.
FULL set (187), not sampled -- local_min_patch is a cheap 5x5-point probe
regardless of the domain size, so no need to subsample this task.
Output: r1b_origonly_forensics.json
"""
import json
import time

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_origonly_forensics.log')


def main():
    AU0, Ls = C.load_au0()
    match = json.load(open('r1b_match.json', encoding='utf-8'))
    only_orig = match['only_orig']
    log('ORIG-ONLY zeros to analyze: %d' % len(only_orig))

    step24_meta = json.load(open('H6_step24.json', encoding='utf-8'))['meta']['scan']
    shift_re = step24_meta['shift_re']
    shift_im = step24_meta['shift_im']
    step = step24_meta['step']
    log('STEP24 grid: step=%.6f shift=(%.6f,%.6f) [from H6_step24.json meta]' % (step, shift_re, shift_im))

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
                 'is_local_min_on_step24_grid': patch['is_local_min']}

        if not patch['center_valid']:
            entry['classification'] = 'missed_invalid_node'
        elif not patch['is_local_min']:
            entry['classification'] = 'missed_not_local_min'
        else:
            node_u0 = complex(patch['node_re'], patch['node_im'])
            pipe = C.full_candidate_pipeline(node_u0, AU0, m=8)
            entry['pipeline_from_node'] = {
                'newton_absH': pipe.get('newton_absH'),
                'domain_ok': pipe.get('domain', {}).get('ok'),
                'polish_absH': pipe.get('polish_absH'),
                'stage_failed': pipe.get('stage_failed'),
            }
            if pipe.get('stage_failed') is not None:
                entry['classification'] = 'missed_newton_diverged_or_rejected'
            else:
                pu = pipe['polish_u']
                if abs(pu - u0) < 1e-4:
                    entry['classification'] = 'recovered_but_different_node(contradiction-check)'
                else:
                    d = np.abs(orig_arr - np.array([pu.real, pu.imag])).max(axis=1)
                    nearest_j = int(d.argmin())
                    entry['converged_to_orig_idx'] = nearest_j
                    entry['converged_to_dist'] = float(d.min())
                    entry['classification'] = 'missed_converged_elsewhere'
        results.append(entry)
        if (k + 1) % 30 == 0:
            log('  %d/%d (%.0fs)' % (k + 1, len(only_orig), time.time() - t0))

    from collections import Counter
    cls_counts = Counter(r['classification'] for r in results)
    log('classification counts: %s' % dict(cls_counts))

    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'step24_shift': [shift_re, shift_im], 'step24_step': step},
           'n': len(results), 'classification_counts': dict(cls_counts), 'zeros': results}
    C.dump_json(out, 'r1b_origonly_forensics.json')
    log('wrote r1b_origonly_forensics.json')


if __name__ == '__main__':
    main()

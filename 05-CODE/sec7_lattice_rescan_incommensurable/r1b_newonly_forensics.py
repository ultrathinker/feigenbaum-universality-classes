#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_newonly_forensics.py -- Task 3 (d=6): for the SAMPLED 60 NEW-ONLY zeros
analyzed in r1b_authenticity.json (NEW_ONLY_sample), determine why absent
from ORIG: not genuine, beyond p2's own Newton rmax=70.0, no level-7 partner,
not a local minimum on the ORIGINAL unshifted 1/16 grid, or basin-captured by
a different catalogued zero.  Same protocol as r1a_newonly_forensics.py.
Requires r1b_authenticity.json, r1b_match.json.
Output: r1b_newonly_forensics.json
"""
import json
import time
from collections import Counter

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_newonly_forensics.log')


def main():
    AU0, Ls = C.load_au0()
    auth = json.load(open('r1b_authenticity.json', encoding='utf-8'))
    match = json.load(open('r1b_match.json', encoding='utf-8'))
    sample_idx = auth['NEW_ONLY_sample_idx']
    sample_recs = auth['NEW_ONLY_sample']
    only_new_by_idx = {d['step24_idx']: d for d in match['only_new']}
    log('NEW-ONLY sampled zeros to analyze: %d' % len(sample_idx))

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    orig_arr = np.array(orig_zs)

    results = []
    t0 = time.time()
    for k, (i, a) in enumerate(zip(sample_idx, sample_recs)):
        mrec = only_new_by_idx[i]
        u0 = complex(a['re'], a['im'])
        beyond_rmax = abs(u0) > C.P2_DEFAULT_RMAX
        entry = {'step24_idx': i, 're': a['re'], 'im': a['im'], 'abs_u': abs(u0),
                 'beyond_orig_rmax': beyond_rmax,
                 'nearest_orig_dist': mrec['nearest_orig_dist'],
                 'genuine': a['genuine'], 'level7_partner': a['level7_partner'],
                 'level7_nearest_dist': a['level7_nearest_dist'],
                 'n_return_of_4': a['n_return_of_4'],
                 'absH8_mp_independent': a['absH8_mp_independent']}

        if beyond_rmax:
            entry['classification'] = ('not_genuine_spurious' if not a['genuine']
                                       else 'genuine_but_beyond_orig_rmax70_unreachable')
            results.append(entry)
            continue
        if not a['genuine']:
            entry['classification'] = 'not_genuine_spurious'
            results.append(entry)
            continue

        d_orig = mrec['nearest_orig_dist']
        entry['dedup_boundary_candidate'] = bool(d_orig is not None and 1e-6 <= d_orig < 1e-4)

        patch = C.local_min_patch(u0, AU0, C.STEP16, 0.0, 0.0, m=8, half=2)
        entry['orig16_node_valid'] = patch['center_valid']
        entry['orig16_node_Habs'] = patch['center_Habs']
        entry['orig16_is_local_min'] = patch['is_local_min']

        if not a['level7_partner']:
            reason = 'genuine_but_no_level7_partner'
        elif not patch['center_valid']:
            reason = 'genuine_but_orig16_node_invalid'
        elif not patch['is_local_min']:
            reason = 'genuine_but_not_local_min_on_orig16_grid'
        else:
            node_u0 = complex(patch['node_re'], patch['node_im'])
            pipe = C.full_candidate_pipeline(node_u0, AU0, m=8)
            entry['pipeline_from_orig16_node'] = {'stage_failed': pipe.get('stage_failed')}
            if pipe.get('stage_failed') is not None:
                reason = 'genuine_local_min_but_pipeline_rejects_from_orig16_node'
            else:
                pu = pipe['polish_u']
                if abs(pu - u0) < 1e-4:
                    reason = 'CONTRADICTION_should_be_in_orig'
                else:
                    d = np.abs(orig_arr - np.array([pu.real, pu.imag])).max(axis=1)
                    entry['converged_to_orig_idx'] = int(d.argmin())
                    entry['converged_to_dist'] = float(d.min())
                    reason = 'genuine_but_orig16_node_converges_to_different_zero'
        entry['classification'] = reason
        results.append(entry)
        if (k + 1) % 15 == 0:
            log('  %d/%d (%.0fs)' % (k + 1, len(sample_idx), time.time() - t0))

    cls_counts = Counter(r['classification'] for r in results)
    log('classification counts (of %d SAMPLED NEW-ONLY): %s' % (len(results), dict(cls_counts)))

    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'sample_n': len(results),
                    'sample_of_total': 486},
           'n': len(results), 'classification_counts': dict(cls_counts), 'zeros': results}
    C.dump_json(out, 'r1b_newonly_forensics.json')
    log('wrote r1b_newonly_forensics.json')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_newonly_forensics.py -- Task 3: for each of the 115 NEW-ONLY zeros (found
by the 1/24+shift rescan, absent from ORIG at |du|<1e-6), determine:
  (a) genuineness as a zero of H_8 (from r1a_authenticity.json: independent
      high-precision |H_8|, Newton-return test);
  (b) if genuine, why the ORIGINAL 1/16 unshifted pipeline never produced it:
        - was the nearest 1/16 node NOT a local minimum (<=8 neighbours, <0.2)?
        - if it WAS a local min: does the full original-pipeline run from that
          node fail (Newton doesn't return / domain_check fails / mp_polish
          fails), or converge to a DIFFERENT catalogued zero (basin capture)?
        - does it lack a level-7 partner (pairing filter would have dropped it
          even if the m=8 stage had found it)?
        - is there an ORIG zero within [1e-6, 1e-4) of it (the coarser 1e-4
          ORIG dedup could have merged it into a neighbouring catalogued zero)?
  (c) if not genuine: flag as spurious (an artifact of the loose 1e-8 dedup /
      candidate seeding in the R1 driver, not a real root of H_8).

Requires r1a_authenticity.json and r1a_H4_level7.json.
Output: r1a_newonly_forensics.json with per-zero diagnosis + reason breakdown.
"""
import json
import time
from collections import Counter

import numpy as np

import r1a_common as C

log = C.mk_logger('r1a_newonly_forensics.log')


def main():
    AU0, Ls = C.load_au0()
    match = json.load(open('r1a_match.json', encoding='utf-8'))
    only_new = match['only_new']
    log('NEW-ONLY zeros to analyze: %d' % len(only_new))

    auth = json.load(open('r1a_authenticity.json', encoding='utf-8'))
    step24_auth = auth['STEP24']  # indexed same as H4_step24.json zeros_level8

    orig_zs, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    orig_arr = np.array(orig_zs)

    results = []
    t0 = time.time()
    for k, rec in enumerate(only_new):
        i = rec['step24_idx']
        u0 = complex(rec['re'], rec['im'])
        a = step24_auth[i]
        beyond_rmax = abs(u0) > C.P2_DEFAULT_RMAX
        entry = {'step24_idx': i, 're': rec['re'], 'im': rec['im'],
                 'abs_u': abs(u0),
                 'beyond_orig_rmax': beyond_rmax,
                 'absH_stored': rec['absH_stored'],
                 'nearest_orig_dist': rec['nearest_orig_dist'],
                 'genuine': a['genuine'],
                 'level7_partner': a['level7_partner'],
                 'level7_nearest_dist': a['level7_nearest_dist'],
                 'n_return_of_4': a['n_return_of_4'],
                 'absH8_mp_independent': a['absH8_mp_independent']}

        if beyond_rmax:
            # |u| exceeds p2_H4_zeros.newton_refine's own default rmax=22.0;
            # r1_d4_rescan.py used rmax=25.0 (undeclared deviation) so STEP24
            # can reach candidates the ORIGINAL pipeline's Newton damping
            # could never converge to, regardless of any lattice effect.
            entry['classification'] = ('not_genuine_spurious' if not a['genuine']
                                       else 'genuine_but_beyond_orig_rmax22_unreachable')
            results.append(entry)
            if (k + 1) % 20 == 0:
                log('  %d/%d (%.0fs)' % (k + 1, len(only_new), time.time() - t0))
            continue

        if not a['genuine']:
            entry['classification'] = 'not_genuine_spurious'
            results.append(entry)
            if (k + 1) % 20 == 0:
                log('  %d/%d (%.0fs)' % (k + 1, len(only_new), time.time() - t0))
            continue

        # genuine zero of H_8 -- why absent from ORIG?
        # 1) dedup-boundary check against ORIG (1e-4 vs 1e-6 tolerance gap)
        d_orig = rec['nearest_orig_dist']
        if d_orig is not None and 1e-6 <= d_orig < 1e-4:
            entry['dedup_boundary_candidate'] = True
        else:
            entry['dedup_boundary_candidate'] = False

        # 2) local-min status on the ORIGINAL unshifted 1/16 grid
        patch = C.local_min_patch(u0, AU0, C.STEP16, 0.0, 0.0, m=8, half=2)
        entry['orig16_node'] = {'re': patch['node_re'], 'im': patch['node_im']}
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
            entry['pipeline_from_orig16_node'] = {
                'stage_failed': pipe.get('stage_failed'),
                'polish_u': ([pipe['polish_u'].real, pipe['polish_u'].imag]
                             if pipe.get('polish_u') else None),
                'polish_absH': pipe.get('polish_absH'),
            }
            if pipe.get('stage_failed') is not None:
                reason = 'genuine_local_min_but_pipeline_rejects_from_orig16_node'
            else:
                pu = pipe['polish_u']
                if abs(pu - u0) < 1e-4:
                    reason = 'CONTRADICTION_should_be_in_orig(check_pairing_or_dedup)'
                else:
                    d = np.abs(orig_arr - np.array([pu.real, pu.imag])).max(axis=1)
                    nearest_j = int(d.argmin())
                    entry['converged_to_orig_idx'] = nearest_j
                    entry['converged_to_dist'] = float(d.min())
                    reason = 'genuine_but_orig16_node_converges_to_different_zero'
        entry['classification'] = reason
        results.append(entry)
        if (k + 1) % 20 == 0:
            log('  %d/%d (%.0fs)' % (k + 1, len(only_new), time.time() - t0))

    cls_counts = Counter(r['classification'] for r in results)
    n_dedup_boundary = sum(1 for r in results if r.get('dedup_boundary_candidate'))
    log('classification counts: %s' % dict(cls_counts))
    log('of which dedup-boundary candidates (1e-6<=d<1e-4 to some ORIG zero): %d' % n_dedup_boundary)

    out = {'meta': {'date': time.strftime('%Y-%m-%d %H:%M:%S')},
           'n': len(results), 'classification_counts': dict(cls_counts),
           'n_dedup_boundary_candidates': n_dedup_boundary,
           'zeros': results}
    C.dump_json(out, 'r1a_newonly_forensics.json')
    log('wrote r1a_newonly_forensics.json')


if __name__ == '__main__':
    main()

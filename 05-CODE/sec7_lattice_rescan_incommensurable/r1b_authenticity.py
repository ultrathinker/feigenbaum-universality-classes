#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_authenticity.py -- Task 1 (d=6): genuineness check, SAMPLED per
coordinator instruction (full 342+641=983 would take too long): all 155 BOTH
(STEP24-side values), 60 of 187 ORIG-ONLY (ORIG-side values, seeded sample),
60 of 486 NEW-ONLY (STEP24-side values, seeded sample).  Same 3 checks as
r1a_authenticity.py: independent high-precision |H_8| (mpmath dps=50),
Newton-return test (eps=1e-3, 4 directions, tol 1e-8), level-7 partner
(nearest dist to r1b_H6_level7.json's domain_ok catalog, tol 5e-3).

Requires r1b_H6_level7.json, r1b_match.json.
Output: r1b_authenticity.json
"""
import json
import time

import numpy as np

import r1b_common as C

log = C.mk_logger('r1b_authenticity.log')

EPS = 1e-3
DIRECTIONS = [1.0, -1.0, 1.0j, -1.0j]
NEWTON_RETURN_TOL = 1e-8
PARTNER_TOL = 5e-3
SAMPLE_SEED = 20260816
SAMPLE_N_ORIG_ONLY = 60
SAMPLE_N_NEW_ONLY = 60


def analyze_list(name, zeros_level8, AU0, level7_arr):
    log('--- %s: %d zeros ---' % (name, len(zeros_level8)))
    recs = []
    t0 = time.time()
    for i, z in enumerate(zeros_level8):
        u_float = complex(z['re'], z['im'])
        try:
            u_mp = C.parse_u_mp(z)
            absH8_mp = C.H8_mp_at(u_mp, AU0, m=8)
        except Exception as exc:  # noqa: BLE001
            absH8_mp = None
            log('  [%d] mp eval failed: %s' % (i, exc))
        returns = []
        for d in DIRECTIONS:
            u0 = u_float + EPS * d
            u_ret, res = C.newton_from(u0, 8, AU0)
            if u_ret is None:
                returns.append({'returned': False})
                continue
            dist = abs(u_ret - u_float)
            returns.append({'returned': bool(dist < NEWTON_RETURN_TOL), 'dist': float(dist)})
        n_return = sum(1 for r in returns if r['returned'])
        if len(level7_arr) > 0:
            d7 = np.abs(level7_arr - np.array([u_float.real, u_float.imag])).max(axis=1)
            nearest7 = float(d7.min())
        else:
            nearest7 = None
        has_partner = bool(nearest7 is not None and nearest7 < PARTNER_TOL)
        genuine = bool((absH8_mp is not None and absH8_mp < 1e-15) and n_return >= 3)
        survive = bool(genuine and has_partner)
        recs.append({'re': z['re'], 'im': z['im'], 'absH8_stored': z.get('absH'),
                     'absH8_mp_independent': absH8_mp, 'n_return_of_4': n_return,
                     'level7_nearest_dist': nearest7, 'level7_partner': has_partner,
                     'genuine': genuine, 'would_survive_orig_pipeline': survive})
        if (i + 1) % 30 == 0:
            log('  %d/%d (%.0fs)' % (i + 1, len(zeros_level8), time.time() - t0))
    log('%s done in %.1fs' % (name, time.time() - t0))
    return recs


def summarize(recs):
    n = len(recs)
    return dict(n=n, n_genuine=sum(1 for r in recs if r['genuine']),
                n_highprec_absH_lt_1e15=sum(1 for r in recs if r['absH8_mp_independent'] is not None and r['absH8_mp_independent'] < 1e-15),
                n_newton_return_4of4=sum(1 for r in recs if r['n_return_of_4'] == 4),
                n_newton_return_3plus=sum(1 for r in recs if r['n_return_of_4'] >= 3),
                n_level7_partner=sum(1 for r in recs if r['level7_partner']),
                n_would_survive_orig_pipeline=sum(1 for r in recs if r['would_survive_orig_pipeline']))


def main():
    AU0, Ls = C.load_au0()
    lvl7 = json.load(open('r1b_H6_level7.json', encoding='utf-8'))
    level7_pts = [(z['re'], z['im']) for z in lvl7['zeros_level7_domain_ok']]
    level7_arr = np.array(level7_pts) if level7_pts else np.empty((0, 2))
    log('loaded level-7 catalog: %d domain_ok points' % len(level7_pts))

    _, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    _, step24_zl, _ = C.load_catalog('H6_step24.json')
    match = json.load(open('r1b_match.json', encoding='utf-8'))

    rng = np.random.default_rng(SAMPLE_SEED)
    only_orig_idx = [d['orig_idx'] for d in match['only_orig']]
    only_new_idx = [d['step24_idx'] for d in match['only_new']]
    both_step24_idx = [d['step24_idx'] for d in match['both']]

    sample_orig_only = list(rng.choice(only_orig_idx, size=min(SAMPLE_N_ORIG_ONLY, len(only_orig_idx)), replace=False))
    sample_new_only = list(rng.choice(only_new_idx, size=min(SAMPLE_N_NEW_ONLY, len(only_new_idx)), replace=False))
    log('sample sizes: ORIG-ONLY %d/%d, NEW-ONLY %d/%d, BOTH(full, STEP24-side) %d'
        % (len(sample_orig_only), len(only_orig_idx), len(sample_new_only), len(only_new_idx), len(both_step24_idx)))

    out = {'meta': {'eps': EPS, 'newton_return_tol': NEWTON_RETURN_TOL, 'partner_tol': PARTNER_TOL,
                     'sample_seed': SAMPLE_SEED,
                     'sample_note': 'ORIG-ONLY and NEW-ONLY are SAMPLED (60 each); BOTH is full (155, STEP24-side values)',
                     'genuine_def': 'absH8_mp<1e-15 AND n_return>=3/4 (level-7 partner tracked separately)',
                     'date': time.strftime('%Y-%m-%d %H:%M:%S')}}

    out['ORIG_ONLY_sample_idx'] = [int(i) for i in sample_orig_only]
    out['NEW_ONLY_sample_idx'] = [int(i) for i in sample_new_only]

    out['ORIG_ONLY_sample'] = analyze_list('ORIG-ONLY sample(%d)' % len(sample_orig_only),
                                            [orig_zl[i] for i in sample_orig_only], AU0, level7_arr)
    out['NEW_ONLY_sample'] = analyze_list('NEW-ONLY sample(%d)' % len(sample_new_only),
                                          [step24_zl[i] for i in sample_new_only], AU0, level7_arr)
    out['BOTH_full'] = analyze_list('BOTH full(%d, STEP24-side)' % len(both_step24_idx),
                                    [step24_zl[i] for i in both_step24_idx], AU0, level7_arr)

    for nm in ('ORIG_ONLY_sample', 'NEW_ONLY_sample', 'BOTH_full'):
        out[nm + '_summary'] = summarize(out[nm])
        log('%s summary: %s' % (nm, out[nm + '_summary']))

    C.dump_json(out, 'r1b_authenticity.json')
    log('wrote r1b_authenticity.json')


if __name__ == '__main__':
    main()

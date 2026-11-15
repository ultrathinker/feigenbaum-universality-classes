#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1a_authenticity.py -- Task 1: for each of the 62 ORIG and 146 STEP24 zeros,
check genuineness as a zero of H_8:
  (a) independent |H_8(u)| recompute via mpmath orbit_H_mpmath (dps=50),
      starting from the stored high-precision string (not just the float re/im);
  (b) Newton return test: perturb u by |eps|=1e-3 in 4 directions
      (+1, -1, +i, -i), run p2.newton_refine(tol=1e-9) from u+eps, check the
      result lands back within 1e-8 of u;
  (c) level-7 partner: nearest distance to r1a_H4_level7.json's domain_ok
      catalog, partner iff distance < 5e-3 (same threshold as p2_H4_zeros
      pairing step).

Requires r1a_H4_level7.json (produced by r1a_level7_catalog.py) to exist.
Output: r1a_authenticity.json (per-catalog list of per-zero records + summary
counts of the classification).
"""
import json
import time

import numpy as np
from mpmath import mpc, mpf

import r1a_common as C

# nstr(u_mp, 25) produces e.g. "(-11.95073329739872527775069 + 16.24571531861244522678393j)"
# -- mpc()/mpf() can't parse this compound form directly (mpf.from_str chokes
# on '+'/'j'), so split on the SPACE-SURROUNDED '+'/'-' separator (the only
# unspaced signs are the real part's own leading '-' or an exponent's 'e-NN',
# neither of which is surrounded by spaces) into (real, imag) full-precision
# decimal strings, then build mpc(mpf(re_str), mpf(im_str)).

log = C.mk_logger('r1a_authenticity.log')

EPS = 1e-3
DIRECTIONS = [1.0, -1.0, 1.0j, -1.0j]
NEWTON_RETURN_TOL = 1e-8
PARTNER_TOL = 5e-3


def parse_u_mp(zrec):
    """Parse the high-precision complex string ('(-11.95... + 16.24...j)')
    stored in a catalog record into a full-precision mpc."""
    s = zrec.get('str')
    if s:
        body = s.strip()
        if body.startswith('(') and body.endswith(')'):
            body = body[1:-1]
        body = body.rstrip()
        assert body.endswith('j'), body
        body = body[:-1]  # drop trailing 'j'
        if ' + ' in body:
            re_s, im_s = body.split(' + ', 1)
        elif ' - ' in body:
            re_s, im_s0 = body.split(' - ', 1)
            im_s = '-' + im_s0
        else:
            raise ValueError('unrecognized mpc string: %r' % s)
        return mpc(mpf(re_s.strip()), mpf(im_s.strip()))
    return mpc(mpf(zrec['re']), mpf(zrec['im']))


def analyze_catalog(name, zeros_level8, AU0, level7_pts):
    log('--- %s: %d zeros ---' % (name, len(zeros_level8)))
    level7_arr = np.array(level7_pts) if level7_pts else np.empty((0, 2))
    recs = []
    t0 = time.time()
    for i, z in enumerate(zeros_level8):
        u_float = complex(z['re'], z['im'])
        # (a) independent high-precision |H_8|
        try:
            u_mp = parse_u_mp(z)
            absH8_mp = C.H8_mp_at(u_mp, AU0, m=8)
        except Exception as exc:  # noqa: BLE001
            absH8_mp = None
            log('  [%d] mp eval failed: %s' % (i, exc))
        # (b) Newton-return test, 4 perturbed directions
        returns = []
        for d in DIRECTIONS:
            u0 = u_float + EPS * d
            u_ret, res = C.newton_from(u0, 8, AU0)
            if u_ret is None:
                returns.append({'dir': str(d), 'returned': False, 'dist': None, 'absH': None})
                continue
            dist = abs(u_ret - u_float)
            returns.append({'dir': str(d), 'returned': bool(dist < NEWTON_RETURN_TOL),
                            'dist': float(dist), 'absH': float(res)})
        n_return = sum(1 for r in returns if r['returned'])
        # (c) level-7 partner
        if len(level7_arr) > 0:
            d7 = np.abs(level7_arr - np.array([u_float.real, u_float.imag])).max(axis=1)
            nearest7 = float(d7.min())
        else:
            nearest7 = None
        has_partner = bool(nearest7 is not None and nearest7 < PARTNER_TOL)

        # 'genuine' = it IS a true zero of H_8 (tiny |H|, Newton is locally
        # attracting to it).  This is deliberately independent of the level-7
        # pairing filter, which is a SEPARATE (survival, not genuineness)
        # criterion -- see 'has_partner' / 'would_survive_orig_pipeline'.
        genuine = bool((absH8_mp is not None and absH8_mp < 1e-15) and n_return >= 3)
        would_survive_orig_pipeline = bool(genuine and has_partner)
        recs.append({
            're': z['re'], 'im': z['im'],
            'absH8_stored': z.get('absH'),
            'absH8_mp_independent': absH8_mp,
            'newton_return': returns,
            'n_return_of_4': n_return,
            'level7_nearest_dist': nearest7,
            'level7_partner': has_partner,
            'genuine': genuine,
            'would_survive_orig_pipeline': would_survive_orig_pipeline,
        })
        if (i + 1) % 40 == 0:
            log('  %d/%d (%.0fs)' % (i + 1, len(zeros_level8), time.time() - t0))
    log('%s done in %.1fs' % (name, time.time() - t0))
    return recs


def main():
    AU0, Ls = C.load_au0()
    lvl7 = json.load(open('r1a_H4_level7.json', encoding='utf-8'))
    level7_pts = [(z['re'], z['im']) for z in lvl7['zeros_level7_domain_ok']]
    log('loaded level-7 catalog: %d domain_ok points' % len(level7_pts))

    _, orig_zl, _ = C.load_catalog(C.ORIG_PATH)
    _, step24_zl, _ = C.load_catalog('H4_step24.json')

    out = {'meta': {'eps': EPS, 'newton_return_tol': NEWTON_RETURN_TOL,
                     'partner_tol': PARTNER_TOL,
                     'genuine_def': 'absH8_mp<1e-15 AND n_return>=3/4 (level-7 partner tracked separately)',
                     'survive_def': 'genuine AND level7_partner (== would pass p2_H4_zeros pairing stage)',
                     'date': time.strftime('%Y-%m-%d %H:%M:%S')}}

    out['ORIG'] = analyze_catalog('ORIG(62)', orig_zl, AU0, level7_pts)
    out['STEP24'] = analyze_catalog('STEP24(146)', step24_zl, AU0, level7_pts)

    for nm in ('ORIG', 'STEP24'):
        recs = out[nm]
        n = len(recs)
        n_genuine = sum(1 for r in recs if r['genuine'])
        n_hp_ok = sum(1 for r in recs if r['absH8_mp_independent'] is not None and r['absH8_mp_independent'] < 1e-15)
        n_newton4 = sum(1 for r in recs if r['n_return_of_4'] == 4)
        n_newton3p = sum(1 for r in recs if r['n_return_of_4'] >= 3)
        n_partner = sum(1 for r in recs if r['level7_partner'])
        n_survive = sum(1 for r in recs if r['would_survive_orig_pipeline'])
        out[nm + '_summary'] = dict(n=n, n_genuine=n_genuine,
                                     n_highprec_absH_lt_1e15=n_hp_ok,
                                     n_newton_return_4of4=n_newton4,
                                     n_newton_return_3plus=n_newton3p,
                                     n_level7_partner=n_partner,
                                     n_would_survive_orig_pipeline=n_survive)
        log('%s summary: %s' % (nm, out[nm + '_summary']))

    C.dump_json(out, 'r1a_authenticity.json')
    log('wrote r1a_authenticity.json')


if __name__ == '__main__':
    main()

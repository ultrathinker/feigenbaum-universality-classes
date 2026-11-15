#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fix_pos_n40.py -- clean N=40 partners for the pos branch at high d (the
truncated-series shortcut degrades there), plus N=56 spot checks at d=5,6,7
(tail calibration, degree-6 lesson of agent B: series convergence slows with d).
[oid:spectraF] NOT RIGOROUS.
"""
import json
import os
import time

from mpmath import mp, mpf, nstr

import spectra_grid as sg

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'spectra_grid.json')
FIX_AT = [r['d'] for r in json.load(open(OUTJ, encoding='utf-8'))['points']]
N56_AT = ['5.0', '6.0', '7.0']


def main():
    t0 = time.time()
    data = sg.load_points()
    pts = {r['d']: r for r in data['points']}
    for ds in FIX_AT:
        rec = pts[ds]
        d = mpf(ds)
        a48 = [mpf(s) for s in rec['_a48']]
        al48 = mpf(rec['alpha48'])
        # neighbour-d seed (different d -> no rounding-floor stall)
        nb = min((r for r in data['points'] if r['d'] != ds),
                 key=lambda r: abs(float(r['d']) - float(ds)))
        a_nb = [mpf(s) for s in nb['_a48']]
        al_nb = mpf(nb['alpha48'])
        a40, al40, _it, _r = sg.solve_fp(d, sg.N1, a_nb, al_nb)
        assert abs(al40 - mpf(rec['alpha40'])) < mpf('1e-25'), 'alpha40 mismatch'
        A40 = sg.build_DT_matrix(a40, al40, d, sg.N1)
        sigma = mpf(rec['pos']['N48']['lam'])
        lam = v = mres = None
        for jit in (0, '3e-18', '-3e-18', '1e-12', '-1e-12'):
            sig_try = sigma if jit == 0 else sigma * (1 + mpf(jit))
            lam, v, mres = sg.inverse_iter(A40, sg.N1, sig_try)
            if lam is not None:
                break
        assert lam is not None, 'inverse iteration failed at N40 for d=%s' % ds
        c = [mpf(0)] + [x / sum(v) for x in v]
        dres, phmax = sg.dense_eigen_residual(a40, al40, d, c, lam)
        rec['pos']['N40'] = {'lam': nstr(lam, 30), 'lam_f': float(lam),
                             'mat_resid': float(mres), 'dense_resid': float(dres),
                             'dense_resid_rel': float(dres / max(phmax, mpf('1e-30')))}
        rec['pos']['N40_via'] = 'solved N=40 fixed point (neighbour-d seeded)'
        lam48 = mpf(rec['pos']['N48']['lam'])
        rec['pos']['drift'] = float(abs(lam - lam48))
        rec['pos']['drift_rel'] = float(abs(lam - lam48) / abs(lam48))
        rec['pos']['ok'] = bool(float(dres) < 1e-18 and rec['pos']['drift'] < 1e-6)
        if ds in N56_AT:
            a56, al56, _i2, _r2 = sg.solve_fp(d, 56, a_nb, al_nb)
            A56 = sg.build_DT_matrix(a56, al56, d, 56)
            lam56, v56, mres56 = sg.inverse_iter(A56, 56, sigma)
            if lam56 is not None:
                c56 = [mpf(0)] + [x / sum(v56) for x in v56]
                dres56, phmax56 = sg.dense_eigen_residual(a56, al56, d, c56, lam56)
                rec['pos']['N56'] = {'lam': nstr(lam56, 30), 'lam_f': float(lam56),
                                     'mat_resid': float(mres56),
                                     'dense_resid': float(dres56),
                                     'drift_vs_N48': float(abs(lam56 - lam48))}
        sg.save_points(data)
        sg.plog('FIXPOS d=%s: N40=%s drift=%.2e dense=%.1e  %s'
                % (ds, nstr(lam, 14), rec['pos']['drift'], float(dres),
                   ('N56=%s drift56=%.1e' % (nstr(lam56, 12), rec['pos']['N56']['drift_vs_N48']))
                   if ds in N56_AT and lam56 is not None else ''))
    sg.plog('FIXPOS DONE in %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectrum_window.py -- [oid:vecM] Step 0: independent full-spectrum window
|lam| <= 0.42 at every d in [1.9, 4.0] already fixed-pointed by wave1
(grid+flank+bridge roles only, bisection-only points skipped).

For each d: load the already-certified fixed point (a48, alpha48) from
wave1_spectra_d/spectra_grid.json (NOT re-solved -- Newton is the expensive
part and wave1 already certified it to dense residual <=1e-30; only the
LINEAR spectral step is redone here, independently, via full mpmath eig
(right AND left eigenvectors) instead of wave1's targeted inverse iteration
to 4 named branches -- this is what lets us see ALL modes in the window,
not just the ones wave1 chose to track by continuation).

N-stability ladder (pseudospectral-pollution filter, same principle as
wave1/agentB_deg6): compute at N_LO=28 and N_HI=40; keep only eigenvalues
that agree between the two truncations to rel. tol 1e-6 (physical modes are
N-stable; truncation-edge pollution N-wanders -- verified on this exact
problem: a N-wandering ghost at d=4.0 moves -0.415(N20) -> -0.366(24) ->
-0.324(28) -> -0.288(32) -> gone by N=40, while all 5 genuine modes in that
window are already stable from N=24).

Output: spectrum_window.json, keyed by d string, each holding the STABLE
modes only: lam (N_HI value), rvec/lvec at N_HI (length-N_HI mpf lists, same
abstract {t,...,t^N} monomial basis at every d -> directly comparable across
d without re-basing). Written atomically after every d. NOT RIGOROUS.
CPU: 1 process.
"""
import json
import os
import time

from mpmath import mp, mpf, nstr

import engine as E

mp.dps = 35

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'spectrum_window.json')
PROG = os.path.join(HERE, 'progress.log')
T0 = time.time()

N_LO = 28
N_HI = 40
WINDOW = 0.42
REL_TOL = mpf('1e-6')

DGRID = ['1.9', '1.95', '1.99', '2.0', '2.001', '2.01', '2.05', '2.1', '2.2',
         '2.3', '2.4', '2.5', '2.6', '2.7', '2.8', '2.9', '3.0', '3.1', '3.2',
         '3.3', '3.4', '3.5', '3.6', '3.7', '3.8', '3.9', '4.0']


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def load_out():
    if os.path.exists(OUTJ):
        with open(OUTJ, encoding='utf-8') as f:
            return json.load(f)
    return {'meta': {'tag': '[oid:vecM] BRIEF-M spectrum window, NOT RIGOROUS',
                      'N_LO': N_LO, 'N_HI': N_HI, 'window': WINDOW,
                      'rel_tol_stability': str(REL_TOL), 'dps': mp.dps,
                      'grid': DGRID,
                      'source': 'fixed points reused verbatim from '
                                'wave1_spectra_d/spectra_grid.json (_a48, alpha48)'},
            'points': {}}


def save_out(data):
    tmp = OUTJ + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, OUTJ)


def stable_modes(a, al, d):
    lo = E.full_eig_window(a, al, d, N_LO, window=WINDOW)
    hi = E.full_eig_window(a, al, d, N_HI, window=WINDOW)
    stable = []
    pollution = []
    for m in hi:
        if not lo:
            pollution.append(m['lam'])
            continue
        near = min(lo, key=lambda x: abs(x['lam'] - m['lam']))
        if abs(near['lam'] - m['lam']) < REL_TOL * max(abs(m['lam']), mpf('1e-8')):
            stable.append(m)
        else:
            pollution.append(m['lam'])
    stable.sort(key=lambda r: -abs(r['lam']))
    return stable, pollution


def main():
    plog('SPECTRUM WINDOW START: %d d-points, N_LO=%d N_HI=%d, window=%.2f, dps=%d'
         % (len(DGRID), N_LO, N_HI, WINDOW, mp.dps))
    wave1 = E.load_wave1_points()
    data = load_out()
    times = []
    for i, ds in enumerate(DGRID):
        if ds in data['points']:
            plog('d=%s: already done, skip' % ds)
            continue
        t0 = time.time()
        rec = wave1[ds]
        a48, al48 = E.fixedpoint_from_wave1(rec)
        d = mpf(ds)
        stable, pollution = stable_modes(a48, al48, d)
        pt = {'d': ds,
              'n_stable': len(stable),
              'lams': [nstr(m['lam'], 25) for m in stable],
              'rvecs': [[nstr(x, 22) for x in m['rvec']] for m in stable],
              'lvecs': [[nstr(x, 22) for x in m['lvec']] for m in stable],
              'n_pollution_rejected': len(pollution),
              'pollution_lams': [nstr(p, 10) for p in pollution],
              'elapsed_s': round(time.time() - t0, 1)}
        data['points'][ds] = pt
        save_out(data)
        times.append(pt['elapsed_s'])
        eta = sum(times[-5:]) / len(times[-5:]) * (len(DGRID) - i - 1)
        plog('d=%s DONE (%.1fs): %d stable modes %s | %d pollution rejected %s '
             '[ETA %.0fs]'
             % (ds, pt['elapsed_s'], len(stable),
                '[' + ', '.join(x[:12] for x in pt['lams'][:8]) + (' ...]' if len(stable) > 8 else ']'),
                len(pollution), pt['pollution_lams'][:3], eta))
    plog('SPECTRUM WINDOW COMPLETE in %.0fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

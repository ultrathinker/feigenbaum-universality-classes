#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
refine_bridge.py -- [oid:vecM] BRIEF-M Part 1.3: fine-step refinement.

pairwise_map.py found the pos/m3 raw-cosine confusion is SYSTEMATIC across
the whole d in [2,4] range (not special to d=2) -- but ALSO found a genuine
near-tie right in the disputed bridge: neg(2.01) best-matches pos(2.05) with
cosine 0.8923, gap only 0.020 over the wave1 grid's coarse step (0.04). This
script tests whether that ambiguity is a coarse-step artifact (should vanish,
gap growing, as the step shrinks -- same-branch overlap -> 1 while cross-
branch overlap stays pinned near the ~0.8-0.9 background non-orthogonality
floor measured elsewhere) or a genuine phenomenon (ambiguity persists/worsens
at fine steps).

Solves 3 NEW fixed points (Newton, warm-started continuation, own code in
engine.py -- nothing written to wave1's folder) at d = 2.02, 2.03, 2.04,
bridging 2.01 -> 2.05 at 0.01 steps instead of 0.04. For each, computes the
same N_LO/N_HI=28/40 N-stability-filtered window as spectrum_window.py, and
appends to spectrum_window.json (same schema, so pairwise_map.py can be
re-run unmodified over the enlarged grid).

progress.log with flush + ETA. NOT RIGOROUS. CPU: 1 process.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, nstr

import engine as E

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

mp.dps = 35

HERE = os.path.dirname(os.path.abspath(__file__))
WINJ = os.path.join(HERE, 'spectrum_window.json')
PROG = os.path.join(HERE, 'progress.log')
T0 = time.time()

N_LO = 28
N_HI = 40
WINDOW = 0.42
REL_TOL = mpf('1e-6')

NEW_DS = ['2.02', '2.03', '2.04']


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def stable_modes(a, al, d):
    lo = E.full_eig_window(a, al, d, N_LO, window=WINDOW)
    hi = E.full_eig_window(a, al, d, N_HI, window=WINDOW)
    stable = []
    for m in hi:
        if not lo:
            continue
        near = min(lo, key=lambda x: abs(x['lam'] - m['lam']))
        if abs(near['lam'] - m['lam']) < REL_TOL * max(abs(m['lam']), mpf('1e-8')):
            stable.append(m)
    stable.sort(key=lambda r: -abs(r['lam']))
    return stable


def main():
    plog('REFINE BRIDGE START: solving %s (Newton, warm-started) + N-stable window' % NEW_DS)
    with open(WINJ, encoding='utf-8') as f:
        data = json.load(f)

    # seed fixed point: nearest already-solved wave1 point (2.01, N=48 coeffs)
    wave1 = E.load_wave1_points()
    seed_rec = wave1['2.01']
    a_seed, al_seed = E.fixedpoint_from_wave1(seed_rec)

    times = []
    a_cur, al_cur = a_seed, al_seed
    for i, ds in enumerate(NEW_DS):
        if ds in data['points']:
            plog('d=%s: already done, skip' % ds)
            a_cur, al_cur = [mpf(s) for s in data['points'][ds]['_a48']], mpf(data['points'][ds]['_al48'])
            continue
        t0 = time.time()
        a_new, al_new, it, r = E.solve_fp_continuation(ds, 48, a_cur, al_cur)
        dense_res = E.fp_residual_dense(a_new, al_new, mpf(ds))
        plog('d=%s: fixed point solved (Newton it=%d, res=%.2e, dense_res=%.2e)'
             % (ds, it, float(r), float(dense_res)))
        stable = stable_modes(a_new, al_new, mpf(ds))
        pt = {'d': ds,
              'n_stable': len(stable),
              'lams': [nstr(m['lam'], 25) for m in stable],
              'rvecs': [[nstr(x, 22) for x in m['rvec']] for m in stable],
              'lvecs': [[nstr(x, 22) for x in m['lvec']] for m in stable],
              '_a48': [nstr(x, 40) for x in a_new],
              '_al48': nstr(al_new, 40),
              'elapsed_s': round(time.time() - t0, 1)}
        data['points'][ds] = pt
        tmp = WINJ + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=1, ensure_ascii=False)
        os.replace(tmp, WINJ)
        a_cur, al_cur = a_new, al_new
        times.append(pt['elapsed_s'])
        eta = sum(times) / len(times) * (len(NEW_DS) - i - 1)
        plog('d=%s DONE (%.1fs): %d stable modes %s [ETA %.0fs]'
             % (ds, pt['elapsed_s'], len(stable),
                '[' + ', '.join(x[:12] for x in pt['lams'][:6]) + ']', eta))

    # update meta.grid to include new points, in order
    grid = data['meta']['grid']
    for ds in NEW_DS:
        if ds not in grid:
            grid.append(ds)
    grid.sort(key=float)
    data['meta']['grid'] = grid
    tmp = WINJ + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, WINJ)
    plog('REFINE BRIDGE COMPLETE in %.0fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

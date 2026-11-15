#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
multiplicity_d2.py -- [oid:vecM] BRIEF-M Part 2: multiplicity structure of the
eigenvalue alpha^-2 = 0.15962844038269977... at d=2.0, and how the cluster
disperses at the flanks d=1.99 / d=2.01 (fixed points reused from wave1,
identical convention as spectrum_window.py: same converged a48 fixed point
used for all N truncations, only the LINEAR operator matrix is resized).

For d in {1.99, 2.0, 2.01} and N in {24, 32, 40, 48}:
  - full mpmath eig (ALL eigenvalues, right+left, complex, UNFILTERED -- this
    is deliberately the raw spectrum, not the |lam|<=0.42 window used in
    spectrum_window.py, so a close complex-conjugate pair near the target
    would still show up);
  - target = alpha(2)^-2 = 0.159628440383 FIXED for all three d (not
    re-centred at each d's own tower value) so the window literally tracks
    where the d=2 cluster's members have moved to;
  - nearest real eigenvalue e1, second-nearest (real or from a complex pair)
    e2, gap = |e1-e2|; nearest non-negligible-imaginary eigenvalue to target
    (Jordan/near-degenerate-complex-pair signature check);
  - at N=48: dense functional residual of e1 and e2 (genuineness check),
    cosine angle between their right eigenvectors, and the condition number
    of the 2-column matrix [r1|r2] (via eigenvalues of the 2x2 Gram matrix
    r_i . r_j) -- large condition number / cos near 1 => near-parallel =>
    effectively rank-1 (Jordan-like); cos far from 1 and condition number
    O(1) => genuinely 2-dimensional invariant subspace.

Interpretation (per brief): simple (isolated e1, e2 far / doesn't exist) =>
favours B; (near-)double or Jordan-like (e2 close, ill-conditioned pair) =>
favours A.

Output: multiplicity_d2.json, incremental after every (d,N). progress.log
with flush + ETA. NOT RIGOROUS. CPU: 1 process.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, matrix, nstr
from mpmath import eig as mp_eig

import engine as E

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

mp.dps = 35

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'multiplicity_d2.json')
PROG = os.path.join(HERE, 'progress.log')
T0 = time.time()

TARGET = mpf('0.159628440382699770039480839795')  # alpha(2)^-2, exact wave1 value
DS = ['1.99', '2.0', '2.01']
NS = [24, 32, 40, 48]
WIN_REPORT = mpf('0.02')  # cluster window per brief


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
    return {'meta': {'tag': '[oid:vecM] BRIEF-M Part 2: multiplicity at d=2, NOT RIGOROUS',
                      'target_alpha2_pow_-2': str(TARGET), 'window_report': str(WIN_REPORT),
                      'Ns': NS, 'ds': DS, 'dps': mp.dps},
            'runs': {}}


def save_out(data):
    tmp = OUTJ + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=1, ensure_ascii=False)
    os.replace(tmp, OUTJ)


def raw_full_eig(a, alpha, d, N):
    A = E.build_DT_matrix(a, alpha, d, N)
    Evals, EL, ER = mp_eig(A, left=True, right=True)
    out = []
    for k in range(N):
        lam = Evals[k]
        rvec = [ER[i, k] for i in range(N)]
        lvec = [EL[i, k] for i in range(N)]
        out.append({'lam': lam, 'rvec': rvec, 'lvec': lvec})
    return out


def cosine_c(v1, v2):
    """cosine on the real parts (both should be ~real for genuine near-target modes)."""
    r1 = [x.real for x in v1]
    r2 = [x.real for x in v2]
    n1 = sum(x * x for x in r1) ** mpf('0.5')
    n2 = sum(x * x for x in r2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(r1[i] * r2[i] for i in range(len(r1)))
    return abs(dot / (n1 * n2))


def main():
    plog('MULTIPLICITY-D2 START: ds=%s Ns=%s target=%s window=%s dps=%d'
         % (DS, NS, nstr(TARGET, 18), WIN_REPORT, mp.dps))
    wave1 = E.load_wave1_points()
    data = load_out()
    todo = [(ds, N) for ds in DS for N in NS]
    times = []
    for idx, (ds, N) in enumerate(todo):
        key = '%s_N%d' % (ds, N)
        if key in data['runs']:
            plog('%s: already done, skip' % key)
            continue
        t0 = time.time()
        rec = wave1[ds]
        a48, al48 = E.fixedpoint_from_wave1(rec)
        d = mpf(ds)
        full = raw_full_eig(a48, al48, d, N)

        # sort all by distance to TARGET
        full_sorted = sorted(full, key=lambda m: abs(m['lam'] - TARGET))
        # split real (small imag) vs complex
        def is_real(m, tol=mpf('1e-15')):
            return abs(m['lam'].imag) < tol * max(abs(m['lam']), mpf('1e-10')) + mpf('1e-40')

        reals = [m for m in full_sorted if is_real(m)]
        complexes = [m for m in full_sorted if not is_real(m)]

        e1 = reals[0] if reals else None
        e2 = reals[1] if len(reals) > 1 else None
        nearest_complex = complexes[0] if complexes else None

        in_window = [m for m in reals if abs(m['lam'] - TARGET) < WIN_REPORT]

        run = {'d': ds, 'N': N,
               'n_in_window_0.02': len(in_window),
               'in_window_lams': [nstr(m['lam'].real, 20) for m in in_window],
               'e1_lam': nstr(e1['lam'].real, 22) if e1 else None,
               'e1_dist_to_target': float(abs(e1['lam'] - TARGET)) if e1 else None,
               'e2_lam': nstr(e2['lam'].real, 22) if e2 else None,
               'e2_dist_to_target': float(abs(e2['lam'] - TARGET)) if e2 else None,
               'gap_e1_e2': float(abs(e1['lam'] - e2['lam'])) if (e1 and e2) else None,
               'nearest_complex_lam': nstr(nearest_complex['lam'], 20) if nearest_complex else None,
               'nearest_complex_dist_to_target': float(abs(nearest_complex['lam'] - TARGET)) if nearest_complex else None,
               'elapsed_s': round(time.time() - t0, 1)}

        if N == max(NS) and e1 is not None and e2 is not None:
            # dense functional residuals (genuineness) + pair geometry at largest N
            def hcoef_of(vec):
                return [mpf(0)] + [x.real for x in vec]
            res1, ph1 = E.dense_eigen_residual(a48, al48, d, hcoef_of(e1['rvec']), e1['lam'].real, npts=121)
            res2, ph2 = E.dense_eigen_residual(a48, al48, d, hcoef_of(e2['rvec']), e2['lam'].real, npts=121)
            cos12 = cosine_c(e1['rvec'], e2['rvec'])
            # 2x2 Gram matrix of [r1|r2] (real parts, raw not normalised) -> condition number via eigenvalues
            r1 = [x.real for x in e1['rvec']]
            r2 = [x.real for x in e2['rvec']]
            g11 = sum(x * x for x in r1)
            g22 = sum(x * x for x in r2)
            g12 = sum(r1[i] * r2[i] for i in range(len(r1)))
            Gm = matrix([[g11, g12], [g12, g22]])
            Geig = mp_eig(Gm, left=False, right=False)
            Geig = sorted([abs(x) for x in Geig])
            cond = float((Geig[1] / Geig[0]) ** mpf('0.5')) if Geig[0] != 0 else float('inf')
            run['pair_geometry_at_maxN'] = {
                'e1_dense_resid': float(res1), 'e1_phmax': float(ph1),
                'e2_dense_resid': float(res2), 'e2_phmax': float(ph2),
                'cos_angle_e1_e2': float(cos12),
                'gram_eigs': [float(x) for x in Geig],
                'condition_number_sqrt_ratio': cond}

        data['runs'][key] = run
        save_out(data)
        times.append(run['elapsed_s'])
        eta = sum(times[-5:]) / len(times[-5:]) * (len(todo) - idx - 1)
        plog('%s DONE (%.1fs): n_in_window=%d e1=%s (d=%.2e) e2=%s (d=%.2e) gap=%s '
             'nearest_complex=%s(d=%.2e) [ETA %.0fs]'
             % (key, run['elapsed_s'], run['n_in_window_0.02'],
                (run['e1_lam'] or '-')[:14], run['e1_dist_to_target'] or -1,
                (run['e2_lam'] or '-')[:14], run['e2_dist_to_target'] or -1,
                ('%.3e' % run['gap_e1_e2']) if run['gap_e1_e2'] is not None else '-',
                (run['nearest_complex_lam'] or '-')[:20], run['nearest_complex_dist_to_target'] or -1,
                eta))
        if 'pair_geometry_at_maxN' in run:
            pg = run['pair_geometry_at_maxN']
            plog('   pair geometry @N=%d: cos(e1,e2)=%.6f cond=%.3e e1_resid=%.2e e2_resid=%.2e'
                 % (N, pg['cos_angle_e1_e2'], pg['condition_number_sqrt_ratio'],
                    pg['e1_dense_resid'], pg['e2_dense_resid']))
    plog('MULTIPLICITY-D2 COMPLETE in %.0fs' % (time.time() - T0))


if __name__ == '__main__':
    main()

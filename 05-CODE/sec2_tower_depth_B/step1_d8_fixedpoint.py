#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step1_d8_fixedpoint.py -- BRIEF-S3-F5 task 1, part A.
Newton continuation d=7.0 (_a48 from wave1_spectra_d/spectra_grid.json) -> d=8.0,
at N=48 AND N=56, dps=40. Own re-typed engine (engine_own.py), independent of
MiniMax's b1_4_F5_depth/engine.py (dps=30, N=48 only for the fixed point).

Acceptance per BRIEF-S3-F5: dense residual <= 1e-25 at both N; alpha drift
between N=48 and N=56 solutions <= 1e-15.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, 'step1_d8_fixedpoint.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import load_wave1_points, solve_fp_continuation, fp_residual_dense, newton_fp  # noqa


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def main():
    plog('=== STEP1 d8 fixedpoint: START (dps=40) ===')
    pts = load_wave1_points()
    rec7 = pts['7.0']
    a_seed = [mpf(s) for s in rec7['_a48']]
    al_seed = mpf(rec7['alpha48'])
    plog(f'seed d=7.0: alpha48(wave1,dps30)={rec7["alpha48"]}')
    plog(f'MiniMax d=8 alpha (dps30, N=48 continuation): 1.35801727913805035487376 (for comparison)')

    out = {'meta': {'dps': 40, 'tag': 'S3-F5 step1 d8 fixed point, N=48 and N=56',
                     'created': time.strftime('%Y-%m-%d %H:%M:%S')},
           'points': {}}

    # --- N=48 ---
    t0 = time.time()
    plog('--- solving d=8.0 at N=48 (continuation ladder 10/18/28/40/48 from d=7 seed) ---')
    a48, al48, it48, r48 = solve_fp_continuation(mpf('8.0'), 48, a_seed, al_seed)
    dr48 = fp_residual_dense(a48, al48, mpf('8.0'), npts=801)
    plog(f'N=48: newton_it={it48} newton_resid={nstr(r48,6)} dense_resid={nstr(dr48,6)} '
         f'alpha={nstr(al48,30)} elapsed={time.time()-t0:.1f}s')

    # --- N=56, warm-started from N=48 solution ---
    t1 = time.time()
    plog('--- extending to N=56 (warm start from N=48 solution, Newton refine) ---')
    a56_seed = list(a48) + [mpf(0)] * (57 - len(a48))
    a56, al56, it56, r56 = newton_fp(mpf('8.0'), 56, a56_seed, al48)
    dr56 = fp_residual_dense(a56, al56, mpf('8.0'), npts=801)
    plog(f'N=56: newton_it={it56} newton_resid={nstr(r56,6)} dense_resid={nstr(dr56,6)} '
         f'alpha={nstr(al56,30)} elapsed={time.time()-t1:.1f}s')

    drift = abs(al56 - al48)
    plog(f'alpha drift |alpha56 - alpha48| = {nstr(drift, 6)}  '
         f'(criterion <= 1e-15: {"PASS" if drift <= mpf("1e-15") else "FAIL"})')
    plog(f'dense residual N=48: {nstr(dr48,6)} (criterion <=1e-25: '
         f'{"PASS" if dr48 <= mpf("1e-25") else "FAIL"})')
    plog(f'dense residual N=56: {nstr(dr56,6)} (criterion <=1e-25: '
         f'{"PASS" if dr56 <= mpf("1e-25") else "FAIL"})')

    minimax_alpha = mpf('1.35801727913805035487376')
    plog(f'compare to MiniMax alpha(8)=1.35801727913805035487376: '
         f'diff(N48)={nstr(abs(al48-minimax_alpha),6)} diff(N56)={nstr(abs(al56-minimax_alpha),6)}')

    out['points']['N48'] = {
        'N': 48, 'alpha': nstr(al48, 35), 'newton_iters': it48,
        'newton_resid': float(r48), 'dense_resid': float(dr48),
        'a': [nstr(x, 35) for x in a48],
    }
    out['points']['N56'] = {
        'N': 56, 'alpha': nstr(al56, 35), 'newton_iters': it56,
        'newton_resid': float(r56), 'dense_resid': float(dr56),
        'a': [nstr(x, 35) for x in a56],
    }
    out['checks'] = {
        'alpha_drift_N48_N56': float(drift),
        'alpha_drift_pass_1e-15': bool(drift <= mpf('1e-15')),
        'dense_resid_N48_pass_1e-25': bool(dr48 <= mpf('1e-25')),
        'dense_resid_N56_pass_1e-25': bool(dr56 <= mpf('1e-25')),
        'minimax_alpha8': '1.35801727913805035487376',
        'diff_vs_minimax_N48': float(abs(al48 - minimax_alpha)),
        'diff_vs_minimax_N56': float(abs(al56 - minimax_alpha)),
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP1 d8 fixedpoint: DONE ===')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2b_n64_confirm.py -- third truncation N=64 to (a) confirm the |lam|~0.4295
mode found stable at N=48/56 is not itself an N=56 artifact (three-truncation
ladder 48/56/64), and (b) check whether the N=56-only mode near 0.098 (not
matched at N=48 in step2) persists at N=64 -- would raise n_>(8) to 8 if real.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
FP_JSON = os.path.join(HERE, 'step1_d8_fixedpoint.json')
OUT_JSON = os.path.join(HERE, 'step2b_n64_confirm.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import full_eig_window, newton_fp, fp_residual_dense  # noqa

WINDOW = 0.5
LAM_MIN = mpf('0.03')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def main():
    plog('=== STEP2b N=64 confirmation: START ===')
    with open(FP_JSON, encoding='utf-8') as f:
        fp = json.load(f)
    a56 = [mpf(s) for s in fp['points']['N56']['a']]
    al = mpf(fp['points']['N56']['alpha'])
    d = mpf('8.0')

    t0 = time.time()
    a64_seed = list(a56) + [mpf(0)] * (65 - len(a56))
    a64, al64, it64, r64 = newton_fp(d, 64, a64_seed, al)
    dr64 = fp_residual_dense(a64, al64, d, npts=801)
    plog(f'N=64 fixed point: it={it64} resid={nstr(r64,6)} dense_resid={nstr(dr64,6)} '
         f'alpha={nstr(al64,30)} drift_vs_N56={nstr(abs(al64-al),6)} elapsed={time.time()-t0:.1f}s')

    t1 = time.time()
    modes64 = full_eig_window(a64, al64, d, 64, window=WINDOW)
    modes64 = [m for m in modes64 if abs(m['lam']) > LAM_MIN]
    plog(f'N=64: full eig done ({time.time()-t1:.1f}s), {len(modes64)} modes |lam|>{LAM_MIN}')
    lam_list = [nstr(m['lam'], 12) for m in modes64]
    plog(f'N=64 lam list: {lam_list}')

    target_high = mpf('0.4294647751')
    target_amb = mpf('0.09819507302')
    tower1 = al64 ** (-d)

    def nearest(target):
        best = min(modes64, key=lambda m: abs(m['lam'] - target) if not m['complex'] else mpf('inf'))
        return best

    m_high = nearest(target_high)
    m_amb = nearest(target_amb)
    plog(f'nearest to 0.4295 candidate: lam={nstr(m_high["lam"],15)} '
         f'diff={nstr(abs(m_high["lam"]-target_high),6)}')
    plog(f'nearest to 0.0982 candidate: lam={nstr(m_amb["lam"],15)} '
         f'diff={nstr(abs(m_amb["lam"]-target_amb),6)}')

    rel_high = abs(m_high['lam'] - target_high) / abs(target_high)
    rel_amb = abs(m_amb['lam'] - target_amb) / abs(target_amb)
    plog(f'0.4295 mode: rel diff N56->N64 = {nstr(rel_high,6)} '
         f'({"CONFIRMED stable" if rel_high < mpf("1e-6") else "NOT confirmed at this tol"})')
    plog(f'0.0982 mode: rel diff N56->N64 = {nstr(rel_amb,6)} '
         f'({"present, close" if rel_amb < mpf("0.05") else "NOT present / moved a lot -> pollution"})')
    plog(f'tower1(N64) = {nstr(tower1,15)}, is 0.0982 candidate above tower1? '
         f'{abs(m_amb["lam"]) > tower1}')

    out = {
        'N': 64, 'alpha': nstr(al64, 30), 'dense_resid': float(dr64),
        'alpha_drift_vs_N56': float(abs(al64 - al)),
        'lam_list': lam_list,
        'high_mode_check': {'target': float(target_high), 'found': nstr(m_high['lam'], 15),
                             'rel_diff': float(rel_high), 'confirmed': bool(rel_high < mpf('1e-6'))},
        'ambiguous_mode_check': {'target': float(target_amb), 'found': nstr(m_amb['lam'], 15),
                                  'rel_diff': float(rel_amb),
                                  'plausibly_real': bool(rel_amb < mpf('0.05'))},
        'tower1': float(tower1),
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP2b N=64 confirmation: DONE ===')


if __name__ == '__main__':
    main()

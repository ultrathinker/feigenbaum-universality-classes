#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2c_full_confirm.py -- FINAL, careful pass for task 1 (BRIEF-S3-F5).
Lesson learned from step2/step2b: a fixed window cutoff (MiniMax used 0.42;
my own step2 used 0.5) can itself silently exclude a genuine, stable
essential mode -- exactly the kind of bug this brief asked to hunt for.
So here: window = 1.6 (well above the 0.69 candidate seen in the exploratory
window=2.0 probe, well below delta(8)~10.9), THREE truncations N=48/56/64,
full mpmath eig each, 3-way ladder-stability match (rel<1e-8 pairwise),
dense functional residual as ground truth for every surviving real mode.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
FP_JSON = os.path.join(HERE, 'step1_d8_fixedpoint.json')
OUT_JSON = os.path.join(HERE, 'step2c_full_confirm.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import full_eig_window, newton_fp, eigenpair_dense_residual  # noqa

WINDOW = 1.6
LAM_MIN = mpf('0.02')
TOWER_TOL = mpf('1e-8')
STABLE_REL_TOL = mpf('1e-8')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def get_modes(a, al, d, N):
    modes = full_eig_window(a, al, d, N, window=WINDOW)
    return [m for m in modes if abs(m['lam']) > LAM_MIN]


def match3(m48, m56, m64):
    """3-way stability match: for each N=56 mode, find nearest at N=48 and at
    N=64 (same real/complex type); keep only if BOTH pairwise rel diffs < tol."""
    stable = []
    for m56i in m56:
        cands48 = [m for m in m48 if m['complex'] == m56i['complex']]
        cands64 = [m for m in m64 if m['complex'] == m56i['complex']]
        if not cands48 or not cands64:
            continue
        best48 = min(cands48, key=lambda m: abs(m['lam'] - m56i['lam']))
        best64 = min(cands64, key=lambda m: abs(m['lam'] - m56i['lam']))
        denom48 = max(abs(m56i['lam']), abs(best48['lam']), mpf('1e-30'))
        denom64 = max(abs(m56i['lam']), abs(best64['lam']), mpf('1e-30'))
        rel48 = abs(m56i['lam'] - best48['lam']) / denom48
        rel64 = abs(m56i['lam'] - best64['lam']) / denom64
        if rel48 < STABLE_REL_TOL and rel64 < STABLE_REL_TOL:
            stable.append({'lam': m56i['lam'], 'complex': m56i['complex'],
                            'rel48': rel48, 'rel64': rel64,
                            'rvec': m56i['rvec'], 'lvec': m56i['lvec']})
    stable.sort(key=lambda z: -abs(z['lam']))
    return stable


def main():
    plog('=== STEP2c FULL CONFIRM: START (window=1.6, N=48/56/64, dps=40) ===')
    with open(FP_JSON, encoding='utf-8') as f:
        fp = json.load(f)
    a56 = [mpf(s) for s in fp['points']['N56']['a']]
    al = mpf(fp['points']['N56']['alpha'])
    d = mpf('8.0')

    t0 = time.time()
    a64_seed = list(a56) + [mpf(0)] * (65 - len(a56))
    a64, al64, it64, r64 = newton_fp(d, 64, a64_seed, al)
    plog(f'N=64 fixed point re-solved: it={it64} resid={nstr(r64,4)} '
         f'alpha_drift_vs_N56={nstr(abs(al64-al),4)} ({time.time()-t0:.1f}s)')

    t1 = time.time()
    m48 = get_modes(a56[:49], al, d, 48)
    plog(f'N=48: {len(m48)} modes |lam|>{LAM_MIN}, window<={WINDOW} ({time.time()-t1:.1f}s)')
    t1 = time.time()
    m56 = get_modes(a56, al, d, 56)
    plog(f'N=56: {len(m56)} modes ({time.time()-t1:.1f}s)')
    t1 = time.time()
    m64 = get_modes(a64, al64, d, 64)
    plog(f'N=64: {len(m64)} modes ({time.time()-t1:.1f}s)')

    plog('N=48 list: ' + ', '.join(nstr(m['lam'], 10) + ('C' if m['complex'] else '') for m in m48))
    plog('N=56 list: ' + ', '.join(nstr(m['lam'], 10) + ('C' if m['complex'] else '') for m in m56))
    plog('N=64 list: ' + ', '.join(nstr(m['lam'], 10) + ('C' if m['complex'] else '') for m in m64))

    stable = match3(m48, m56, m64)
    plog(f'3-way ladder-stable (rel<{STABLE_REL_TOL} both N48 and N64 vs N56): {len(stable)} modes')
    for s in stable:
        plog(f'  lam={nstr(s["lam"],16)} complex={s["complex"]} '
             f'rel48={nstr(s["rel48"],3)} rel64={nstr(s["rel64"],3)}')

    tower1 = al64 ** (-d)
    delta_guess = None  # delta not in this window by construction (window<10)
    plog(f'tower1 = alpha8^-8 = {nstr(tower1,18)}')

    classified = []
    for s in stable:
        absval = abs(s['lam'])
        is_tower = False
        k_d = None
        if not s['complex']:
            for k in range(1, 4):
                tk = al64 ** (-d * k)
                if tk > mpf('1e-30') and abs(absval - tk) / tk < TOWER_TOL:
                    is_tower = True
                    k_d = k
                    break
        is_delta = absval > mpf('1.0')
        is_ess = (not is_tower) and (not is_delta)
        dres = None
        if not s['complex']:
            try:
                dres = float(eigenpair_dense_residual(a64, al64, d, s['lam'], s['rvec'], npts=401))
            except Exception:
                dres = None
        classified.append({'lam': s['lam'], 'complex': s['complex'], 'is_tower': is_tower,
                            'k_d': k_d, 'is_delta': is_delta, 'is_essential': is_ess,
                            'dense_resid_rel': dres, 'rvec': s['rvec']})

    ess = sorted([c for c in classified if c['is_essential']], key=lambda z: -abs(z['lam']))
    n_ess_total = len(ess)
    n_gt = sum(1 for c in ess if abs(c['lam']) > tower1)
    plog(f'n_essential(total, stable, window<=1.6) = {n_ess_total}')
    plog('essential list (desc |lam|): ' +
         ', '.join(f'{nstr(c["lam"],14)}(dres={c["dense_resid_rel"]:.2e})' if c['dense_resid_rel'] is not None
                    else nstr(c['lam'], 14) for c in ess))
    plog(f'*** n_>(8) [3-truncation ladder-stable, window=1.6, dps=40] = {n_gt} ***')
    plog('for comparison: MiniMax reported n_>(8)=6 (window=0.42, ladder N=28/40, dps=30)')
    plog(f'MOST-LIKELY-MISSING-MODE audit: 0.4295 mode (>0.42=MiniMax window) '
         f'{"IS" if any(abs(abs(c["lam"])-mpf("0.4295"))<mpf("0.001") for c in ess) else "is NOT"} in our essential list')

    out = {
        'meta': {'window': WINDOW, 'lam_min': float(LAM_MIN), 'N': [48, 56, 64], 'dps': 40,
                 'stable_rel_tol': float(STABLE_REL_TOL), 'created': time.strftime('%Y-%m-%d %H:%M:%S')},
        'alpha_N64': nstr(al64, 30),
        'tower1': float(tower1),
        'n48_raw': [nstr(m['lam'], 14) + ('C' if m['complex'] else '') for m in m48],
        'n56_raw': [nstr(m['lam'], 14) + ('C' if m['complex'] else '') for m in m56],
        'n64_raw': [nstr(m['lam'], 14) + ('C' if m['complex'] else '') for m in m64],
        'stable_modes': [{'lam': nstr(c['lam'], 18), 'complex': c['complex'], 'is_tower': c['is_tower'],
                           'k_d': c['k_d'], 'is_delta': c['is_delta'], 'is_essential': c['is_essential'],
                           'dense_resid_rel': c['dense_resid_rel']} for c in classified],
        'n_essential_total': n_ess_total,
        'n_gt_8': n_gt,
        'essential_sorted': [nstr(c['lam'], 16) for c in ess],
        'minimax_n_gt_8_reported': 6,
        'minimax_window': 0.42,
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP2c FULL CONFIRM: DONE ===')


if __name__ == '__main__':
    main()

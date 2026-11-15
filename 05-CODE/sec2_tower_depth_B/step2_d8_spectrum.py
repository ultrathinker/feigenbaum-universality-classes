#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step2_d8_spectrum.py -- BRIEF-S3-F5 task 1, part B.
Full spectrum of DT at d=8 (fixed point from step1_d8_fixedpoint.json),
TWO truncations N=48 and N=56, dps=40, full mpmath eig (complex-capable),
window |lam| > 0.03 (lower than MiniMax's 0.05, to not miss borderline
modes), stability check rel 1e-8 between N=48/N=56 (tighter than MiniMax's
1e-6 between N=28/40). Classifies tower alpha^{-8k} (k=1,2,3), delta,
essential; computes n_>(8) with and without delta.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
IN_JSON = os.path.join(HERE, 'step1_d8_fixedpoint.json')
OUT_JSON = os.path.join(HERE, 'step2_d8_spectrum.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import full_eig_window, eigenpair_dense_residual  # noqa

WINDOW = 0.5
LAM_MIN = mpf('0.03')
TOWER_TOL = mpf('1e-8')
STABLE_REL_TOL = mpf('1e-8')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def classify_one(lam_val, alpha, d, is_cplx):
    out = {'is_tower_dk': False, 'k_d': None, 'is_delta': False}
    absval = abs(lam_val)
    if not is_cplx:
        for k in range(1, 4):
            t = alpha ** (-d * k)
            if t > mpf('1e-30') and abs(absval - t) / t < TOWER_TOL:
                out['is_tower_dk'] = True
                out['k_d'] = k
                break
    out['is_delta'] = absval > mpf('1.0')
    return out


def main():
    plog('=== STEP2 d8 spectrum: START (dps=40, N=48 & N=56, window>0.03) ===')
    with open(IN_JSON, encoding='utf-8') as f:
        fp = json.load(f)
    a56 = [mpf(s) for s in fp['points']['N56']['a']]
    al = mpf(fp['points']['N56']['alpha'])
    d = mpf('8.0')
    plog(f'loaded fixed point: alpha={nstr(al,30)}')

    results = {}
    for N in (48, 56):
        t0 = time.time()
        modes = full_eig_window(a56[:N + 1], al, d, N, window=WINDOW)
        # filter by |lam|>LAM_MIN (keep complex too if |lam|>LAM_MIN)
        modes = [m for m in modes if abs(m['lam']) > LAM_MIN]
        plog(f'N={N}: full eig done ({time.time()-t0:.1f}s), '
             f'{len(modes)} modes with |lam|>{LAM_MIN}')
        results[N] = modes

    m48, m56 = results[48], results[56]
    plog(f'N=48 lam list: {[nstr(m["lam"],10) for m in m48]}')
    plog(f'N=56 lam list: {[nstr(m["lam"],10) for m in m56]}')

    # stability match N48 <-> N56, rel tol 1e-8, real-vs-real / complex-vs-complex
    stable = []
    used56 = set()
    for m48i in m48:
        best_j, best_d = None, None
        for j, m56j in enumerate(m56):
            if j in used56:
                continue
            if m48i['complex'] != m56j['complex']:
                continue
            diff = abs(m48i['lam'] - m56j['lam'])
            if best_d is None or diff < best_d:
                best_d, best_j = diff, j
        if best_j is None:
            continue
        m56j = m56[best_j]
        denom = max(abs(m48i['lam']), abs(m56j['lam']), mpf('1e-30'))
        rel = abs(m48i['lam'] - m56j['lam']) / denom
        if rel < STABLE_REL_TOL:
            used56.add(best_j)
            stable.append({'lam48': m48i['lam'], 'lam56': m56j['lam'], 'rel_diff': rel,
                            'complex': m48i['complex'], 'rvec56': m56j['rvec'],
                            'lvec56': m56j['lvec']})

    stable.sort(key=lambda z: -abs(z['lam56']))
    plog(f'stable (rel<{STABLE_REL_TOL}) modes N48<->N56: {len(stable)}')
    for s in stable:
        plog(f'  lam56={nstr(s["lam56"],15)} rel_diff={nstr(s["rel_diff"],4)} '
             f'complex={s["complex"]}')

    # unstable (present at N=48 or N=56 but not matched) -- report for audit
    unmatched48 = [m for m in m48 if not any(abs(m['lam'] - s['lam48']) < mpf('1e-20') for s in stable)]
    unmatched56 = [m for i, m in enumerate(m56) if i not in used56]
    plog(f'N=48 modes NOT matched to a stable pair: {[nstr(m["lam"],10) for m in unmatched48]}')
    plog(f'N=56 modes NOT matched to a stable pair: {[nstr(m["lam"],10) for m in unmatched56]}')

    # classify stable modes
    tower1 = al ** (-d)
    plog(f'tower1 = alpha^-8 = {nstr(tower1, 15)}')
    classified = []
    for s in stable:
        cls = classify_one(s['lam56'], al, d, s['complex'])
        cls['lam'] = s['lam56']
        cls['complex'] = s['complex']
        cls['rvec'] = s['rvec56']
        cls['lvec'] = s['lvec56']
        cls['is_essential'] = (not cls['is_tower_dk']) and (not cls['is_delta'])
        classified.append(cls)

    # dense functional residual check for each classified real mode (ground truth,
    # independent of the N=56 truncation that produced lam/vec)
    for cls in classified:
        if not cls['complex']:
            try:
                dres = eigenpair_dense_residual(a56, al, d, cls['lam'], cls['rvec'], npts=401)
                cls['dense_resid_rel'] = float(dres)
            except Exception as e:
                cls['dense_resid_rel'] = None
        else:
            cls['dense_resid_rel'] = None

    n_ess = sum(1 for c in classified if c['is_essential'])
    n_tower = sum(1 for c in classified if c['is_tower_dk'])
    n_delta = sum(1 for c in classified if c['is_delta'])
    n_gt_ess = sum(1 for c in classified if c['is_essential'] and abs(c['lam']) > tower1)
    n_gt_with_delta = sum(1 for c in classified if abs(c['lam']) > tower1)
    ess_sorted = sorted([c for c in classified if c['is_essential']], key=lambda z: -abs(z['lam']))

    plog(f'n_essential={n_ess} n_tower={n_tower} n_delta={n_delta}')
    plog(f'n_>(8) essential (no delta, no tower) = {n_gt_ess}')
    plog(f'n_>(8) with delta counted = {n_gt_with_delta}')
    plog('essential modes sorted by |lam| desc: ' +
         ', '.join(f'{nstr(c["lam"],10)}(dres={c["dense_resid_rel"]})' for c in ess_sorted))

    # position of tower1 among essentials + among all non-delta
    pos_ess = -1
    for i, c in enumerate(ess_sorted):
        if abs(c['lam']) > tower1:
            pos_ess = i + 1
    # tower1 slots in right after the last essential mode above it (1-indexed rank
    # it would occupy if inserted into the essential-only list)
    pos_ess = sum(1 for c in ess_sorted if abs(c['lam']) > tower1)

    plog(f'tower1 position among essentials (1-indexed, modes above it) = {pos_ess}')

    minimax_essentials = [-0.40617372393628219, 0.275446127163221061, -0.254095508051430718,
                           -0.201100184613817362, 0.157655033619081369, -0.118406994118896572]
    plog(f'MiniMax essentials (N=28/40,dps30): {minimax_essentials}')
    plog(f'MiniMax n_essential(reported)=6, all above tower1')

    out = {
        'meta': {'dps': 40, 'window': WINDOW, 'lam_min': float(LAM_MIN),
                 'stable_rel_tol': float(STABLE_REL_TOL), 'tower_tol': float(TOWER_TOL),
                 'created': time.strftime('%Y-%m-%d %H:%M:%S')},
        'alpha': nstr(al, 30),
        'tower1': float(tower1),
        'tower1_str': nstr(tower1, 20),
        'n48_raw': [nstr(m['lam'], 15) for m in m48],
        'n56_raw': [nstr(m['lam'], 15) for m in m56],
        'unmatched_N48': [nstr(m['lam'], 10) for m in unmatched48],
        'unmatched_N56': [nstr(m['lam'], 10) for m in unmatched56],
        'stable_modes': [{'lam': nstr(c['lam'], 18), 'complex': c['complex'],
                           'is_tower_dk': c['is_tower_dk'], 'k_d': c['k_d'],
                           'is_delta': c['is_delta'], 'is_essential': c['is_essential'],
                           'dense_resid_rel': c['dense_resid_rel']} for c in classified],
        'n_essential': n_ess,
        'n_tower': n_tower,
        'n_delta': n_delta,
        'n_gt_essential': n_gt_ess,
        'n_gt_with_delta': n_gt_with_delta,
        'tower1_position_in_essentials': pos_ess,
        'essential_list_sorted': [nstr(c['lam'], 18) for c in ess_sorted],
        'minimax_comparison': {
            'minimax_essentials': minimax_essentials,
            'minimax_n_essential': 6,
            'minimax_tower1': 0.08644853405611026,
        },
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP2 d8 spectrum: DONE ===')


if __name__ == '__main__':
    main()

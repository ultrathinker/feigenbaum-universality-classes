#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_neighborhood.py -- BRIEF-S3-F5 task 2.
Neighborhood of d=4 (3.98,3.99,4.0,4.01,4.02) and d=6 (5.98,6.0,6.02):
own build (independent 2nd code, per BRIEF-S3-F5 wording "..\\..\\wave2_check_d205
\\core.py OR своя сборка" -- own assembly used here: same equations, own
LU-cached matrix build, dps=40, N=40/48 two truncations, full eig with
eigenvectors for the cosine/"coordinate mode shadow" test).

For each d: Newton-continue fixed point from nearest wave1 grid anchor
(3.9->3.98/3.99/4.0, 4.0->4.01/4.02, 6.0->5.98/6.02), N=48 (truncate to 40
for the second truncation of the LINEAR step only, same convention as
b1_4_F5_depth/engine.py's fixedpoint_from_wave1).
Find modes within |lam - alpha^-d| < 0.02 at N=40 and N=48 (multiplicity,
stability). Sign of B(d)-alpha^-d on both sides of d=4 and d=6 (slope).
Coordinate-mode-shadow test: cosine(eigvec of branch B(d) at d=4.01,
eigvec of tower mode alpha_4^-4 at d=4.0) in the shared {t,...,t^N} basis.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, 'step3_neighborhood.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import (load_wave1_points, solve_fp_continuation, fp_residual_dense,
                         full_eig_window, cosine_overlap)  # noqa

WINDOW = 0.5
NEIGH = mpf('0.02')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def solve_at(d_val, seed_a, seed_al):
    a, al, it, r = solve_fp_continuation(mpf(d_val), 48, seed_a, seed_al)
    dr = fp_residual_dense(a, al, mpf(d_val), npts=601)
    return a, al, it, r, dr


def analyze(d_val, a, al):
    d = mpf(d_val)
    tower1 = al ** (-d)
    entry = {'d': d_val, 'alpha': nstr(al, 20), 'tower1': float(tower1)}
    modes_by_N = {}
    for N in (40, 48):
        modes = full_eig_window(a[:N + 1], al, d, N, window=WINDOW)
        near = [m for m in modes if abs(abs(m['lam']) - tower1) < NEIGH]
        near.sort(key=lambda m: abs(abs(m['lam']) - tower1))
        modes_by_N[N] = near
        entry[f'near_N{N}'] = [{'lam': nstr(m['lam'], 16), 'delta_to_tower1': float(abs(m['lam']) - tower1)}
                                for m in near]
    entry['multiplicity_N40'] = len(modes_by_N[40])
    entry['multiplicity_N48'] = len(modes_by_N[48])
    # nearest branch at N=48 (finest truncation) -- excluding an exact tower hit if d integer
    return entry, modes_by_N


def main():
    plog('=== STEP3 neighborhood d=4, d=6: START ===')
    w1 = load_wave1_points()
    seed39 = w1['3.9']
    seed40 = w1['4.0']
    seed60 = w1['6.0']
    a39 = [mpf(s) for s in seed39['_a48']]
    al39 = mpf(seed39['alpha48'])
    a40 = [mpf(s) for s in seed40['_a48']]
    al40 = mpf(seed40['alpha48'])
    a60 = [mpf(s) for s in seed60['_a48']]
    al60 = mpf(seed60['alpha48'])

    results = {}
    mode_store = {}

    # d=4 neighborhood: 3.98,3.99 from 3.9 seed; 4.0 exact from wave1 (refine anyway);
    # 4.01,4.02 from 4.0 seed
    plan4 = [('3.98', a39, al39), ('3.99', a39, al39), ('4.0', a40, al40),
             ('4.01', a40, al40), ('4.02', a40, al40)]
    for d_val, sa, sal in plan4:
        t0 = time.time()
        a, al, it, r, dr = solve_at(d_val, sa, sal)
        plog(f'd={d_val}: newton it={it} dense_resid={nstr(dr,4)} alpha={nstr(al,18)} '
             f'({time.time()-t0:.1f}s)')
        entry, modes_by_N = analyze(d_val, a, al)
        results[d_val] = entry
        mode_store[d_val] = modes_by_N
        near48 = entry['near_N48']
        plog(f'  d={d_val}: tower1={entry["tower1"]:.8f} near(N48, window<0.02): ' +
             ', '.join(f"{m['lam']}(Δ={m['delta_to_tower1']:.6f})" for m in near48))

    # d=6 neighborhood: 5.98,6.02 from 6.0 seed; 6.0 exact refine
    plan6 = [('5.98', a60, al60), ('6.0', a60, al60), ('6.02', a60, al60)]
    for d_val, sa, sal in plan6:
        t0 = time.time()
        a, al, it, r, dr = solve_at(d_val, sa, sal)
        plog(f'd={d_val}: newton it={it} dense_resid={nstr(dr,4)} alpha={nstr(al,18)} '
             f'({time.time()-t0:.1f}s)')
        entry, modes_by_N = analyze(d_val, a, al)
        results[d_val] = entry
        mode_store[d_val] = modes_by_N
        near48 = entry['near_N48']
        plog(f'  d={d_val}: tower1={entry["tower1"]:.8f} near(N48, window<0.02): ' +
             ', '.join(f"{m['lam']}(Δ={m['delta_to_tower1']:.6f})" for m in near48))

    # --- slope of B(d)-alpha^-d around d=4 and d=6 ---
    def branch_delta(d_val):
        near = mode_store[d_val][48]
        if not near:
            return None
        # take the one with smallest |delta| that is NOT the exact tower hit if d=int
        return near[0]

    plog('--- linearity / slope check around d=4 ---')
    for d_val in ['3.98', '3.99', '4.0', '4.01', '4.02']:
        b = branch_delta(d_val)
        tower1 = results[d_val]['tower1']
        if b is not None:
            delta = abs(b['lam']) if isinstance(b['lam'], (int, float)) else None
        near = results[d_val]['near_N48']
        plog(f'  d={d_val}: nearest-to-tower entries and signed delta: {near}')

    # explicit slope estimate using the OFF-tower branch (2nd nearest if 1st is the
    # exact tower hit) at d=3.99 and d=4.01 (symmetric pair)
    def offtower_delta(d_val):
        near = results[d_val]['near_N48']
        # exclude entries with |delta| < 1e-6 (those ARE the tower coincidence)
        cand = [m for m in near if abs(m['delta_to_tower1']) > 1e-6]
        return cand[0] if cand else (near[0] if near else None)

    d399 = offtower_delta('3.99')
    d401 = offtower_delta('4.01')
    d398 = offtower_delta('3.98')
    d402 = offtower_delta('4.02')
    plog(f'off-tower branch delta: d=3.98:{d398}, d=3.99:{d399}, d=4.01:{d401}, d=4.02:{d402}')
    if d399 and d401:
        slope = (d401['delta_to_tower1'] - d399['delta_to_tower1']) / (4.01 - 3.99)
        plog(f'slope estimate (3.99->4.01) of [ |branch|-tower1 ] vs d: {slope:.6f} per unit d')

    plog('--- linearity / slope check around d=6 ---')
    for d_val in ['5.98', '6.0', '6.02']:
        near = results[d_val]['near_N48']
        plog(f'  d={d_val}: {near}')
    d598 = offtower_delta('5.98')
    d602 = offtower_delta('6.02')
    if d598 and d602:
        slope6 = (d602['delta_to_tower1'] - d598['delta_to_tower1']) / (6.02 - 5.98)
        plog(f'slope estimate (5.98->6.02) of [|branch|-tower1] vs d: {slope6:.6f} per unit d')

    # --- coordinate-mode-shadow cosine test ---
    plog('--- coordinate-mode-shadow test: cosine(branch@4.01, tower-mode@4.0) ---')
    # tower mode eigenvector at d=4.0, N=48: the entry with |delta|<1e-6 (exact tower hit)
    near40_48 = mode_store['4.0'][48]
    tower_mode_40 = None
    for m in near40_48:
        if abs(abs(m['lam']) - results['4.0']['tower1']) < 1e-6:
            tower_mode_40 = m
            break
    if tower_mode_40 is None and near40_48:
        tower_mode_40 = near40_48[0]

    branch_401_48 = mode_store['4.01'][48]
    branch_mode_401 = branch_401_48[0] if branch_401_48 else None

    cos_val = None
    if tower_mode_40 is not None and branch_mode_401 is not None:
        cos_val = cosine_overlap(tower_mode_40['rvec'], branch_mode_401['rvec'])
        plog(f'tower mode @4.0: lam={nstr(tower_mode_40["lam"],12)}; '
             f'branch mode @4.01: lam={nstr(branch_mode_401["lam"],12)}')
        plog(f'cosine overlap = {nstr(cos_val, 10)} '
             f'({"~1 -> SAME mode (shadow confirmed)" if cos_val > mpf("0.99") else "NOT the same mode"})')
    else:
        plog('cosine test: could not find both modes')

    # also d=6 analogous cosine test (branch@6.02 vs tower@6.0)
    near60_48 = mode_store['6.0'][48]
    tower_mode_60 = None
    for m in near60_48:
        if abs(abs(m['lam']) - results['6.0']['tower1']) < 1e-6:
            tower_mode_60 = m
            break
    if tower_mode_60 is None and near60_48:
        tower_mode_60 = near60_48[0]
    branch_602_48 = mode_store['6.02'][48]
    branch_mode_602 = branch_602_48[0] if branch_602_48 else None
    cos_val6 = None
    if tower_mode_60 is not None and branch_mode_602 is not None:
        cos_val6 = cosine_overlap(tower_mode_60['rvec'], branch_mode_602['rvec'])
        plog(f'[d=6] tower mode @6.0: lam={nstr(tower_mode_60["lam"],12)}; '
             f'branch mode @6.02: lam={nstr(branch_mode_602["lam"],12)}')
        plog(f'[d=6] cosine overlap = {nstr(cos_val6, 10)}')

    out = {
        'meta': {'window_neigh': float(NEIGH), 'eig_window': WINDOW, 'N': [40, 48], 'dps': 40,
                 'created': time.strftime('%Y-%m-%d %H:%M:%S')},
        'points': results,
        'cosine_test_d4': {'cosine': float(cos_val) if cos_val is not None else None,
                            'tower_mode_d4_lam': nstr(tower_mode_40['lam'], 16) if tower_mode_40 else None,
                            'branch_mode_d401_lam': nstr(branch_mode_401['lam'], 16) if branch_mode_401 else None},
        'cosine_test_d6': {'cosine': float(cos_val6) if cos_val6 is not None else None,
                            'tower_mode_d6_lam': nstr(tower_mode_60['lam'], 16) if tower_mode_60 else None,
                            'branch_mode_d602_lam': nstr(branch_mode_602['lam'], 16) if branch_mode_602 else None},
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP3 neighborhood: DONE ===')


if __name__ == '__main__':
    main()

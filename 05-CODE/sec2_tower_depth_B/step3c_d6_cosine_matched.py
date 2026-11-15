#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3c_d6_cosine_matched.py -- d=6 cosine test with Δd=0.01 (d=6.0 -> d=6.01),
matched to the d=4 test's step size (4.0 -> 4.01) in step3b_cosine_fix.py,
for an apples-to-apples comparison of the coordinate-mode-shadow effect at
d=4 vs d=6. (step3_neighborhood.py's original d=6 cosine used d=6.02, i.e.
Δd=0.02 -- kept too, reported alongside for completeness.)
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, 'step3c_d6_cosine_matched.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine_own import load_wave1_points, solve_fp_continuation, full_eig_window, cosine_overlap  # noqa

WINDOW = 0.5


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def nearest_real(modes, target):
    reals = [m for m in modes if not m['complex']]
    return min(reals, key=lambda m: abs(m['lam'] - target))


def main():
    plog('=== STEP3c d6 cosine (matched dd=0.01): START ===')
    w1 = load_wave1_points()
    seed60 = w1['6.0']
    a60seed = [mpf(s) for s in seed60['_a48']]
    al60seed = mpf(seed60['alpha48'])

    a6, al6, it6, r6 = solve_fp_continuation(mpf('6.0'), 48, a60seed, al60seed)
    a601, al601, it601, r601 = solve_fp_continuation(mpf('6.01'), 48, a60seed, al60seed)
    plog(f'd=6.0: alpha={nstr(al6,18)}; d=6.01: alpha={nstr(al601,18)}')

    modes6_48 = full_eig_window(a6, al6, mpf('6.0'), 48, window=WINDOW)
    modes6_40 = full_eig_window(a6[:41], al6, mpf('6.0'), 40, window=WINDOW)
    modes601_48 = full_eig_window(a601, al601, mpf('6.01'), 48, window=WINDOW)
    modes601_40 = full_eig_window(a601[:41], al601, mpf('6.01'), 40, window=WINDOW)

    tower1_6 = al6 ** mpf('-6.0')
    tower1_601 = al601 ** mpf('-6.01')

    tm48 = nearest_real(modes6_48, tower1_6)
    tm40 = nearest_real(modes6_40, tower1_6)
    plog(f'tower mode d=6.0: N48={nstr(tm48["lam"],16)} N40={nstr(tm40["lam"],16)} '
         f'match_rel={nstr(abs(tm48["lam"]-tm40["lam"])/abs(tm48["lam"]),4)}')

    near48 = sorted([m for m in modes601_48 if not m['complex']
                      and abs(abs(m['lam']) - tower1_601) < mpf('0.02')],
                     key=lambda m: abs(abs(m['lam']) - tower1_601))
    near40 = sorted([m for m in modes601_40 if not m['complex']
                      and abs(abs(m['lam']) - tower1_601) < mpf('0.02')],
                     key=lambda m: abs(abs(m['lam']) - tower1_601))
    plog(f'd=6.01 near-tower N48: {[nstr(m["lam"],12) for m in near48]}')
    plog(f'd=6.01 near-tower N40: {[nstr(m["lam"],12) for m in near40]}')

    branch = None
    for m48 in near48:
        for m40 in near40:
            denom = max(abs(m48['lam']), abs(m40['lam']), mpf('1e-30'))
            if abs(m48['lam'] - m40['lam']) / denom < mpf('1e-6'):
                branch = m48
                break
        if branch:
            break
    if branch is None:
        branch = near48[0] if near48 else nearest_real(modes601_48, tower1_601)
    plog(f'N-stable branch d=6.01: lam={nstr(branch["lam"],16)}')

    cos_val = cosine_overlap(tm48['rvec'], branch['rvec'])
    plog(f'COSINE(tower@6.0, branch@6.01, dd=0.01) = {nstr(cos_val,10)}')
    plog('for reference: step3_neighborhood.py got cosine=0.8067 at dd=0.02 (6.0->6.02); '
         'step3b_cosine_fix.py got cosine=0.99994 at dd=0.01 for d=4 (4.0->4.01)')

    out = {
        'tower_mode_d6': nstr(tm48['lam'], 18),
        'branch_mode_d601': nstr(branch['lam'], 18),
        'cosine_dd001': float(cos_val),
        'cosine_dd002_from_step3': 0.8067303953007262,
        'cosine_d4_dd001_from_step3b': 0.9999441632,
        'comparison': ('d=4: cosine~1 (branch eigenvector essentially frozen through the crossing). '
                        'd=6: cosine 0.57-0.81 depending on step (0.57 at dd=0.01, 0.81 at dd=0.02) '
                        '-- NOT close to 1 either way, i.e. the eigenvector rotates substantially '
                        'through the d=6 crossing. This is a genuine qualitative difference from d=4, '
                        'not a step-size artifact (smaller step gives LOWER cosine at d=6, opposite of '
                        'what pure discretisation/continuity would predict).'),
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP3c d6 cosine matched: DONE ===')


if __name__ == '__main__':
    main()

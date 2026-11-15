#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3b_cosine_fix.py -- correction of the cosine ('coordinate-mode-shadow')
test from step3_neighborhood.py.

step3 found that at d=3.99, 4.01, 4.02 a SECOND entry appears inside the
|lam-tower1|<0.02 window at N=48 but is ABSENT at N=40 -- i.e. it fails the
basic N-stability check used everywhere else in this project, so it is
pseudospectral pollution, not a physical branch. The genuine, N-stable
physical branch is the OTHER entry (present at both N=40 and N=48), and it
is exactly the one MiniMax reported (symmetric +-0.000622 at d=4.0+-0.01).
step3's automatic branch_mode_401 picker took near_N48[0] (nearest to
tower1), which at d=4.01 is the SPURIOUS N48-only entry (lam=0.1227441885),
not the real branch (lam=0.1229760097448126, present at N=40 AND N=48).
This script redoes ONLY the d=4 cosine test with the correct eigenvector.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, 'step3b_cosine_fix.json')
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


def main():
    plog('=== STEP3b cosine fix: START ===')
    w1 = load_wave1_points()
    seed40 = w1['4.0']
    a40seed = [mpf(s) for s in seed40['_a48']]
    al40seed = mpf(seed40['alpha48'])

    a4, al4, it4, r4 = solve_fp_continuation(mpf('4.0'), 48, a40seed, al40seed)
    a401, al401, it401, r401 = solve_fp_continuation(mpf('4.01'), 48, a40seed, al40seed)
    plog(f'd=4.0: alpha={nstr(al4,18)}; d=4.01: alpha={nstr(al401,18)}')

    modes4_48 = full_eig_window(a4, al4, mpf('4.0'), 48, window=WINDOW)
    modes4_40 = full_eig_window(a4[:41], al4, mpf('4.0'), 40, window=WINDOW)
    modes401_48 = full_eig_window(a401, al401, mpf('4.01'), 48, window=WINDOW)
    modes401_40 = full_eig_window(a401[:41], al401, mpf('4.01'), 40, window=WINDOW)

    tower1_4 = al4 ** mpf('-4.0')
    tower1_401 = al401 ** mpf('-4.01')

    # tower mode at d=4.0: nearest to tower1_4, must match at both N
    def nearest(modes, target):
        reals = [m for m in modes if not m['complex']]
        return min(reals, key=lambda m: abs(m['lam'] - target))

    tm48 = nearest(modes4_48, tower1_4)
    tm40 = nearest(modes4_40, tower1_4)
    plog(f'tower mode @4.0: N48 lam={nstr(tm48["lam"],16)}, N40 lam={nstr(tm40["lam"],16)}, '
         f'match_rel={nstr(abs(tm48["lam"]-tm40["lam"])/abs(tm48["lam"]),4)}')

    # branch mode candidates @4.01 within 0.02 window, both N
    near401_48 = sorted([m for m in modes401_48 if not m['complex'] and abs(abs(m['lam']) - tower1_401) < mpf('0.02')],
                         key=lambda m: abs(abs(m['lam']) - tower1_401))
    near401_40 = sorted([m for m in modes401_40 if not m['complex'] and abs(abs(m['lam']) - tower1_401) < mpf('0.02')],
                         key=lambda m: abs(abs(m['lam']) - tower1_401))
    plog(f'@4.01 near-tower entries N48: {[nstr(m["lam"],12) for m in near401_48]}')
    plog(f'@4.01 near-tower entries N40: {[nstr(m["lam"],12) for m in near401_40]}')

    # the physical branch = the N48 entry that matches (rel<1e-6) an N40 entry
    branch48 = None
    for m48 in near401_48:
        for m40 in near401_40:
            denom = max(abs(m48['lam']), abs(m40['lam']), mpf('1e-30'))
            if abs(m48['lam'] - m40['lam']) / denom < mpf('1e-6'):
                branch48 = m48
                break
        if branch48:
            break
    if branch48 is None:
        plog('WARNING: no N-stable branch found near tower at 4.01 within 0.02 -- using nearest N48')
        branch48 = near401_48[0] if near401_48 else nearest(modes401_48, tower1_401)

    plog(f'CONFIRMED physical branch mode @4.01 (N-stable): lam={nstr(branch48["lam"],16)}')

    cos_val = cosine_overlap(tm48['rvec'], branch48['rvec'])
    plog(f'CORRECTED cosine overlap(tower-mode@4.0, branch-mode@4.01) = {nstr(cos_val,10)} '
         f'({"~1 -> SAME mode (coordinate-mode shadow CONFIRMED)" if cos_val > mpf("0.9") else "NOT close to 1 -> different mode / hypothesis REJECTED"})')

    # also report cosine using the (wrong) N48-only pollution entry, for the record
    pollution_entry = None
    for m in near401_48:
        if abs(m['lam'] - branch48['lam']) > mpf('1e-6'):
            pollution_entry = m
            break
    cos_pollution = None
    if pollution_entry:
        cos_pollution = cosine_overlap(tm48['rvec'], pollution_entry['rvec'])
        plog(f'(for reference) cosine vs the N48-only POLLUTION entry lam={nstr(pollution_entry["lam"],12)}: '
             f'{nstr(cos_pollution,6)} (this is what step3.py mistakenly used)')

    out = {
        'tower_mode_d4': {'lam': nstr(tm48['lam'], 18), 'N40_match': nstr(tm40['lam'], 18)},
        'branch_mode_d401_physical': {'lam': nstr(branch48['lam'], 18)},
        'branch_mode_d401_pollution_N48only': (nstr(pollution_entry['lam'], 18) if pollution_entry else None),
        'cosine_corrected': float(cos_val),
        'cosine_using_pollution_entry_forref': (float(cos_pollution) if cos_pollution is not None else None),
        'verdict': ('coordinate-mode shadow CONFIRMED (cos~1)' if cos_val > mpf('0.9')
                    else 'hypothesis REJECTED (cos far from 1) -- branch is a genuinely different, independent essential mode that merely crosses tower1 in value, not the same eigenvector'),
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== STEP3b cosine fix: DONE ===')


if __name__ == '__main__':
    main()

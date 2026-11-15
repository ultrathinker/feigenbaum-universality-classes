#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r9_main.py -- BRIEF-R9-SHADOW-D6: cosine test of "coordinate-mode shadow" at
d=4 and d=6 for SMALLER steps Delta d (down to 0.0005), plus Phi_2
sanity check (k=2 tower mode) and N sensitivity.

Method (UNCHANGED from s3_F5_d8/engine_own.py):
  - engine.py (mpmath float, dps=40): same DT-linearization, same Newton
    continuation, same full eig (right+left eigenvectors), same cosine on
    shared {t,...,t^N} monomial basis.
  - WINDOW = 0.5 (eig window), NEIGH = 0.02 (window around |lam-tower1|).
  - branch B^+ = nearest entry to |tower1| inside NEIGH that is N-stable
    (present at both N=40 and N=48, rel<1e-6).
  - Phi_1 eigenvector = tower1 mode at the integer d=4.0 (resp. 6.0).
  - Phi_2 eigenvector = tower-2 mode at the integer d=4.0 (resp. 6.0),
    eigenvalue alpha^{-2d} -- to confirm B^+ is NOT confused with Phi_2.

For each d in the plan: list 3-4 nearest real modes to |tower1| within
NEIGH, with |lambda|, delta, cos(Phi_1), cos(Phi_2), and flag the
N-stable branch B^+ explicitly.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr

mp.dps = 40

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, 'r9_results.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine import (load_wave1_points, solve_fp_continuation, full_eig_window,
                    cosine_overlap)  # noqa

WINDOW = 0.5
NEIGH = mpf('0.02')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def rvec_real(m):
    """Take real part of rvec (mpmath eig can return mpc for real eigenvalues
    since DT matrix is non-symmetric); for non-defective real pairs it IS real."""
    rv = m['rvec']
    if hasattr(rv[0], 'real'):
        return [mpf(x.real) for x in rv]
    return [mpf(x) for x in rv]


def analyze_point(d_val, a, al, phi1_rvec, phi2_rvec, phi1_lam_at_int, phi2_lam_at_int,
                  Ns=(40, 48), neigh=NEIGH, top_k=4):
    """For point d_val: tower1, list of up to top_k nearest real modes to |tower1|,
    each with cos(Phi_1) and cos(Phi_2). Branch B^+ flagged as N-stable."""
    d = mpf(d_val)
    tower1 = al ** (-d)
    by_N = {}
    for N in Ns:
        aN = list(a[:N + 1]) + [mpf(0)] * max(0, N + 1 - len(a))
        by_N[N] = full_eig_window(aN, al, d, N, window=WINDOW)
    # N-stable branch: nearest entry at Ns[-1] that matches Ns[0] within rel<1e-6
    near_fine = sorted([m for m in by_N[Ns[-1]] if not m['complex']
                        and abs(abs(m['lam']) - tower1) < neigh],
                       key=lambda m: abs(abs(m['lam']) - tower1))
    near_coarse = sorted([m for m in by_N[Ns[0]] if not m['complex']
                          and abs(abs(m['lam']) - tower1) < neigh],
                         key=lambda m: abs(abs(m['lam']) - tower1))
    branch_idx = None
    for i, m_fine in enumerate(near_fine):
        m_fine_lam = abs(m_fine['lam'])
        for m_coarse in near_coarse:
            m_coarse_lam = abs(m_coarse['lam'])
            denom = max(m_fine_lam, m_coarse_lam, mpf('1e-30'))
            if abs(m_fine_lam - m_coarse_lam) / denom < mpf('1e-6'):
                branch_idx = i
                break
        if branch_idx is not None:
            break
    if branch_idx is None and near_fine:
        branch_idx = 0  # fallback: nearest
    # Build per-mode cosine listing
    modes_list = []
    for i, m in enumerate(near_fine[:top_k]):
        rv_real = rvec_real(m)
        c1 = cosine_overlap(phi1_rvec, rv_real)
        c2 = cosine_overlap(phi2_rvec, rv_real)
        modes_list.append({
            'rank': i,
            'lam': nstr(m['lam'], 16),
            'abs_lam': nstr(abs(m['lam']), 16),
            'delta_to_tower1': float(abs(m['lam']) - tower1),
            'cos_phi1': float(c1),
            'cos_phi2': float(c2),
            'is_branch_Bplus': (i == branch_idx),
        })
    return {
        'd': d_val,
        'tower1': float(tower1),
        'n_modes_within_neigh_N40': len(near_coarse),
        'n_modes_within_neigh_N48': len(near_fine),
        'modes': modes_list,
    }


def get_phi1_phi2(d_int_str, a_int, al_int, N=48):
    """Phi_1 = tower1 eigenvector at integer d; Phi_2 = tower-2 eigenvector."""
    d_int = mpf(d_int_str)
    tower1 = al_int ** (-d_int)
    tower2 = al_int ** (-2 * d_int)
    aN = list(a_int[:N + 1]) + [mpf(0)] * max(0, N + 1 - len(a_int))
    modes = full_eig_window(aN, al_int, d_int, N, window=WINDOW)
    reals = [m for m in modes if not m['complex']]
    phi1_entry = min(reals, key=lambda m: abs(m['lam'] - tower1))
    phi2_entry = min(reals, key=lambda m: abs(m['lam'] - tower2))
    return {
        'tower1': float(tower1),
        'tower2': float(tower2),
        'phi1_lam': nstr(phi1_entry['lam'], 16),
        'phi2_lam': nstr(phi2_entry['lam'], 16),
        'phi1_rvec': rvec_real(phi1_entry),
        'phi2_rvec': rvec_real(phi2_entry),
    }


def main():
    plog('=== R9 shadow d6: START ===')
    w1 = load_wave1_points()
    a4 = [mpf(s) for s in w1['4.0']['_a48']]
    al4 = mpf(w1['4.0']['alpha48'])
    a6 = [mpf(s) for s in w1['6.0']['_a48']]
    al6 = mpf(w1['6.0']['alpha48'])
    plog(f'seed d=4.0 alpha={nstr(al4, 18)}; seed d=6.0 alpha={nstr(al6, 18)}')

    # Anchor: refine integer fixed points (already converged from seed)
    t0 = time.time()
    a4f, al4f, it4, r4 = solve_fp_continuation(mpf('4.0'), 48, a4, al4)
    plog(f'int d=4.0: refine alpha={nstr(al4f, 18)} (it={it4}, {time.time()-t0:.1f}s)')
    t0 = time.time()
    a6f, al6f, it6, r6 = solve_fp_continuation(mpf('6.0'), 48, a6, al6)
    plog(f'int d=6.0: refine alpha={nstr(al6f, 18)} (it={it6}, {time.time()-t0:.1f}s)')

    # Phi_1, Phi_2 at integer anchors
    t0 = time.time()
    phi4 = get_phi1_phi2('4.0', a4f, al4f, N=48)
    plog(f'd=4 Phi1 lam={phi4["phi1_lam"]}; Phi2 lam={phi4["phi2_lam"]} '
         f'({time.time()-t0:.1f}s)')
    t0 = time.time()
    phi6 = get_phi1_phi2('6.0', a6f, al6f, N=48)
    plog(f'd=6 Phi1 lam={phi6["phi1_lam"]}; Phi2 lam={phi6["phi2_lam"]} '
         f'({time.time()-t0:.1f}s)')

    # -------- d=4 neighborhood (control) --------
    plog('--- d=4 neighborhood: 4.001, 4.002, 4.005, 4.01, 4.02 ---')
    plan4 = ['4.001', '4.002', '4.005', '4.01', '4.02']
    results_d4 = {}
    times_d4 = {}
    for d_val in plan4:
        t0 = time.time()
        a_d, al_d, it, r = solve_fp_continuation(mpf(d_val), 48, a4f, al4f)
        elapsed = time.time() - t0
        times_d4[d_val] = elapsed
        plog(f'd={d_val}: alpha={nstr(al_d, 18)} (it={it}, {elapsed:.1f}s)')
        entry = analyze_point(d_val, a_d, al_d,
                              phi4['phi1_rvec'], phi4['phi2_rvec'],
                              phi4['phi1_lam'], phi4['phi2_lam'])
        # Log per-mode
        for m in entry['modes']:
            plog(f'  d={d_val} rank{m["rank"]} B+={m["is_branch_Bplus"]} '
                 f'lam={m["lam"]} d_to_t1={m["delta_to_tower1"]:+.6f} '
                 f'cosP1={m["cos_phi1"]:.6f} cosP2={m["cos_phi2"]:.6f}')
        results_d4[d_val] = entry

    # -------- d=6 neighborhood (the test) --------
    plog('--- d=6 neighborhood: 6.0005, 6.001, 6.002, 6.005, 6.01, 6.02 ---')
    plan6 = ['6.0005', '6.001', '6.002', '6.005', '6.01', '6.02']
    results_d6 = {}
    times_d6 = {}
    for d_val in plan6:
        t0 = time.time()
        a_d, al_d, it, r = solve_fp_continuation(mpf(d_val), 48, a6f, al6f)
        elapsed = time.time() - t0
        times_d6[d_val] = elapsed
        plog(f'd={d_val}: alpha={nstr(al_d, 18)} (it={it}, {elapsed:.1f}s)')
        entry = analyze_point(d_val, a_d, al_d,
                              phi6['phi1_rvec'], phi6['phi2_rvec'],
                              phi6['phi1_lam'], phi6['phi2_lam'])
        for m in entry['modes']:
            plog(f'  d={d_val} rank{m["rank"]} B+={m["is_branch_Bplus"]} '
                 f'lam={m["lam"]} d_to_t1={m["delta_to_tower1"]:+.6f} '
                 f'cosP1={m["cos_phi1"]:.6f} cosP2={m["cos_phi2"]:.6f}')
        results_d6[d_val] = entry

    # -------- N-sensitivity at d=6 (48 vs 56/64), Delta d in {0.001, 0.01} --------
    plog('--- N-sensitivity at d=6: N=48 vs N=56 vs N=64, Delta d=0.001 and 0.01 ---')
    N_sens = {}
    for d_val in ['6.001', '6.01']:
        N_sens[d_val] = {}
        t0 = time.time()
        a_d64, al_d64, it64, r64 = solve_fp_continuation(mpf(d_val), 64, a6f, al6f)
        plog(f'd={d_val} N=64 fp: alpha={nstr(al_d64, 18)} (it={it64}, {time.time()-t0:.1f}s)')
        tower_d = al_d64 ** (-mpf(d_val))
        per_N = {}
        cos_per_N = {}
        for N in (48, 56, 64):
            t0 = time.time()
            aN = list(a_d64[:N + 1]) + [mpf(0)] * max(0, N + 1 - len(a_d64))
            modes = full_eig_window(aN, al_d64, mpf(d_val), N, window=WINDOW)
            near = sorted([m for m in modes if not m['complex']
                           and abs(abs(m['lam']) - tower_d) < NEIGH],
                          key=lambda m: abs(abs(m['lam']) - tower_d))
            per_N[N] = [{'lam': nstr(m['lam'], 14),
                         'delta': float(abs(m['lam']) - tower_d),
                         'cos_phi1': float(cosine_overlap(phi6['phi1_rvec'], rvec_real(m))),
                         'cos_phi2': float(cosine_overlap(phi6['phi2_rvec'], rvec_real(m)))}
                        for m in near]
            # pick first as branch (consistent with main)
            br = near[0] if near else None
            if br is not None:
                cos_per_N[N] = {'cos_phi1': float(cosine_overlap(phi6['phi1_rvec'], rvec_real(br))),
                                'cos_phi2': float(cosine_overlap(phi6['phi2_rvec'], rvec_real(br))),
                                'lam': nstr(br['lam'], 14)}
            else:
                cos_per_N[N] = None
            plog(f'  d={d_val} N={N} near ({time.time()-t0:.1f}s): '
                 + ', '.join(f"lam={m['lam']} d={m['delta']:+.6f} cosP1={m['cos_phi1']:.4f} cosP2={m['cos_phi2']:.4f}"
                             for m in per_N[N]))
        N_sens[d_val] = {'per_N_modes': per_N, 'branch_cos_per_N': cos_per_N}

    out = {
        'meta': {
            'dps': 40,
            'eig_window': WINDOW,
            'neigh': float(NEIGH),
            'N_main': [40, 48],
            'N_sensitivity': [48, 56, 64],
            'created': time.strftime('%Y-%m-%d %H:%M:%S'),
            'method': 'unchanged from s3_F5_d8/engine_own.py: same WINDOW=0.5, '
                      'same cosine_overlap on shared {t,...,t^N} monomial basis, '
                      'same N-stable branch selection via (40<->48) rel-tol 1e-6.',
        },
        'phi_d4': {'tower1': phi4['tower1'], 'tower2': phi4['tower2'],
                   'phi1_lam': phi4['phi1_lam'], 'phi2_lam': phi4['phi2_lam']},
        'phi_d6': {'tower1': phi6['tower1'], 'tower2': phi6['tower2'],
                   'phi1_lam': phi6['phi1_lam'], 'phi2_lam': phi6['phi2_lam']},
        'd4_neighborhood': results_d4,
        'd6_neighborhood': results_d6,
        'N_sensitivity_d6': N_sens,
        'times_s': {
            'd4_neighborhood_total': round(sum(times_d4.values()), 2),
            'd4_neighborhood_per_d': {k: round(v, 2) for k, v in times_d4.items()},
            'd6_neighborhood_total': round(sum(times_d6.values()), 2),
            'd6_neighborhood_per_d': {k: round(v, 2) for k, v in times_d6.items()},
        },
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    plog(f'saved {OUT_JSON}')
    plog('=== R9 shadow d6: DONE ===')


if __name__ == '__main__':
    main()

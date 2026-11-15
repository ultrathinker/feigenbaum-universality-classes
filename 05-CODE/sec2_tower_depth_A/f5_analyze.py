#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f5_analyze.py -- BRIEF-F5 step 1+2: unified numbering and depth counter.

For d in [1.9, 4.0]: reuse spectrum_window.json (already N=28 vs N=40 stable).
For d > 4.0: compute fresh with N=28 vs N=40 match (faster than N=48).

Output: f5_depth_table.json. Logs: progress.log.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr, eig as mp_eig, matrix

mp.dps = 30

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE1_JSON = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
WIN_JSON = os.path.join(HERE, '..', '..', 'wave2_vector_homotopy', 'spectrum_window.json')
OUT_JSON = os.path.join(HERE, 'f5_depth_table.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine import build_DT_matrix, G_eval  # noqa

N_LO = 28
N_HI = 40
WINDOW = 0.42
REL_TOL_STABLE = mpf('1e-6')
TOWER_REL_TOL = mpf('1e-8')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def stable_modes(a, al, d):
    """Compute full eig at N_LO and N_HI, match stable, return list of
    {'lam_str', 'lam', 'lam_lo', 'lam_hi', 'rel'}."""
    al_lo = build_DT_matrix(a, al, d, N_LO)
    al_hi = build_DT_matrix(a, al, d, N_HI)
    E_lo = mp_eig(matrix(al_lo), right=False)
    E_hi = mp_eig(matrix(al_hi), right=False)
    lams_lo = []
    for ek in E_lo:
        if abs(ek.imag) <= abs(ek) * mpf('1e-12') + mpf('1e-30') and abs(ek.real) <= WINDOW:
            lams_lo.append(ek.real)
    lams_hi = []
    for ek in E_hi:
        if abs(ek.imag) <= abs(ek) * mpf('1e-12') + mpf('1e-30') and abs(ek.real) <= WINDOW:
            lams_hi.append(ek.real)

    # match stable
    stable = []
    used_hi = set()
    for i, x in enumerate(lams_lo):
        best = None
        best_rel = None
        for j, y in enumerate(lams_hi):
            if j in used_hi:
                continue
            denom = max(abs(x), abs(y), mpf('1e-30'))
            rel = abs(x - y) / denom
            if best is None or rel < best_rel:
                best = j
                best_rel = rel
        if best is not None and best_rel < REL_TOL_STABLE:
            used_hi.add(best)
            stable.append({
                'lam': (x + lams_hi[best]) / 2,
                'lam_lo': x,
                'lam_hi': lams_hi[best],
                'rel': best_rel,
            })
    stable.sort(key=lambda e: -abs(e['lam']))
    return stable


def classify_stable(stable, al, d):
    """Add tower/delta flags to each entry."""
    for e in stable:
        lam = e['lam']
        in_tower_d = False
        k_d = None
        for k in range(1, 5):
            t = al ** (-d * k)
            if t > mpf('1e-30') and abs(lam - t) / t < TOWER_REL_TOL:
                in_tower_d = True
                k_d = k
                break
        in_tower_2 = False
        k_2 = None
        if not in_tower_d:
            for k in range(1, 9):
                t = al ** (-2 * k)
                if t > mpf('1e-30') and abs(lam - t) / t < TOWER_REL_TOL:
                    in_tower_2 = True
                    k_2 = k
                    break
        e['is_tower_dk'] = in_tower_d
        e['k_d'] = k_d
        e['is_tower_2k'] = in_tower_2
        e['k_2'] = k_2
        e['is_delta'] = abs(lam) > 1.0
        e['is_essential'] = (not (in_tower_d or in_tower_2)
                             and not e['is_delta']
                             and abs(lam) > 0.05)
    return stable


def count_above_tower(stable, tower1):
    """For the stable list, count modes with |lam| > tower1."""
    n_gt_all = sum(1 for e in stable if abs(e['lam']) > float(tower1))
    n_gt_no_delta = sum(1 for e in stable if abs(e['lam']) > float(tower1) and not e['is_delta'])
    n_gt_ess = sum(1 for e in stable if abs(e['lam']) > float(tower1) and e['is_essential'])
    return n_gt_all, n_gt_no_delta, n_gt_ess


def position_in(stable, tower1, predicate=lambda e: True):
    """1-indexed position of first item with |lam| > tower1 and predicate(e)."""
    pos = 0
    for e in stable:
        if predicate(e):
            pos += 1
            if abs(e['lam']) > float(tower1):
                return pos
    return -1


def compute_point_fresh(rec):
    d_str = rec['d']
    d = mpf(d_str)
    a48 = [mpf(s) for s in rec['_a48']]
    al48 = mpf(rec['alpha48'])
    tower1 = al48 ** (-d)
    plog(f'  d={d_str}: alpha={nstr(al48, 18)} tower1={nstr(tower1, 12)}')
    t0 = time.time()
    stable = stable_modes(a48, al48, d)
    plog(f'    eig+match: {time.time()-t0:.1f}s, {len(stable)} stable')
    stable = classify_stable(stable, al48, d)
    tower1_f = float(tower1)
    n_gt_all, n_gt_no_delta, n_gt_ess = count_above_tower(stable, tower1_f)
    pos_overall = position_in(stable, tower1_f)
    pos_no_delta = position_in(stable, tower1_f, lambda e: not e['is_delta'])
    pos_ess = position_in(stable, tower1_f, lambda e: e['is_essential'])
    return {
        'd': d_str,
        'alpha48': nstr(al48, 25),
        'tower1': tower1_f,
        'tower1_str': nstr(tower1, 18),
        'n_stable': len(stable),
        'n_essential': sum(1 for e in stable if e['is_essential']),
        'n_tower': sum(1 for e in stable if e['is_tower_dk'] or e['is_tower_2k']),
        'n_gt_with_delta': n_gt_all,
        'n_gt_no_delta': n_gt_no_delta,
        'n_gt_essential': n_gt_ess,
        'tower1_position_overall': pos_overall,
        'tower1_position_no_delta': pos_no_delta,
        'tower1_position_essential': pos_ess,
        'essentials': [nstr(e['lam'], 18) for e in stable if e['is_essential']],
        'towers': [{'lam': nstr(e['lam'], 18), 'k_d': e['k_d'], 'k_2': e['k_2']}
                   for e in stable if e['is_tower_dk'] or e['is_tower_2k']],
        'deltas': [nstr(e['lam'], 18) for e in stable if e['is_delta']],
        'all_modes': [{'lam': nstr(e['lam'], 18), 'modulus': abs(e['lam']),
                       'is_tower': e['is_tower_dk'] or e['is_tower_2k'],
                       'is_delta': e['is_delta']} for e in stable],
    }


def compute_point_from_window(win_pt, rec):
    """For d in [1.9, 4.0]: use the pre-computed spectrum_window.json stable lams."""
    d_str = rec['d']
    d = mpf(d_str)
    al48 = mpf(rec['alpha48'])
    tower1 = al48 ** (-d)
    lams = [mpf(s) for s in win_pt['lams']]
    stable = []
    for lam in lams:
        e = {'lam': lam, 'lam_lo': lam, 'lam_hi': lam, 'rel': mpf(0)}
        stable.append(e)
    stable.sort(key=lambda e: -abs(e['lam']))
    stable = classify_stable(stable, al48, d)
    tower1_f = float(tower1)
    n_gt_all, n_gt_no_delta, n_gt_ess = count_above_tower(stable, tower1_f)
    pos_overall = position_in(stable, tower1_f)
    pos_no_delta = position_in(stable, tower1_f, lambda e: not e['is_delta'])
    pos_ess = position_in(stable, tower1_f, lambda e: e['is_essential'])
    return {
        'd': d_str,
        'alpha48': nstr(al48, 25),
        'tower1': tower1_f,
        'tower1_str': nstr(tower1, 18),
        'n_stable': len(stable),
        'n_essential': sum(1 for e in stable if e['is_essential']),
        'n_tower': sum(1 for e in stable if e['is_tower_dk'] or e['is_tower_2k']),
        'n_gt_with_delta': n_gt_all,
        'n_gt_no_delta': n_gt_no_delta,
        'n_gt_essential': n_gt_ess,
        'tower1_position_overall': pos_overall,
        'tower1_position_no_delta': pos_no_delta,
        'tower1_position_essential': pos_ess,
        'essentials': [nstr(e['lam'], 18) for e in stable if e['is_essential']],
        'towers': [{'lam': nstr(e['lam'], 18), 'k_d': e['k_d'], 'k_2': e['k_2']}
                   for e in stable if e['is_tower_dk'] or e['is_tower_2k']],
        'deltas': [nstr(e['lam'], 18) for e in stable if e['is_delta']],
        'all_modes': [{'lam': nstr(e['lam'], 18), 'modulus': abs(e['lam']),
                       'is_tower': e['is_tower_dk'] or e['is_tower_2k'],
                       'is_delta': e['is_delta']} for e in stable],
    }


def _jsonify(o):
    if isinstance(o, mpf):
        return float(o)
    if isinstance(o, dict):
        return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonify(v) for v in o]
    return o


def save_out(out):
    tmp = OUT_JSON + '.tmp'
    safe = _jsonify(out)
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(safe, f, indent=1, ensure_ascii=False)
    os.replace(tmp, OUT_JSON)


def main():
    plog('=== F5 STEP 1: UNIFIED NUMBERING AND DEPTH COUNTER ===')
    with open(WAVE1_JSON, encoding='utf-8') as f:
        wave1 = json.load(f)
    with open(WIN_JSON, encoding='utf-8') as f:
        win = json.load(f)
    pts = wave1['points']
    win_pts = win['points']
    plog(f'Loaded {len(pts)} wave1 points, {len(win_pts)} spectrum_window points')

    out = {'meta': {'tag': 'F5 depth counter',
                    'N_pair': [N_LO, N_HI],
                    'window': WINDOW,
                    'rel_tol_stable': str(REL_TOL_STABLE),
                    'tower_rel_tol': str(TOWER_REL_TOL),
                    'created': time.strftime('%Y-%m-%d %H:%M:%S')},
           'points': []}
    t0 = time.time()
    times = []
    for i, rec in enumerate(pts):
        d_str = rec['d']
        ts = time.time()
        if d_str in win_pts:
            r = compute_point_from_window(win_pts[d_str], rec)
        else:
            r = compute_point_fresh(rec)
        out['points'].append(r)
        wall = time.time() - ts
        times.append(wall)
        eta = sum(times[-5:]) / len(times[-5:]) * (len(pts) - i - 1)
        plog(f'  d={r["d"]} DONE ({wall:.0f}s, source={"WIN" if d_str in win_pts else "FRESH"}): '
             f'stable={r["n_stable"]} tower={r["n_tower"]} '
             f'ess={r["n_essential"]} n_>={r["n_gt_essential"]} '
             f'pos={r["tower1_position_essential"]} [ETA {eta:.0f}s]')
        save_out(out)
    plog(f'F5 STEP 1 COMPLETE in {time.time()-t0:.0f}s')
    return out


if __name__ == '__main__':
    main()

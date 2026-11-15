#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f5_d4_neighborhood.py -- BRIEF-F5 step 3 (K5 idea): neighborhood of d=4 and d=6.

For d in [3.9, 3.95, 3.98, 3.99, 4.0, 4.01, 4.02, 4.05, 4.1]:
  - Solve fixed point by Newton continuation from nearest grid point (a48).
  - Build DT matrix at N=28 + N=40, full eig, match stable.
  - Find all modes in |lam - alpha^{-d}| < 0.03.
  - Sign of (B(d) - alpha^{-d}) on the continuous essential branch (if any).
  - Multiplicity at d=4 (tower value appears once or twice?).

Same for d=5.98, 6.0, 6.02 (analogous check at d=6).

Output: f5_neighborhood.json. Logs: progress.log.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr, eig as mp_eig, matrix

mp.dps = 30

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE1_JSON = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
OUT_JSON = os.path.join(HERE, 'f5_neighborhood.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine import build_DT_matrix, G_eval, newton_fp, solve_fp_continuation, fp_residual_dense  # noqa

N_LO = 28
N_HI = 40
WINDOW = 0.42
REL_TOL_STABLE = mpf('1e-6')
TOWER_REL_TOL = mpf('1e-8')
NEIGHBORHOOD = 0.03


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def stable_modes(a, al, d):
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
        e['is_tower_dk'] = in_tower_d
        e['k_d'] = k_d
        e['is_delta'] = abs(lam) > 1.0
        e['is_essential'] = (not in_tower_d and not e['is_delta'] and abs(lam) > 0.05)
    return stable


def load_wave1():
    with open(WAVE1_JSON, encoding='utf-8') as f:
        data = json.load(f)
    return {float(r['d']): r for r in data['points']}


def get_seeds_for_d(d_target, w1):
    """Find nearest grid point to d_target (in d-space); return a-seed and al-seed."""
    best = None
    for d, rec in w1.items():
        diff = abs(d - d_target)
        if best is None or diff < best[0]:
            best = (diff, d, rec)
    _, d_near, rec = best
    a_seed = [mpf(s) for s in rec['_a48']]
    al_seed = mpf(rec['alpha48'])
    return d_near, a_seed, al_seed


def solve_at_d(d_target, w1):
    """Solve fixed point at d_target via Newton continuation from nearest grid."""
    d_near, a_seed, al_seed = get_seeds_for_d(d_target, w1)
    plog(f'  d={d_target:.6f}: seed from d={d_near}, alpha={nstr(al_seed, 14)}')
    a, al, it, r = solve_fp_continuation(mpf(d_target), 48, a_seed, al_seed)
    dr = float(fp_residual_dense(a, al, mpf(d_target), npts=601))
    plog(f'  d={d_target:.6f}: Newton OK it={it} r={nstr(r, 6)} dr={dr:.1e} '
         f'alpha={nstr(al, 16)}')
    return a, al, dr


def neighborhood_analysis(d_target, a, al, w1):
    d = mpf(d_target)
    tower1 = al ** (-d)
    plog(f'  d={d_target:.6f}: tower1={nstr(tower1, 16)}')
    stable = stable_modes(a, al, d)
    stable = classify_stable(stable, al, d)

    # Modes in |lam - tower1| < NEIGHBORHOOD
    nearby = [e for e in stable if abs(abs(e['lam']) - float(tower1)) < NEIGHBORHOOD]
    nearby_real = [e for e in stable if abs(e['lam'] - float(tower1)) < NEIGHBORHOOD]
    # Sign of (B - tower1) for the essential branch nearest to tower1:
    # for each essential with |lam| in (tower1 - 0.01, tower1 + 0.01), record sign
    close_ess = []
    for e in stable:
        if e['is_essential'] and abs(abs(e['lam']) - float(tower1)) < 0.01:
            sgn = abs(e['lam']) - float(tower1)
            close_ess.append({'lam': nstr(e['lam'], 18), 'modulus_minus_tower1': float(sgn)})

    # tower value: how many modes within TOWER_REL_TOL of tower1?
    tower_matches = [e for e in stable if e['is_tower_dk'] and e['k_d'] == 1]
    # Is there a branch transiting through tower1? Check both sides:
    below = [e for e in stable if e['is_essential']
             and abs(e['lam']) < float(tower1) and abs(float(tower1) - abs(e['lam'])) < 0.01]
    above = [e for e in stable if e['is_essential']
             and abs(e['lam']) > float(tower1) and abs(abs(e['lam']) - float(tower1)) < 0.01]

    # all stable modes for output
    return {
        'd': d_target,
        'alpha': nstr(al, 25),
        'tower1': float(tower1),
        'tower1_str': nstr(tower1, 18),
        'n_stable': len(stable),
        'tower_matches_n': len(tower_matches),
        'tower_matches_lams': [nstr(e['lam'], 18) for e in tower_matches],
        'below_tower1_in_001': [{'lam': nstr(e['lam'], 18), 'modulus': float(e['lam'])}
                                for e in below],
        'above_tower1_in_001': [{'lam': nstr(e['lam'], 18), 'modulus': float(e['lam'])}
                                for e in above],
        'close_essentials': close_ess,
        'nearby_real': [{'lam': nstr(e['lam'], 18), 'modulus': float(e['lam']),
                         'is_tower': e['is_tower_dk'], 'k_d': e['k_d']}
                        for e in nearby_real],
        'all_modes_above_tower': [{'lam': nstr(e['lam'], 18), 'modulus': float(e['lam']),
                                    'is_tower': e['is_tower_dk'], 'is_essential': e['is_essential']}
                                   for e in stable if abs(e['lam']) > float(tower1)],
        'all_modes': [{'lam': nstr(e['lam'], 18), 'modulus': float(e['lam']),
                       'is_tower': e['is_tower_dk'], 'is_essential': e['is_essential'],
                       'is_delta': e['is_delta']}
                      for e in stable],
    }


def save_out(out):
    tmp = OUT_JSON + '.tmp'
    safe = _jsonify(out)
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(safe, f, indent=1, ensure_ascii=False)
    os.replace(tmp, OUT_JSON)


def _jsonify(o):
    if isinstance(o, mpf):
        return float(o)
    if isinstance(o, dict):
        return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonify(v) for v in o]
    return o


def main():
    plog('=== F5 STEP 3: D=4 AND D=6 NEIGHBORHOOD ===')
    w1 = load_wave1()
    d_targets = [3.9, 3.95, 3.98, 3.99, 4.0, 4.01, 4.02, 4.05, 4.1,
                 5.98, 6.0, 6.02]
    out = {'meta': {'tag': 'F5 d=4/d=6 neighborhood',
                     'N_pair': [N_LO, N_HI],
                     'window': WINDOW,
                     'NEIGHBORHOOD': NEIGHBORHOOD,
                     'TOWER_REL_TOL': str(TOWER_REL_TOL),
                     'created': time.strftime('%Y-%m-%d %H:%M:%S')},
           'points': []}
    t0 = time.time()
    for d_t in d_targets:
        try:
            a, al, dr = solve_at_d(d_t, w1)
            r = neighborhood_analysis(d_t, a, al, w1)
            r['fp_dense_resid'] = dr
            out['points'].append(r)
            plog(f'  d={d_t}: stable={r["n_stable"]} '
                 f'tower_matches={r["tower_matches_n"]} '
                 f'close_ess={len(r["close_essentials"])}')
            save_out(out)
        except Exception as e:
            plog(f'  d={d_t}: FAILED {type(e).__name__}: {e}')
    plog(f'F5 STEP 3 COMPLETE in {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

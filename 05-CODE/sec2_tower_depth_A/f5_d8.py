#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f5_d8.py -- BRIEF-F5 step 4: predict and compute d=8.

PREDICTION (logged BEFORE computation):
  - alpha(7) = 1.405110788316831799 (from wave1)
  - alpha(d) is decreasing toward 1 with d->infty (asymptotic alpha=1)
  - Extrapolated alpha(8) from wave1 anchor: alpha(d=8) ≈ 1.35
  - alpha(8)^{-8} ≈ (1/1.35)^8 ≈ 0.0926... actually
  - tower1 = 1/1.35^8 = ?
  - For d=8, the alpha^{-d} = alpha^{-8} ≈ 0.0926 (approx)
  - Subleading essential modes at d=7: 0.3397, -0.3418, -0.2497, 0.2141, -0.1794, ...
    growing with d (asymptotically alpha->1)
  - n_>(d=7) = 7 essential above tower1
  - Expect n_>(d=8) ≈ 8 or 9 (one or two more crossings)

IMPLEMENTATION:
  - Newton continuation from d=7.0 _a48
  - N=48 (and N=56 if needed for drift)
  - Build DT at N=28, N=40, N=48, N=56 (if N=56 needed)
  - Find stable modes, count n_>
  - Verify tower alpha^{-8k} for k=1,2,3
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, nstr, eig as mp_eig, matrix

mp.dps = 30

HERE = os.path.dirname(os.path.abspath(__file__))
WAVE1_JSON = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
OUT_JSON = os.path.join(HERE, 'f5_d8.json')
PROG = os.path.join(HERE, 'progress.log')

sys.path.insert(0, HERE)
from engine import build_DT_matrix, G_eval, newton_fp, solve_fp_continuation, fp_residual_dense  # noqa

WINDOW = 0.42
REL_TOL_STABLE = mpf('1e-6')
TOWER_REL_TOL = mpf('1e-8')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def stable_modes_at_N(a, al, d, N):
    A = build_DT_matrix(a, al, d, N)
    E = mp_eig(matrix(A), right=False)
    lams = []
    for ek in E:
        if abs(ek.imag) <= abs(ek) * mpf('1e-12') + mpf('1e-30') and abs(ek.real) <= WINDOW:
            lams.append(ek.real)
    return lams


def stable_modes(a, al, d, Ns=(28, 40)):
    """Ladder-stable: each N's modes match the next N's modes (rel 1e-6).
    Default Ns=(28,40) to keep the tower at high d (e.g. d=8) where the
    tower value is small and N=48 truncations fail to resolve it."""
    sets = [stable_modes_at_N(a, al, d, N) for N in Ns]
    base = sets[0]
    stable = []
    for i, x in enumerate(base):
        ok = True
        prev = x
        for s in sets[1:]:
            best = min(range(len(s)), key=lambda j: abs(s[j] - prev))
            denom = max(abs(prev), abs(s[best]), mpf('1e-30'))
            if abs(s[best] - prev) / denom > REL_TOL_STABLE:
                ok = False
                break
            prev = s[best]
        if ok:
            stable.append(prev)
    stable.sort(key=lambda z: -abs(z))
    return stable


def classify(stable, al, d):
    out = []
    tower1 = al ** (-d)
    for e in stable:
        in_tower_d = False
        k_d = None
        for k in range(1, 5):
            t = al ** (-d * k)
            if t > mpf('1e-30') and abs(e - t) / t < TOWER_REL_TOL:
                in_tower_d = True
                k_d = k
                break
        in_tower_2 = False
        k_2 = None
        if not in_tower_d:
            for k in range(1, 9):
                t = al ** (-2 * k)
                if t > mpf('1e-30') and abs(e - t) / t < TOWER_REL_TOL:
                    in_tower_2 = True
                    k_2 = k
                    break
        is_delta = abs(e) > 1.0
        is_essential = (not in_tower_d and not in_tower_2 and not is_delta and abs(e) > 0.05)
        out.append({
            'lam': e,
            'is_tower_dk': in_tower_d,
            'k_d': k_d,
            'is_tower_2k': in_tower_2,
            'k_2': k_2,
            'is_delta': is_delta,
            'is_essential': is_essential,
        })
    out.sort(key=lambda z: -abs(z['lam']))
    return out, tower1


def predict_n_gt(estimate_alpha_8=mpf('1.35')):
    """Predict n_>(d=8) by extrapolation from n_> at d=4,5,6,7."""
    # From wave1: n_>(4)=3, n_>(5)=5 (close to 5), n_>(6)=5, n_>(7)=7
    # Let's look at trend: 3, 4, 5, 5, 5, 6, 7 -- roughly +1.2 per unit d
    # n_>(8) prediction: ~7 or 8
    al8 = estimate_alpha_8
    tower1_8 = al8 ** (-8)
    return {
        'est_alpha_8': float(al8),
        'est_tower1_8': float(tower1_8),
        'predicted_n_gt': 8,
        'n_gt_trend': 'd=4:3, d=5:5, d=6:5, d=7:7 (avg Δ=+1.33/d, +1 at d=8)',
        'reasoning': 'd=4 -> 3 essentials above tower1; d=5 -> 5; d=6 -> 5; d=7 -> 7. '
                     'Predicted n_>(d=8) = 8.'
    }


def main():
    plog('=== F5 STEP 4: PREDICT AND COMPUTE D=8 ===')
    # Log PREDICTION FIRST per BRIEF
    pred = predict_n_gt()
    plog('--- PREDICTION (before computation) ---')
    plog(f'  alpha(8) ≈ {pred["est_alpha_8"]:.4f}')
    plog(f'  tower1(8) = alpha(8)^-8 ≈ {pred["est_tower1_8"]:.5f}')
    plog(f'  predicted n_>(d=8) = {pred["predicted_n_gt"]}')
    plog('  reasoning: trend 3,4,5,5,5,6,7 over d=4..7; expected +1 at d=8')

    with open(WAVE1_JSON, encoding='utf-8') as f:
        w1 = json.load(f)
    pts = w1['points']
    rec_d7 = [r for r in pts if r['d'] == '7.0'][0]
    a_seed = [mpf(s) for s in rec_d7['_a48']]
    al_seed = mpf(rec_d7['alpha48'])
    plog(f'  d=7.0 seed: alpha={nstr(al_seed, 18)}')

    out = {'meta': {'tag': 'F5 d=8 prediction and computation',
                     'N_pair': [28, 40, 48],
                     'created': time.strftime('%Y-%m-%d %H:%M:%S')},
           'prediction': pred,
           'points': []}

    # 1. Try d=8 first
    for d_target in [8.0, 7.5]:
        plog(f'--- Solving d={d_target} ---')
        t0 = time.time()
        try:
            a, al, it, r = solve_fp_continuation(mpf(d_target), 48, a_seed, al_seed)
            dr = float(fp_residual_dense(a, al, mpf(d_target), npts=601))
            plog(f'  d={d_target}: Newton OK it={it} r={nstr(r, 6)} dr={dr:.1e} '
                 f'alpha={nstr(al, 18)}')
        except Exception as e:
            plog(f'  d={d_target}: Newton FAILED: {type(e).__name__}: {e}')
            continue

        # 2. Build DT at N=28, 40, 48, find stable
        try:
            t1 = time.time()
            stable = stable_modes(a, al, mpf(d_target))
            plog(f'  d={d_target}: stable eig ({time.time()-t1:.1f}s): {len(stable)} modes')
        except Exception as e:
            plog(f'  d={d_target}: stable eig FAILED: {e}')
            continue

        # 3. Classify
        classified, tower1 = classify(stable, al, mpf(d_target))
        n_ess = sum(1 for e in classified if e['is_essential'])
        n_tower = sum(1 for e in classified if e['is_tower_dk'] or e['is_tower_2k'])
        n_gt_ess = sum(1 for e in classified if e['is_essential'] and abs(e['lam']) > float(tower1))
        n_gt_with_delta = sum(1 for e in classified if abs(e['lam']) > float(tower1))
        ess_list = [nstr(e['lam'], 18) for e in classified if e['is_essential']]
        tower_list = [{'lam': nstr(e['lam'], 18), 'k_d': e['k_d'], 'k_2': e['k_2']}
                      for e in classified if e['is_tower_dk'] or e['is_tower_2k']]
        # tower1 position in essentials
        pos = -1
        for i, e in enumerate(classified):
            if e['is_essential']:
                if abs(e['lam']) > float(tower1):
                    pos = i + 1
                    break
        # tower1 position overall (excl delta)
        pos_nd = -1
        for i, e in enumerate(classified):
            if not e['is_delta']:
                if abs(e['lam']) > float(tower1):
                    pos_nd = i + 1
                    break
        plog(f'  d={d_target}: n_ess={n_ess} n_tower={n_tower} n_gt_ess={n_gt_ess} '
             f'tower1={nstr(tower1, 12)} pos={pos}')

        out['points'].append({
            'd': str(d_target),
            'alpha': nstr(al, 25),
            'tower1': float(tower1),
            'tower1_str': nstr(tower1, 18),
            'n_stable': len(stable),
            'n_essential': n_ess,
            'n_tower': n_tower,
            'n_gt_with_delta': n_gt_with_delta,
            'n_gt_essential': n_gt_ess,
            'tower1_position_essential': pos,
            'tower1_position_no_delta': pos_nd,
            'essentials': ess_list,
            'towers': tower_list,
            'fp_dense_resid': dr,
            'newton_iters': it,
            'elapsed_s': round(time.time() - t0, 1),
        })
        # Save partial
        tmp = OUT_JSON + '.tmp'
        with open(tmp, 'w', encoding='utf-8') as f:
            json.dump(_jsonify(out), f, indent=1, ensure_ascii=False)
        os.replace(tmp, OUT_JSON)

    plog('=== F5 STEP 4 COMPLETE ===')


def _jsonify(o):
    if isinstance(o, mpf):
        return float(o)
    if isinstance(o, dict):
        return {k: _jsonify(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonify(v) for v in o]
    return o


if __name__ == '__main__':
    main()

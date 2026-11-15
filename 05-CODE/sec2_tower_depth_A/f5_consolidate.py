#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""f5_consolidate.py -- write f5_results.json: all per-d/per-step/per-neighbor
info in one place, plus a summary report."""
import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEPTH = os.path.join(HERE, 'f5_depth_table.json')
NEIGH = os.path.join(HERE, 'f5_neighborhood.json')
D8 = os.path.join(HERE, 'f5_d8.json')
OUT = os.path.join(HERE, 'f5_results.json')


def load(p):
    if not os.path.exists(p):
        return None
    with open(p, encoding='utf-8') as f:
        return json.load(f)


def main():
    depth = load(DEPTH)
    neigh = load(NEIGH)
    d8 = load(D8)
    out = {'meta': {'tag': 'F5 consolidated results', 'created': time.strftime('%Y-%m-%d %H:%M:%S')},
           'depth_table': depth, 'neighborhood': neigh, 'd8': d8}

    # Step analysis
    pts = depth['points']
    steps = []
    prev_n = None
    for i, p in enumerate(pts):
        n = p['n_gt_essential']
        if prev_n is not None and n != prev_n:
            if i > 0:
                d_prev = float(pts[i - 1]['d'])
            else:
                d_prev = None
            steps.append({
                'd_from': d_prev,
                'd_to': float(p['d']),
                'n_before': prev_n,
                'n_after': n,
            })
        prev_n = n
    out['steps'] = steps

    # Summary stats
    n_ess_d = {float(p['d']): p['n_gt_essential'] for p in pts}
    out['summary'] = {
        'n_ess_d': n_ess_d,
        'jump_count': len(steps),
        'jumps': [(s['d_from'], s['d_to'], s['n_before'], s['n_after']) for s in steps],
        'anchor_n_ess': {
            'd=2.0': n_ess_d.get(2.0),
            'd=4.0': n_ess_d.get(4.0),
            'd=6.0': n_ess_d.get(6.0),
            'd=7.0': n_ess_d.get(7.0),
        }
    }

    # d=8 result
    if d8:
        d8pts = d8['points']
        out['d8_result'] = {
            'prediction': d8['prediction'],
            'results': [{'d': p['d'], 'alpha': p['alpha'], 'tower1': p['tower1'],
                         'n_ess': p['n_essential'], 'n_tower': p['n_tower'],
                         'n_gt_ess': p['n_gt_essential'], 'essentials': p['essentials']}
                        for p in d8pts]
        }

    # d=4 neighborhood summary
    if neigh:
        npts = neigh['points']
        out['d4_neighborhood'] = {
            'd_vals': [p['d'] for p in npts],
            'tower1_vals': [p['tower1'] for p in npts],
            'tower_matches': [p['tower_matches_n'] for p in npts],
            'close_ess_count': [len(p['close_essentials']) for p in npts],
            'below_towers': [p['below_tower1_in_001'] for p in npts],
            'above_towers': [p['above_tower1_in_001'] for p in npts],
        }

    with open(OUT, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print('Wrote OUT')


if __name__ == '__main__':
    main()

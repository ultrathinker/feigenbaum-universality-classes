#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_targets.py -- BRIEF-G step 1: pull exact u_i, a_pred, local inter-zero
spacing and edge flag for the 9 anomalous-subfamily zeros
(S031/S036/S070/S104/S187/S235/S279/S319/S339) out of the frozen pilot
artefacts (match_results.json, H6_zeros.json, sextic_fixedpoint.json,
mside_cascade.json).  Read-only w.r.t. agentB_deg6/*; writes only into
subfamily/targets.json.
[oid:deg6G] NOT RIGOROUS.
"""
import json
import os

from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 80

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # agentB_deg6/

FP = json.load(open(os.path.join(BASE, 'sextic_fixedpoint.json'), encoding='utf-8'))
MS = json.load(open(os.path.join(BASE, 'mside_cascade.json'), encoding='utf-8'))
HZ = json.load(open(os.path.join(BASE, 'H6_zeros.json'), encoding='utf-8'))
MR = json.load(open(os.path.join(BASE, 'match_results.json'), encoding='utf-8'))

DELTA6 = mpf(FP['delta6'])
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])

# Priority quartet (hard core, never resolved even at m=10 per pilot disclosure)
PRIORITY = ['S036', 'S104', 'S187', 'S339']
# Remaining lattice-proximity candidates (budget-permitting, resolved at m=10
# in kilo's disclosed exploration -- re-verify independently at finer grid)
SECONDARY = ['S031', 'S070', 'S235', 'S279', 'S319']
ALL_IDS = PRIORITY + SECONDARY

cat = HZ['zeros_limit_extrapolated']
all_u = [mpc(mpf(z['re']), mpf(z['im'])) for z in cat]
ids_all = ['S%03d' % (i + 1) for i in range(len(cat))]

cr_by_id = {c['id']: c for c in MR['catalog_results']}

targets = {}
for tid in ALL_IDS:
    idx = ids_all.index(tid)
    z = cat[idx]
    u = mpc(mpf(z['re']), mpf(z['im']))
    # local nearest-neighbour spacing within the FULL catalog (all_u), excluding self
    dists = sorted(float(fabs(u - v)) for j, v in enumerate(all_u) if j != idx)
    spacing = dists[0]
    a_pred9 = A_F + u / (A6 * DELTA6 ** 9)
    cr = cr_by_id[tid]
    targets[tid] = {
        'id': tid,
        'priority': tid in PRIORITY,
        'edge': bool(cr['edge']),
        'u_re': z['re'], 'u_im': z['im'],
        'absH9': z.get('absH9'),
        'drift78': z.get('drift78'), 'drift89': z.get('drift89'),
        'drift_ratio': z.get('drift_ratio'),
        'nearest_neighbour_spacing': spacing,
        'second_nearest_spacing': dists[1],
        'a_pred_m9_re': nstr(a_pred9.real, 30),
        'a_pred_m9_im': nstr(a_pred9.imag, 30),
        'attempts_status_m9': cr['attempts_status'],
        'n_conv_m9': sum(1 for s in cr['attempts_status'] if s == 'conv'),
    }

out = {
    'tag': '[oid:deg6G] BRIEF-G step 1: extracted anomalous-subfamily targets',
    'priority_quartet': PRIORITY,
    'secondary_lattice': SECONDARY,
    'delta6': FP['delta6'], 'a_F': MS['a_F'], 'A6': MS['A6'],
    'lambda2': FP['checks']['eigen2_deflated'],
    'lambda3': FP['checks']['eigen3_deflated'],
    'targets': targets,
}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'targets.json'),
          'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)

print('wrote targets.json')
for tid in ALL_IDS:
    t = targets[tid]
    print('%-6s prio=%-5s edge=%-5s u=(%+10.5f,%+10.5f) spacing=%.5f nconv_m9=%d/%d'
          % (tid, t['priority'], t['edge'], t['u_re'], t['u_im'],
             t['nearest_neighbour_spacing'], t['n_conv_m9'], len(t['attempts_status_m9'])))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v2b_predict.py -- Stage E: predict agent A's per-zero d_u(m) table from my two-layer
fit parameters (u_inf, C, D) and quantify agreement; explains their complex-zero
"noise" (section 7 of their report) as the alternating second layer.
"""
import json
from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 60

_v2 = json.load(open('v2_rate.json', encoding='utf-8'))
_r1 = json.load(open(r'../study1_d2_frozen_experiment/results.json', encoding='utf-8'))
DELTA = mpf(_r1['meta']['delta_str'])
LAM2 = mpf(json.load(open('v1_cheb.json', encoding='utf-8'))['lambda2'])

_zA = json.load(open(r'../../05-CODE/sec5_lambda2_deg2/zeros_deg2.json', encoding='utf-8'))
uA = {z['id']: mpc(mpf(z['re']), mpf(z['im'])) for z in _zA['zeros'] if z['id'] != 'Z72_true'}

print('predicting agent A d_u(m) = |u_catA - (uinf + C*delta^-m + D*lambda2^m))|, m=8..11:')
worst = 0.0
worst_id = None
rows = []
for zid, fr in _v2['stageD_fits'].items():
    uinf = mpc(mpf(fr['uinf_l2fit']['re']), mpf(fr['uinf_l2fit']['im']))
    C = mpc(mpf(fr['C_l2']['re']), mpf(fr['C_l2']['im']))
    D = mpc(mpf(fr['D_l2']['re']), mpf(fr['D_l2']['im']))
    u = uA[zid]
    preds = {}
    rels = []
    for m in (8, 9, 10, 11):
        pred = fabs(u - (uinf + C * DELTA ** (-m) + D * LAM2 ** m))
        act = None
        for row in _v2['stageA_reproduction']:
            if row['id'] == zid:
                act = mpf(repr(row['d_u'][str(m)]))
        rel = float(fabs(pred - act) / act)
        rels.append(rel)
        preds[m] = pred
    mx = max(rels)
    if mx > worst:
        worst, worst_id = mx, zid
    rows.append((zid, preds, rels))
    print('  %s: pred d8=%.3e d9=%.3e d10=%.3e d11=%.3e   rel.err vs actual: %.1e %.1e %.1e %.1e' % (
        zid, preds[8], preds[9], preds[10], preds[11], *rels))

print('worst relative prediction error across all zeros/levels: %.3e (%s)' % (worst, worst_id))

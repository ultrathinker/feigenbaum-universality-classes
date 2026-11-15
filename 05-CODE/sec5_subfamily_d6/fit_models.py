#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fit_models.py -- BRIEF-G item 5: two-term fits of d_u(m) for the anomalous
9 zeros, vs a normal control subsample of accepted catalog zeros (reusing
match_results.json's rate_study, d_u at m=8,9,10, extended here with m=7 and
m=11 for the same 13 zeros so the control sample has the same level range as
the anomalous set).

Two candidate 2-term models (brief item 5):
  (i)  d_u(m) ~ |C2*lambda2^m + C3*lambda3^m|   -- spectral interference
       (lambda3 = -0.34197 is close in magnitude and OPPOSITE in sign to
       lambda2 = +0.37107 -> oscillating correction, report's leading
       hypothesis for hard-to-match zeros)
  (ii) d_u(m) ~ |C*delta6^(-m) + D*lambda2^m|   -- fast-decaying delta6^-m
       term would signal a level-dependent artifact (grid/extrapolation)
       riding on top of the normal lambda2 term; already disqualified for
       these 9 by hside_scan (hypothesis b rejected), kept here as a
       quantitative cross-check.

Fit method (exploratory, NOT rigorous nonlinear regression): for each zero
with >=3 valid (m, d_u) pairs, solve the 2x2 LINEAR system for (C2,C3) [or
(C,D)] from the two SMALLEST-m pairs (treating d_u as if it carries the sign
of the dominant term at that m -- both sign assignments are tried, the
better out-of-sample predictor is kept), then predict the remaining pair(s)
and report relative error. This is a coarse but disclosed and honest probe,
not a real weighted nonlinear least-squares fit.
[oid:deg6G] NOT RIGOROUS.
"""
import itertools
import json
import os

from mpmath import mp, mpf, fabs, nstr

mp.dps = 50

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)

LAM2 = C.LAM2
LAM3 = C.LAM3
DELTA6 = C.DELTA6

MR = json.load(open(os.path.join(BASE, 'match_results.json'), encoding='utf-8'))
rate_full_path = os.path.join(HERE, 'rate_full.json')
rate_full = json.load(open(rate_full_path, encoding='utf-8')) if os.path.exists(rate_full_path) else {}
rate_priority_path = os.path.join(HERE, 'rate_priority.json')
rate_secondary_path = os.path.join(HERE, 'rate_secondary.json')


def load_du_table(tid):
    """Prefer rate_full.json (full-tier search); fall back to the cheaper
    2-seed rate_priority/rate_secondary if rate_full incomplete/missing."""
    tab = {}
    for src_path in (rate_priority_path, rate_secondary_path):
        if os.path.exists(src_path):
            src = json.load(open(src_path, encoding='utf-8'))
            if tid in src:
                for k, v in src[tid]['d_u'].items():
                    if v is not None:
                        tab[int(k)] = v
    if tid in rate_full:
        for k, v in rate_full[tid]['d_u'].items():
            if v is not None:
                tab[int(k)] = v  # full-tier overrides cheap 2-seed
    return tab


def solve_2term(m_a, d_a, m_b, d_b, base1, base2, sign_a, sign_b):
    """Solve C*base1^m + D*base2^m = sign*d for two (m,d) pairs (linear in
    C,D); returns (C,D) or None if singular."""
    b1a, b2a = base1 ** m_a, base2 ** m_a
    b1b, b2b = base1 ** m_b, base2 ** m_b
    det = b1a * b2b - b2a * b1b
    if fabs(det) < mpf('1e-60'):
        return None
    rhs_a = sign_a * mpf(d_a)
    rhs_b = sign_b * mpf(d_b)
    C_ = (rhs_a * b2b - rhs_b * b2a) / det
    D_ = (b1a * rhs_b - b1b * rhs_a) / det
    return C_, D_


def try_model(du_table, base1, base2, name):
    ms = sorted(du_table.keys())
    if len(ms) < 3:
        return {'model': name, 'status': 'insufficient_points (%d)' % len(ms), 'n_points': len(ms)}
    m_a, m_b = ms[0], ms[1]
    d_a, d_b = du_table[m_a], du_table[m_b]
    best = None
    for sign_a, sign_b in itertools.product((1, -1), repeat=2):
        sol = solve_2term(m_a, d_a, m_b, d_b, base1, base2, sign_a, sign_b)
        if sol is None:
            continue
        C_, D_ = sol
        # predict remaining points
        errs = []
        preds = {}
        for m in ms[2:]:
            pred = float(fabs(C_ * base1 ** m + D_ * base2 ** m))
            obs = du_table[m]
            rel = abs(pred - obs) / obs if obs else float('inf')
            errs.append(rel)
            preds[m] = {'pred': pred, 'obs': obs, 'rel_err': rel}
        score = max(errs) if errs else 0.0
        if best is None or score < best['score']:
            best = {'C': nstr(C_, 10), 'D': nstr(D_, 10),
                    'fit_m': [m_a, m_b], 'signs': [sign_a, sign_b],
                    'predictions': preds, 'score': score}
    return {'model': name, 'status': 'ok', 'n_points': len(ms), 'points': du_table, **best}


# --- control sample: extend match_results.json rate_study to m=7,11 --------
rate_study = MR['rate_study']  # list of {'id','absu','d_u':{'8':...,'9':...,'10':...}}
control_ids = [r['id'] for r in rate_study]
control_extra_path = os.path.join(HERE, 'control_extra.json')
control_extra = json.load(open(control_extra_path, encoding='utf-8')) if os.path.exists(control_extra_path) else {}

report = {'targets': {}, 'controls': {}}

TARGET_IDS = ['S036', 'S104', 'S187', 'S339', 'S031', 'S070', 'S235', 'S279', 'S319']
for tid in TARGET_IDS:
    du = load_du_table(tid)
    du = {m: v for m, v in du.items() if v is not None and v < 0.5}  # drop wrong-branch outliers
    entry = {'d_u_table': du}
    entry['model_interference_l2l3'] = try_model(du, LAM2, LAM3, 'C2*l2^m+C3*l3^m')
    entry['model_delta_l2'] = try_model(du, 1 / DELTA6, LAM2, 'C*delta6^-m+D*l2^m')
    report['targets'][tid] = entry

for r in rate_study:
    tid = r['id']
    du = dict(r['d_u'])
    if tid in control_extra:
        du.update(control_extra[tid])
    du = {int(k): v for k, v in du.items() if v is not None}
    entry = {'d_u_table': du}
    entry['model_interference_l2l3'] = try_model(du, LAM2, LAM3, 'C2*l2^m+C3*l3^m')
    entry['model_delta_l2'] = try_model(du, 1 / DELTA6, LAM2, 'C*delta6^-m+D*l2^m')
    report['controls'][tid] = entry

with open(os.path.join(HERE, 'fit_results.json'), 'w', encoding='utf-8') as f:
    json.dump(report, f, indent=2, ensure_ascii=False, default=str)

print('=== TARGETS ===')
for tid, e in report['targets'].items():
    m1 = e['model_interference_l2l3']
    print('%-6s pts=%s  interference: %s  (score=%s)'
          % (tid, e['d_u_table'], m1.get('status'), m1.get('score')))
print('=== CONTROLS (subset) ===')
for tid, e in list(report['controls'].items())[:5]:
    m1 = e['model_interference_l2l3']
    print('%-6s pts=%s  interference: %s  (score=%s)'
          % (tid, e['d_u_table'], m1.get('status'), m1.get('score')))
print('wrote fit_results.json')

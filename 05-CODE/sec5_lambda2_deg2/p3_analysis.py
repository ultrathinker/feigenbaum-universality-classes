#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3_analysis.py -- Phase 3 of agent A's lambda2-degree2 retro-check.
[oid:lam2A]  NOT RIGOROUS.

Loads p1_deg2.json (spectrum) and rate_deg2.json (d_u(m) per zero, m=8..11),
computes per-zero ratios d(m)/d(m+1), medians per level-pair, and compares
against the candidate table required by BRIEF Phase 3:
    1/lambda2_essential (main hypothesis), 1/lambda3_essential, delta/alpha^2,
    alpha^2, delta^0.5.
Z72_true is EXCLUDED from ratio statistics (see progress.log / build script
docstring: its own construction is m=10-circular) but reported separately.

Output: p3_analysis.json (machine-readable) + prints the tables used to write
REPORT_LAMBDA2_DEG2.md by hand.
"""
import json
import time
from statistics import median

from mpmath import mp, mpf, nstr

mp.dps = 60

TAG = "[oid:lam2A] NOT RIGOROUS float pilot (agent A), phase 3 analysis"
T0 = time.time()
PROGRESS_LOG = 'progress.log'


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    with open(PROGRESS_LOG, 'a', encoding='utf-8') as f:
        f.write(time.strftime('[%Y-%m-%d %H:%M:%S] ') + '[phase3] ' + msg + '\n')


def main():
    log(TAG)
    p1 = json.load(open('p1_deg2.json', encoding='utf-8'))
    rate = json.load(open('rate_deg2.json', encoding='utf-8'))

    alpha = mpf(p1['alpha'])
    delta = mpf(p1['delta'])
    lam2 = mpf(p1['checks']['lambda2_essential'])
    lam3 = mpf(p1['checks']['lambda3_essential'])

    candidates = {
        '1/lambda2_essential (MAIN)': abs(1 / lam2),
        '1/lambda3_essential': abs(1 / lam3),
        'delta/alpha^2': delta / alpha ** 2,
        'alpha^2': alpha ** 2,
        'delta^0.5': delta ** mpf('0.5'),
        'delta (leading eigenvalue)': delta,
    }
    log('candidate table (magnitudes):')
    for name, val in candidates.items():
        log('  %-28s = %s' % (name, nstr(val, 15)))

    rows = rate['rows']
    per_zero_ratios = []
    for row in rows:
        if row.get('is_diagnostic'):
            continue
        du = {int(k): v for k, v in row['d_u'].items()}
        if not all(m in du for m in (8, 9, 10, 11)):
            log('  SKIP %s: incomplete d_u levels %s' % (row['id'], sorted(du)))
            continue
        r89 = du[8] / du[9]
        r910 = du[9] / du[10]
        r1011 = du[10] / du[11]
        per_zero_ratios.append({'id': row['id'], 'abs_u': row['abs_u'],
                                 'd8': du[8], 'd9': du[9], 'd10': du[10], 'd11': du[11],
                                 'r89': r89, 'r910': r910, 'r1011': r1011})

    log('=== per-zero table: d_u(m) and ratios ===')
    for r in per_zero_ratios:
        log('  %-6s |u|=%6.3f  d8=%.3e d9=%.3e d10=%.3e d11=%.3e   r89=%.3f r910=%.3f r1011=%.3f'
            % (r['id'], r['abs_u'], r['d8'], r['d9'], r['d10'], r['d11'],
               r['r89'], r['r910'], r['r1011']))

    med89 = median(r['r89'] for r in per_zero_ratios)
    med910 = median(r['r910'] for r in per_zero_ratios)
    med1011 = median(r['r1011'] for r in per_zero_ratios)
    log('MEDIANS: r89=%.4f  r910=%.4f  r1011=%.4f' % (med89, med910, med1011))

    # calibration zero: Z01 is the |u| closest to purely-real-axis-like "clean" one used
    # in quartic as the calibration anchor (here: smallest |u|, most benign geometry)
    cal = min(per_zero_ratios, key=lambda r: r['abs_u'])
    log('Calibration zero (smallest |u|, cleanest): %s  r89=%.4f r910=%.4f r1011=%.4f'
        % (cal['id'], cal['r89'], cal['r910'], cal['r1011']))

    log('=== comparison of MEDIAN ratios (all 3 level-pairs) against candidates ===')
    for label, med in (('r89', med89), ('r910', med910), ('r1011', med1011)):
        log('  %s median = %.4f' % (label, med))
        for name, val in candidates.items():
            log('      vs %-28s = %-10.4f  ratio(median/cand)=%.4f  diff%%=%.1f%%'
                % (name, float(val), med / float(val), 100 * abs(med - float(val)) / float(val)))

    # Z72 diagnostic (excluded from stats, reported qualitatively)
    z72row = next((r for r in rows if r['id'] == 'Z72_true'), None)
    z72_info = None
    if z72row:
        du72 = {int(k): v for k, v in z72row['d_u'].items()}
        z72_info = {'d8': du72.get(8), 'd9': du72.get(9), 'd11': du72.get(11)}
        if 8 in du72 and 9 in du72:
            log('Z72_true diagnostic: d8=%.3e d9=%.3e r89=%.3f (m=10 excluded by construction, d11=%.3e)'
                % (du72[8], du72[9], du72[8] / du72[9], du72.get(11, float('nan'))))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'not_rigorous': True,
                 'elapsed_sec': time.time() - T0},
        'lambda2_essential': nstr(lam2, 30), 'lambda3_essential': nstr(lam3, 30),
        'alpha': nstr(alpha, 30), 'delta': nstr(delta, 30),
        'candidates': {k: float(v) for k, v in candidates.items()},
        'per_zero_ratios': per_zero_ratios,
        'medians': {'r89': med89, 'r910': med910, 'r1011': med1011},
        'calibration_zero': cal['id'],
        'z72_diagnostic': z72_info,
    }
    with open('p3_analysis.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote p3_analysis.json')


if __name__ == '__main__':
    main()

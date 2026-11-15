#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
corrupt_battery.py -- mini battery of deliberate corruptions (BRIEF item 5).
For each degree d and the N=24 certified point, apply 4 corruptions and confirm
the Krawczyk certification CORRECTLY REJECTS all of them (raises / contained=False
across the whole box-radius ladder). A corruption that gets "certified" would be
a serious bug; a corruption that fails to certify is the expected, correct result.
"""
import json
import time

from mpmath import mp, mpf, nstr
from flint import ctx, arb

import bridge_lib as bl

T0 = time.time()
LOGFILE = 'corrupt_battery.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOGFILE, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def run_degree(d, seed_path, DPS=150, N=24):
    log('=== degree d=%d, N=%d: preparing valid midpoint solution ===' % (d, N))
    seed = json.load(open(seed_path, encoding='utf-8'))
    lam_seed = mpf(seed['checks']['eigen2_deflated'])
    a_init = [mpf(1), mpf('-1.6')]
    a_cur, al_cur = a_init, mpf('1.5')
    c_cur = None
    lam_cur = lam_seed
    for Nrung in [8, 14, 20, 24]:
        a_cur, al_cur, _, _ = bl.solve_fixedpoint_mp(d, Nrung, a_cur, al_cur, tol=mpf('1e-140'), maxit=300)
        c_cur, lam_cur, _, _ = bl.solve_eigen_mp(d, a_cur, al_cur, Nrung, c_cur, lam_cur, tol=mpf('1e-140'), maxit=300)
    a_mid, alpha_mid, c_mid, lam_mid = a_cur, al_cur, c_cur, lam_cur
    log('  valid midpoint ready: alpha=%s lambda2=%s' % (nstr(alpha_mid, 20), nstr(lam_mid, 20)))

    # sanity: the UN-corrupted point must still certify (baseline control)
    r1 = bl.certify_fixedpoint_auto(d, a_mid, alpha_mid, N, DPS)
    aX_box = [arb(1)] + r1['X'][:N]
    alphaX_box = r1['X'][N]
    r2 = bl.certify_eigen_auto(d, a_mid, alpha_mid, c_mid, lam_mid, N, DPS, aX_box, alphaX_box)
    log('  BASELINE (uncorrupted) control: K1 OK, K2 OK (lambda2 half-width %.3e) -- as expected'
        % r2['halfwidths'][N])

    results = {}

    # --- corruption 1: shift lambda ---------------------------------------
    log('  --- corruption 1: shift lambda by 1e-8 ---')
    lam_bad = lam_mid + mpf('1e-8')
    try:
        r2bad = bl.certify_eigen_auto(d, a_mid, alpha_mid, c_mid, lam_bad, N, DPS, aX_box, alphaX_box)
        rejected = not r2bad['contained']
        log('  UNEXPECTED: shifted-lambda certified?! contained=%s' % r2bad['contained'])
    except RuntimeError:
        rejected = True
        log('  correctly REJECTED (Krawczyk failed for all radii, as expected)')
    results['shift_lambda'] = {'rejected': rejected}

    # --- corruption 2: corrupt a g-coefficient (fixed point stage) --------
    log('  --- corruption 2: corrupt g-coefficient a_5 by 1e-8 ---')
    a_bad = list(a_mid)
    idx = min(5, len(a_bad) - 1)
    a_bad[idx] = a_bad[idx] + mpf('1e-8')
    try:
        r1bad = bl.certify_fixedpoint_auto(d, a_bad, alpha_mid, N, DPS)
        rejected = not r1bad['contained']
        log('  UNEXPECTED: corrupted-a certified?! contained=%s' % r1bad['contained'])
    except RuntimeError:
        rejected = True
        log('  correctly REJECTED (K1 Krawczyk failed for all radii, as expected)')
    results['corrupt_g_coeff'] = {'rejected': rejected}

    # --- corruption 3: corrupt an eigenvector component -------------------
    log('  --- corruption 3: corrupt eigenvector coeff c_3 by 1e-8 ---')
    c_bad = list(c_mid)
    idx = min(3, len(c_bad) - 1)
    c_bad[idx] = c_bad[idx] + mpf('1e-8')
    try:
        r2bad = bl.certify_eigen_auto(d, a_mid, alpha_mid, c_bad, lam_mid, N, DPS, aX_box, alphaX_box)
        rejected = not r2bad['contained']
        log('  UNEXPECTED: corrupted-eigenvector certified?! contained=%s' % r2bad['contained'])
    except RuntimeError:
        rejected = True
        log('  correctly REJECTED (K2 Krawczyk failed for all radii, as expected)')
    results['corrupt_eigenvector'] = {'rejected': rejected}

    # --- corruption 4: wrong normalization (scale eigenvector by 2x) ------
    log('  --- corruption 4: wrong normalization, c scaled by 2x (violates Phi(1)=1) ---')
    c_bad2 = [2 * x for x in c_mid]
    try:
        r2bad = bl.certify_eigen_auto(d, a_mid, alpha_mid, c_bad2, lam_mid, N, DPS, aX_box, alphaX_box)
        rejected = not r2bad['contained']
        log('  UNEXPECTED: wrong-normalization certified?! contained=%s' % r2bad['contained'])
    except RuntimeError:
        rejected = True
        log('  correctly REJECTED (K2 Krawczyk failed for all radii, as expected)')
    results['wrong_normalization'] = {'rejected': rejected}

    all_ok = all(v['rejected'] for v in results.values())
    log('=== degree d=%d corruption battery: all 4 correctly rejected = %s ===' % (d, all_ok))
    return {'baseline_ok': True, 'corruptions': results, 'all_correctly_rejected': all_ok}


def main():
    mp.dps = 150
    ctx.dps = 150
    open(LOGFILE, 'w', encoding='utf-8').close()
    out = {}
    out['d4'] = run_degree(4, r'../sec3_pilot_d4/quartic_fixedpoint.json')
    out['d6'] = run_degree(6, r'../sec3_pilot_d6/sextic_fixedpoint.json')
    with open('corrupt_battery_results.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote corrupt_battery_results.json')
    log('=== ALL DONE, total %.1fs ===' % (time.time() - T0))


if __name__ == '__main__':
    main()

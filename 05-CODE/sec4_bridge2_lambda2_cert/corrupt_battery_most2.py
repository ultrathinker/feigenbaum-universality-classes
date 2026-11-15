#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
corrupt_battery_most2.py -- MOST-2 (Chebyshev-basis) corruption battery.

For each corruption we test Krawczyk containment at a FIXED box radius equal
to the radius that was actually certified for the true (uncorrupted) point at
that (d, N) -- we do NOT sweep the full DEFAULT_RADII ladder for corrupted
points, because sweeping up to a much larger radius would trivially "recover"
the true root (which lies inside a big enough ball around a small
perturbation) and prove nothing about the corruption. Testing at the ORIGINAL
certified radius is the sharp, meaningful test: does the corrupted candidate
still admit a Krawczyk box of the same width that worked for the true one?

Every corruption is expected to make containment FAIL (contained=False, or
K1/K2 chain raising because the box is not even well formed). A positive
control (the true, uncorrupted point) is run first/alongside and MUST
succeed, proving the harness itself is not just always failing.

Runs entirely inside . Reads centers
from most2_run_d4.json / most2_run_d6.json (already-certified midpoints from
run_most2.py). Uses <=1 process (well inside the 4-core budget).
"""
import json
import time

from mpmath import mp, mpf, nstr
from flint import ctx, arb

import most2_lib as ml

T0 = time.time()
LOGFILE = 'TAMPER_MOST2.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOGFILE, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def log_progress(msg):
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    line = '%s %s' % (ts, msg)
    with open('progress.log', 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


# Decoy eigenvalues / calibration constants (brief §2, литералы из БРИФА-агента)
LAMBDA3 = {4: mpf('-0.25566455847429740633'), 6: mpf('-0.34196832254780210510')}


def load_record(d, N):
    fn = 'most2_run_d%d.json' % d
    data = json.load(open(fn, encoding='utf-8'))
    rec = data[str(N)]
    a_full = [mpf(v) for v in rec['b_mid']]
    alpha_mid = mpf(rec['alpha_mid'])
    c_mid = [mpf(v) for v in rec['c_mid']]
    lam_mid = mpf(rec['lambda2_mid'])
    k1_r = rec['k1_box_r']
    k2_r = rec['k2_box_r']
    return a_full, alpha_mid, c_mid, lam_mid, k1_r, k2_r


def k1_test(d, N, a_full, alpha_mid, r, dps):
    """Direct single-radius K1 Krawczyk test (no ladder sweep)."""
    res = ml.krawczyk_generic(
        list(a_full[1:]) + [alpha_mid],
        ml.fp_residual_jacobian_mp(list(a_full[1:]), alpha_mid, d, ml.cheb_vals_mp(N + 1))[1],
        lambda v: ml.fp_F_J_arb(d, N, [ml.b0_of(v[:N], N)] + v[:N], v[N],
                                ml.gauss_nodes_arb(N + 1), want_J=False)[0],
        lambda v: ml.fp_F_J_arb(d, N, [ml.b0_of(v[:N], N)] + v[:N], v[N],
                                ml.gauss_nodes_arb(N + 1), want_J=True),
        dps, r)
    return res


def k2_test(d, N, a_full, alpha_mid, aX_box, alphaX_box, c_mid, lam_mid, r, dps):
    """Direct single-radius K2 Krawczyk test (no ladder sweep), chained on a
    valid (uncorrupted) K1 box aX_box/alphaX_box unless the corruption targets
    K1 itself."""
    res = ml.krawczyk_generic(
        list(c_mid) + [lam_mid],
        ml.eigen_F_J_mp(a_full, alpha_mid, c_mid, lam_mid, d, ml.cheb_vals_mp(N))[1],
        lambda v: ml.eigen_F_J_arb(d, N, aX_box, alphaX_box, v[:N], v[N], ml.gauss_nodes_arb(N))[0],
        lambda v: ml.eigen_F_J_arb(d, N, aX_box, alphaX_box, v[:N], v[N], ml.gauss_nodes_arb(N)),
        dps, r)
    return res


def valid_k1_box(d, N, a_full, alpha_mid, r, dps):
    res1 = k1_test(d, N, a_full, alpha_mid, r, dps)
    aX_box = [ml.b0_of(res1['X'][:N], N)] + res1['X'][:N]
    alphaX_box = res1['X'][N]
    return res1, aX_box, alphaX_box


def run_case(name, fn):
    """fn() -> dict with 'contained' (bool) or raises. Returns result dict."""
    try:
        res = fn()
        contained = bool(res['contained'])
        log('  [%s] contained=%s%s' % (name, contained,
            '' if contained else ' fail_idx=%s' % res.get('fail_idx')))
        return {'ran': True, 'contained': contained, 'error': None}
    except Exception as e:  # noqa: BLE001 -- battery must never crash on expected failures
        log('  [%s] EXCEPTION (%s): %s' % (name, type(e).__name__, e))
        return {'ran': True, 'contained': False, 'error': '%s: %s' % (type(e).__name__, e)}


def battery_for(d, N, dps, extra_tag=''):
    log('=== battery d=%d N=%d%s ===' % (d, N, (' [' + extra_tag + ']') if extra_tag else ''))
    a_full, alpha_mid, c_mid, lam_mid, k1_r, k2_r = load_record(d, N)
    results = {}

    # --- positive control: uncorrupted point MUST certify at its own radius --
    res1_ctrl, aX_box, alphaX_box = valid_k1_box(d, N, a_full, alpha_mid, k1_r, dps)
    res2_ctrl = k2_test(d, N, a_full, alpha_mid, aX_box, alphaX_box, c_mid, lam_mid, k2_r, dps)
    ctrl_ok = res1_ctrl['contained'] and res2_ctrl['contained']
    log('  [control] K1 contained=%s (r=%s)  K2 contained=%s (r=%s)  -> %s'
        % (res1_ctrl['contained'], k1_r, res2_ctrl['contained'], k2_r,
           'OK (positive control passes, as required)' if ctrl_ok else 'FAILED (harness bug!)'))
    results['control_uncorrupted'] = {'k1_contained': res1_ctrl['contained'],
                                       'k2_contained': res2_ctrl['contained'],
                                       'expected': 'contained=True (must certify)',
                                       'pass': ctrl_ok}

    def k2_at(lam_bad=None, c_bad=None, r=None):
        lam = lam_bad if lam_bad is not None else lam_mid
        c = c_bad if c_bad is not None else c_mid
        return k2_test(d, N, a_full, alpha_mid, aX_box, alphaX_box, c, lam, r or k2_r, dps)

    # 1. shift lambda center by 10x the certified K2 radius -------------------
    shift10 = 10 * mpf(k2_r)
    r1 = run_case('shift_lambda_10r', lambda: k2_at(lam_bad=lam_mid + shift10))
    results['shift_lambda_10r'] = dict(r1, expected='contained=False (broken)',
                                        detail='lambda center shifted by +10*box_r=%s at r=%s' % (shift10, k2_r))

    # 2. shift lambda center by 1e-30 (>> box radius 1e-40 or 1e-38) ----------
    shift_tiny = mpf('1e-30')
    r2 = run_case('shift_lambda_1e-30', lambda: k2_at(lam_bad=lam_mid + shift_tiny))
    results['shift_lambda_1e-30'] = dict(r2, expected='contained=False (broken)',
                                          detail='lambda center shifted by +1e-30 at r=%s (1e-30 >> box radius)' % k2_r)

    # 3. substitute lambda2 with lambda3 (wrong branch) ------------------------
    r3 = run_case('swap_to_lambda3', lambda: k2_at(lam_bad=LAMBDA3[d]))
    results['swap_to_lambda3'] = dict(r3, expected='contained=False (broken)',
                                       detail='lambda2 replaced by lambda3=%s' % LAMBDA3[d])

    # 4. substitute lambda2 with calibration alpha^{-d} ------------------------
    calib = alpha_mid ** (-d)
    r4 = run_case('swap_to_calibration', lambda: k2_at(lam_bad=calib))
    results['swap_to_calibration'] = dict(r4, expected='contained=False (broken)',
                                           detail='lambda2 replaced by alpha^-%d=%s' % (d, nstr(calib, 15)))

    # 5. corrupt one Chebyshev coefficient b_k of the fixed point (K1) --------
    kbad = min(5, N - 1)
    a_bad = list(a_full)
    a_bad[kbad + 1] = a_bad[kbad + 1] + mpf('1e-20')  # a_bad index 0 is b0 (recomputed anyway)
    r5 = run_case('corrupt_b_coeff', lambda: k1_test(d, N, a_bad, alpha_mid, k1_r, dps))
    results['corrupt_b_coeff'] = dict(r5, expected='contained=False (broken)',
                                       detail='b_%d += 1e-20 at r=%s' % (kbad + 1, k1_r))

    # 6. corrupt one eigenvector coefficient c_k (K2) --------------------------
    kbad2 = min(3, N - 1)
    c_bad = list(c_mid)
    c_bad[kbad2] = c_bad[kbad2] + mpf('1e-20')
    r6 = run_case('corrupt_c_coeff', lambda: k2_at(c_bad=c_bad))
    results['corrupt_c_coeff'] = dict(r6, expected='contained=False (broken)',
                                       detail='c_%d += 1e-20 at r=%s' % (kbad2 + 1, k2_r))

    # 7. corrupt alpha by 1e-25 (K1) --------------------------------------------
    alpha_bad = alpha_mid + mpf('1e-25')
    r7 = run_case('corrupt_alpha', lambda: k1_test(d, N, a_full, alpha_bad, k1_r, dps))
    results['corrupt_alpha'] = dict(r7, expected='contained=False (broken)',
                                     detail='alpha += 1e-25 at r=%s' % k1_r)

    # 8. corrupt normalization: scale eigenvector by 2x (breaks Phi(1)=1) ------
    c_scaled = [2 * x for x in c_mid]
    r8 = run_case('wrong_normalization', lambda: k2_at(c_bad=c_scaled))
    results['wrong_normalization'] = dict(r8, expected='contained=False (broken)',
                                           detail='c scaled by 2x (Phi(1)=2 instead of 1) at r=%s' % k2_r)

    # 9. inflate box radius to 1e-3 on the TRUE (uncorrupted) center -----------
    # Not a corrupted CENTER -- a stress test of an absurdly loose radius.
    # Report honestly what happens (brief explicitly asks for this).
    try:
        r9k1 = k1_test(d, N, a_full, alpha_mid, '1e-3', dps)
        k1_note = 'K1 contained=%s at r=1e-3' % r9k1['contained']
        k1_ok = r9k1['contained']
    except Exception as e:  # noqa: BLE001
        k1_note = 'K1 EXCEPTION %s: %s' % (type(e).__name__, e)
        k1_ok = False
    try:
        # chain on the (now huge) K1 box if it succeeded, else reuse the tight one
        if k1_ok:
            aX_box9 = [ml.b0_of(r9k1['X'][:N], N)] + r9k1['X'][:N]
            alphaX_box9 = r9k1['X'][N]
        else:
            aX_box9, alphaX_box9 = aX_box, alphaX_box
        r9k2 = k2_test(d, N, a_full, alpha_mid, aX_box9, alphaX_box9, c_mid, lam_mid, '1e-3', dps)
        k2_note = 'K2 contained=%s at r=1e-3' % r9k2['contained']
        k2_ok = r9k2['contained']
    except Exception as e:  # noqa: BLE001
        k2_note = 'K2 EXCEPTION %s: %s' % (type(e).__name__, e)
        k2_ok = False
    log('  [inflate_radius_1e-3] %s | %s' % (k1_note, k2_note))
    results['inflate_radius_1e-3'] = {
        'ran': True, 'k1_contained': k1_ok, 'k2_contained': k2_ok,
        'note_k1': k1_note, 'note_k2': k2_note,
        'expected': 'documented honestly (K1 may still trivially "succeed" since alpha is extremely '
                    'well-determined; K2 typically degenerates/fails -- not a hard pass/fail requirement)',
        'pass': None,  # judged qualitatively in the report, not pass/fail
    }

    n_meaningful = [k for k in results if k not in ('control_uncorrupted', 'inflate_radius_1e-3')]
    n_broken = sum(1 for k in n_meaningful if not results[k]['contained'])
    log('=== d=%d N=%d: control OK=%s, corruptions correctly broken = %d/%d ==='
        % (d, N, ctrl_ok, n_broken, len(n_meaningful)))
    return results, ctrl_ok, n_broken, len(n_meaningful)


def main():
    mp.dps = 170
    ctx.dps = 160
    open(LOGFILE, 'w', encoding='utf-8').close()
    log_progress('STEP9 Батарея порч most2: старт (N=24 d=4/d=6, plus N=40 d=4).')

    out = {}
    total_broken = 0
    total_meaningful = 0
    all_controls_ok = True

    for d, N in [(4, 24), (6, 24), (4, 40)]:
        dps = ml.ctx.dps - 30
        res, ctrl_ok, nb, nm = battery_for(d, N, dps)
        out['d%d_N%d' % (d, N)] = res
        total_broken += nb
        total_meaningful += nm
        all_controls_ok = all_controls_ok and ctrl_ok
        log_progress('STEP9 d=%d N=%d: control_ok=%s, corruptions broken=%d/%d'
                      % (d, N, ctrl_ok, nb, nm))

    out['_summary'] = {
        'all_controls_ok': all_controls_ok,
        'total_corruptions_meaningful': total_meaningful,
        'total_corruptions_correctly_broken': total_broken,
        'all_broken': total_broken == total_meaningful,
    }
    with open('corrupt_battery_most2_results.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False, default=str)
    log('wrote corrupt_battery_most2_results.json')
    log('=== FINAL: controls_ok=%s, broken=%d/%d, ALL_BROKEN=%s ===' %
        (all_controls_ok, total_broken, total_meaningful, total_broken == total_meaningful))
    log_progress('STEP9 Батарея порч most2 ЗАВЕРШЕНА: controls_ok=%s, broken=%d/%d.'
                  % (all_controls_ok, total_broken, total_meaningful))


if __name__ == '__main__':
    main()

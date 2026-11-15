#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_full.py -- [K-D205-CHECK] NOT RIGOROUS. Blind independent check, degree
family 1-a|x|^d, target d=2.05, calibration d=2.0. Writes progress.log
(flushed after every step) and d205_result.json at the end (and after each
major stage, so partial results survive a crash).
"""
import json
import time
import sys
from mpmath import mp, mpf, nstr

mp.dps = 60
import core

T0 = time.time()
LOG = open('progress.log', 'a', encoding='utf-8')


def log(msg):
    core.log(LOG, T0, msg)


def eig_stack(a, al, d, N, ndeflate=5):
    """Return list of (eigenvalue, eigenvector) for the top ndeflate modes of
    the collocation matrix, in decreasing |eigenvalue| order."""
    A = core.build_DT_matrix(a, al, d, N)
    out = []
    M = A
    for k in range(ndeflate):
        lam, v, it = core.power_iter(M, N)
        out.append((lam, v, it))
        M = core.deflate(M, lam, v, N)
    return out


def classify(a, al, d, N, eigs):
    """Identify the alpha^-d gauge mode among the deflated eigenvalues; return
    dict with delta, gauge_mode, lambda2, lambda3, lambda4_extra (list), plus
    diagnostics."""
    alpha_neg_d = al ** (-d)
    delta = eigs[0][0]
    rest = eigs[1:]
    # find whichever of the 'rest' eigenvalues is closest (relative) to alpha^-d
    best_idx = None
    best_rel = None
    for i, (lam, v, it) in enumerate(rest):
        rel = abs(lam - alpha_neg_d) / abs(alpha_neg_d)
        if best_rel is None or rel < best_rel:
            best_rel = rel
            best_idx = i
    gauge = rest[best_idx]
    physical = [e for i, e in enumerate(rest) if i != best_idx]
    return {
        'delta': delta,
        'gauge_eig': gauge[0],
        'gauge_rel_diff_from_alpha_negd': best_rel,
        'alpha_neg_d': alpha_neg_d,
        'lambda2': physical[0][0] if len(physical) > 0 else None,
        'lambda2_vec': physical[0][1] if len(physical) > 0 else None,
        'lambda3': physical[1][0] if len(physical) > 1 else None,
        'lambda3_vec': physical[1][1] if len(physical) > 1 else None,
        'lambda4': physical[2][0] if len(physical) > 2 else None,
    }


def full_pipeline_at_N(N, d0=mpf(2), dtarget=mpf('2.05'), dsteps=10,
                        a_start=None, al_start=None, N_ladder=None):
    """Solve fixed point at d0 (with N-ladder continuation to N), then continue
    in d up to dtarget in dsteps, returns dict with everything needed."""
    result = {}
    if N_ladder is None:
        N_ladder = [4, 8, 12, 16, 20, 24, 28, 32, 36, 40]
    N_ladder = [n for n in N_ladder if n <= N]
    if N_ladder[-1] != N:
        N_ladder.append(N)

    log('=== N=%d: fixed point at d=2.0, N-ladder %s ===' % (N, N_ladder))
    a_cur, al_cur = a_start, al_start
    for Nl in N_ladder:
        a_cur, al_cur, it, res = core.solve_fixedpoint(Nl, d0, a_cur, al_cur)
        log('  N=%2d d=2.0  alpha=%s  iters=%d  colloc_res=%.3e'
            % (Nl, nstr(al_cur, 25), it, float(res)))
    a20, al20 = a_cur, al_cur
    dense20 = core.fp_residual_dense(a20, al20, d0)
    log('  d=2.0 N=%d dense residual = %.3e' % (N, float(dense20)))

    log('=== N=%d: spectrum at d=2.0 ===' % N)
    eigs20 = eig_stack(a20, al20, d0, N, ndeflate=5)
    for i, (lam, v, it) in enumerate(eigs20):
        log('  eig#%d = %s  (power-iter %d its)' % (i + 1, nstr(lam, 20), it))
    cls20 = classify(a20, al20, d0, N, eigs20)
    log('  identified gauge mode ~ alpha^-d: eig=%s alpha^-d=%s rel_diff=%.3e'
        % (nstr(cls20['gauge_eig'], 20), nstr(cls20['alpha_neg_d'], 20),
           float(cls20['gauge_rel_diff_from_alpha_negd'])))
    log('  => delta=%s  lambda2=%s  lambda3=%s'
        % (nstr(cls20['delta'], 20), nstr(cls20['lambda2'], 20), nstr(cls20['lambda3'], 20)))

    # eigenpair functional residuals (delta, lambda2, lambda3) -- build coeff
    # vectors (c_0=0, c_1..c_N = eigenvector components) and check DT h = lam h
    # on a dense grid.
    def hvec(v):
        return [mpf(0)] + [v[i] for i in range(N)]

    res_delta = core.eigenpair_residual(a20, al20, d0, hvec(eigs20[0][1]), cls20['delta'])
    res_l2 = core.eigenpair_residual(a20, al20, d0, hvec(cls20['lambda2_vec']), cls20['lambda2'])
    res_l3 = core.eigenpair_residual(a20, al20, d0, hvec(cls20['lambda3_vec']), cls20['lambda3'])
    log('  eigenpair residuals (dense grid, d=2.0): delta=%.3e lambda2=%.3e lambda3=%.3e'
        % (float(res_delta), float(res_l2), float(res_l3)))

    result['d2.00'] = {
        'N': N, 'alpha': nstr(al20, 40), 'dense_fp_residual': float(dense20),
        'delta': nstr(cls20['delta'], 30), 'delta_eigresidual': float(res_delta),
        'gauge_eig_alpha_negd': nstr(cls20['gauge_eig'], 30),
        'gauge_rel_diff': float(cls20['gauge_rel_diff_from_alpha_negd']),
        'lambda2': nstr(cls20['lambda2'], 30), 'lambda2_eigresidual': float(res_l2),
        'lambda3': nstr(cls20['lambda3'], 30), 'lambda3_eigresidual': float(res_l3),
        'lambda4_extra': nstr(cls20['lambda4'], 20) if cls20['lambda4'] is not None else None,
    }

    # ---- continuation in d up to dtarget --------------------------------------
    log('=== N=%d: continuation d: 2.00 -> %s in %d steps ===' % (N, nstr(dtarget, 6), dsteps))
    a_cur, al_cur = a20, al20
    d_cur = d0
    step = (dtarget - d0) / dsteps
    for k in range(1, dsteps + 1):
        d_new = d0 + step * k
        a_cur, al_cur, it, res = core.solve_fixedpoint(N, d_new, a_cur, al_cur)
        log('  d=%s  alpha=%s  iters=%d  colloc_res=%.3e' % (nstr(d_new, 8), nstr(al_cur, 20), it, float(res)))
        d_cur = d_new
    a_t, al_t = a_cur, al_cur
    dense_t = core.fp_residual_dense(a_t, al_t, dtarget)
    log('  d=%s N=%d dense residual = %.3e' % (nstr(dtarget, 6), N, float(dense_t)))

    log('=== N=%d: spectrum at d=%s ===' % (N, nstr(dtarget, 6)))
    eigs_t = eig_stack(a_t, al_t, dtarget, N, ndeflate=5)
    for i, (lam, v, it) in enumerate(eigs_t):
        log('  eig#%d = %s  (power-iter %d its)' % (i + 1, nstr(lam, 20), it))
    cls_t = classify(a_t, al_t, dtarget, N, eigs_t)
    log('  identified gauge mode ~ alpha^-d: eig=%s alpha^-d=%s rel_diff=%.3e'
        % (nstr(cls_t['gauge_eig'], 20), nstr(cls_t['alpha_neg_d'], 20),
           float(cls_t['gauge_rel_diff_from_alpha_negd'])))
    log('  => delta=%s  lambda2=%s  lambda3=%s'
        % (nstr(cls_t['delta'], 20), nstr(cls_t['lambda2'], 20), nstr(cls_t['lambda3'], 20)))

    res_delta_t = core.eigenpair_residual(a_t, al_t, dtarget, hvec(eigs_t[0][1]), cls_t['delta'])
    res_l2_t = core.eigenpair_residual(a_t, al_t, dtarget, hvec(cls_t['lambda2_vec']), cls_t['lambda2'])
    res_l3_t = core.eigenpair_residual(a_t, al_t, dtarget, hvec(cls_t['lambda3_vec']), cls_t['lambda3'])
    log('  eigenpair residuals (dense grid, d=%s): delta=%.3e lambda2=%.3e lambda3=%.3e'
        % (nstr(dtarget, 6), float(res_delta_t), float(res_l2_t), float(res_l3_t)))

    # ---- probe #3: eigenvalue near alpha(d)^-d ? -------------------------------
    alpha_negd_t = al_t ** (-dtarget)
    all_eigs_t = [cls_t['delta'], cls_t['gauge_eig'], cls_t['lambda2'], cls_t['lambda3']]
    if cls_t['lambda4'] is not None:
        all_eigs_t.append(cls_t['lambda4'])
    rel_diffs = [(abs(e - alpha_negd_t) / abs(alpha_negd_t), e) for e in all_eigs_t]
    rel_diffs.sort(key=lambda p: p[0])
    log('  PROBE#3 at d=%s: alpha^-d=%s ; closest eigenvalue=%s rel_diff=%.3e'
        % (nstr(dtarget, 6), nstr(alpha_negd_t, 20), nstr(rel_diffs[0][1], 20), float(rel_diffs[0][0])))

    result['d%s' % nstr(dtarget, 4)] = {
        'N': N, 'alpha': nstr(al_t, 40), 'dense_fp_residual': float(dense_t),
        'delta': nstr(cls_t['delta'], 30), 'delta_eigresidual': float(res_delta_t),
        'gauge_eig_alpha_negd': nstr(cls_t['gauge_eig'], 30),
        'gauge_rel_diff': float(cls_t['gauge_rel_diff_from_alpha_negd']),
        'lambda2': nstr(cls_t['lambda2'], 30), 'lambda2_eigresidual': float(res_l2_t),
        'lambda3': nstr(cls_t['lambda3'], 30), 'lambda3_eigresidual': float(res_l3_t),
        'lambda4_extra': nstr(cls_t['lambda4'], 20) if cls_t['lambda4'] is not None else None,
        'probe3_alpha_neg_d': nstr(alpha_negd_t, 30),
        'probe3_closest_eig': nstr(rel_diffs[0][1], 30),
        'probe3_rel_diff': float(rel_diffs[0][0]),
    }

    return result, (a20, al20), (a_t, al_t)


def main():
    log('=== START run_full.py [K-D205-CHECK] mp.dps=%d ===' % mp.dps)
    all_results = {}

    # N=32 pipeline
    res32, fp20_32, fpT_32 = full_pipeline_at_N(32)
    all_results['N32'] = res32
    with open('d205_result.json', 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    log('wrote partial d205_result.json (N=32 done)')

    # N=40 pipeline (fresh N-ladder continuation, independent of N=32 run except
    # sharing the same N-ladder start)
    res40, fp20_40, fpT_40 = full_pipeline_at_N(40)
    all_results['N40'] = res40
    with open('d205_result.json', 'w', encoding='utf-8') as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    log('wrote d205_result.json (N=40 done)')

    # N-drift summary
    log('=== N-drift summary (N32 vs N40) ===')
    for dkey in res32:
        d32 = res32[dkey]
        d40 = res40[dkey]
        for key in ['alpha', 'delta', 'lambda2', 'lambda3']:
            v32 = mpf(d32[key])
            v40 = mpf(d40[key])
            log('  %s %s: N32=%s N40=%s |diff|=%.3e' % (dkey, key, nstr(v32, 20), nstr(v40, 20), float(abs(v32 - v40))))

    log('=== DONE ===')
    LOG.close()


if __name__ == '__main__':
    main()

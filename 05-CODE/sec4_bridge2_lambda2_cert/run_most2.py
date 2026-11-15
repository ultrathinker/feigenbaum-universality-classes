#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_most2.py -- production run: Chebyshev-basis rigorous certificates for
lambda2(4) and lambda2(6) at N = 16,24,32,40 (+48 if easy), with K1/K2
Krawczyk + independent Gauss-Seidel + anchor-coverage + truncation ladder."""
import json
import sys
import time

from mpmath import mp, mpf, mpmathify, nstr
from most2_lib import (solve_fixedpoint_mp, solve_eigen_mp, certify_fixedpoint,
                       certify_eigen, gauss_seidel_verify, cond_of, cheb_vals_mp,
                       b0_of, eigen_F_J_mp, fp_residual_jacobian_mp)

LOG_PATH = 'progress.log'
T0 = time.time()


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def log_mid(msg):
    """intermediate progress line (no timestamp bloat)."""
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def main():
    d = int(sys.argv[1])
    ladder = [int(x) for x in sys.argv[2].split(',')]
    anchor_str = sys.argv[3]
    anchor = mpf(anchor_str)
    dps = 130
    log('=== run_most2 d=%d ladder=%s anchor=%s ===' % (d, ladder, anchor_str))
    out = {}
    x_cur, al_cur = None, None
    prev_lam = None
    prev_lam_str = None
    for N in ladder:
        log('-- N=%d --' % N)
        x_cur, al_cur, it, res = solve_fixedpoint_mp(d, N, x_cur, al_cur, log=log_mid)
        a_full = [mpf(1) - sum(((-1) ** k) * x_cur[k - 1] for k in range(1, N + 1))] + list(x_cur)
        log('  FP N=%d it=%d res=%.2e alpha=%s' % (N, it, float(res), nstr(al_cur, 40)))
        c, lam, ite, rese = solve_eigen_mp(d, a_full, al_cur, N, anchor,
                                           lam_init=prev_lam, log=log_mid)
        log('  EIG N=%d it=%d res=%.2e lambda2=%s' % (N, ite, float(rese), nstr(lam, 45)))
        # K1
        log_mid('  K1 cert...')
        res1 = certify_fixedpoint(d, N, a_full, al_cur, dps, log=log_mid)
        aX_box = [b0_of(res1['X'][:N], N)] + res1['X'][:N]
        alphaX_box = res1['X'][N]
        # GS verify K1
        pts_mp = cheb_vals_mp(N + 1)
        Fmp, Jmp = fp_residual_jacobian_mp(list(a_full[1:]), al_cur, d, pts_mp)
        xhat1 = list(a_full[1:]) + [al_cur]

        def F1_at(v):
            b = [b0_of(v[:N], N)] + v[:N]
            from most2_lib import gauss_nodes_arb, fp_F_J_arb
            return fp_F_J_arb(d, N, b, v[N], gauss_nodes_arb(N + 1), want_J=False)[0]

        def F1_box(v):
            b = [b0_of(v[:N], N)] + v[:N]
            from most2_lib import gauss_nodes_arb, fp_F_J_arb
            return fp_F_J_arb(d, N, b, v[N], gauss_nodes_arb(N + 1), want_J=True)

        gs1 = gauss_seidel_verify(xhat1, Jmp, F1_at, F1_box, dps, res1['box_r_used'])
        log('  K1 OK contained=%s box_r=%s alpha_hw=%.3e GS=%s'
            % (res1['contained'], res1['box_r_used'], res1['halfwidths'][N], gs1['contained']))
        # K2
        log_mid('  K2 cert...')
        res2 = certify_eigen(d, N, a_full, al_cur, c, lam, dps, aX_box, alphaX_box, log=log_mid)
        # GS verify K2
        from most2_lib import eigen_F_J_arb, gauss_nodes_arb
        Fmp2, Jmp2 = eigen_F_J_mp([b0_of(a_full[1:], N)] + list(a_full[1:]), al_cur,
                                  c, lam, d, cheb_vals_mp(N))
        xhat2 = list(c) + [lam]

        def F2_at(v):
            return eigen_F_J_arb(d, N, aX_box, alphaX_box, v[:N], v[N], gauss_nodes_arb(N))[0]

        def F2_box(v):
            return eigen_F_J_arb(d, N, aX_box, alphaX_box, v[:N], v[N], gauss_nodes_arb(N))

        gs2 = gauss_seidel_verify(xhat2, Jmp2, F2_at, F2_box, dps, res2['box_r_used'])
        condJ2 = float(cond_of(Jmp2, N + 1))
        lam_rad = res2['halfwidths'][N]
        lam_low = lam - mpf('%.3e' % lam_rad)
        lam_high = lam + mpf('%.3e' % lam_rad)
        anchor_covered = lam_low <= anchor <= lam_high
        # truncation ladder vs previous N
        trunc = None
        if prev_lam_str:
            trunc = abs(lam - mpf(prev_lam_str))
        log('  K2 OK contained=%s box_r=%s lambda_hw=%.3e GS=%s cond_J2=%.3e'
            % (res2['contained'], res2['box_r_used'], lam_rad, gs2['contained'], condJ2))
        log('  anchor_covered=%s  anchor_diff=%.3e  trunc_prevN=%.3e'
            % (anchor_covered, abs(lam - anchor), (trunc if trunc is not None else float('nan'))))
        rec = {
            'd': d, 'N': N,
            'alpha_mid': nstr(al_cur, 70), 'alpha_rad': '%.3e' % res1['halfwidths'][N],
            'k1_box_r': res1['box_r_used'], 'k1_gs': gs1['contained'],
            'lambda2_mid': nstr(lam, 90), 'lambda2_rad': '%.3e' % lam_rad,
            'k2_box_r': res2['box_r_used'], 'k2_gs': gs2['contained'],
            'cond_J2': '%.6e' % condJ2,
            'anchor_covered': bool(anchor_covered),
            'anchor_diff': '%.3e' % abs(lam - anchor),
            'trunc_vs_prev': ('%.3e' % trunc) if trunc is not None else None,
            'b_mid': [nstr(v, 70) for v in a_full],
            'c_mid': [nstr(v, 70) for v in c],
        }
        out[str(N)] = rec
        with open('most2_run_d%d.json' % d, 'w', encoding='utf-8') as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
        prev_lam = lam
        prev_lam_str = nstr(lam, 90)
        log('  saved most2_run_d%d.json (%d entries)' % (d, len(out)))
    log('=== run_most2 d=%d COMPLETE ===' % d)


if __name__ == '__main__':
    main()

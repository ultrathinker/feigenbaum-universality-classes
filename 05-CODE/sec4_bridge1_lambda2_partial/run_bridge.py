#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_bridge.py -- production run: bridge lambda2 certificate for one degree d,
at two collocation orders (N, N+8). Writes incremental results + progress log.

Usage: python run_bridge.py --d 4 --N 48 56 --seed <path> --out <out.json> --log <log.txt>
"""
import argparse
import json
import time

from mpmath import mp, mpf, nstr
from flint import ctx, arb

import bridge_lib as bl

T0 = time.time()
LOGFILE = None


def log(msg):
    line = '[%8.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    if LOGFILE:
        with open(LOGFILE, 'a', encoding='utf-8') as f:
            f.write(line + '\n')
            f.flush()


def main():
    global LOGFILE
    ap = argparse.ArgumentParser()
    ap.add_argument('--d', type=int, required=True)
    ap.add_argument('--N', type=int, nargs='+', required=True)
    ap.add_argument('--seed', type=str, required=True)
    ap.add_argument('--out', type=str, required=True)
    ap.add_argument('--log', type=str, required=True)
    ap.add_argument('--dps', type=int, default=110)
    args = ap.parse_args()

    LOGFILE = args.log
    open(LOGFILE, 'w', encoding='utf-8').close()
    mp.dps = args.dps
    ctx.dps = args.dps
    DPS = args.dps

    d = args.d
    targets = sorted(args.N)
    log('=== bridge cert d=%d targets N=%s dps=%d ===' % (d, targets, DPS))

    seed = json.load(open(args.seed, encoding='utf-8'))
    lam2_seed = mpf(seed['checks']['eigen2_deflated'])
    alpha_key = 'alpha%d' % d
    alpha_lit = mpf(seed[alpha_key])
    log('seed loaded: alpha_lit=%s lambda2_seed=%s' % (nstr(alpha_lit, 25), nstr(lam2_seed, 25)))

    out = {
        'meta': {'d': d, 'targets': targets, 'dps': DPS, 'seed_path': args.seed,
                  'method': 'bridge: collocation fixed point + eigenpair (Phi,lambda) '
                            '+ normalization Phi(1)=1, chained Krawczyk K1(fixed pt)->K2(eigen), '
                            'python-flint arb balls; independent verifier = preconditioned '
                            'interval Gauss-Seidel (Hansen-Sengupta), different code path/criterion',
                  'started': time.strftime('%Y-%m-%d %H:%M:%S')},
        'results': {},
    }

    def flush_out():
        with open(args.out, 'w', encoding='utf-8') as f:
            json.dump(out, f, indent=2, ensure_ascii=False)

    flush_out()

    log('--- fixed point continuation ladder -> targets %s ---' % targets)
    fp_by_N = bl.continue_fixedpoint_mp(d, targets, tol=mpf('1e-90'), log=log)
    log('--- eigenpair continuation ladder (chained on fixed pt) -> targets %s ---' % targets)
    a_by_N = {}
    al_by_N = {}
    a_cur = None
    al_cur = None
    # need a,alpha at EVERY ladder rung too (continue_eigen_mp needs a_by_N dict for
    # ladder AND target Ns); recompute the same ladder path for fixed point
    # (cheap; solve_fixedpoint_mp is fast) so we have a/alpha at every rung.
    ladder_all = sorted(set(bl.DEFAULT_LADDER + targets))
    a_cur = [mpf(1), mpf('-1.6')]
    al_cur = mpf('1.5')
    for N in ladder_all:
        a_cur, al_cur, it, res = bl.solve_fixedpoint_mp(d, N, a_cur, al_cur, tol=mpf('1e-90'), maxit=300)
        a_by_N[N] = list(a_cur)
        al_by_N[N] = al_cur
        log('  FP ladder N=%d it=%d resid=%.3e alpha=%s' % (N, it, float(res), nstr(al_cur, 25)))

    eig_by_N = bl.continue_eigen_mp(d, a_by_N, al_by_N, targets, lam2_seed, tol=mpf('1e-90'), log=log)

    for N in targets:
        t_stage = time.time()
        log('=== CERTIFY N=%d ===' % N)
        a_mid, alpha_mid, fp_it, fp_res = fp_by_N[N]
        c_mid, lam_mid, eig_it, eig_res = eig_by_N[N]
        log('  mid values: alpha=%s lambda2=%s (Newton resid fp=%.2e eig=%.2e)'
            % (nstr(alpha_mid, 25), nstr(lam_mid, 25), float(fp_res), float(eig_res)))
        dense_fp = bl.fp_residual_dense_mp(a_mid, alpha_mid, d, npts=1001)
        dense_eig = bl.eig_residual_dense_mp(a_mid, alpha_mid, c_mid, lam_mid, d, npts=1001)
        log('  dense scalar residual: fixed pt=%.3e  eigenpair=%.3e' % (float(dense_fp), float(dense_eig)))

        log('  K1 (fixed point) Krawczyk, auto box radius...')
        r1 = bl.certify_fixedpoint_auto(d, a_mid, alpha_mid, N, DPS, log=log)
        alpha_cert_mid = r1['xhat_arb'][N]
        alpha_cert_rad = r1['halfwidths'][N]
        log('  K1 OK: box_r=%s alpha certified +-%.4e' % (r1['box_r_used'], alpha_cert_rad))

        aX_box = [arb(1)] + r1['X'][:N]
        alphaX_box = r1['X'][N]

        log('  K2 (eigenpair, chained on K1 box) Krawczyk, auto box radius...')
        r2 = bl.certify_eigen_auto(d, a_mid, alpha_mid, c_mid, lam_mid, N, DPS, aX_box, alphaX_box, log=log)
        lam2_cert_mid = r2['xhat_arb'][N]
        lam2_cert_rad = r2['halfwidths'][N]
        log('  K2 OK: box_r=%s lambda2 certified +-%.4e' % (r2['box_r_used'], lam2_cert_rad))
        log('  lambda2(d=%d,N=%d) = %s +- %.4e' % (d, N, lam2_cert_mid.mid().str(30, radius=False), lam2_cert_rad))

        log('  independent verifier: interval Gauss-Seidel re-check (different code path)...')
        pts1_mp = bl.colloc_points_mp(N + 1)
        pts1_arb = bl.colloc_points_arb(N + 1, DPS)
        Fmp1, Jmp1 = bl.fp_residual_jacobian_mp(a_mid, alpha_mid, pts1_mp, d)

        def F_at_point1(v, N=N, pts1_arb=pts1_arb, d=d):
            aa = [arb(1)] + v[:N]
            al = v[N]
            F, _ = bl.fp_F_J_arb(aa, al, pts1_arb, d)
            return F

        def F_J_on_box1(v, N=N, pts1_arb=pts1_arb, d=d):
            aa = [arb(1)] + v[:N]
            al = v[N]
            return bl.fp_F_J_arb(aa, al, pts1_arb, d)

        xhat1 = list(a_mid[1:]) + [alpha_mid]
        gs1 = bl.gauss_seidel_verify(xhat1, Jmp1, F_at_point1, F_J_on_box1, DPS, r1['box_r_used'], sweeps=3)
        log('  GS independent K1: contained=%s' % gs1['contained'])

        pts2_mp = bl.colloc_points_mp(N)
        pts2_arb = bl.colloc_points_arb(N, DPS)
        Fmp2, Jmp2 = bl.eig_residual_jacobian_mp(a_mid, alpha_mid, c_mid, lam_mid, pts2_mp, d)
        G1_box = bl.G_eval_arb(aX_box, arb(1))

        def F_at_point2(v, aX_box=aX_box, alphaX_box=alphaX_box, N=N, pts2_arb=pts2_arb, d=d, G1_box=G1_box):
            F, _ = bl.eig_F_J_arb(aX_box, alphaX_box, v[:N], v[N], pts2_arb, d, G1_box)
            return F

        def F_J_on_box2(v, aX_box=aX_box, alphaX_box=alphaX_box, N=N, pts2_arb=pts2_arb, d=d, G1_box=G1_box):
            return bl.eig_F_J_arb(aX_box, alphaX_box, v[:N], v[N], pts2_arb, d, G1_box)

        xhat2 = list(c_mid) + [lam_mid]
        gs2 = bl.gauss_seidel_verify(xhat2, Jmp2, F_at_point2, F_J_on_box2, DPS, r2['box_r_used'], sweeps=3)
        log('  GS independent K2: contained=%s' % gs2['contained'])

        stage_elapsed = time.time() - t_stage
        log('  N=%d stage done in %.1fs' % (N, stage_elapsed))

        out['results'][str(N)] = {
            'N': N, 'd': d,
            'alpha_mid': nstr(alpha_mid, 40),
            'alpha_cert_halfwidth': alpha_cert_rad,
            'alpha_box_r_used': r1['box_r_used'],
            'k1_contained': r1['contained'],
            'gs1_contained': gs1['contained'],
            'lambda2_mid': lam2_cert_mid.mid().str(40, radius=False),
            'lambda2_cert_halfwidth': lam2_cert_rad,
            'lambda2_box_r_used': r2['box_r_used'],
            'k2_contained': r2['contained'],
            'gs2_contained': gs2['contained'],
            'newton_resid_fp': float(fp_res),
            'newton_resid_eig': float(eig_res),
            'dense_resid_fp': float(dense_fp),
            'dense_resid_eig': float(dense_eig),
            'anchor_alpha_diff': float(abs(alpha_mid - alpha_lit)),
            'anchor_lambda2_diff': float(abs(lam_mid - lam2_seed)),
            'elapsed_sec': stage_elapsed,
        }
        flush_out()
        log('  wrote incremental results to %s' % args.out)

    # cross-N stability check (two orders N, N+8): compare certified balls
    if len(targets) >= 2:
        N1, N2 = targets[0], targets[-1]
        r_lo = out['results'][str(N1)]
        r_hi = out['results'][str(N2)]
        lam1 = mpf(r_lo['lambda2_mid'])
        lam2 = mpf(r_hi['lambda2_mid'])
        diff = abs(lam1 - lam2)
        combined_rad = r_lo['lambda2_cert_halfwidth'] + r_hi['lambda2_cert_halfwidth']
        overlap = float(diff) <= combined_rad
        log('=== cross-N stability: N=%d vs N=%d ===' % (N1, N2))
        log('  |lambda2(N=%d) - lambda2(N=%d)| = %.4e ; sum of half-widths = %.4e ; overlap=%s'
            % (N1, N2, float(diff), combined_rad, overlap))
        out['cross_N_stability'] = {
            'N_lo': N1, 'N_hi': N2, 'diff': float(diff),
            'combined_halfwidth': combined_rad, 'balls_overlap': overlap,
        }
        flush_out()

    out['meta']['finished'] = time.strftime('%Y-%m-%d %H:%M:%S')
    out['meta']['total_elapsed_sec'] = time.time() - T0
    flush_out()
    log('=== ALL DONE d=%d, total %.1fs ===' % (d, time.time() - T0))


if __name__ == '__main__':
    main()

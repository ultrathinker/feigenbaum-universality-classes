#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""End-to-end small-N test: mpmath solve -> Krawczyk K1 (fixed pt) -> K2 (eigen,
chained) -> Gauss-Seidel independent verify. Validate machinery before scaling."""
import time
from mpmath import mp, mpf, nstr
from flint import ctx, arb

mp.dps = 90
ctx.dps = 90
DPS = 90

import bridge_lib as bl

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


def run(d, N, lam_seed):
    log('=== d=%d N=%d ===' % (d, N))
    a_init = [mpf(1), mpf('-1.6')]
    a, alpha, it, res = bl.solve_fixedpoint_mp(d, N, a_init, mpf('1.5'), tol=mpf('1e-75'), maxit=300)
    log('fixed pt: it=%d resid=%.3e alpha=%s' % (it, float(res), nstr(alpha, 25)))
    c, lam, it2, res2 = bl.solve_eigen_mp(d, a, alpha, N, lam_init=mpf(lam_seed), tol=mpf('1e-75'), maxit=300)
    log('eigen: it=%d resid=%.3e lambda2=%s' % (it2, float(res2), nstr(lam, 25)))

    r1 = bl.certify_fixedpoint_auto(d, a, alpha, N, DPS, log=log)
    log('K1 OK box_r=%s. alpha half-width=%.3e' % (r1['box_r_used'], r1['halfwidths'][N]))

    aX_box = [arb(1)] + r1['X'][:N]
    alphaX_box = r1['X'][N]
    r2 = bl.certify_eigen_auto(d, a, alpha, c, lam, N, DPS, aX_box, alphaX_box, log=log)
    log('K2 OK box_r=%s. lambda2 half-width=%.3e' % (r2['box_r_used'], r2['halfwidths'][N]))
    lam2_mid = r2['xhat_arb'][N]
    lam2_ball = lam2_mid + bl.ball('%.6e' % r2['halfwidths'][N])
    log('lambda2 certified ball mid=%s halfwidth=%.3e' % (lam2_mid.mid().str(20), r2['halfwidths'][N]))

    # independent Gauss-Seidel re-verify of BOTH stages
    pts1_mp = bl.colloc_points_mp(N + 1)
    pts1_arb = bl.colloc_points_arb(N + 1, DPS)
    Fmp1, Jmp1 = bl.fp_residual_jacobian_mp(a, alpha, pts1_mp, d)

    def F_at_point1(v):
        aa = [arb(1)] + v[:N]
        al = v[N]
        F, _ = bl.fp_F_J_arb(aa, al, pts1_arb, d)
        return F

    def F_J_on_box1(v):
        aa = [arb(1)] + v[:N]
        al = v[N]
        return bl.fp_F_J_arb(aa, al, pts1_arb, d)

    xhat1 = list(a[1:]) + [alpha]
    gs1 = bl.gauss_seidel_verify(xhat1, Jmp1, F_at_point1, F_J_on_box1, DPS, r1['box_r_used'], sweeps=3)
    log('GS independent verify K1: contained=%s' % gs1['contained'])
    assert gs1['contained'], 'GS K1 FAILED'

    pts2_mp = bl.colloc_points_mp(N)
    pts2_arb = bl.colloc_points_arb(N, DPS)
    Fmp2, Jmp2 = bl.eig_residual_jacobian_mp(a, alpha, c, lam, pts2_mp, d)
    G1_box = bl.G_eval_arb(aX_box, arb(1))

    def F_at_point2(v):
        F, _ = bl.eig_F_J_arb(aX_box, alphaX_box, v[:N], v[N], pts2_arb, d, G1_box)
        return F

    def F_J_on_box2(v):
        return bl.eig_F_J_arb(aX_box, alphaX_box, v[:N], v[N], pts2_arb, d, G1_box)

    xhat2 = list(c) + [lam]
    gs2 = bl.gauss_seidel_verify(xhat2, Jmp2, F_at_point2, F_J_on_box2, DPS, r2['box_r_used'], sweeps=3)
    log('GS independent verify K2: contained=%s' % gs2['contained'])
    assert gs2['contained'], 'GS K2 FAILED'

    log('=== d=%d N=%d ALL OK ===\n' % (d, N))


run(4, 12, '0.29183840898354510453')
run(6, 12, '0.3710714676819301324709834')

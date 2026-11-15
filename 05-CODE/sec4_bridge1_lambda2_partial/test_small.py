#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quick sanity check of bridge_lib mpmath solvers at small N, both degrees,
against the known float-pilot values (alpha4, alpha6, eigen2_deflated)."""
import json
from mpmath import mp, mpf, nstr

mp.dps = 60
import bridge_lib as bl


def check_degree(d, seed_path, N):
    seed = json.load(open(seed_path, encoding='utf-8'))
    alpha_key = 'alpha%d' % d
    alpha_lit = mpf(seed[alpha_key])
    print('--- d=%d  N=%d ---' % (d, N))
    a_init = [mpf(1), mpf('-1.6')]
    alpha_init = mpf('1.5')
    a, alpha, it, res = bl.solve_fixedpoint_mp(d, N, a_init=a_init, alpha_init=alpha_init,
                                                tol=mpf('1e-55'), maxit=300)
    print('  fixed point converged it=%d resid=%.3e' % (it, float(res)))
    print('  alpha  = %s' % nstr(alpha, 30))
    print('  alpha* = %s (literature/seed, full N=%s)' % (nstr(alpha_lit, 30), seed['N']))
    print('  |alpha - alpha_lit| = %.3e  (expect small-N truncation error, not machine-eps)'
          % float(abs(alpha - alpha_lit)))
    dense = bl.fp_residual_dense_mp(a, alpha, d, npts=1001)
    print('  dense residual = %.3e' % float(dense))

    lam_seed = mpf(seed['checks']['eigen2_deflated'])
    c, lam, it2, res2 = bl.solve_eigen_mp(d, a, alpha, N, lam_init=lam_seed, tol=mpf('1e-55'), maxit=300)
    print('  eigen converged it=%d resid=%.3e' % (it2, float(res2)))
    print('  lambda2       = %s' % nstr(lam, 30))
    print('  lambda2(seed) = %s' % nstr(lam_seed, 30))
    print('  |lambda2 - lambda2(seed)| = %.3e' % float(abs(lam - lam_seed)))
    dres = bl.eig_residual_dense_mp(a, alpha, c, lam, d, npts=1001)
    print('  eigen dense residual = %.3e' % float(dres))
    print()


check_degree(4, r'../sec3_pilot_d4/quartic_fixedpoint.json', 20)
check_degree(6, r'../sec3_pilot_d6/sextic_fixedpoint.json', 20)

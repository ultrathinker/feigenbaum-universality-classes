#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""smoke.py -- quick validation of K1+K2 arb certification on d=4, N=16."""
import time
from mpmath import mpf, nstr
from most2_lib import (solve_fixedpoint_mp, solve_eigen_mp, certify_fixedpoint,
                       certify_eigen, run_certification)

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


if __name__ == '__main__':
    d = 4
    anchor = mpf('0.29183840898354510453')
    x, al, it, res = solve_fixedpoint_mp(d, 16)
    log('FP N=16 alpha=%s res=%.2e' % (nstr(al, 30), float(res)))
    a_full = [mpf(1) - sum(((-1) ** k) * x[k - 1] for k in range(1, 17))] + list(x)
    c, lam, ite, rese = solve_eigen_mp(d, a_full, al, 16, anchor)
    log('EIG N=16 lambda=%s res=%.2e' % (nstr(lam, 40), float(rese)))
    log('=== K1 ===')
    res1 = certify_fixedpoint(d, 16, a_full, al, 130)
    log('K1 contained=%s box_r=%s alpha_hw=%.3e' % (res1['contained'], res1['box_r_used'], res1['halfwidths'][16]))
    aX_box = [__import__('most2_lib').b0_of(res1['X'][:16], 16)] + res1['X'][:16]
    alphaX_box = res1['X'][16]
    log('=== K2 ===')
    res2 = certify_eigen(d, 16, a_full, al, c, lam, 130, aX_box, alphaX_box)
    log('K2 contained=%s box_r=%s lambda_hw=%.3e' % (res2['contained'], res2['box_r_used'], res2['halfwidths'][16]))
    log('DONE')

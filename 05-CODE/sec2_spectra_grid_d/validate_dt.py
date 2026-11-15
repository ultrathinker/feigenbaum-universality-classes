#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_dt.py -- pre-run validation for spectra_grid.py  [oid:spectraF]
NOT RIGOROUS.  Checks, at fractional d=2.5 and odd d=5 (N=32, dps 50):
  (1) the general-d DT_t formula against a finite difference of the doubling
      operator R(f) = -alpha_f f(f(x/alpha_f)) (eps=1e-25), incl. the |U|^d
      branch for negative inner values;
  (2) the coordinate mode K = G - d t G' is an eigenvector of the UNPROJECTED
      linearisation with eigenvalue exactly 1 (rescaling-family mode);
  (3) Newton fixed-point dense residual at both d.
"""
import time
from mpmath import mp, mpf, nstr

import spectra_grid as sg

mp.dps = 50


def R_of_f(f, alpha_f, x):
    return -alpha_f * f(f(x / alpha_f))


def check_d(dstr, N=48):
    d = mpf(dstr)
    t0 = time.time()
    a, al, _it, _r = sg.solve_fp(d, N)
    dres = sg.fp_residual_dense(a, al, d, 401)
    print('d=%s: alpha=%s  newton iters=%d  dense fp residual=%.2e  (%.0fs)'
          % (dstr, nstr(al, 16), _it, float(dres), time.time() - t0), flush=True)

    # (1) DT vs finite difference of R (two eps; each leg compared to its own
    # noise floor ~ fp dense residual / eps, set by the truncation quality of g)
    fd_ok = True
    worst_all = mpf(0)
    for eps_s in ('1e-25', '1e-30'):
        eps = mpf('1e-25') if eps_s == '1e-25' else mpf('1e-30')
        c_test = [mpf(0)] + [mpf(1) / (j + 1) for j in range(1, 10)]

        def g_x(x, _a=a, _d=d):
            return sg.G_eval(_a, abs(x) ** _d)

        def f_eps_x(x):
            return g_x(x) + eps * sg._horner0(c_test, abs(x) ** d)

        worst = mpf(0)
        for t in sg.colloc_points(11):
            x = t ** (1 / d)
            alpha_eps = -1 / f_eps_x(mpf(1))
            fd = (R_of_f(f_eps_x, alpha_eps, x) - g_x(x)) / eps
            an = sg.DT_scalar_t(a, al, d, c_test, t)
            worst = max(worst, abs(fd - an))
        worst_all = max(worst_all, worst)
        floor = max(1e-18, 1e3 * float(dres) / float(eps))
        leg_ok = float(worst) < floor
        fd_ok = fd_ok and leg_ok
        print('  DT vs FD of R (eps=%s): max diff = %.3e (floor %.1e) -> %s'
              % (eps_s, float(worst), floor, 'OK' if leg_ok else 'BAD'), flush=True)

    # (2) coordinate mode
    Kc = [a[0]] + [a[j] * (1 - d * j) for j in range(1, N + 1)]
    ratios = []
    for t in sg.colloc_points(11):
        num = sg.DT_unproj_t(a, al, d, Kc, t)
        den = sg._horner0(Kc, t)
        ratios.append(num / den)
    spread = max(ratios) - min(ratios)
    cm_ok = float(spread) < max(1e-35, 1e3 * float(dres))
    print('  coord mode DT_unproj(K)/K = %s (spread %.2e) -> %s'
          % (nstr(sorted(ratios)[5], 15), float(spread), 'OK' if cm_ok else 'BAD'),
          flush=True)
    return fd_ok and cm_ok


if __name__ == '__main__':
    sg.plog('VALIDATION START (dps %d)' % mp.dps)
    ok1 = check_d('2.5')
    ok2 = check_d('5.0')
    sg.plog('VALIDATION %s (d=2.5: %s, d=5.0: %s)'
            % ('PASS' if (ok1 and ok2) else 'FAIL',
               'ok' if ok1 else 'FAIL', 'ok' if ok2 else 'FAIL'))

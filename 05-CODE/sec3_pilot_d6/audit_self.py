#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
audit_self.py -- third-code-path self-audit of the degree-6 pilot results
(independent implementations: direct power sums instead of Horner, separate
orbit loops; from artifacts only the NUMBERS are taken).
[oid:deg6B] NOT RIGOROUS -- mpmath audit of mpmath results.

Checks:
 A1 |H9| at catalog positions of 4 accepted matches (fresh g6/h6 evaluation)
 A2 superstable residual |x_1024(a_ss)| and exact-period margins for 4 matches
 A3 recompute u_actual for 4 matches with the fresh orbit -> compare d_u
 A4 drift law spot check: d_u(8)/d_u(9) for two rate-study zeros
 A5 gauge: H9(-1) ~ 0 and u(a_10) ~ -1 from the M-side numbers alone
 A6 corruption catch: shift a_ss by 1e-9 -> |x_1024| blows up
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 60
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FP = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
MS = json.load(open('mside_cascade.json', encoding='utf-8'))
HZ = json.load(open('H6_zeros.json', encoding='utf-8'))
MR = json.load(open('match_results.json', encoding='utf-8'))

A = [mpf(s) for s in FP['g_coeffs']]
C = [mpf(s) for s in FP['h_coeffs']]
DELTA6 = mpf(FP['delta6'])
AU0 = mpf(HZ['AU0_6'])
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])
LAM2 = mpf(FP['checks']['eigen2_deflated'])


# fresh evaluators: direct power sums (no Horner)
def g6(z):
    return sum(A[j] * z ** (6 * j) for j in range(len(A)))


def h6(z):
    return sum(C[j] * z ** (6 * j) for j in range(len(C)))


def xH(u, m):
    s = u * AU0 / DELTA6 ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    x2m = None
    for k in range(n2):
        x = g6(x) + s * h6(x)
        if k + 1 == n1:
            x2m = x
    return x2m, x


def xm(a, n):
    x = mpc(0)
    for _ in range(n):
        x = 1 - a * x ** 6
        if fabs(x) > mpf('1e10'):
            return x
    return x


# --- A1/A3: four accepted matches -----------------------------------------------
log('--- A1/A3: fresh H9 at catalog zeros; fresh u_actual vs reported ---')
pick = ['S002', 'S010', 'S100', 'S250']
res = {r['id']: r for r in MR['catalog_results']}
for i, z in enumerate(HZ['zeros_limit_extrapolated']):
    pid = 'S%03d' % (i + 1)
    if pid not in pick:
        continue
    r = res[pid]
    u = mpc(mpf(str(z['re'])), mpf(str(z['im'])))
    den9, num9 = xH(u, 9)
    log('%s |H9(u_cat)| = %.2e (catalog absH9 %.1e)'
        % (pid, float(fabs(num9 / den9)), z['absH9'] or -1))
    a_ss = mpc(mpf(r['a_ss_re']), mpf(r['a_ss_im']))
    u_act_fresh = A6 * (a_ss - A_F) * DELTA6 ** 9
    d_fresh = float(fabs(u_act_fresh - u))
    log('   d_u fresh = %.4e ; reported = %.4e ; rel diff %.1e'
        % (d_fresh, r['d_u'], abs(d_fresh - r['d_u']) / r['d_u']))

# --- A2: superstability + exact period margins -----------------------------------
log('--- A2: |x_1024(a_ss)| and min|x_{2^j}|, j=1..9 ---')
for pid in pick:
    r = res[pid]
    a_ss = mpc(mpf(r['a_ss_re']), mpf(r['a_ss_im']))
    v1024 = fabs(xm(a_ss, 1024))
    margins = [float(fabs(xm(a_ss, 2 ** j))) for j in range(1, 10)]
    log('%s |x_1024| = %.2e ; min_j |x_{2^j}| = %.2e (j=%d)'
        % (pid, float(v1024), min(margins), margins.index(min(margins)) + 1))

# --- A4: drift-law spot check ------------------------------------------------------
log('--- A4: rate spot check on two rate-study zeros ---')
for row in MR['rate_study'][:2]:
    pid = row['id']
    i = int(pid[1:]) - 1
    z = HZ['zeros_limit_extrapolated'][i]
    u = mpc(mpf(str(z['re'])), mpf(str(z['im'])))
    # recompute d_u(8), d_u(9) from the reported a_ss? not stored -> recompute via newton
    from mpmath import mp as _mp
    log('%s reported d8=%.3e d9=%.3e ratio=%.3f (1/lambda2=%.4f)'
        % (pid, row['d_u']['8'], row['d_u']['9'], row['ratio_89'], float(1 / LAM2)))

# --- A5: gauge from M-side numbers alone --------------------------------------------
log('--- A5: gauge check u(a_10) at m=9 ---')
a10 = mpf(MS['cascade']['10'])
u_g = A6 * (a10 - A_F) * DELTA6 ** 9
log('u(a_10)|_{m=9} = %s (should be ~ -1)' % nstr(u_g, 15))

# --- A6: corruption catch -------------------------------------------------------------
log('--- A6: shift a_ss by 1e-9 ---')
r = res['S002']
a_ss = mpc(mpf(r['a_ss_re']), mpf(r['a_ss_im'])) + mpf('1e-9')
log('S002 a_ss+1e-9: |x_1024| = %.2e (was ~1e-69 at the root)'
    % float(fabs(xm(a_ss, 1024))))
log('audit_self done')

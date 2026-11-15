#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_orbit.py -- measure sextic orbit amplitudes (max|x|^6) for the universal
perturbed family and the usable t-range of the N=48 vs N=56 series; compare with
the quartic (maxt4 stats of its zero catalog).
[oid:deg6B] NOT RIGOROUS -- mpmath diagnostic.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 50
T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FP6 = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A48 = [mpf(s) for s in FP6['g_coeffs']]
C48 = [mpf(s) for s in FP6['h_coeffs']]
AL6 = mpf(FP6['alpha6'])
D6 = mpf(FP6['delta6'])
AU0 = mpf('-0.667969539632466139489385700663')  # from p2 calibration log

FP4 = json.load(open('../sec3_pilot_d4/quartic_fixedpoint.json', encoding='utf-8'))
A4 = [mpf(s) for s in FP4['g_coeffs']]
AL4 = mpf(FP4['alpha4'])
D4 = mpf(FP4['delta4'])
HZ4 = json.load(open('../sec3_pilot_d4/H4_zeros.json', encoding='utf-8'))
AU04 = mpf(HZ4['AU0_4'])


def G(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def g6_48(x):
    return G(A48, x ** 6)


def h6_48(x):
    return G(C48, x ** 6)


def g4(x):
    return G(A4, x ** 4)


def orbit_max_t(f, power, n, s=None, hfun=None, cap=6.0):
    x = mpc(0)
    mt = mpf(0)
    for _ in range(n):
        if s is None:
            x = f(x)
        else:
            x = f(x) + s * hfun(x)
        t = abs(x) ** power
        if t > mt:
            mt = t
        if t > cap:
            return mt, False
    return mt, True


log('--- sextic g6 (s=0) orbit amplitude ---')
mt, ok = orbit_max_t(g6_48, 6, 1024)
log('max|x|^6 over 1024 steps (s=0): %.6f (bounded=%s)' % (float(mt), ok))
for u in (-1, 1, 2, 4, 8, 12, 20, 1 + 4j, 3 + 8j, 5 + 15j, 8 + 8j):
    for m in (7, 8):
        s = u * AU0 / D6 ** m
        mt, ok = orbit_max_t(g6_48, 6, 2 ** (m + 1), s=s, hfun=h6_48)
        log('u=%8.3f m=%d: max|x|^6 = %.4f (bounded=%s)'
            % (float(u.real if hasattr(u, 'real') else u), m, float(mt), ok))

log('--- quartic comparison ---')
mt4, ok4 = orbit_max_t(g4, 4, 1024)
log('quartic max|x|^4 (s=0): %.6f' % float(mt4))
mx4 = [z['maxt4'] for z in HZ4['zeros_level8']]
mx4s = sorted(mx4)
log('quartic zero catalog maxt4: min %.3f median %.3f max %.3f'
    % (mx4s[0], mx4s[len(mx4s) // 2], mx4s[-1]))

log('--- sextic series accuracy |G48 - G56|(t) and |h48..| tail proxy ---')
A56 = A48  # placeholder: N=56 coefficients are NOT in the JSON; use ladder tail
# instead estimate tail by |G48(t) - G40(t)|
log('t: |G48-G40|  (proxy for truncation error at N=48)')
for t in (mpf('0.8'), mpf('1.0'), mpf('1.1'), mpf('1.2'), mpf('1.3'),
          mpf('1.4'), mpf('1.5'), mpf('1.6'), mpf('1.8'), mpf('2.0')):
    d = fabs(G(A48, t) - G(A48[:41], t))
    dh = fabs(G(C48, t) - G(C48[:41], t))
    log('t=%.2f  |dG|=%.3e  |dh|=%.3e' % (float(t), float(d), float(dh)))

log('--- how large must N be? tail estimate (a_j ~ rho^-j, rho=1.677) ---')
rho = mpf('1.677')
for tmax, tol in ((mpf('1.2'), mpf('1e-11')), (mpf('1.3'), mpf('1e-11')),
                  (mpf('1.4'), mpf('1e-11')), (mpf('1.25'), mpf('1e-11'))):
    r = tmax / rho
    import mpmath
    n_need = int(mpmath.ceil(mpmath.log(tol * (1 - r)) / mpmath.log(r))) + 1
    log('tmax=%.2f tol=1e-11: need N ~ %d' % (float(tmax), n_need))
log('diag_orbit done')

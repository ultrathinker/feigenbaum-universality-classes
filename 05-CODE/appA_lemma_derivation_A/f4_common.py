# -*- coding: utf-8 -*-
"""
f4_common.py -- общие функции F4: загрузка неподвижных точек, оператор DT
(точная формула из spectra_grid.py), кандидат башенной моды.

[oid:f4GLM] NOT RIGOROUS -- mpmath float exploration, no interval certificates.
"""
import json
import os
import time

from mpmath import mp, mpf, cos, pi

ROOT = r'D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип'
HERE = os.path.dirname(os.path.abspath(__file__))
PROG = os.path.join(HERE, 'progress.log')


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


# ---------------------------------------------------------------- series ----
def horner(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def load_fp(tag):
    """-> (a, alpha, d, meta). tag in {'2.0','4.0','6.0'}"""
    if tag == '2.0':
        j = json.load(open(os.path.join(ROOT, 'wave1_spectra_d', 'spectra_grid.json'),
                           encoding='utf-8'))
        rec = next(r for r in j['points'] if r['d'] == '2.0')
        a = [mpf(s) for s in rec['_a48']]
        return a, mpf(rec['alpha48']), mpf(2), {'src': 'spectra_grid _a48 N=48',
                                                'fp_dense_resid': rec['fp_dense_resid_N48']}
    if tag == '4.0':
        j = json.load(open(os.path.join(ROOT, 'quartic_pilot', 'quartic_fixedpoint.json'),
                           encoding='utf-8'))
        return ([mpf(s) for s in j['g_coeffs']], mpf(j['alpha4']), mpf(4),
                {'src': 'quartic_fixedpoint N=%s dps=%s' % (j['N'], j['meta']['mp_dps'])})
    if tag == '6.0':
        j = json.load(open(os.path.join(ROOT, 'agentB_deg6', 'sextic_fixedpoint.json'),
                           encoding='utf-8'))
        return ([mpf(s) for s in j['g_coeffs']], mpf(j['alpha6']), mpf(6),
                {'src': 'sextic_fixedpoint N=%s dps=%s' % (j['N'], j['meta']['mp_dps'])})
    raise KeyError(tag)


# ---------------------------------------------------------------- DT --------
def DT_t(a, alpha, d, Phi, t):
    """(DT Phi)(t) -- точная формула из spectra_grid.DT_scalar_t, Phi callable."""
    ad = alpha ** d
    G1 = horner(a, mpf(1))
    ta = t / ad
    U = horner(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    Kt = horner(a, t) - d * t * Gp(a, t)
    return (-(Phi(mpf(1)) / G1) * Kt - alpha * Phi(W)
            - alpha * d * sU * aU ** (d - 1) * Gp(a, W) * Phi(ta))


def DT_unproj_t(a, alpha, d, Phi, t):
    """Bare linearisation (fixed alpha_f) -- coordinate-mode diagnostic."""
    ad = alpha ** d
    ta = t / ad
    U = horner(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    return (-alpha * Phi(W)
            - alpha * d * sU * aU ** (d - 1) * Gp(a, W) * Phi(ta))


# ------------------------------------------------------- tower candidates --
# C1: касание чистого сопряжения аргумента  g(psi_eps(x)), psi=x(1+eps x^{dk})
#     -> Phi = d t^{k+1} G'(t)
# C2: сопряжение + нормировка f(0)=1:
#     f_eps = psi_eps(g(psi_eps^{-1}(x)))/psi_eps(1) -> Phi = G^{dk+1} - d t^{k+1} G' - G
def cand_C1(a, d, k):
    return lambda t: d * t ** (k + 1) * Gp(a, t)


def cand_C2(a, d, k):
    return lambda t: (horner(a, t) ** (d * k + 1) - d * t ** (k + 1) * Gp(a, t)
                      - horner(a, t))


def K_fun(a, d):
    return lambda t: horner(a, t) - d * t * Gp(a, t)


def grid_t(npts, tmax):
    # плотнее к 0 и к 1: чебышёвские узлы, как в кодах
    return [tmax * (1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def resid_table(a, alpha, d, Phi, lam, npts=200, tmax=mpf('1.2')):
    """max|DT Phi - lam Phi| / max|Phi| на сетке."""
    rmax = mpf(0)
    pmax = mpf(0)
    worst_t = None
    for t in grid_t(npts, tmax):
        r = abs(DT_t(a, alpha, d, Phi, t) - lam * Phi(t))
        p = abs(Phi(t))
        if r > rmax:
            rmax, worst_t = r, t
        pmax = max(pmax, p)
    return rmax / pmax, worst_t, pmax

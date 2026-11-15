#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_spectrum.py -- diagnostic: float spectrum of the DT collocation matrix
for the sextic class across N, to locate the physical subleading eigenvalues
(lambda2, lambda3 < 1) versus pseudo-spectral pollution.
[oid:deg6B] NOT RIGOROUS -- float diagnostic only.
"""
import json
import time

import numpy as np

T0 = time.time()


def log(msg):
    print("[%6.1fs] %s" % (time.time() - T0, msg), flush=True)


FPJ = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
A_CO = np.array([float(s) for s in FPJ['g_coeffs']])
ALPHA = float(FPJ['alpha6'])


def cheb_nodes(npts):
    i = np.arange(npts) + 0.5
    return (1.0 - np.cos(np.pi * i / npts)) / 2.0


def G_eval(a, t):
    r = np.full_like(t, a[-1])
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    n = len(a) - 1
    s = np.zeros_like(t)
    for j in range(1, n + 1):
        s += j * a[j] * t ** (j - 1)
    return s


def DT_of_monomial(a, alpha, j, x):
    """DT(phi)(x) for phi(x) = x^{6j}, scalar-exact, float, vectorised."""
    t = x ** 6
    G1 = G_eval(a, np.array([1.0]))[0]
    ta = t / alpha ** 6
    U = G_eval(a, ta)
    W = U ** 6
    gx = G_eval(a, t)
    gpx = 6.0 * x ** 5 * Gp_eval(a, t)
    u = U
    gpu = 6.0 * u ** 5 * Gp_eval(a, W)
    phi1 = 1.0
    phix = x ** (6 * j)
    phiu = u ** (6 * j)
    phita = ta ** j  # phi(x/alpha) = (x/alpha)^{6j} = ta^j
    return (-(phi1 / G1) * (gx - x * gpx)
            - alpha * phiu - alpha * gpu * phita)


def build_A(N):
    pts = cheb_nodes(N)
    x = pts ** (1.0 / 6.0)
    V = np.empty((N, N))
    F = np.empty((N, N))
    for i in range(N):
        for j in range(1, N + 1):
            V[i, j - 1] = pts[i] ** j
            F[i, j - 1] = DT_of_monomial(A_CO, ALPHA, j, np.array([x[i]]))[0]
    return np.linalg.solve(V, F)


log('building float DT collocation matrices for N = 16,24,32,40,48,56')
spectra = {}
for N in (16, 24, 32, 40, 48, 56):
    A = build_A(N)
    ev = np.linalg.eigvals(A)
    ev = ev[np.argsort(-np.abs(ev))]
    spectra[N] = ev
    log('N=%d: top |eig|: %s' % (N, ' '.join('%.4f%+.4fi|%.3f' % (e.real, e.imag, abs(e)) for e in ev[:8])))

# look for N-stable eigenvalues: track eigenvalues of N=56 present in all ladders
ref = spectra[56]
log('--- stability tracking (nearest eigenvalue across N to each N=56 eigenvalue with |eig|<1.5) ---')
for e in ref:
    if abs(e) > 1.5 or abs(e) < 1e-3:
        continue
    track = []
    for N in (16, 24, 32, 40, 48):
        d = spectra[N] - e
        k = np.argmin(np.abs(d))
        track.append(abs(d[k]))
    log('cand %+.8f%+.8fi (|.|=%.6f): drifts N16..N48 = %s'
        % (e.real, e.imag, abs(e), ' '.join('%.1e' % t for t in track)))

log('--- full list |eig|<1.5 at N=56 ---')
for e in ref:
    if abs(e) < 1.5:
        log('   %+.10f %+.10f i   |.|=%.10f' % (e.real, e.imag, abs(e)))
log('done')

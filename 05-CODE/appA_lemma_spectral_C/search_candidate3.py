#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_candidate3.py -- [oid:F4-M3]

Look at the actual numerical eigenvectors and try to fit a known form.
"""
import json
import os
from mpmath import mp, mpf, pi, cos, fsum

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))


def G_eval(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def colloc_points(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def DT_scalar_t(a, alpha, d, c, t):
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    ta = t / ad
    U = G_eval(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    Kt = G_eval(a, t) - d * t * Gp_eval(a, t)

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    Phi1 = sum(c)
    return (-(Phi1 / G1) * Kt - alpha * horner0(c, W)
            - alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W) * horner0(c, ta))


def main():
    # Load d=4 fixed point
    with open(os.path.join(HERE, '..', '..', 'quartic_pilot', 'quartic_fixedpoint.json'), encoding='utf-8') as f:
        qd = json.load(f)
    a_g = [mpf(s) for s in qd['g_coeffs']]
    # Convert g_coeffs (in x) to G_coeffs (in t=|x|^4) by taking every 4th
    d = 4
    N_g = len(a_g)
    N_G = N_g // d
    a = [mpf(0)] * (N_G + 1)
    for j in range(N_g):
        m, r = divmod(j, d)
        if r == 0 and m <= N_G:
            a[m] = a_g[j]
    alpha = mpf(qd['alpha4'])

    print(f'alpha_4 = {float(alpha):.10f}')
    print(f'alpha^-4 = {float(alpha**-4):.10f}')
    print(f'N_G = {N_G}')
    print(f'G coefficients (first 10): {[float(a[i]) for i in range(10)]}')

    # Get numerical rvec for tower1 at d=4 from wave2
    with open(os.path.join(HERE, '..', '..', 'wave2_vector_homotopy', 'spectrum_window.json'), encoding='utf-8') as f:
        sw = json.load(f)
    rec = sw['points']['4.0']
    lams = rec['lams']
    rvecs = rec['rvecs']
    # alpha^-4 is at index 3
    tv = [mpf(s) for s in rvecs[3]]
    print(f'tower1 rvec (lam={float(lams[3])}), first 10:')
    print(f'  {[float(tv[i]) for i in range(10)]}')
    print(f'normalized to tv[1]=1:')
    rn = [tv[i]/tv[1] for i in range(len(tv))]
    print(f'  {[float(rn[i]) for i in range(10)]}')

    # Now try various candidates at d=4 for eigenvalue alpha^{-4}
    N = 48

    def series_diff_test(candidate, name):
        c = [mpf(0)] * N
        for j in range(1, N):
            jj = min(j, len(a) - 1)
            c[j] = candidate(jj)
        target = alpha ** -4
        tpts = colloc_points(201)
        maxr = mpf(0); maxp = mpf(0)
        for t in tpts:
            lhs = DT_scalar_t(a, alpha, d, c, t)
            v = c[0]
            for k in range(1, len(c)):
                v = v * t + c[k]
            rhs = target * v
            ar = abs(lhs - rhs)
            if ar > maxr: maxr = ar
            if abs(v) > maxp: maxp = abs(v)
        rel = maxr / maxp if maxp > 0 else mpf(0)
        # Cosine with tv
        L = min(len(c), len(tv))
        n1 = fsum(x*x for x in c[:L]) ** mpf('0.5')
        n2 = fsum(x*x for x in tv[:L]) ** mpf('0.5')
        cos = abs(fsum(c[i]*tv[i] for i in range(L)) / (n1*n2)) if n1 > 0 and n2 > 0 else mpf(0)
        print(f'  {name:40s} | rel={float(rel):.4e} | cos={float(cos):.4e}')

    print('\n=== Candidates at d=4, target=alpha^-4 ===')
    # d·t·G'(t)
    series_diff_test(lambda j: d * j * a[j], 'd*t*Gprime(t) [dilation]')
    # d·t²·G'(t)
    series_diff_test(lambda j: d * (j - 1) * a[j - 1] if j - 1 >= 1 else mpf(0), 'd*t^2*Gprime(t) [Fable]')
    # d·t^k·G'(t) for various k
    for kk in [3, 4, 5]:
        series_diff_test(lambda j, k=kk: d * (j - k + 1) * a[j - k + 1] if j - k + 1 >= 1 else mpf(0), f'd*t^{kk}*Gprime(t)')
    # d·t·G(t)
    series_diff_test(lambda j: d * (a[j - 1] if j - 1 >= 0 else mpf(1)), 'd*t*G(t)')
    # d·t²·G(t)
    series_diff_test(lambda j: d * (a[j - 2] if j - 2 >= 0 else mpf(0)), 'd*t^2*G(t)')
    # (G-1)
    series_diff_test(lambda j: a[j], '(G-1)')
    # t·(G-1)
    series_diff_test(lambda j: a[j - 1] if j - 1 >= 0 else mpf(0), 't*(G-1)')
    # K(t) = G - d·t·G'
    series_diff_test(lambda j: a[j] - d * j * a[j], 'K(t) gauge')
    # Combination
    series_diff_test(lambda j: a[j] + a[j - 1] if j - 1 >= 0 else a[j], 'G(t) + t·G(t)')


if __name__ == '__main__':
    main()

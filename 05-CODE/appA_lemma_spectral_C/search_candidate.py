#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_candidate.py -- [oid:F4-M3]

Systematic numerical search for the right Phi_k(t) formula.
For d=2, target eigenvalue alpha^{-2}, look at the actual rvec and try
many candidate forms. We minimize |residual| / |Phi| over (a, b, c, ...)
for candidates of the form Phi(t) = a * A(t) + b * B(t) + ... where
A, B, ... are basis functions derived from G, G', t, alpha, d.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, pi, cos, matrix, lu_solve

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
    """Scalar-exact (DT Phi)(t)."""
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


def dense_residual_target(a, alpha, d, c, target, npts=201):
    tpts = colloc_points(npts)

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    maxr = mpf(0)
    maxp = mpf(0)
    for t in tpts:
        lhs = DT_scalar_t(a, alpha, d, c, t)
        rhs = target * horner0(c, t)
        ar = abs(lhs - rhs)
        if ar > maxr:
            maxr = ar
        ph = abs(horner0(c, t))
        if ph > maxp:
            maxp = ph
    rel = maxr / maxp if maxp > 0 else mpf(0)
    return maxr, maxp, rel


def series_diff(a, b):
    """Element-wise difference of two coefficient lists."""
    n = max(len(a), len(b))
    out = [mpf(0)] * n
    for i in range(len(a)):
        out[i] = a[i]
    for i in range(len(b)):
        out[i] = out[i] - b[i]
    return out


def series_sum(a, b, c_a=mpf(1), c_b=mpf(1)):
    n = max(len(a), len(b))
    out = [mpf(0)] * n
    for i in range(len(a)):
        out[i] += c_a * a[i]
    for i in range(len(b)):
        out[i] += c_b * b[i]
    return out


def main():
    # Load d=2 fixed point
    with open(os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json'), encoding='utf-8') as f:
        sg = json.load(f)
    for r in sg['points']:
        if r['d'] == '2.0':
            a = [mpf(s) for s in r['_a48']]
            alpha = mpf(r['alpha48'])
            d = mpf(2)
            tv = [mpf(s) for s in r['tower1']['N48_vec']]
            target = alpha ** (-2)
            break

    print(f'alpha = {float(alpha):.10f}')
    print(f'target = alpha^-2 = {float(target):.10f}')
    print(f'tower1 rvec (first 8): {[float(tv[i]) for i in range(8)]}')

    # Try many candidates
    N = 48
    candidates = {}

    # 1) d·t·G'(t) (dilation, my first try)
    c = [mpf(0)] * N
    for j in range(1, N):
        if j - 1 < N:
            c[j] = d * j * a[j]
    candidates['d*t*Gprime(t)'] = c

    # 2) d·t²·G'(t) (Fable x^{d+1})
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] = d * j * a[j]
    candidates['d*t^2*Gprime(t) [Fable]'] = c

    # 3) G(t) - 1
    c = [mpf(0)] * N
    for j in range(1, N):
        c[j] = a[j]
    candidates['G(t) - 1'] = c

    # 4) t·G(t) - (G(t) - 1)
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] += a[j]
    for j in range(1, N):
        c[j] -= a[j]
    candidates['t*G(t) - (G(t)-1)'] = c

    # 5) t²·G(t)
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] += a[j - 1] if j - 1 >= 1 else a[j - 1]
    candidates['t^2*G(t)'] = c

    # 6) d*t*(G(t)-1)
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] = d * a[j]
    candidates['d*t*(G(t)-1)'] = c

    # 7) (1-1/alpha)*t*G(t) - (G(t)-1)  -- some combination
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] += (1 - 1 / alpha) * (a[j - 1] if j >= 2 else 1)
    for j in range(1, N):
        c[j] -= a[j]
    candidates['(1-1/alpha)*t*G(t) - (G(t)-1)'] = c

    # 8) d*G(t) - some derivative combination
    # Actually let me check t*G'(t) without d factor
    c = [mpf(0)] * N
    for j in range(1, N):
        c[j] = j * a[j]
    candidates['t*Gprime(t)'] = c

    # 9) G(t) - alpha * d * t * G'(t)
    c = [mpf(0)] * N
    for j in range(1, N):
        c[j] = a[j] - alpha * d * j * a[j]
    candidates['G - alpha*d*t*Gprime'] = c

    # 10) K(t) = G(t) - d*t*G'(t) -- gauge
    c = [mpf(0)] * N
    for j in range(1, N):
        c[j] = a[j] - d * j * a[j]
    candidates['K(t) [gauge]'] = c

    # 11) Try: d*t*G(t) (NOT G')
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        c[j + 1] = d * (a[j - 1] if j - 1 >= 0 else 1)
    candidates['d*t*G(t)'] = c

    # 12) d*t²*G(t)
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        if j >= 1:
            c[j + 1] = d * (a[j - 1] if j - 1 >= 0 else 1)
    # Actually t²·G(t): c[j+2] = a[j], so c[i] = a[i-2]
    c = [mpf(0)] * N
    for j in range(1, N - 1):
        if j + 1 < N:
            c[j + 1] = d * a[j - 1] if j - 1 >= 0 else d * 1
    candidates['d*t^2*G(t)'] = c

    # 13) Combination: G(t) + c*d*t*G'(t) where c is fitted
    for c_factor in [mpf('-2'), mpf('-1'), mpf('1'), mpf('2'), mpf('-0.5'), mpf('0.5')]:
        c = [mpf(0)] * N
        for j in range(1, N):
            c[j] = a[j] + c_factor * d * j * a[j]
        candidates[f'G(t) + {float(c_factor)}*d*t*Gprime(t)'] = c

    # Test all candidates
    print('\n=== Candidate residuals (target = alpha^-2) ===')
    print(f'{"name":40s} | {"max|r|":12s} | {"max|Phi|":12s} | {"rel":12s} | {"cos":10s}')
    for name, c in candidates.items():
        maxr, maxp, rel = dense_residual_target(a, alpha, d, c, target, npts=201)
        # cosine with tv
        L = min(len(c), len(tv))
        if maxp > 0:
            n1 = sum(x * x for x in c[:L]) ** mpf('0.5')
            n2 = sum(x * x for x in tv[:L]) ** mpf('0.5')
            if n1 > 0 and n2 > 0:
                dot = sum(c[i] * tv[i] for i in range(L))
                cos = abs(dot / (n1 * n2))
            else:
                cos = mpf(0)
        else:
            cos = mpf(0)
        print(f'{name:40s} | {float(maxr):.4e} | {float(maxp):.4e} | {float(rel):.4e} | {float(cos):.4e}')


if __name__ == '__main__':
    main()

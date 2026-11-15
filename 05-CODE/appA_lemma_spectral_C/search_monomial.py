#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_monomial.py -- [oid:F4-M3]

Test the monomial candidates Phi(t) = t^k for various k.
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


def load_fp_d2():
    with open(os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json'), encoding='utf-8') as f:
        sg = json.load(f)
    for r in sg['points']:
        if r['d'] == '2.0':
            return [mpf(s) for s in r['_a48']], mpf(r['alpha48'])


def main():
    a, alpha = load_fp_d2()
    d = mpf(2)
    N = 48

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    # Phi(t) = t^k for k = 1, 2, 3, 4
    # Coefficients: c[k] = 1, all others 0
    print(f'd=2, alpha = {float(alpha):.6f}')
    print(f'alpha^-2 = {float(alpha**-2):.6f}')
    print(f'alpha^-4 = {float(alpha**-4):.6f}')
    print(f'alpha^-6 = {float(alpha**-6):.6f}')
    print()

    for k in range(1, 5):
        c = [mpf(0)] * N
        c[k] = mpf(1)
        # Test against various eigenvalues
        print(f'Candidate Phi(t) = t^{k}:')
        for n in range(1, 5):
            target = alpha ** (-d * n)
            tpts = colloc_points(201)
            maxr = mpf(0); maxp = mpf(0)
            for t in tpts:
                lhs = DT_scalar_t(a, alpha, d, c, t)
                Phi = horner0(c, t)
                rhs = target * Phi
                ar = abs(lhs - rhs)
                if ar > maxr: maxr = ar
                if abs(Phi) > maxp: maxp = abs(Phi)
            rel = maxr / maxp if maxp > 0 else mpf(0)
            print(f'  target = alpha^(-{n}d) = {float(target):.4e}: rel = {float(rel):.4e}')
        # Compute "best" eigenvalue via least-squares
        Phi_Phi = mpf(0)
        DTPhi_Phi = mpf(0)
        for t in tpts:
            Phi = horner0(c, t)
            DT = DT_scalar_t(a, alpha, d, c, t)
            Phi_Phi += Phi * Phi
            DTPhi_Phi += DT * Phi
        if Phi_Phi > 0:
            best_lam = DTPhi_Phi / Phi_Phi
            print(f'  Best (LSQ) eigenvalue: {float(best_lam):.4f}')
        print()


if __name__ == '__main__':
    main()

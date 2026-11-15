#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_eigenvalue.py -- [oid:F4-M3]

Compute the actual eigenvalue of the candidate Phi(t) = d*t^k*G'(t) at d=2.
The candidate might be an eigenvector with a specific eigenvalue that we
need to discover.
"""
import json
import os
from mpmath import mp, mpf, pi, cos, fsum, matrix, lu_solve, eigsy

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

    # For Phi(t) = d*t^k*G'(t), find the best eigenvalue.
    # Use least-squares: minimize sum_t ((DT Phi)(t) - lambda*Phi(t))^2
    # Then eigenvalue = (DT Phi, Phi) / (Phi, Phi)

    def make_phi(k, N):
        c = [mpf(0)] * N
        for j in range(1, N):
            jj = min(j, len(a) - 1)
            c[j] = d * (jj - k + 1) * a[jj - k + 1] if (jj - k + 1) >= 1 else mpf(0)
        return c

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    N = 48
    tpts = colloc_points(101)
    print('Best eigenvalue for candidate Phi_k(t) = d*t^k*Gprime(t) at d=2:')
    for k in range(0, 8):
        c = make_phi(k, N)
        # Compute <Phi, Phi> and <DT Phi, Phi>
        Phi_Phi = mpf(0)
        DTPhi_Phi = mpf(0)
        for t in tpts:
            Phi = horner0(c, t)
            DT = DT_scalar_t(a, alpha, d, c, t)
            Phi_Phi += Phi * Phi
            DTPhi_Phi += DT * Phi
        if Phi_Phi > 0:
            best_lam = DTPhi_Phi / Phi_Phi
            print(f'  k={k}: best eigenvalue = {float(best_lam):.6f}')
            # Also compute alpha powers
            print(f'    alpha^-d·{k+1} = alpha^{-2*(k+1)} = {float(alpha**(-2*(k+1))):.6f}')
            print(f'    alpha^-d·{k} = alpha^{-2*k} = {float(alpha**(-2*k)):.6f}')
        else:
            print(f'  k={k}: Phi = 0')

    # Also try the Fable candidate without d factor: t^{k+1}*G'(t)
    print('\nWithout d factor: Phi_k(t) = t^k*Gprime(t) at d=2:')
    for k in range(0, 8):
        c = [mpf(0)] * N
        for j in range(1, N):
            jj = min(j, len(a) - 1)
            c[j] = (jj - k + 1) * a[jj - k + 1] if (jj - k + 1) >= 1 else mpf(0)
        Phi_Phi = mpf(0)
        DTPhi_Phi = mpf(0)
        for t in tpts:
            Phi = horner0(c, t)
            DT = DT_scalar_t(a, alpha, d, c, t)
            Phi_Phi += Phi * Phi
            DTPhi_Phi += DT * Phi
        if Phi_Phi > 0:
            best_lam = DTPhi_Phi / Phi_Phi
            print(f'  k={k}: best eigenvalue = {float(best_lam):.6f}')


if __name__ == '__main__':
    main()

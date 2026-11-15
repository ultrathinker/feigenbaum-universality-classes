#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eig_candidate.py -- [oid:F4-M3]

Compute eigenvectors of DT at d in {2, 4, 6} and look at the structure
of the tower eigenvectors (eigenvalue alpha^{-dk}).
Then try to identify a closed-form candidate.
"""
import json
import os
import sys
from mpmath import mp, mpf, pi, cos, matrix, lu_solve, eig, fsum

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


def build_DT(a, alpha, d, N):
    pts = colloc_points(N)
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    V = matrix(N, N)
    F = matrix(N, N)
    for i, t in enumerate(pts):
        ta = t / ad
        U = G_eval(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        W = aU ** d
        Kt = G_eval(a, t) - d * t * Gp_eval(a, t)
        gauge = -(1 / G1) * Kt
        c3 = -alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W)
        pow_t = mpf(1); pow_W = mpf(1); pow_ta = mpf(1)
        for j in range(1, N + 1):
            pow_t *= t; pow_W *= W; pow_ta *= ta
            V[i, j - 1] = pow_t
            F[i, j - 1] = gauge - alpha * pow_W + c3 * pow_ta
    A = matrix(N, N)
    for jcol in range(N):
        b = matrix(N, 1)
        for i in range(N):
            b[i] = F[i, jcol]
        x = lu_solve(V, b)
        for i in range(N):
            A[i, jcol] = x[i]
    return A


def load_fp(d_target):
    if float(d_target) == 2.0:
        with open(os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json'), encoding='utf-8') as f:
            sg = json.load(f)
        for r in sg['points']:
            if r['d'] == '2.0':
                return [mpf(s) for s in r['_a48']], mpf(r['alpha48'])
    elif float(d_target) == 4.0:
        with open(os.path.join(HERE, '..', '..', 'quartic_pilot', 'quartic_fixedpoint.json'), encoding='utf-8') as f:
            qd = json.load(f)
        return [mpf(s) for s in qd['g_coeffs']], mpf(qd['alpha4'])
    elif float(d_target) == 6.0:
        with open(os.path.join(HERE, '..', '..', 'agentB_deg6', 'sextic_fixedpoint.json'), encoding='utf-8') as f:
            sd = json.load(f)
        return [mpf(s) for s in sd['g_coeffs']], mpf(sd['alpha6'])


def main():
    N = 40

    for d_target in [2, 4]:
        print(f'\n========== d = {d_target} ==========')
        a, alpha = load_fp(d_target)
        d = mpf(d_target)
        A = build_DT(a, alpha, d, N)
        E, EL, ER = eig(A, left=True, right=True)
        # Find tower modes
        print(f'alpha^-d = {float(alpha**-d):.6f}')
        print(f'alpha^-2d = {float(alpha**-2*d):.6f}')
        for k_target in [1, 2, 3]:
            target = float(alpha ** (-d * k_target))
            # Find eigenvalue closest to target
            best_idx = None
            best_diff = None
            for i in range(N):
                lam = float(E[i].real)
                diff = abs(lam - target)
                if best_diff is None or diff < best_diff:
                    best_diff = diff
                    best_idx = i
            if best_idx is not None:
                lam = float(E[best_idx].real)
                print(f'\nTower k={k_target} eigenvalue={lam:.8f} (target={target:.8f}, diff={best_diff:.2e})')
                rvec = [ER[i, best_idx].real for i in range(N)]
                print(f'  rvec (first 12):')
                for j in range(min(12, N)):
                    print(f'    t^{j+1}: {float(rvec[j]):.6f}')
                # Normalize to rvec[1]=1 (t^2 = 1)
                if abs(rvec[1]) > 1e-10:
                    rn = [rvec[i] / rvec[1] for i in range(N)]
                    print(f'  normalized to t^2=1:')
                    for j in range(min(12, N)):
                        print(f'    t^{j+1}: {float(rn[j]):.6f}')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_candidate4.py -- [oid:F4-M3]

Try: x·g'(x) with normalization correction.
Phi(t) = d·t·G'(t) - (d·G'(1)/K(1))·K(t)
where K(t) = G(t) - d·t·G'(t).
This makes Phi(1) = 0 (no gauge component).

Also try: variants with different k (d·t^k·G'(t) with correction).
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


def load_fp(d_target):
    """Load (a, alpha) for given d. NOTE: 'g_coeffs' in the data files are
    actually G(t) coefficients (in t = |x|^d), NOT g(x) coefficients in x.
    The naming is misleading but verified by G(1) = -1/alpha check."""
    if float(d_target) == 2.0:
        with open(os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json'), encoding='utf-8') as f:
            sg = json.load(f)
        for r in sg['points']:
            if r['d'] == '2.0':
                return [mpf(s) for s in r['_a48']], mpf(r['alpha48']), 'wave1'
    elif float(d_target) == 4.0:
        with open(os.path.join(HERE, '..', '..', 'quartic_pilot', 'quartic_fixedpoint.json'), encoding='utf-8') as f:
            qd = json.load(f)
        return [mpf(s) for s in qd['g_coeffs']], mpf(qd['alpha4']), 'quartic'
    elif float(d_target) == 6.0:
        with open(os.path.join(HERE, '..', '..', 'agentB_deg6', 'sextic_fixedpoint.json'), encoding='utf-8') as f:
            sd = json.load(f)
        return [mpf(s) for s in sd['g_coeffs']], mpf(sd['alpha6']), 'sextic'
    raise ValueError(d_target)


def main():
    for d_target in [2, 4, 6]:
        print(f'\n========== d = {d_target} ==========')
        a, alpha, src = load_fp(mpf(d_target))
        d = mpf(d_target)
        # Compute G(1), G'(1), K(1)
        G1 = G_eval(a, mpf(1))
        Gp1 = Gp_eval(a, mpf(1))
        K1 = G1 - d * Gp1
        print(f'alpha = {float(alpha):.10f}')
        print(f'G(1) = {float(G1):.10f} (should be -1/alpha = {float(-1/alpha):.10f})')
        print(f'Gprime(1) = {float(Gp1):.10f}')
        print(f'K(1) = G(1) - d*Gprime(1) = {float(K1):.10f}')
        print(f'alpha^-d = {float(alpha**-d):.10f}')
        print(f'alpha^-2d = {float(alpha**-2*d):.10f}')
        print(f'alpha^-3d = {float(alpha**-3*d):.10f}')

        # Test candidate: Phi_k(t) = d·t^k·G'(t) - (d·G'(1)/K(1))·K(t) [correction]
        # But wait, we want Phi_k(1) = 0. Let me compute:
        # d·t^k·G'(t) at t=1: d·G'(1)
        # K(t) at t=1: K(1)
        # So corrected: d·t^k·G'(t) - (d·G'(1)/K(1))·K(t) -> value at t=1 is d·G'(1) - d·G'(1) = 0. ✓
        #
        # But which k corresponds to eigenvalue alpha^{-dk}?
        # Let me try k=1, 2, 3 and see which gives the right eigenvalue.

        N = 48

        def make_candidate(k):
            """Phi_k(t) = d·t^k·G'(t) - (d·G'(1)/K(1))·K(t)
            where K(t) = G(t) - d·t·G'(t)."""
            c = [mpf(0)] * N
            for j in range(1, N):
                jj = min(j, len(a) - 1)
                c[j] = d * (jj - k + 1) * a[jj - k + 1] if (jj - k + 1) >= 1 else mpf(0)
            # Subtract correction
            correction = d * Gp1 / K1
            for j in range(1, N):
                jj = min(j, len(a) - 1)
                c[j] -= correction * (a[jj] - d * jj * a[jj])
            return c

        def test_residual(c, target):
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
            return maxr, maxp, rel

        # Try k = 1, 2, 3 against eigenvalues alpha^{-d}, alpha^{-2d}, alpha^{-3d}
        print('Candidate Phi_k(t) = d·t^k·Gprime(t) - (d·Gprime(1)/K(1))·K(t):')
        for k in [1, 2, 3]:
            c = make_candidate(k)
            for n_lam in [1, 2, 3]:
                target = alpha ** (-d * n_lam)
                maxr, maxp, rel = test_residual(c, target)
                print(f'  k={k}, target=alpha^(-{n_lam}d)={float(target):.4e}: rel = {float(rel):.4e}')

        # Also try k=1, 2, 3 with NO correction
        print('\nWithout correction:')
        for k in [1, 2, 3]:
            c = [mpf(0)] * N
            for j in range(1, N):
                jj = min(j, len(a) - 1)
                c[j] = d * (jj - k + 1) * a[jj - k + 1] if (jj - k + 1) >= 1 else mpf(0)
            for n_lam in [1, 2, 3]:
                target = alpha ** (-d * n_lam)
                maxr, maxp, rel = test_residual(c, target)
                print(f'  k={k}, target=alpha^(-{n_lam}d)={float(target):.4e}: rel = {float(rel):.4e}')


if __name__ == '__main__':
    main()

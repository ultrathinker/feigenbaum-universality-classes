#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_candidate2.py -- [oid:F4-M3]

Try a fit-based approach. Build a dictionary of basis functions
B_1(t), B_2(t), ... derived from G, G', t, alpha, d. Then for each (d, k),
find a linear combination that minimizes the residual ||DT Phi - alpha^{-dk} Phi||.
If the residual can be made <= 1e-25, the candidate is essentially an
identity, and the formula gives us the eigenvector.
"""
import json
import os
from mpmath import mp, mpf, pi, cos, matrix, lu_solve, fsum

mp.dps = 60

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
    print(f'tower1 rvec (first 10): {[float(tv[i]) for i in range(10)]}')

    # Build basis functions: G(t), t·G(t), t²·G(t), t·G'(t), t²·G'(t), (G-1), t·(G-1), t²·(G-1)
    # Each basis is a t-series coefficient list
    N = 48

    def series_from_formula(formula):
        """formula(j): returns coefficient of t^j in the t-series. Returns c[0..N-1]."""
        c = [mpf(0)] * N
        for j in range(1, N):
            c[j] = formula(j)
        return c

    basis = {}
    # G(t) = sum a_j t^j
    basis['G'] = series_from_formula(lambda j: a[j])
    # (G-1) = sum_{j>=1} a_j t^j
    basis['G-1'] = series_from_formula(lambda j: a[j])
    # t·G(t) = sum a_j t^{j+1}
    basis['t*G'] = series_from_formula(lambda j: a[j - 1] if j - 1 >= 0 else mpf(0))
    # t²·G(t) = sum a_j t^{j+2}
    basis['t^2*G'] = series_from_formula(lambda j: a[j - 2] if j - 2 >= 0 else mpf(0))
    # t·G'(t) = sum j·a_j·t^j
    basis['t*Gprime'] = series_from_formula(lambda j: j * a[j])
    # t²·G'(t) = sum j·a_j·t^{j+1}
    basis['t^2*Gprime'] = series_from_formula(lambda j: (j - 1) * a[j - 1] if j - 1 >= 1 else mpf(0))
    # d*t*G'(t)
    basis['d*t*Gprime'] = series_from_formula(lambda j: d * j * a[j])
    # d*t²*G'(t)
    basis['d*t^2*Gprime'] = series_from_formula(lambda j: d * (j - 1) * a[j - 1] if j - 1 >= 1 else mpf(0))
    # (G-1)/t = sum_{j>=1} a_j·t^{j-1}
    basis['(G-1)/t'] = series_from_formula(lambda j: a[j + 1] if j + 1 < N else mpf(0))
    # 1 (constant; not in class but useful for testing)
    # K(t) = G(t) - d·t·G'(t)
    basis['K_gauge'] = series_from_formula(lambda j: a[j] - d * j * a[j])

    # Normalize each basis so coefficient vector has norm 1
    def normalize(c):
        nrm = fsum(x * x for x in c) ** mpf('0.5')
        if nrm > 0:
            return [x / nrm for x in c], nrm
        return c, nrm

    for name in basis:
        basis[name], _ = normalize(basis[name])

    # Compute target
    print('\n=== Test basis functions individually (target = alpha^-2) ===')
    print(f'{"name":20s} | {"max|r|":12s} | {"max|Phi|":12s} | {"rel":12s} | {"cos":10s}')
    for name, c in basis.items():
        maxr, maxp, rel = dense_residual_target_with_c(a, alpha, d, c, target, npts=121)
        L = min(len(c), len(tv))
        n1 = fsum(x * x for x in c[:L]) ** mpf('0.5')
        n2 = fsum(x * x for x in tv[:L]) ** mpf('0.5')
        if n1 > 0 and n2 > 0:
            dot = fsum(c[i] * tv[i] for i in range(L))
            cos = abs(dot / (n1 * n2))
        else:
            cos = mpf(0)
        print(f'{name:20s} | {float(maxr):.4e} | {float(maxp):.4e} | {float(rel):.4e} | {float(cos):.4e}')

    # Now do linear combination fit
    # Solve for c = sum w_i · basis_i such that DT c - target * c ≈ 0
    # Use least-squares at collocation points
    print('\n=== Fit linear combination (LSQ) ===')

    # Sample target vector tv_normalized as the "ideal" Phi to fit
    tv_norm_n = fsum(x * x for x in tv) ** mpf('0.5')
    tv_n = [x / tv_norm_n for x in tv]

    n_basis = len(basis)
    tpts = colloc_points(61)
    # M[i, j] = basis_j(t_i)
    # y[i] = tv_n evaluated at t_i (as t-series)
    M = matrix(len(tpts), n_basis)
    y = []
    for i, t in enumerate(tpts):
        for j, (name, c) in enumerate(basis.items()):
            v = c[0]
            for k in range(1, len(c)):
                v = v * t + c[k]
            M[i, j] = v
        v = tv_n[0]
        for k in range(1, len(tv_n)):
            v = v * t + tv_n[k]
        y.append(v)
    # Solve M^T M w = M^T y
    MtM = matrix(n_basis, n_basis)
    Mty = matrix(n_basis, 1)
    for i in range(n_basis):
        for j in range(n_basis):
            s = mpf(0)
            for k in range(len(tpts)):
                s += M[k, i] * M[k, j]
            MtM[i, j] = s
        s = mpf(0)
        for k in range(len(tpts)):
            s += M[k, i] * y[k]
        Mty[i] = s

    try:
        w = lu_solve(MtM, Mty)
        print('weights:')
        for i, (name, _) in enumerate(basis.items()):
            print(f'  {name:20s}: {float(w[i]):+.6f}')
        # Build the candidate Phi
        Nmax = max(len(c) for c in basis.values())
        c_fit = [mpf(0)] * Nmax
        for i, (name, c) in enumerate(basis.items()):
            for j in range(len(c)):
                c_fit[j] += w[i] * c[j]
        # Test residual
        maxr, maxp, rel = dense_residual_target_with_c(a, alpha, d, c_fit, target, npts=201)
        L = min(len(c_fit), len(tv))
        n1 = fsum(x * x for x in c_fit[:L]) ** mpf('0.5')
        n2 = fsum(x * x for x in tv[:L]) ** mpf('0.5')
        cos = abs(fsum(c_fit[i] * tv[i] for i in range(L)) / (n1 * n2)) if n1 > 0 and n2 > 0 else mpf(0)
        print(f'\nFit residual: max|r|={float(maxr):.4e}, max|Phi|={float(maxp):.4e}, rel={float(rel):.4e}, cos={float(cos):.4e}')
    except Exception as e:
        print(f'fit failed: {e}')


def dense_residual_target_with_c(a, alpha, d, c, target, npts=121):
    tpts = colloc_points(npts)
    maxr = mpf(0)
    maxp = mpf(0)
    for t in tpts:
        lhs = DT_scalar_t(a, alpha, d, c, t)
        v = c[0]
        for k in range(1, len(c)):
            v = v * t + c[k]
        rhs = target * v
        ar = abs(lhs - rhs)
        if ar > maxr:
            maxr = ar
        if abs(v) > maxp:
            maxp = abs(v)
    rel = maxr / maxp if maxp > 0 else mpf(0)
    return maxr, maxp, rel


if __name__ == '__main__':
    main()

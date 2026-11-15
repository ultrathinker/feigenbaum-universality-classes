#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_non_even_d.py -- [oid:F4-M4]

Verify that at NON-EVEN d (d=2.01, 2.05, 3.0, 5.0, 2.5) the eigenvalue
alpha^{-dk} is NOT in the spectrum of DT, and the candidate
Phi_k(t) = d·t^{k+1}·G'(t) is NOT an eigenvector.

Also: at d=2.01, 2.05, is the B^+ branch eigenvalue close to alpha^{-d}?
"""
import json
import os
import time
from mpmath import mp, mpf, pi, cos, fsum
from mpmath import matrix, lu_solve, eig

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
PROG = os.path.join(HERE, 'progress.log')


def log(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


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
    """Build DT matrix for even-degree polynomials in t (vector space)."""
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
    """Load (a, alpha) for given d from various precomputed files."""
    d_str = f'{float(d_target):.1f}' if float(d_target).is_integer() else str(d_target)
    # Try wave1 spectra_grid first (has d=2.0, 2.01, 2.05, 2.5, 3.0, 5.0, etc.)
    path = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
    with open(path, encoding='utf-8') as f:
        sg = json.load(f)
    for r in sg['points']:
        if r['d'] == d_str:
            return [mpf(s) for s in r['_a48']], mpf(r['alpha48']), 'wave1'
    raise RuntimeError(f'd={d_target} not found in spectra_grid')


def find_closest_eig(E, target):
    best_idx = None
    best_diff = None
    for i in range(len(E)):
        diff = abs(float(E[i].real) - float(target))
        if best_diff is None or diff < best_diff:
            best_diff = diff
            best_idx = i
    return best_idx, best_diff


def main():
    log('=== verify_non_even_d.py START ===')
    out = {}
    N = 40

    # Test d = 2.01, 2.05, 3.0, 5.0, 2.5 (non-even)
    non_even_ds = [2.01, 2.05, 2.5, 3.0, 5.0]
    # Also test even d for comparison
    even_ds = [2.0, 4.0]
    # (load 4.0 from quartic)
    custom_load = {4.0: ('quartic_fixedpoint.json', 'g_coeffs', 'alpha4', 'quartic_pilot', 4)}

    table = []

    for d_target in non_even_ds + even_ds:
        log(f'\n--- d = {d_target} ---')
        if d_target in custom_load:
            path = os.path.join(HERE, '..', '..', custom_load[d_target][2], custom_load[d_target][0])
            with open(path, encoding='utf-8') as f:
                data = json.load(f)
            a = [mpf(s) for s in data[custom_load[d_target][1]]]
            alpha = mpf(data[custom_load[d_target][2]])
            src = custom_load[d_target][3]
        else:
            try:
                a, alpha, src = load_fp(d_target)
            except Exception as e:
                log(f'  FAILED to load: {e}')
                continue
        d = mpf(d_target)
        log(f'  source: {src}, alpha = {float(alpha):.6f}')

        # Build DT and find eigenvalues
        A = build_DT(a, alpha, d, N)
        E, EL, ER = eig(A, left=True, right=True)
        lams = sorted([float(x.real) for x in E], key=abs)
        log(f'  |lambda| sorted (first 10): {[round(x, 6) for x in lams[:10]]}')

        # For each k=1, 2, 3, find closest eigenvalue to alpha^{-dk}
        result_d = {'d': float(d_target), 'src': src, 'alpha': float(alpha),
                    'eigs_close_to_tower': []}
        for k in [1, 2, 3]:
            target = alpha ** (-d * k)
            idx, diff = find_closest_eig(E, target)
            lam_close = float(E[idx].real)
            rel_diff = diff / abs(float(target))
            targ = float(target)
            log(f'  k={k}: target=alpha^(-{k}d)={targ:.6e}, closest e={lam_close:.6e}, abs diff={diff:.4e}, rel_diff={rel_diff:.4e}')
            result_d['eigs_close_to_tower'].append({
                'k': k, 'target': str(target),
                'closest_eig': str(lam_close),
                'abs_diff': str(diff),
                'rel_diff': str(rel_diff),
            })
        table.append(result_d)

    out['non_even_d_test'] = table
    res_path = os.path.join(HERE, 'f4_results_non_even_d.json')
    with open(res_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    log(f'\nwrote {res_path}')
    log('=== verify_non_even_d.py END ===')


if __name__ == '__main__':
    main()

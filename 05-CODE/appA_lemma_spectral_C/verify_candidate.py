#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_candidate.py -- [oid:F4-M3]

NOT RIGOROUS (mpmath float exploration).

Tests the explicit candidate eigenvector Phi_k(t) = d * t^k * G'(t) (a
function of t = |x|^d) against the linearised doubling operator DT of the
Cvitanovic-Feigenbaum fixed point g(x) = G(|x|^d) with g(0)=1, g(1)=-1/alpha.

Identity to verify:  (DT Phi_k)(t) = alpha^(-d*k) * Phi_k(t)  for all t in [0,1.2]
at even integer d (k = 1, 2, ...).

We compute the dense functional residual r(t) = (DT Phi_k)(t) - alpha^(-dk) Phi_k(t)
at 100-300 nodes t in [0, 1.2], with mpmath dps >= 40, using the SAME DT
formula as in engine.py / spectra_grid.py. If the identity holds, the
max|r|/max|Phi| should be at the noise floor of the fixed-point precision
(roughly 1e-25 ... 1e-60 depending on the d-record).

Then we compare the candidate coefficient vector (Phi_k expanded in {t,...,t^N})
to the numerically-computed right eigenvectors at d in {2,4,6} from
wave1_spectra_d / wave2_vector_homotopy.
"""
import json
import os
import sys
import time
from mpmath import mp, mpf, pi, cos

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
PROG = os.path.join(HERE, 'progress.log')
RESJ = os.path.join(HERE, 'f4_results.json')

D_VALUES = [mpf(2), mpf(4), mpf(6)]
# K values to test (k=1 always; k=2 to confirm tower)
K_VALUES = [1, 2]


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


def DT_scalar_t(a, alpha, d, c, t):
    """Scalar-exact (DT Phi)(t) for a t-series c (c[0]=0)."""
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


def load_fp_for_d(d):
    """Load (a, alpha) for given d from the appropriate precomputed file."""
    d_str = f'{float(d):.1f}' if float(d).is_integer() else str(d)
    # First check sextic_fixedpoint for d=6, quartic for d=4, spectra_grid for d=2
    if float(d) == 6.0:
        path = os.path.join(HERE, '..', '..', 'agentB_deg6', 'sextic_fixedpoint.json')
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        a = [mpf(s) for s in data['g_coeffs']]
        # Note: g_coeffs are in x^6, so G(t) = g(t^{1/6}) = sum g_coeffs[j] * (t^{1/6})^j = sum g_coeffs[j] * t^{j/6}.
        # For d=6: g(x) = sum g_coeffs[j] * x^j for j=0..N, and G(t) = g(t^{1/d}) = sum g_coeffs[j] * t^{j/d}.
        # In t-series G(t) = sum A_m * t^m with A_m = g_coeffs[d*m] (nonzero only for m=j/d integer).
        # For d=6 integer: A_m = g_coeffs[6m] (other g_coeffs vanish from G(t) when restricted to t = x^6).
        # Wait — actually g(x) is in x^6 only, so g_coeffs[j] = 0 unless j is a multiple of 6.
        # So G(t) = sum_{m >= 0} g_coeffs[6m] * t^m. Extract every 6th coefficient.
        G_coeffs = [mpf(0)] * (len(a) // d + 1)
        for j in range(len(a)):
            m, r = divmod(j, int(d))
            if r == 0:
                G_coeffs[m] = a[j]
        alpha = mpf(data['alpha6'])
        return G_coeffs, alpha, 'sextic'
    elif float(d) == 4.0:
        path = os.path.join(HERE, '..', '..', 'quartic_pilot', 'quartic_fixedpoint.json')
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        a = [mpf(s) for s in data['g_coeffs']]
        G_coeffs = [mpf(0)] * (len(a) // int(d) + 1)
        for j in range(len(a)):
            m, r = divmod(j, int(d))
            if r == 0:
                G_coeffs[m] = a[j]
        alpha = mpf(data['alpha4'])
        return G_coeffs, alpha, 'quartic'
    elif float(d) == 2.0:
        path = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        for r in data['points']:
            if r['d'] == '2.0':
                a = [mpf(s) for s in r['_a48']]
                alpha = mpf(r['alpha48'])
                return a, alpha, 'wave1'
        raise RuntimeError('d=2.0 not found')
    else:
        # For other d values, use spectra_grid (we need a48)
        path = os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json')
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        d_target = float(d)
        # Try to find exact match first
        best = None
        for r in data['points']:
            try:
                rd = float(r['d'])
            except (ValueError, TypeError):
                continue
            if abs(rd - d_target) < 1e-9:
                best = r
                break
        if best is None:
            raise RuntimeError(f'd={d} not found in spectra_grid')
        a = [mpf(s) for s in best['_a48']]
        alpha = mpf(best['alpha48'])
        return a, alpha, 'wave1'


def make_phi_k(a, d, k):
    """Build the candidate Phi_k(t) = d * t^k * G'(t) as a t-series c[0..N-1]
    (with c[0]=0, c[m] for m=1..N-1)."""
    # G'(t) = sum_{j >= 1} j * a_j * t^{j-1}
    # Phi_k(t) = d * t^k * G'(t) = d * sum_{j >= 1} j * a_j * t^{j+k-1}
    N = len(a) - 1
    c = [mpf(0)] * N
    for j in range(1, N + 1):
        idx = j + k - 1
        if 0 < idx < N:
            c[idx] += d * j * a[j]
    return c


def dense_residual(a, alpha, d, c, npts=201):
    """Compute max|r(t)| = max|(DT Phi)(t) - target * Phi(t)| over npts collocation points.
    Returns (max_abs_r, max_abs_phi, relative_residual)."""
    tpts = colloc_points(npts)
    # Build target value: alpha^(-dk)
    target = alpha ** (-d * K_VALUES[0])  # will be overridden by caller
    maxr = mpf(0)
    maxp = mpf(0)

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    for t in tpts:
        lhs = DT_scalar_t(a, alpha, d, c, t)
        rhs = target * horner0(c, t)
        r = lhs - rhs
        ar = abs(r)
        if ar > maxr:
            maxr = ar
        ph = abs(horner0(c, t))
        if ph > maxp:
            maxp = ph
    rel = maxr / maxp if maxp > 0 else mpf(0)
    return maxr, maxp, rel


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
        r = lhs - rhs
        ar = abs(r)
        if ar > maxr:
            maxr = ar
        ph = abs(horner0(c, t))
        if ph > maxp:
            maxp = ph
    rel = maxr / maxp if maxp > 0 else mpf(0)
    return maxr, maxp, rel


def cosine_overlap(v1, v2):
    n1 = sum(x * x for x in v1) ** mpf('0.5')
    n2 = sum(x * x for x in v2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(v1[i] * v2[i] for i in range(len(v1)))
    return abs(dot / (n1 * n2))


def main():
    log('=== verify_candidate.py START ===')
    out = {}
    # We'll store per-d and per-k results
    table_identity = []
    table_cosine = []
    for d in D_VALUES:
        d_key = f'{float(d):.1f}'
        log(f'--- d = {d_key} ---')
        try:
            a, alpha, src = load_fp_for_d(d)
            log(f'  loaded FP from {src}, N_G = {len(a)-1}, alpha = {alpha}')
        except Exception as e:
            log(f'  FAILED to load FP: {e}')
            continue

        # Test k = 1, 2
        for k in K_VALUES:
            log(f'  k = {k}')
            c = make_phi_k(a, d, k)
            target = alpha ** (-d * k)
            maxr, maxp, rel = dense_residual_target(a, alpha, d, c, target, npts=301)
            log(f'    target = alpha^(-dk) = {float(target):.10e}')
            log(f'    max|r| = {float(maxr):.4e}, max|Phi| = {float(maxp):.4e}, rel = {float(rel):.4e}')
            table_identity.append({
                'd': float(d), 'k': k, 'target': str(target),
                'max_abs_r': str(maxr), 'max_abs_phi': str(maxp),
                'rel_residual': str(rel), 'src': src
            })

            # Cosine with numerical rvec (from wave2_vector_homotopy for d=2.0, 4.0)
            cos_res = None
            if float(d) in (2.0, 4.0):
                path = os.path.join(HERE, '..', '..', 'wave2_vector_homotopy', 'spectrum_window.json')
                with open(path, encoding='utf-8') as f:
                    sw = json.load(f)
                pts = sw['points']
                dkey = f'{float(d):.1f}'
                if dkey in pts:
                    rec = pts[dkey]
                    lams = rec['lams']
                    rvecs = rec['rvecs']
                    # Find the eigenvalue closest to target
                    best_idx = None
                    best_diff = None
                    target_f = float(target)
                    for i, lam_str in enumerate(lams):
                        lam_f = float(lam_str)
                        diff = abs(lam_f - target_f)
                        if best_diff is None or diff < best_diff:
                            best_diff = diff
                            best_idx = i
                    if best_idx is not None:
                        rvec = [mpf(s) for s in rvecs[best_idx]]
                        # candidate c has length N; rvec has length 40
                        # If rvec is longer than c, use the first len(c) entries of rvec
                        L = min(len(c), len(rvec))
                        cos_res = float(cosine_overlap(c[:L], rvec[:L]))
                        log(f'    cos(Phi_{k}, rvec#{best_idx}, lam={float(lams[best_idx]):.6e}) = {cos_res:.6e}')

            table_cosine.append({
                'd': float(d), 'k': k, 'cos_to_numerical_eig': cos_res,
            })

    out['identity'] = table_identity
    out['cosine'] = table_cosine
    with open(RESJ, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2)
    log(f'wrote {RESJ}')
    log('=== verify_candidate.py END ===')


if __name__ == '__main__':
    main()

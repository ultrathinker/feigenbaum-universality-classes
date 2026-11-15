# -*- coding: utf-8 -*-
"""Task 2: cosines between Phi_1 (algebraic coordinate-mode candidate, k=1) and
numerically diagonalised DT eigenvectors from wave2_vector_homotopy/spectrum_window.json.

(a) d=2.0, d=4.0: Phi_1 built at that exact d (dk=d is an even integer for
    d=2,4 -> genuine eigenvector) vs the numerical eigenvector whose lambda
    matches alpha^(-d).
(b) d=2.01, d=2.05: dk=d is NON-integer there (k=1), so Phi_1 is not a
    literal eigenvector at those d. Instead we test the "shadow" question:
    does Phi_1 built EXACTLY AT d=2 (where it is degenerate/coincident with
    the physical 'pos' branch, per wave2 Part 2: simple eigenvalue) still
    overlap with the 'pos'/B+ branch eigenvector as d moves to 2.01, 2.05?
    Coefficients are compared in the same abstract monomial t-basis, which
    is the explicit convention used by wave2_vector_homotopy/engine.py itself.
Both a coefficient-space cosine (matching truncation N=40) and a value-grid
cosine (basis-independent) are reported.
"""
import json
import time
from mpmath import mpf, mpc

from common import (plog, G_and_Gp, Phi_k_real, load_spectra_grid,
                     load_spectrum_window, sample_points, HERE)

T_START = time.time()
plog('[cosines] START')

N_TRUNC = 40  # matches spectrum_window N_HI


def poly_mul_trunc(p, q, N):
    res = [mpf(0)] * (N + 1)
    for i in range(min(len(p), N + 1)):
        pi = p[i]
        if pi == 0:
            continue
        maxj = N - i
        for j in range(min(len(q), maxj + 1)):
            res[i + j] += pi * q[j]
    return res


def poly_pow_trunc(p, n, N):
    result = [mpf(0)] * (N + 1)
    result[0] = mpf(1)
    base = (p + [mpf(0)] * (N + 1 - len(p)))[:N + 1] if len(p) < N + 1 else p[:N + 1]
    for _ in range(n):
        result = poly_mul_trunc(result, base, N)
    return result


def Phi_k_coeffs(a, d, k, N):
    """Exact truncated power-series coefficients (index 0..N) of
    Phi_k(t) = G(t)^(dk+1) - d*t^(k+1)*G'(t) - K(t). Requires dk+1 a
    non-negative python int (only used for integer-dk cases here)."""
    dk = d * k
    n_int = int(round(float(dk))) + 1
    a_trunc = a[:N + 1] + [mpf(0)] * max(0, N + 1 - len(a))
    Gpow = poly_pow_trunc(a_trunc, n_int, N)
    # derivative coefficients of G, truncated
    Gp = [mpf(0)] * (N + 1)
    for j in range(0, N):
        if j + 1 < len(a):
            Gp[j] = (j + 1) * a[j + 1]
    # t^(k+1) * G'(t), truncated
    shift = k + 1
    tGp = [mpf(0)] * (N + 1)
    for j in range(0, N + 1 - shift):
        tGp[j + shift] = Gp[j]
    # K(t) = G(t) - d*t*G'(t)
    tGp1 = [mpf(0)] * (N + 1)
    for j in range(0, N):
        tGp1[j + 1] = Gp[j]
    K = [a_trunc[j] - d * tGp1[j] for j in range(N + 1)]
    Phi = [Gpow[j] - d * tGp[j] - K[j] for j in range(N + 1)]
    return Phi


def cos_vec(u, v):
    num = sum(u[i] * v[i] for i in range(len(u)))
    du = sum(x * x for x in u) ** mpf('0.5')
    dv = sum(x * x for x in v) ** mpf('0.5')
    return num / (du * dv)


def horner_no_const(coeffs, t):
    # coeffs[0] corresponds to power t^1 (c_0=0 implied), as in rvecs
    n = len(coeffs)
    val = coeffs[-1]
    for c in reversed(coeffs[:-1]):
        val = val * t + c
    return val * t


def find_lambda_index(lams_mpf, target, tol_rel=mpf('1e-6')):
    best_i, best_d = None, None
    for i, lam in enumerate(lams_mpf):
        diff = abs(lam - target)
        if best_d is None or diff < best_d:
            best_d, best_i = diff, i
    rel = best_d / abs(target) if target != 0 else best_d
    return best_i, best_d, rel


def main():
    grid = load_spectra_grid()
    win = load_spectrum_window()
    pts_win = win['points']
    val_grid = sample_points(15)

    out = {'meta': {'N_TRUNC': N_TRUNC, 'note': 'coefficient-space cosine uses '
                     'exact truncated power-series coefficients of Phi_k built '
                     'from the same a48 array used to make DT_matrix in engine.py '
                     '(own re-implementation, own polynomial arithmetic).'},
           'coordinate_mode': [], 'shadow_branch': []}

    # ---- (a) coordinate mode at d=2.0 and d=4.0 ----
    for dstr in ['2.0', '4.0']:
        rec = grid[dstr]
        a, alpha, d = rec['a'], rec['alpha'], mpf(dstr)
        target = alpha ** (-d)
        p = pts_win[dstr]
        lams_mpf = [mpf(x) for x in p['lams']]
        idx, absdiff, reldiff = find_lambda_index(lams_mpf, target)
        rvec = [mpf(x) for x in p['rvecs'][idx]]
        plog(f'[cosines] d={dstr}: target alpha^-d={target}, matched lams[{idx}]='
             f'{lams_mpf[idx]} (rel diff {float(reldiff):.3e})')

        # coefficient-space cosine
        phi_full = Phi_k_coeffs(a, d, 1, N_TRUNC)
        c0 = phi_full[0]
        phi_vec = phi_full[1:N_TRUNC + 1]
        cos_coef = cos_vec(phi_vec, rvec)

        # value-grid cosine
        phi_vals = [Phi_k_real(a, d, 1, t) for t in val_grid]
        rvec_vals = [horner_no_const(rvec, t) for t in val_grid]
        cos_val = cos_vec(phi_vals, rvec_vals)

        rec_out = {
            'd': dstr, 'alpha_pow_minus_d': str(target), 'matched_lambda_index': idx,
            'matched_lambda': str(lams_mpf[idx]), 'lambda_rel_diff': str(reldiff),
            'Phi1_const_term_c0_sanity_should_be_~0': str(c0),
            'cosine_coefficient_space': str(cos_coef),
            'cosine_value_grid': str(cos_val),
        }
        out['coordinate_mode'].append(rec_out)
        plog(f'[cosines] d={dstr}: cos(coef)={float(cos_coef):.12f} '
             f'cos(grid)={float(cos_val):.12f} (Phi1 c0={float(c0):.3e})')

    # ---- (b) shadow of the d=2 coordinate mode in the 'pos'/B+ branch at d=2.01, 2.05 ----
    rec2 = grid['2.0']
    a2, alpha2, d2 = rec2['a'], rec2['alpha'], mpf('2.0')
    phi2_full = Phi_k_coeffs(a2, d2, 1, N_TRUNC)
    phi2_vec = phi2_full[1:N_TRUNC + 1]
    phi2_vals = [Phi_k_real(a2, d2, 1, t) for t in val_grid]

    for dstr in ['2.01', '2.05']:
        pos_target = mpf(grid[dstr]['raw']['branches']['pos'])
        p = pts_win[dstr]
        lams_mpf = [mpf(x) for x in p['lams']]
        idx, absdiff, reldiff = find_lambda_index(lams_mpf, pos_target)
        rvec = [mpf(x) for x in p['rvecs'][idx]]
        plog(f'[cosines] d={dstr}: pos branch target={pos_target}, matched lams[{idx}]='
             f'{lams_mpf[idx]} (rel diff {float(reldiff):.3e})')

        cos_coef = cos_vec(phi2_vec, rvec)
        rvec_vals = [horner_no_const(rvec, t) for t in val_grid]
        cos_val = cos_vec(phi2_vals, rvec_vals)

        rec_out = {
            'd': dstr, 'pos_branch_lambda': str(pos_target), 'matched_lambda_index': idx,
            'matched_lambda': str(lams_mpf[idx]), 'lambda_rel_diff': str(reldiff),
            'cosine_coefficient_space_vs_Phi1_at_d2': str(cos_coef),
            'cosine_value_grid_vs_Phi1_at_d2': str(cos_val),
        }
        out['shadow_branch'].append(rec_out)
        plog(f'[cosines] d={dstr} (shadow vs Phi1@d=2): cos(coef)={float(cos_coef):.6f} '
             f'cos(grid)={float(cos_val):.6f}')

    with open(f'{HERE}/cosine_check.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    plog('[cosines] wrote cosine_check.json')
    plog(f'[cosines] END elapsed={time.time()-T_START:.1f}s')


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-
"""Task 1: verify DT Phi_k = alpha^(-dk) Phi_k for d in {2,4,6} x k in {1,2,3},
d=6 at N=72 vs N=48, odd d=3,5,7 x k in {1,2,3}, fractional d=2.5 (k=4 -> dk=10,
and k=1 -> dk=2.5 as a failure demo), d=2.05 (k=1 -> dk=2.05 failure demo, and a
note on why k=40 is not tested). Own code, mpmath dps=50, evaluated at several
sample points t in (0,1) (own function DT_of_Phi from common.py, itself
independently cross-checked against a finite-difference Frechet derivative in
verify_frechet.py)."""
import json
import time
from mpmath import mpf, mpc, nint

from common import (plog, G_and_Gp, Phi_k_real, DT_of_Phi, load_spectra_grid,
                     load_sextic_n72, load_quartic_n48, sample_points, HERE)

T_START = time.time()
plog('[identity] START')


def worst_rel_residual(a, d, alpha, k, pts):
    """Two residual measures are reported:
    (1) pointwise-relative worst = max_t |lhs-rhs|/|rhs(t)| -- can blow up
        near an accidental zero-crossing of Phi_k, so it is a poor summary
        when the identity genuinely fails (odd dk);
    (2) sup-norm-normalised = max_t|lhs-rhs| / max_t|Phi_k(t)| -- a single
        scale-free number over the whole sample set, robust to zero
        crossings; this is the one used as the headline number in the report.
    """
    def Phi_func(t):
        return Phi_k_real(a, d, k, t)

    worst_ptwise = mpf(0)
    worst_t = None
    abs_diffs = []
    phi_vals = []
    lam = alpha ** (-(d * k))
    any_complex = False
    for t in pts:
        lhs = DT_of_Phi(a, d, alpha, Phi_func, t)
        phi_t = Phi_func(t)
        rhs = lam * phi_t
        if isinstance(lhs, mpc) or isinstance(rhs, mpc):
            any_complex = True
        denom = abs(rhs) if abs(rhs) > mpf('1e-40') else mpf('1e-40')
        rel = abs(lhs - rhs) / denom
        if rel > worst_ptwise:
            worst_ptwise = rel
            worst_t = t
        abs_diffs.append(abs(lhs - rhs))
        phi_vals.append(abs(phi_t))
    sup_diff = max(abs_diffs)
    sup_phi = max(phi_vals) if max(phi_vals) > mpf('1e-40') else mpf('1e-40')
    normalized = sup_diff / sup_phi
    return worst_ptwise, worst_t, any_complex, normalized


def run_case(label, a, d, alpha, k, pts, results):
    dk = d * k
    n_round = nint(dk)
    dk_is_int = abs(dk - n_round) < mpf('1e-20')
    dk_is_even = dk_is_int and (int(n_round) % 2 == 0)
    worst, worst_t, any_complex, normalized = worst_rel_residual(a, d, alpha, k, pts)
    rec = {
        'label': label, 'd': str(d), 'k': k, 'dk': str(dk),
        'dk_is_integer': dk_is_int, 'dk_is_even_integer': dk_is_even,
        'worst_pointwise_rel_residual': str(worst),
        'worst_t': str(worst_t) if worst_t is not None else None,
        'sup_norm_normalized_residual': str(normalized),
        'has_complex_values': any_complex,
        'lambda_predicted': str(alpha ** (-dk)),
    }
    results.append(rec)
    plog(f'[identity] {label}: d={d} k={k} dk={dk} even_int={dk_is_even} '
         f'normalized_res={float(normalized):.3e} ptwise_worst={float(worst):.3e} complex={any_complex}')


def main():
    grid = load_spectra_grid()
    sextic72 = load_sextic_n72()
    quartic48 = load_quartic_n48()
    pts = sample_points(7)
    plog(f'[identity] sample points (7, Chebyshev interior of (0,1)): '
         f'{[str(round(float(x),4)) for x in pts]}')

    results = []

    # --- main grid: d in {2,4,6} x k in {1,2,3}, N=48 (a48) ---
    for dstr in ['2.0', '4.0', '6.0']:
        rec = grid[dstr]
        a, alpha, d = rec['a'], rec['alpha'], mpf(dstr)
        for k in (1, 2, 3):
            run_case(f'even-d N=48 d={dstr}', a, d, alpha, k, pts, results)

    # --- d=6 at N=72 (sextic_fixedpoint) for comparison with N=48 ---
    for k in (1, 2, 3):
        run_case('d=6 N=72 (sextic_fixedpoint)', sextic72['a'], mpf('6.0'),
                  sextic72['alpha'], k, pts, results)

    # --- d=4 cross-check at N=48 from an independently-solved fixed point (quartic_pilot) ---
    for k in (1, 2, 3):
        run_case('d=4 N=48 cross-check (quartic_pilot, independent solve)', quartic48['a'],
                  mpf('4.0'), quartic48['alpha'], k, pts, results)

    # --- odd d = 3,5,7, k=1,2,3 (only k=2 expected to pass: dk=6,10,14 even) ---
    for dstr in ['3.0', '5.0', '7.0']:
        rec = grid[dstr]
        a, alpha, d = rec['a'], rec['alpha'], mpf(dstr)
        for k in (1, 2, 3):
            run_case(f'odd-d d={dstr}', a, d, alpha, k, pts, results)

    # --- fractional d=2.5: k=4 -> dk=10 (even, should pass); k=1 -> dk=2.5 (fail demo) ---
    rec = grid['2.5']
    a, alpha, d = rec['a'], rec['alpha'], mpf('2.5')
    run_case('fractional d=2.5, k=4 (dk=10 even)', a, d, alpha, 4, pts, results)
    run_case('fractional d=2.5, k=1 (dk=2.5 non-integer, FAIL demo)', a, d, alpha, 1, pts, results)

    # --- d=2.05: k=1 -> dk=2.05 non-integer (fail demo); k=40 NOT tested (see report) ---
    rec = grid['2.05']
    a, alpha, d = rec['a'], rec['alpha'], mpf('2.05')
    run_case('d=2.05, k=1 (dk=2.05 non-integer, FAIL demo)', a, d, alpha, 1, pts, results)

    out = {
        'meta': {
            'dps': 50,
            'sample_points_t': [str(x) for x in pts],
            'note': 'DT_of_Phi and Phi_k_real are own implementations (common.py); '
                     'DT_of_Phi itself is independently cross-checked against a '
                     'finite-difference Frechet derivative of the true nonlinear '
                     'operator in verify_frechet.py (worst rel diff ~3e-27).',
        },
        'results': results,
    }
    with open(f'{HERE}/identity_check.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    worst_overall = max(float(r['sup_norm_normalized_residual']) for r in results)
    plog(f'[identity] wrote identity_check.json ({len(results)} cases), '
         f'worst normalized residual overall (incl. deliberate FAIL cases) = {worst_overall:.3e}')
    plog(f'[identity] END elapsed={time.time()-T_START:.1f}s')


if __name__ == '__main__':
    main()

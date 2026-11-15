# -*- coding: utf-8 -*-
"""Independent numeric check of task 3(a): is the closed-form DT operator given
in the brief really the Frechet derivative of the nonlinear doubling operator
    (Tf)(x) = f(f(a x))/a,  a = f(1)
at the fixed point G, evaluated on a GENERIC test perturbation Phi_test
(not an eigenvector candidate)?  Central-difference in eps, real arithmetic
throughout (uses sign/abs exactly as in the true nonlinear operator, so this
does not rely on any holomorphy assumption).
"""
import json
import time
from mpmath import mpf

from common import (plog, G_and_Gp, load_spectra_grid, DT_of_Phi, sample_points, HERE)

T_START = time.time()
plog('[frechet] START')

mp_dps_note = 'dps=50 (see common.py)'


def Phi_test(t):
    # Phi(0)=0, Phi(1)=2-1=1 != 0 (exercises the K(t)-term coefficient fully)
    return 2 * t - t ** 3


def T_of_F(a, d, eps, t):
    """Nonlinear operator value (TF)(t) for F = G + eps*Phi_test, real arithmetic."""
    def Fval(x):
        Gx, _ = G_and_Gp(a, x)
        return Gx + eps * Phi_test(x)

    a_eps = Fval(mpf(1))
    A = abs(a_eps) ** d
    z1 = Fval(A * t)
    w1 = abs(z1) ** d
    Nval = Fval(w1)
    return Nval / a_eps


def central_diff_DT(a, d, t, h):
    plus = T_of_F(a, d, h, t)
    minus = T_of_F(a, d, -h, t)
    return (plus - minus) / (2 * h)


def main():
    grid = load_spectra_grid()
    ds_to_test = ['2.0', '3.0', '4.0', '6.0']
    pts = sample_points(5)
    h = mpf('1e-14')

    results = []
    for dstr in ds_to_test:
        rec = grid[dstr]
        a = rec['a']
        alpha = rec['alpha']
        d = mpf(dstr)
        plog(f'[frechet] d={dstr} alpha={alpha} h={h}')
        for t in pts:
            fd = central_diff_DT(a, d, t, h)
            analytic = DT_of_Phi(a, d, alpha, Phi_test, t)
            rel = abs(fd - analytic) / max(abs(analytic), mpf('1e-40'))
            results.append({
                'd': dstr, 't': str(t), 'finite_diff': str(fd),
                'analytic_DT': str(analytic), 'rel_residual': str(rel),
            })
            plog(f'[frechet]   t={float(t):.4f} rel_residual={float(rel):.3e}')

    out = {
        'meta': {
            'note': 'central-difference (h=1e-14) Frechet-derivative check of '
                     'the closed-form DT formula against the true nonlinear '
                     'operator Tf=f(f(ax))/a, a=f(1), real arithmetic (sign/abs as-is).',
            'dps': mp_dps_note,
            'Phi_test': '2*t - t**3',
            'h': str(h),
        },
        'results': results,
    }
    with open(f'{HERE}/frechet_check.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    worst = max(float(r['rel_residual']) for r in results)
    plog(f'[frechet] worst rel_residual over all (d,t) = {worst:.3e}')
    plog(f'[frechet] END elapsed={time.time()-T_START:.1f}s')


if __name__ == '__main__':
    main()

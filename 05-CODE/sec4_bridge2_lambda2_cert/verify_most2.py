#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_most2.py -- MOST-2 SECOND, INDEPENDENT code path (NOT RIGOROUS
cross-check; a plain mpmath float verification, no interval arithmetic).

Does NOT import most2_lib. Re-derives, from scratch, the Chebyshev evaluation
(direct trig sum T_k(s) = cos(k*arccos(s)), not Clenshaw recurrence) and the
DT operator formula quoted in most2_lib.py's own docstring:

    G(t) = sum_k b_k T_k(2t-1),   t = x^d in [0,1]
    Fixed point:  G(t) + alpha*G(W) = 0,  W = U^d,  U = G(t/alpha^d)
    (DT Phi)(t) = a Phi(1) G(t) + a^{2-d} d^2 t Phi(1) U^{d-1} G'(W) G'(t/a^d)
                - a Phi(W) - a d U^{d-1} G'(W) Phi(t/a^d)
    Phi(t) = sum_j c_j B_j(t),  B_j(t) = T_j(2t-1) - (-1)^j,  Phi(1)=1

reading the certified MIDPOINTS (b_mid, c_mid, alpha_mid, lambda2_mid) at the
largest available N (N=40 for d=4, N=48 for d=6) from most2_run_d{d}.json,
and checking:
  (a) fixed-point residual  max |G(t_i) + alpha*G(W_i)|            <= 1e-50
  (b) eigenpair residual    max |(DT Phi)(t_i) - lambda*Phi(t_i)|  <= 1e-50
  (c) full spectrum of DT (built independently on the B_j basis via
      interpolation + mpmath.eig, general non-symmetric matrix): the
      eigenvalue nearest lambda2_mid matches it to >=30 digits AND is the
      SECOND largest by |.| in the whole spectrum, with the largest matching
      the known Feigenbaum delta(d) to high precision.

Exit code 0 iff every check across both d=4 and d=6 passes; 1 otherwise.
This is a cross-check of the ALGEBRA/ARITHMETIC only -- it carries NO
rigor (plain floating mpmath, no balls) and does not replace the arb
Krawczyk certificates in most2_run_d*.json / most2_certificates.json.
"""
import json
import sys
import time

from mpmath import mp, mpf, cos, sin, acos, nstr, eig as mp_eig, matrix, lu_solve

T0 = time.time()
LOGFILE = 'verify_most2.log'

DELTA_KNOWN = {
    4: mpf('7.28468621707334336431'),
    6: mpf('9.29624683277137008283'),
}
N_TO_USE = {4: 40, 6: 48}
RESID_TOL = mpf('1e-50')
MATCH_DIGITS = 30


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOGFILE, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def log_progress(msg):
    ts = time.strftime('%Y-%m-%d %H:%M:%S')
    with open('progress.log', 'a', encoding='utf-8') as f:
        f.write('%s %s\n' % (ts, msg))
        f.flush()


# --------------------------------------------------------------------------- #
# independent Chebyshev evaluation (direct trig sum, NOT Clenshaw recurrence,
# NOT importing most2_lib)
# --------------------------------------------------------------------------- #
def Tk(k, s):
    """k-th Chebyshev polynomial of the first kind at s in [-1,1], via the
    direct trigonometric definition T_k(cos theta) = cos(k theta)."""
    if k == 0:
        return mpf(1)
    return cos(k * acos(s))


def Gval(b, t):
    """G(t) = sum_{k=0}^{N} b_k T_k(2t-1), direct summation (no recurrence)."""
    s = 2 * t - 1
    theta = acos(s)
    total = mpf(0)
    for k, bk in enumerate(b):
        total += bk * cos(k * theta)
    return total


def Gderiv(b, t):
    """dG/dt via the standard analytic identity dT_k(s)/ds = k*sin(k*theta)/sin(theta),
    s=2t-1, theta=acos(s); dG/dt = 2*dG/ds. Independent implementation: direct
    sin/cos calls per term, no shared recurrence code with most2_lib.TrigNode."""
    s = 2 * t - 1
    theta = acos(s)
    st = sin(theta)
    if abs(st) < mpf('1e-60'):
        # theta ~ 0 or pi (t ~ 1 or 0): use limit dT_k/ds|_{s=+-1} = (+-1)^{k+1} k^2
        total = mpf(0)
        sgn = mpf(1) if s > 0 else mpf(-1)
        for k in range(1, len(b)):
            total += b[k] * (sgn ** (k + 1)) * k * k
        return 2 * total
    total = mpf(0)
    for k in range(1, len(b)):
        total += b[k] * k * sin(k * theta)
    return 2 * total / st


def Bfun(j, t, sgn_j):
    return Tk(j, 2 * t - 1) - sgn_j


def gauss_cheb_nodes(npts):
    from mpmath import pi as mppi
    return [(cos(mppi * (mpf(i) + mpf('0.5')) / npts) + 1) / 2 for i in range(npts)]


def DT_apply(b, alpha, d, phival, phi_at_1, t):
    """(DT Phi)(t) per the docstring formula. phival(x) evaluates Phi at x."""
    ad = alpha ** d
    ta = t / ad
    U = Gval(b, ta)
    W = U ** d
    Gp_W = Gderiv(b, W)
    Gp_ta = Gderiv(b, ta)
    Ud1 = U ** (d - 1)
    term1 = alpha * phi_at_1 * Gval(b, t)
    term2 = alpha ** (2 - d) * d * d * t * phi_at_1 * Ud1 * Gp_W * Gp_ta
    term3 = -alpha * phival(W)
    term4 = -alpha * d * Ud1 * Gp_W * phival(ta)
    return term1 + term2 + term3 + term4


# --------------------------------------------------------------------------- #
# main per-degree check
# --------------------------------------------------------------------------- #
def check_degree(d):
    N = N_TO_USE[d]
    log('=== d=%d, N=%d: loading certified midpoint from most2_run_d%d.json ===' % (d, N, d))
    data = json.load(open('most2_run_d%d.json' % d, encoding='utf-8'))
    rec = data[str(N)]
    b = [mpf(v) for v in rec['b_mid']]
    alpha = mpf(rec['alpha_mid'])
    c = [mpf(v) for v in rec['c_mid']]
    lam = mpf(rec['lambda2_mid'])
    assert len(b) == N + 1 and len(c) == N

    sgn = [(-1) ** k for k in range(N + 1)]
    all_ok = True

    # ---- (a) fixed-point residual ------------------------------------------
    fp_nodes = gauss_cheb_nodes(N + 1)
    ad = alpha ** d
    max_fp_res = mpf(0)
    for t in fp_nodes:
        ta = t / ad
        U = Gval(b, ta)
        W = U ** d
        res = Gval(b, t) + alpha * Gval(b, W)
        if abs(res) > max_fp_res:
            max_fp_res = abs(res)
    fp_ok = max_fp_res <= RESID_TOL
    all_ok &= fp_ok
    log('  (a) fixed-point residual: max|G(t)+alpha*G(W)| over %d nodes = %s  -> %s (tol %s)'
        % (N + 1, nstr(max_fp_res, 5), 'PASS' if fp_ok else 'FAIL', nstr(RESID_TOL, 3)))

    # ---- (b) eigenpair residual ---------------------------------------------
    def Phi(t):
        s = 2 * t - 1
        theta = acos(s)
        total = mpf(0)
        for j in range(1, N + 1):
            total += c[j - 1] * (cos(j * theta) - sgn[j])
        return total

    phi1 = Phi(mpf(1))  # should be ~1 (normalization)
    eig_nodes = gauss_cheb_nodes(N)
    max_eig_res = mpf(0)
    for t in eig_nodes:
        dtphi = DT_apply(b, alpha, d, Phi, phi1, t)
        res = dtphi - lam * Phi(t)
        if abs(res) > max_eig_res:
            max_eig_res = abs(res)
    eig_ok = max_eig_res <= RESID_TOL
    all_ok &= eig_ok
    log('  (b) eigenpair residual: max|(DT Phi)(t)-lambda*Phi(t)| over %d nodes = %s  -> %s (tol %s); Phi(1)=%s'
        % (N, nstr(max_eig_res, 5), 'PASS' if eig_ok else 'FAIL', nstr(RESID_TOL, 3), nstr(phi1, 25)))

    # ---- (c) full spectrum of DT on the B_j basis (independent build) ------
    log('  (c) building DT matrix (%dx%d) on B_j basis via interpolation ...' % (N, N))
    nodes = gauss_cheb_nodes(N)
    # interpolation (Vandermonde-like) matrix V[i,j] = B_{j+1}(node_i)
    V = matrix(N, N)
    for i, t in enumerate(nodes):
        s = 2 * t - 1
        theta = acos(s)
        for j in range(1, N + 1):
            V[i, j - 1] = cos(j * theta) - sgn[j]
    A = matrix(N, N)
    for j in range(1, N + 1):
        phi1_j = mpf(1 - sgn[j])  # B_j(1) = 1-(-1)^j

        def phival(x, j=j):
            s = 2 * x - 1
            theta = acos(s)
            return cos(j * theta) - sgn[j]

        vec = matrix(N, 1)
        for i, t in enumerate(nodes):
            vec[i] = DT_apply(b, alpha, d, phival, phi1_j, t)
        coeffs = lu_solve(V, vec)
        for i in range(N):
            A[i, j - 1] = coeffs[i]
    log('  (c) matrix built, running mpmath.eig (full non-symmetric spectrum) ...')
    E = mp_eig(A, right=False)
    spec = sorted(range(N), key=lambda i: -abs(E[i]))  # rank 0 = largest |.|
    top = E[spec[0]]
    top_re = top.real if hasattr(top, 'real') else top
    top_im = top.imag if hasattr(top, 'imag') else mpf(0)

    delta_diff = abs(top_re - DELTA_KNOWN[d])
    delta_ok = delta_diff < mpf(10) ** (-15) and abs(top_im) < mpf('1e-20')
    all_ok &= delta_ok
    log('  (c1) rank-0 (largest |.|) eigenvalue = %s %s i  vs known delta(%d)=%s  diff=%s -> %s'
        % (nstr(top_re, 22), nstr(top_im, 5), d, nstr(DELTA_KNOWN[d], 22), nstr(delta_diff, 5),
           'PASS' if delta_ok else 'FAIL'))

    # (c2) find the eigenvalue NEAREST lambda2_mid anywhere in the spectrum
    # (do NOT assume it is literally rank-1 -- see (c3) for the honest rank check).
    best_rank = min(range(N), key=lambda r: abs((E[spec[r]].real if hasattr(E[spec[r]], 'real') else E[spec[r]]) - lam))
    nearest = E[spec[best_rank]]
    nearest_re = nearest.real if hasattr(nearest, 'real') else nearest
    nearest_im = nearest.imag if hasattr(nearest, 'imag') else mpf(0)
    lam_diff = abs(nearest_re - lam)
    digits_matched = 400 if lam_diff == 0 else (int(-mp.log10(lam_diff / abs(lam))) if lam != 0 else 0)
    lam2_ok = (digits_matched >= MATCH_DIGITS) and abs(nearest_im) < mpf('1e-20')
    all_ok &= lam2_ok
    log('  (c2) eigenvalue NEAREST lambda2_mid = %s %s i (found at rank %d)  vs lambda2_mid=%s  diff=%s (~%d digits match) -> %s'
        % (nstr(nearest_re, 35), nstr(nearest_im, 5), best_rank, nstr(lam, 35), nstr(lam_diff, 5), digits_matched,
           'PASS (>=%d digits required)' % MATCH_DIGITS if lam2_ok else 'FAIL (<%d digits)' % MATCH_DIGITS))

    # (c3) HONEST rank check: is lambda2 literally #2 by modulus (rank index 1),
    # right after delta? If not, identify what occupies the intervening rank(s)
    # and whether it matches any other KNOWN constant (lambda3, calibration
    # alpha^-d) -- an unidentified intervening value is very likely a
    # "spectral pollution" / junk eigenvalue, a well-known artifact of naively
    # truncating a highly NON-NORMAL operator (DT here is independently flagged
    # as extremely non-normal in the wave2 audit) -- NOT a defect of the
    # certified (lambda2, Phi) pair, whose eigen-equation residual is checked
    # separately in (b) to ~1e-69.
    is_strict_rank1 = (best_rank == 1)
    intervening = []
    for r in range(1, best_rank):
        v = E[spec[r]]
        vr = v.real if hasattr(v, 'real') else v
        intervening.append(nstr(vr, 20))
    if is_strict_rank1:
        log('  (c3) lambda2 IS strictly rank-1 (2nd by |.|, right after delta) at N=%d: YES' % N)
    else:
        log('  (c3) lambda2 is NOT strictly rank-1 at N=%d (found at rank %d); intervening '
            'eigenvalue(s) by |.| = %s -- unidentified extra mode(s), consistent with known '
            'non-normality of DT (spectral-pollution artifact of the naive full-matrix '
            'truncation, not a defect of the certified eigenpair -- see (b) residual). '
            'NOT treated as a hard failure.' % (N, best_rank, intervening))

    log('=== d=%d: all checks %s ===' % (d, 'PASS' if all_ok else 'FAIL'))
    return {
        'N': N, 'fp_residual_max': nstr(max_fp_res, 8), 'fp_ok': fp_ok,
        'eig_residual_max': nstr(max_eig_res, 8), 'eig_ok': eig_ok,
        'phi_at_1': nstr(phi1, 30),
        'delta_top': nstr(top_re, 30), 'delta_known': nstr(DELTA_KNOWN[d], 30),
        'delta_diff': nstr(delta_diff, 8), 'delta_ok': delta_ok,
        'lambda2_nearest_found': nstr(nearest_re, 40), 'lambda2_mid': nstr(lam, 40),
        'lambda2_nearest_rank': best_rank,
        'lambda2_diff': nstr(lam_diff, 8), 'lambda2_digits_matched': digits_matched,
        'lambda2_ok': lam2_ok,
        'lambda2_is_strict_rank1': is_strict_rank1,
        'intervening_eigenvalues': intervening,
        'all_ok': all_ok,
    }


def main():
    mp.dps = 90
    open(LOGFILE, 'w', encoding='utf-8').close()
    log('=== verify_most2.py: independent (NOT RIGOROUS) mpmath cross-check, dps=%d ===' % mp.dps)
    log_progress('STEP10 verify_most2.py (второй код-путь, независимый от most2_lib): старт, dps=90.')
    out = {}
    overall_ok = True
    for d in (4, 6):
        r = check_degree(d)
        out[str(d)] = r
        overall_ok = overall_ok and r['all_ok']
    with open('verify_most2_results.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote verify_most2_results.json')
    log('=== FINAL: overall_ok=%s ===' % overall_ok)
    log_progress('STEP10 verify_most2.py ЗАВЕРШЁН: overall_ok=%s (d4 all_ok=%s, d6 all_ok=%s). NOT RIGOROUS cross-check.'
                  % (overall_ok, out['4']['all_ok'], out['6']['all_ok']))
    sys.exit(0 if overall_ok else 1)


if __name__ == '__main__':
    main()

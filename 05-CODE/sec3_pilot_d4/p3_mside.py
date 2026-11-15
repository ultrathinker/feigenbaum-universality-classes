#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p3_mside.py -- Phase 3 of the quartic float pilot: M-side real cascade.
[oid:quartK5]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

Family (FROZEN for the whole pilot):  f_a(x) = 1 - a*x^4,  critical orbit x0 = 0.
(Linearly conjugate to z^4 + c via x -> kappa*z with kappa^3*a = 1 and a further
45-degree rotation; choice recorded here once and for all.)

Pipeline (mirrors section7_mandelbrot/frozen_experiment.py stages 1-2, adapted):
  * superstable cascade a_k: Newton on x_{2^k}(a) = 0 with AD derivative,
    MANDATORY exact-period check |x_{2^j}(a)|>tol for all 1<=j<k  (article §7.4:
    without it Newton happily returns roots of neighbouring layers);
  * a_F extraction: Feigenbaum sequence (delta*a_k - a_{k-1})/(delta-1) +
    Richardson r^2 = 1/delta^2;
  * INDEPENDENT delta4 check: ratios (a_{k-1}-a_k)/(a_k-a_{k+1}) -> delta4,
    Richardson-accelerated, compared with the OPERATOR eigenvalue from phase 1;
  * alpha4 check: q_k = |x_{2^{k-1}}(a_k)| (closest return before the exact one),
    q_k/q_{k+1} -> ? compared with alpha4 from phase 1;
  * calibration L4 = lim (a_{k+1} - a_F)*delta4^k (Richardson ratio 1/delta4),
    A4 = -1/L4;  self-consistency u(a_k)|_{m=k-1} = -1.

Output: mside_cascade.json
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr, matrix

mp.dps = 80

TAG = "[oid:quartK5] NOT RIGOROUS float pilot, phase 3"
T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


# ---- frozen constants --------------------------------------------------------
FP = json.load(open('quartic_fixedpoint.json', encoding='utf-8'))
ALPHA4 = mpf(FP['alpha4'])
DELTA4 = mpf(FP['delta4'])
K_MAX = 17
CASCADE_TOL = mpf('1e-60')
CASCADE_MAXITER = 100
EXACT_PERIOD_TOL = mpf('1e-45')


# ---- quartic family dynamics --------------------------------------------------
def orbit_x_dx(a, n):
    """x_0 = 0, x_{i+1} = 1 - a*x_i^4 ; dx_{i+1} = -x_i^4 - 4 a x_i^3 dx_i."""
    x = mpf(0)
    dx = mpf(0)
    for _ in range(n):
        dx = -(x ** 4) - 4 * a * x ** 3 * dx
        x = 1 - a * x ** 4
    return x, dx


def orbit_x(a, n):
    x = mpf(0)
    for _ in range(n):
        x = 1 - a * x ** 4
    return x


def exact_period_ok(a, k):
    """x_{2^j}(a) != 0 for all 1 <= j <= k-1  (minimal period exactly 2^k)."""
    for j in range(1, k):
        if fabs(orbit_x(a, 2 ** j)) < EXACT_PERIOD_TOL:
            return False, j
    return True, None


def newton_superstable(c0, n, tol=CASCADE_TOL, max_iter=CASCADE_MAXITER):
    """Damped Newton on x_n(a) = 0 (real a)."""
    a = c0
    x, dx = orbit_x_dx(a, n)
    if fabs(x) < tol:
        return a, 'conv', fabs(x), 0
    for it in range(max_iter):
        if dx == 0:
            return a, 'singular', fabs(x), it
        step = x / dx
        s = step
        improved = False
        for _h in range(40):
            a_try = a - s
            x_try, dx_try = orbit_x_dx(a_try, n)
            if fabs(x_try) < fabs(x):
                a, x, dx = a_try, x_try, dx_try
                improved = True
                break
            s = s / 2
        if not improved:
            return a, 'stalled', fabs(x), it
        if fabs(x) < tol:
            return a, 'conv', fabs(x), it + 1
    return a, 'maxiter', fabs(x), max_iter


def main():
    log(TAG)
    log('delta4 (operator, phase 1) = %s' % nstr(DELTA4, 30))
    log('alpha4 (operator, phase 1) = %s' % nstr(ALPHA4, 30))

    # ---- stage 1: superstable cascade ------------------------------------------
    cascade = {1: mpf(1)}  # a_1 = 1 exactly: 0 -> 1 -> 0, exact period 2
    # a_2 (period 4): float scan for a sign change of x_4(a), then Newton
    def _x4_float(a):
        x = 0.0
        for _ in range(4):
            x = 1.0 - a * x ** 4
        return x
    lo = None
    prev_a, prev_v = 1.001, _x4_float(1.001)
    aa = 1.01
    while aa < 2.0:
        v = _x4_float(aa)
        if prev_v * v < 0:
            lo, hi = prev_a, aa
            break
        prev_a, prev_v = aa, v
        aa += 0.01
    if lo is None:
        raise RuntimeError('no sign change of x_4(a) found in float scan')
    a2, st, xr, _it = newton_superstable(mpf((lo + hi) / 2), 4)
    ok, fj = exact_period_ok(a2, 2)
    assert st == 'conv' and ok, 'a_2 failed'
    cascade[2] = a2
    log('a_1 = 1 (exact), a_2 = %s (bracket [%.2f,%.2f])' % (nstr(a2, 30), lo, hi))
    for k in range(3, K_MAX + 1):
        n = 2 ** k
        extrap = cascade[k - 1] + (cascade[k - 1] - cascade[k - 2]) / DELTA4
        cands = [extrap] + [extrap + off for off in
                            (mpf('-1e-3'), mpf('1e-3'), mpf('-1e-2'), mpf('1e-2'))]
        chosen = None
        for s0 in cands:
            a, st, xr, _it = newton_superstable(s0, n)
            if st == 'conv':
                ok, fj = exact_period_ok(a, k)
                if ok:
                    chosen = a
                    break
        if chosen is None:
            raise RuntimeError('cascade failed at k=%d' % k)
        cascade[k] = chosen
        log('k=%2d  a_k = %s' % (k, nstr(chosen, 34)))

    # ---- a_F extraction ---------------------------------------------------------
    cF_seq = {k: (DELTA4 * cascade[k] - cascade[k - 1]) / (DELTA4 - 1)
              for k in range(2, K_MAX + 1)}
    r2 = 1 / (DELTA4 * DELTA4)
    cF_rich = {k: (cF_seq[k] - r2 * cF_seq[k - 1]) / (1 - r2)
               for k in range(3, K_MAX + 1)}
    a_F = cF_rich[K_MAX]
    a_F_prev = cF_rich[K_MAX - 1]
    log('a_F (Richardson@%d) = %s' % (K_MAX, nstr(a_F, 40)))
    log('a_F drift last->prev: %.3e (expect ~delta4^{-%d})'
        % (float(fabs(a_F - a_F_prev)), K_MAX))

    # ---- independent delta4 check ------------------------------------------------
    ratio_seq = {k: (cascade[k - 1] - cascade[k]) / (cascade[k] - cascade[k + 1])
                 for k in range(2, K_MAX)}
    rL = 1 / DELTA4
    delta_ratio_rich = {k: (ratio_seq[k] - rL * ratio_seq[k - 1]) / (1 - rL)
                        for k in range(3, K_MAX)}
    kd = sorted(delta_ratio_rich)
    d_from_ratios = delta_ratio_rich[kd[-1]]
    log('delta4 from M-side ratios (Richardson): %s' % nstr(d_from_ratios, 25))
    log('  |ratio-delta4 - operator-delta4| = %.3e  (%d agreeing digits)'
        % (float(fabs(d_from_ratios - DELTA4)),
           int(-mp.log(fabs(d_from_ratios - DELTA4), 10))))

    # ---- alpha4 check --------------------------------------------------------------
    q = {k: fabs(orbit_x(cascade[k], 2 ** (k - 1))) for k in range(3, K_MAX + 1)}
    qratio = {k: q[k] / q[k + 1] for k in range(3, K_MAX)}
    qr_rich = {k: (qratio[k] - rL * qratio[k - 1]) / (1 - rL)
               for k in range(4, K_MAX)}
    alpha_from_q = qr_rich[K_MAX - 1]
    log('q-ratio -> %s   (operator alpha4 = %s, diff %.3e)'
        % (nstr(alpha_from_q, 25), nstr(ALPHA4, 25), float(fabs(alpha_from_q - ALPHA4))))

    # ---- calibration L4, A4 ----------------------------------------------------------
    L_seq = {k: (cascade[k + 1] - a_F) * DELTA4 ** k for k in range(6, K_MAX - 1)}
    ks = sorted(L_seq)
    L = (L_seq[ks[-1]] - rL * L_seq[ks[-2]]) / (1 - rL)
    A4 = -1 / L
    log('L4 = %s' % nstr(L, 40))
    log('A4 = -1/L4 = %s' % nstr(A4, 40))

    # ---- self-consistency: u(a_k) at m = k-1 should -> -1 ------------------------------
    u_check = {}
    for k in (10, 12, 14, 16):
        u_check[str(k)] = A4 * (cascade[k] - a_F) * DELTA4 ** (k - 1)
        log('u(a_%d) at m=%d -> %s' % (k, k - 1, nstr(u_check[str(k)], 20)))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'family': 'f_a(x) = 1 - a*x^4, x0=0 (frozen choice)',
                 'mp_dps': mp.dps, 'k_max': K_MAX,
                 'delta4_source': 'phase 1 operator eigenvalue',
                 'elapsed_sec': time.time() - T0},
        'a_F': nstr(a_F, 50),
        'L4': nstr(L, 50),
        'A4': nstr(A4, 50),
        'delta4': nstr(DELTA4, 50),
        'alpha4': nstr(ALPHA4, 50),
        'cascade': {str(k): nstr(cascade[k], 45) for k in cascade},
        'aF_richardson': {str(k): nstr(cF_rich[k], 40) for k in cF_rich},
        'checks': {
            'delta4_from_ratio_richardson': nstr(d_from_ratios, 30),
            'delta4_ratio_vs_operator_absdiff': float(fabs(d_from_ratios - DELTA4)),
            'alpha4_from_qratio_richardson': nstr(alpha_from_q, 30),
            'alpha4_qratio_vs_operator_absdiff': float(fabs(alpha_from_q - ALPHA4)),
            'q_values': {str(k): nstr(q[k], 25) for k in q},
            'u_selfcheck': {k: nstr(v, 20) for k, v in u_check.items()},
            'aF_last_vs_prev': float(fabs(a_F - a_F_prev)),
        },
    }
    with open('mside_cascade.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote mside_cascade.json')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p2_ratecheck.py -- Phase 2 of agent A's lambda2-degree2 retro-check.
[oid:lam2A]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

For the degree-2 family f_c(z) = z^2 + c (critical orbit z0 = 0), and for each
zero u in zeros_deg2.json, find the superstable parameter c*(m) of EXACT
period 2^(m+1) for m = 8, 9, 10, 11 by damped Newton (companion derivative,
NO finite differences) started from the correspondence-law prediction

    c_pred(m) = c_F + u / (A * delta^m)

with c_F, A, delta taken VERBATIM from Investigation-1's own frozen cascade
(section7_mandelbrot/results.json meta) -- NOT recomputed, so the cross-check
uses the exact same constants as the article's own pipeline (BRIEF requirement).

MANDATORY exact-period check (BRIEF trap, article Sec 7.4): min_{1<=j<=m}
|z_{2^j}(c*)| must be bounded away from zero, i.e. c* is NOT a superstable
point of some smaller period 2^j, j<=m dividing into the Newton target.

Then d_u(m) = |u - A*(c*(m)-c_F)*delta^m|, and the ratios d_u(8)/d_u(9),
d_u(9)/d_u(10), d_u(10)/d_u(11) are the observable this whole retro-check is
about (compared against 1/lambda2_essential(2) from phase 1 in phase 3).

dps=80 (BRIEF requires >=60; period 4096 at m=11 with |dc/dz| ~ delta^11 ~ 3e7
needs headroom -- this matches the quartic pilot's own dps=80 choice for an
analogous period-1024 computation).

Output: written incrementally -- rate_deg2.json rewritten after EVERY zero
(not just at the end), so a session death loses at most the current zero.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 80

TAG = "[oid:lam2A] NOT RIGOROUS float pilot (agent A), phase 2 rate study"
T0 = time.time()
PROGRESS_LOG = 'progress.log'
OUT_JSON = 'rate_deg2.json'

M_LEVELS = [8, 9, 10, 11]
NEWTON_TOL = mpf('1e-65')
EXACT_PERIOD_TOL = mpf('1e-45')
NEWTON_MAXITER = 200
ESCAPE = mpf('1e8')


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    with open(PROGRESS_LOG, 'a', encoding='utf-8') as f:
        f.write(time.strftime('[%Y-%m-%d %H:%M:%S] ') + '[phase2-rate] ' + msg + '\n')


# ---- frozen constants (VERBATIM from Investigation-1, NOT recomputed) --------------
RESULTS_JSON = r'../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json'
RJ = json.load(open(RESULTS_JSON, encoding='utf-8'))
META = RJ['meta']


def parse_complex_str(s):
    # format "(-1.4011... + 0.0j)"
    s = s.strip().lstrip('(').rstrip(')')
    s = s.replace('j', '')
    parts = s.split('+')
    re_s = parts[0].strip()
    im_s = parts[1].strip() if len(parts) > 1 else '0'
    return mpc(mpf(re_s), mpf(im_s))


C_F = parse_complex_str(META['c_F_str'])
A_CAL = parse_complex_str(META['A_str'])
DELTA = mpf(META['delta_str'])


# ---- z^2+c dynamics -----------------------------------------------------------------
def orbit_z_dz(c, n):
    z = mpc(0)
    dz = mpc(0)
    for _ in range(n):
        dz = 2 * z * dz + 1
        z = z * z + c
        if fabs(z) > ESCAPE:
            return z, dz
    return z, dz


def orbit_z(c, n):
    z = mpc(0)
    for _ in range(n):
        z = z * z + c
        if fabs(z) > ESCAPE:
            return z
    return z


def exact_period_ok(c_ss, m):
    for j in range(1, m + 1):
        if fabs(orbit_z(c_ss, 2 ** j)) < EXACT_PERIOD_TOL:
            return False, j
    return True, None


def newton(c0, n, tol=NEWTON_TOL, max_iter=NEWTON_MAXITER):
    c = c0
    z, dz = orbit_z_dz(c, n)
    if fabs(z) < tol:
        return c, 'conv', fabs(z), 0
    if fabs(z) > ESCAPE:
        return c, 'escape', fabs(z), 0
    for it in range(max_iter):
        if dz == 0:
            return c, 'singular', fabs(z), it
        step = z / dz
        if fabs(step) < mpf('1e-78') and fabs(z) > tol:
            return c, 'stalled', fabs(z), it
        s = step
        improved = False
        for _h in range(30):
            c_try = c - s
            z2, dz2 = orbit_z_dz(c_try, n)
            if fabs(z2) < fabs(z):
                c, z, dz = c_try, z2, dz2
                improved = True
                break
            s = s / 2
        if not improved:
            return c, 'stalled', fabs(z), it
        if fabs(z) > ESCAPE:
            return c, 'escape', fabs(z), it
        if fabs(z) < tol:
            return c, 'conv', fabs(z), it + 1
    return c, 'maxiter', fabs(z), max_iter


def c_of_u(u, m):
    return C_F + u / (A_CAL * DELTA ** m)


def u_of_c(c, m):
    return A_CAL * (c - C_F) * DELTA ** m


def solve_one(u, m, seed=None, tag=''):
    """Direct Newton at c_pred(m); if that fails, a couple of fallback probes."""
    c_pred = c_of_u(u, m) if seed is None else seed
    n = 2 ** (m + 1)
    c_ss, st, xr, iters = newton(c_pred, n)
    attempts = [{'tag': 'direct' + tag, 'status': st, 'iters': iters, 'final_abs': float(xr)}]
    if st != 'conv':
        # fallback: small deterministic perturbations (fractions of |c_pred - C_F|)
        scale = fabs(c_pred - C_F) * mpf('0.001')
        for dr, di in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]:
            c0b = c_pred + mpc(dr, di) * scale
            c_ss2, st2, xr2, iters2 = newton(c0b, n)
            attempts.append({'tag': 'fallback_%d_%d' % (dr, di), 'status': st2, 'iters': iters2,
                              'final_abs': float(xr2)})
            if st2 == 'conv':
                c_ss, st, xr, iters = c_ss2, st2, xr2, iters2
                break
    if st != 'conv':
        return None, attempts
    ok, fail_j = exact_period_ok(c_ss, m)
    if not ok:
        attempts[-1]['exact_period_fail_j'] = fail_j
        return None, attempts
    return c_ss, attempts


def process_zero(zid, u, is_diagnostic=False, m_levels=None):
    levels = m_levels if m_levels is not None else M_LEVELS
    row = {'id': zid, 'u_re': float(u.real), 'u_im': float(u.imag), 'abs_u': float(fabs(u)),
           'is_diagnostic': is_diagnostic, 'd_u': {}, 'c_ss': {}, 'exact_period_ok': {},
           'attempts_summary': {}}
    for m in levels:
        t_m0 = time.time()
        c_ss, attempts = solve_one(u, m)
        row['attempts_summary'][str(m)] = [a['status'] for a in attempts]
        if c_ss is None:
            log('  %-12s m=%2d  FAILED (%s)' % (zid, m, attempts[-1]['status']))
            continue
        u_act = u_of_c(c_ss, m)
        d_u = fabs(u_act - u)
        row['d_u'][str(m)] = float(d_u)
        row['c_ss'][str(m)] = nstr(c_ss, 40)
        row['exact_period_ok'][str(m)] = True
        log('  %-12s m=%2d  d_u=%.6e  c*=%s  (%.1fs, %d attempt(s))'
            % (zid, m, float(d_u), nstr(c_ss, 20), time.time() - t_m0, len(attempts)))
    return row


def save(rows, done_ids):
    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'not_rigorous': True,
                 'mp_dps': mp.dps, 'm_levels': M_LEVELS,
                 'c_F': nstr(C_F, 50), 'A': nstr(A_CAL, 50), 'delta': nstr(DELTA, 60),
                 'constants_source': 'Investigation-1/09-CODE/section7_mandelbrot/results.json meta (VERBATIM, not recomputed)',
                 'done_ids': done_ids,
                 'elapsed_sec': time.time() - T0},
        'rows': rows,
    }
    with open(OUT_JSON, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)


def main():
    log(TAG)
    log('mp.dps=%d  c_F=%s  A=%s  delta=%s' % (mp.dps, nstr(C_F, 30), nstr(A_CAL, 30), nstr(DELTA, 30)))
    zcat = json.load(open('zeros_deg2.json', encoding='utf-8'))['zeros']

    rows = []
    done_ids = []
    for z in zcat:
        zid = z['id']
        u = mpc(mpf(z['re']), mpf(z['im']))
        is_diag = (z['kind'] == 'z72_true_diagnostic_only')
        levels = [8, 9, 11] if is_diag else M_LEVELS  # skip m=10 for Z72 (trivial/circular there)
        log('=== zero %s (kind=%s, |u|=%.4f) ===' % (zid, z['kind'], float(fabs(u))))
        row = process_zero(zid, u, is_diagnostic=is_diag, m_levels=levels)
        rows.append(row)
        done_ids.append(zid)
        save(rows, done_ids)  # incremental write after EVERY zero
        log('  -> saved rate_deg2.json (%d/%d zeros done)' % (len(rows), len(zcat)))

    log('ALL DONE. %d zeros processed.' % len(rows))


if __name__ == '__main__':
    main()

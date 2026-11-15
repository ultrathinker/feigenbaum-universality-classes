#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p4_match.py -- Phase 4 of the degree-6 float pilot: the correspondence test.
[oid:deg6B]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

HYPOTHESIS UNDER TEST (sextic analogue of article section 7 / quartic phase 4):
    The zeros of H6 (phase 3 catalog, limit-extrapolated positions u_inf) are
    delta6-scaled limit positions of the superstable parameters of the frozen
    sextic family f_a(x) = 1 - a*x^6 near its Feigenbaum point a_F:

        u = A6 * (a - a_F) * delta6^m ,   period(a) = 2^(m+1) exactly.

A6 and a_F are derived in phase 2 (own cascade, nothing imported from phase 3);
delta6 and lambda2 come from phase 1.  After calibration, every catalog zero is
a pure prediction.

FROZEN BEFORE ANY PER-POINT RESULT (all constants justified by GLOBAL scale
arguments available before this run):
  M_MATCH = 9  (period 1024): expected d_u scale d_u(m) ~ C*lambda2^m with
      lambda2 = 0.37107 -> d_u(9) ~ 1e-4..1e-3 across the catalog (quartic C~7.6);
      acceptance radius 0.01 has 10x margin; orbit cost 1024 steps, multiplier
      |A6|*delta6^9 ~ 1.6e9 -- fine at dps=80 (noise floor ~1e-71 vs
      exact-period tol 1e-45, 26 orders margin).
  R_ACCEPT_U = 0.01 (u-units): smallest spacing between distinct catalog zeros
      is 0.1400 -- 14x above the radius.  Same gap logic as article 7.3.
  EDGE_U = 44: catalog zeros with |u| > 44 live within 4 of the fine-scan
      boundary R=48 (Newton rmax 53); quartic lesson (Q57): edge candidates are
      shaky.  Failures in the edge band are REPORTED but do not count against
      the verdict; the verdict set is |u| <= 44.
  EXACT PERIOD CHECK: |x_{2^j}(a)| > 1e-45 for j=1..m -- mandatory (article
      7.4: without it Newton returns roots of neighbouring layers).
  CONTROLS: 6 deterministic points (seed 20260814), |u| in [2,43], Im(u) >= 0,
      at distance >= 0.5 from every catalog zero (and its conjugate);
      they must NOT be accepted.
  Density probes: 20 zeros (15 smallest |u| + 5 largest |u|<=44), 12 jittered
      Newton starts in a 0.4-wide u-box; specificity = (median distance of
      OTHER found roots) / d_u of the accepted match; the matched root itself
      is excluded.
  Rate study: d_u(m) at m = 8,9,10 for 13 zeros (13 smallest |u| excluding the
      calibration zero u=-1); candidate-ratio table guards against coincidences
      (1/lambda3, alpha6^2, delta6/alpha6^2, ...).

Parallel: <= 6 worker processes; the run log is written by the main process
only (phase-3 lesson: module-level open('w') in children clobbered the log).

Output: match_results.json ; run log p4_run.log.
"""
import json
import random
import time
from concurrent.futures import ProcessPoolExecutor

import cmath
import math as _math

from mpmath import mp, mpf, mpc, fabs, cos, sin, pi, nstr, re as mpre, im as mpim, conj

pi_r = _math.pi

mp.dps = 80

TAG = "[oid:deg6B] NOT RIGOROUS float pilot, phase 4 (match, degree 6)"
T0 = time.time()
LOGF = None  # opened only in the main process (children must not truncate)
MAX_WORKERS = 6


def log(msg):
    line = "[%7.1fs] %s" % (time.time() - T0, msg)
    print(line, flush=True)
    if LOGF is not None:
        LOGF.write(line + "\n")
        LOGF.flush()


# ---- frozen constants ---------------------------------------------------------
FP = json.load(open('sextic_fixedpoint.json', encoding='utf-8'))
MS = json.load(open('mside_cascade.json', encoding='utf-8'))
HZ = json.load(open('H6_zeros.json', encoding='utf-8'))

DELTA6 = mpf(FP['delta6'])
LAM2 = mpf(FP['checks']['eigen2_deflated'])
LAM3 = mpf(FP['checks']['eigen3_deflated'])
CERT = [mpf(s) for s in FP['checks']['eigen_certified_list']]
ALPHA6 = mpf(FP['alpha6'])
A_F = mpf(MS['a_F'])
A6 = mpf(MS['A6'])

M_MATCH = 9
PERIOD = 2 ** (M_MATCH + 1)
R_ACCEPT_U = mpf('0.01')
EDGE_U = mpf('44')
EXACT_PERIOD_TOL = mpf('1e-45')
NEWTON_TOL = mpf('1e-55')
NEWTON_MAXITER = 200
ESCAPE = mpf('1e8')
TIER2_FRAC = mpf('1e-3')
TIER3_FRACS = (mpf('0.02'), mpf('0.05'), mpf('0.10'), mpf('0.20'))
TIER3_ANGLES = (0, 45, 90, 135, 180, 225, 270, 315)
N_CONTROLS = 6
CONTROL_SEED = 20260814
DENS_GRID = 12
DENS_SIDE_U = mpf('0.4')


# ---- M-side dynamics (frozen family f_a(x) = 1 - a*x^6) --------------------------
def orbit_x_dx(a, n):
    """x_0=0, x_{k+1} = 1 - a x_k^6, with d/da; early exit on escape."""
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = -(x ** 6) - 6 * a * x ** 5 * dx
        x = 1 - a * x ** 6
        if fabs(x) > ESCAPE:
            return x, dx
    return x, dx


def orbit_x(a, n):
    x = mpc(0)
    for _ in range(n):
        x = 1 - a * x ** 6
        if fabs(x) > ESCAPE:
            return x
    return x


def exact_period_check(a_ss, m):
    for j in range(1, m + 1):
        if fabs(orbit_x(a_ss, 2 ** j)) < EXACT_PERIOD_TOL:
            return False, j
    return True, None


def newton(a0, n, tol=NEWTON_TOL, max_iter=NEWTON_MAXITER):
    """Damped Newton on x_n(a)=0; statuses: conv/escape/stalled/maxiter/singular."""
    a = a0
    x, dx = orbit_x_dx(a, n)
    if fabs(x) < tol:
        return a, 'conv', fabs(x), 0
    if fabs(x) > ESCAPE:
        return a, 'escape', fabs(x), 0
    for it in range(max_iter):
        if dx == 0:
            return a, 'singular', fabs(x), it
        step = x / dx
        if fabs(step) < mpf('1e-75') and fabs(x) > tol:
            return a, 'stalled', fabs(x), it
        s = step
        improved = False
        for _h in range(24):
            a_try = a - s
            x2, dx2 = orbit_x_dx(a_try, n)
            if fabs(x2) < fabs(x):
                a, x, dx = a_try, x2, dx2
                improved = True
                break
            s = s / 2
        if not improved:
            return a, 'stalled', fabs(x), it
        if fabs(x) > ESCAPE:
            return a, 'escape', fabs(x), it
        if fabs(x) < tol:
            return a, 'conv', fabs(x), it + 1
    return a, 'maxiter', fabs(x), max_iter


def a_of_u(u, m=M_MATCH):
    return A_F + u / (A6 * DELTA6 ** m)


def u_of_a(a, m=M_MATCH):
    return A6 * (a - A_F) * DELTA6 ** m


def attempt_point(pid, u_cat_re, u_cat_im):
    """Tiered deterministic search (worker-process safe: pure function).
    Returns (hit_dict_or_None, attempts_list)."""
    u_cat = mpc(mpf(u_cat_re), mpf(u_cat_im))
    a_pred = a_of_u(u_cat)
    r_acc = R_ACCEPT_U / (fabs(A6) * DELTA6 ** M_MATCH)
    attempts = []

    def try_one(tier, tag, a0, m=M_MATCH):
        a_ss, st, xr, iters = newton(a0, 2 ** (m + 1))
        rec = {'tier': tier, 'tag': tag, 'status': st}
        if st == 'conv':
            ok, fail_j = exact_period_check(a_ss, m)
            dist = fabs(a_ss - a_pred)
            accepted = bool(ok and dist < r_acc and m == M_MATCH)
            u_act = u_of_a(a_ss)
            d_u = fabs(u_act - u_cat)
            rec.update({'exact_period_ok': ok, 'fail_j': fail_j,
                        'dist_a': float(dist), 'accepted': accepted,
                        'a_ss_re': nstr(mpre(a_ss), 30),
                        'a_ss_im': float(mpim(a_ss)),
                        'u_actual_re': float(mpre(u_act)),
                        'u_actual_im': float(mpim(u_act)), 'd_u': float(d_u)})
        attempts.append(rec)
        return rec if rec.get('accepted') else None

    hit = try_one(1, 'direct', a_pred)
    if hit:
        return hit, attempts
    grid = TIER2_FRAC * fabs(a_pred - A_F)
    for dr in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dr == 0 and di == 0:
                continue
            hit = try_one(2, 'grid_%d_%d' % (dr, di), a_pred + mpc(dr, di) * grid)
            if hit:
                return hit, attempts
    for frac in TIER3_FRACS:
        radius = frac * fabs(a_pred - A_F)
        for ang in TIER3_ANGLES:
            th = mpf(ang) * pi / 180
            hit = try_one(3, 'fan%s_%d' % (nstr(frac, 3), ang),
                          a_pred + radius * mpc(cos(th), sin(th)))
            if hit:
                return hit, attempts
    # Tier 4: cascade continuation through level m-1
    m1 = M_MATCH - 1
    a_pred1 = A_F + u_cat / (A6 * DELTA6 ** m1)
    a1, st1, _xr, _it = newton(a_pred1, 2 ** (m1 + 1))
    if st1 == 'conv':
        ok1, _fj = exact_period_check(a1, m1)
        if ok1:
            hit = try_one(4, 'cascade_from_m%d' % m1, a1)
            if hit:
                return hit, attempts
    return None, attempts


def control_job(args):
    pid, u_re, u_im = args
    hit, attempts = attempt_point(pid, u_re, u_im)
    rec = {'id': pid, 'u_re': u_re, 'u_im': u_im, 'accepted': hit is not None,
           'n_attempts': len(attempts), 'statuses': [a['status'] for a in attempts]}
    if not hit:
        conv = [a for a in attempts if a['status'] == 'conv']
        if conv:
            dmin = min(a['dist_a'] for a in conv)
            rec['best_dist_u'] = float(dmin * fabs(A6) * DELTA6 ** M_MATCH)
    return rec


def dens_job(args):
    pid, u_re, u_im = args
    u0 = mpc(mpf(u_re), mpf(u_im))
    rnd2 = random.Random(hash(pid) & 0xffffffff)
    found = []
    for _k in range(DENS_GRID):
        du_re = (rnd2.random() - 0.5) * float(DENS_SIDE_U)
        du_im = (rnd2.random() - 0.5) * float(DENS_SIDE_U)
        u_start = u0 + mpc(mpf(du_re), mpf(du_im))
        a_ss, st, _xr, _it = newton(a_of_u(u_start), PERIOD, tol=mpf('1e-50'),
                                    max_iter=60)
        if st == 'conv':
            ok, _fj = exact_period_check(a_ss, M_MATCH)
            if ok:
                found.append(u_of_a(a_ss))
    ded = []
    for v in found:
        if all(fabs(v - w) > mpf('1e-9') for w in ded):
            ded.append(v)
    d_all = [float(fabs(v - u0)) for v in ded]
    d_others = sorted(d for d in d_all if d > 0.02)
    return (pid, [float(mpre(v)) for v in ded], d_others)


def rate_job(args):
    pid, u_re, u_im = args
    u_cat = mpc(mpf(u_re), mpf(u_im))
    row = {'id': pid, 'absu': float(fabs(u_cat)), 'd_u': {}}
    for m in (8, 9, 10):
        a_pred_m = A_F + u_cat / (A6 * DELTA6 ** m)
        a_ss, st, _xr, _it = newton(a_pred_m, 2 ** (m + 1))
        if st == 'conv':
            ok, _fj = exact_period_check(a_ss, m)
            if ok:
                u_act = A6 * (a_ss - A_F) * DELTA6 ** m
                row['d_u'][str(m)] = float(fabs(u_act - u_cat))
    return row


def main():
    global LOGF
    LOGF = open('p4_run.log', 'w', encoding='utf-8')
    log(TAG)
    log('frozen: M_MATCH=%d period=%d R_ACCEPT_U=%s EDGE_U=%s exact_tol=%s '
        'controls=%d seed=%d workers=%d'
        % (M_MATCH, PERIOD, nstr(R_ACCEPT_U, 3), nstr(EDGE_U, 4),
           nstr(EXACT_PERIOD_TOL, 3), N_CONTROLS, CONTROL_SEED, MAX_WORKERS))
    log('a_F = %s' % nstr(A_F, 30))
    log('A6  = %s' % nstr(A6, 30))
    log('delta6 = %s ; lambda2 = %s (1/lambda2 = %.6f)'
        % (nstr(DELTA6, 30), nstr(LAM2, 20), float(1 / LAM2)))
    log('|A6|*delta6^9 = %.4e  (a-space acceptance radius = %.3e)'
        % (float(fabs(A6) * DELTA6 ** M_MATCH),
           float(R_ACCEPT_U / (fabs(A6) * DELTA6 ** M_MATCH))))

    cat = HZ['zeros_limit_extrapolated']
    catalog = []
    for i, z in enumerate(cat):
        u = mpc(mpf(z['re']), mpf(z['im']))
        catalog.append({'id': 'S%03d' % (i + 1), 'u': u,
                        'u8': complex(z['u8_re'], z['u8_im']),
                        'u9': complex(z['u9_re'], z['u9_im']) if z.get('u9_re') else None,
                        'edge': bool(fabs(u) > EDGE_U)})
    n_int = sum(1 for c in catalog if not c['edge'])
    log('catalog: %d zeros (interior |u|<=%s: %d, edge: %d)'
        % (len(catalog), nstr(EDGE_U, 4), n_int, len(catalog) - n_int))

    # ---- controls -----------------------------------------------------------------
    rnd = random.Random(CONTROL_SEED)
    all_u = [c['u'] for c in catalog] + [conj(c['u']) for c in catalog]
    ctrls = []
    tries = 0
    while len(ctrls) < N_CONTROLS:
        tries += 1
        r = rnd.uniform(2.0, 43.0)
        th = rnd.uniform(0.0, 85.0)
        u = mpc(mpf(r) * cos(mpf(th) * pi / 180), mpf(r) * sin(mpf(th) * pi / 180))
        if min(fabs(u - v) for v in all_u) > mpf('0.5'):
            ctrls.append({'id': 'CTRL%d' % (len(ctrls) + 1), 'u': u})
    log('generated %d controls in %d tries (seed %d)'
        % (N_CONTROLS, tries, CONTROL_SEED))

    # ---- STAGE A: catalog cross-check (parallel) ------------------------------------
    log('=== STAGE A: catalog cross-check (%d zeros, m=%d, period %d) ==='
        % (len(catalog), M_MATCH, PERIOD))
    jobs = [(c['id'], float(mpre(c['u'])), float(mpim(c['u']))) for c in catalog]
    results_by_id = {}
    t_A = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for pid, hit, attempts in _run_attempts(ex, jobs):
            results_by_id[pid] = (hit, attempts)
            done += 1
            if done % 40 == 0 or done == len(jobs):
                el = time.time() - t_A
                eta = el / done * (len(jobs) - done)
                log('  stage A %d/%d (%.0fs, ETA %.0fs)'
                    % (done, len(jobs), el, eta))

    results = []
    n_acc = n_acc_int = 0
    dus = []
    for c in catalog:
        hit, attempts = results_by_id[c['id']]
        rec = {'id': c['id'], 'u_cat_re': float(mpre(c['u'])),
               'u_cat_im': float(mpim(c['u'])), 'edge': c['edge'],
               'n_attempts': len(attempts),
               'attempts_status': [a['status'] for a in attempts],
               'accepted': hit is not None}
        if hit:
            n_acc += 1
            if not c['edge']:
                n_acc_int += 1
                dus.append(hit['d_u'])
            rec.update({'d_u': hit['d_u'], 'tier': hit['tier'],
                        'exact_period_ok': hit['exact_period_ok'],
                        'u_actual_re': hit['u_actual_re'],
                        'u_actual_im': hit['u_actual_im'],
                        'a_ss_re': hit['a_ss_re'], 'a_ss_im': hit['a_ss_im']})
        results.append(rec)
    dus_s = sorted(dus)
    med = dus_s[len(dus_s) // 2] if dus_s else None
    n_fail_int = n_int - n_acc_int
    log('catalog accepted: %d/%d ; INTERIOR accepted: %d/%d (fails: %d)'
        % (n_acc, len(catalog), n_acc_int, n_int, n_fail_int))
    if dus_s:
        log('interior d_u: median %s ; min %s ; max %s'
            % ('%.4e' % med, '%.3e' % dus_s[0], '%.3e' % dus_s[-1]))
    fails = [r for r in results if not r['accepted']]
    if fails:
        log('failed points: %s' % ', '.join(
            '%s(%s, |u|=%.2f)' % (r['id'], 'EDGE' if r['edge'] else 'INTERIOR',
                                  abs(complex(r['u_cat_re'], r['u_cat_im'])))
            for r in fails))

    # bonus self-check: the calibration zero u=-1 must land on the M-side cascade a_10
    for c, r in zip(catalog, results):
        if fabs(c['u'] - mpc(-1)) < mpf('1e-9'):
            if r.get('a_ss_re'):
                a10 = mpf(MS['cascade']['10'])
                log('gauge zero u=-1 -> a_ss = %s vs M-side cascade a_10 = %s '
                    '(diff %.2e)' % (r['a_ss_re'], nstr(a10, 25),
                                     float(fabs(mpf(r['a_ss_re']) - a10))))

    # ---- STAGE B: controls ------------------------------------------------------------
    log('=== STAGE B: %d deterministic controls (must NOT be accepted) ===' % N_CONTROLS)
    cjobs = [(c['id'], float(mpre(c['u'])), float(mpim(c['u']))) for c in ctrls]
    ctrl_results = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for rec in ex.map(control_job, cjobs):
            ctrl_results.append(rec)
            log('  %-7s u=%+8.4f%+8.4fi  %s  (conv attempts: %d, best dist_u %s)'
                % (rec['id'], rec['u_re'], rec['u_im'],
                   'ACCEPTED(!)' if rec['accepted'] else 'rejected',
                   sum(1 for s in rec['statuses'] if s == 'conv'),
                   ('%.3e' % rec['best_dist_u']) if rec.get('best_dist_u') else '-'))
    log('controls accepted: %d/%d (expected 0)'
        % (sum(1 for r in ctrl_results if r['accepted']), N_CONTROLS))

    # ---- STAGE C: density/specificity probes ---------------------------------------------
    log('=== STAGE C: density/specificity probes (12 jittered starts per zero) ===')
    by_u = sorted(catalog, key=lambda c: fabs(c['u']))
    probe_set = by_u[:15] + [c for c in by_u if not c['edge']][-5:]
    djobs = [(c['id'], float(mpre(c['u'])), float(mpim(c['u']))) for c in probe_set]
    dens = {}
    t_C = time.time()
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for k, (pid, ded_re, d_others) in enumerate(ex.map(dens_job, djobs)):
            med_other = d_others[len(d_others) // 2] if d_others else None
            rec = next(r for r in results if r['id'] == pid)
            d_u = rec.get('d_u')
            dens[pid] = {
                'n_found': len(ded_re),
                'n_rival_roots': len(d_others),
                'nearest_rival_u': d_others[0] if d_others else None,
                'median_rival_u': med_other,
                'specificity': (med_other / d_u) if (med_other and d_u) else None,
                'd_u': d_u,
            }
            if (k + 1) % 5 == 0:
                el = time.time() - t_C
                eta = el / (k + 1) * (len(djobs) - k - 1)
                log('  stage C %d/%d (%.0fs, ETA %.0fs)' % (k + 1, len(djobs), el, eta))
    spec = sorted(v['specificity'] for v in dens.values() if v['specificity'])
    if spec:
        log('specificity median = %.0fx (min %.0fx)'
            % (spec[len(spec) // 2], spec[0]))

    # ---- STAGE D: convergence-rate study d_u(m), m=8,9,10 ---------------------------------
    log('=== STAGE D: rate study d_u(m), m=8,9,10 (13 zeros, no gauge u=-1) ===')
    rate_ids = [c['id'] for c in by_u
                if fabs(c['u'] - mpc(-1)) > mpf('1e-6')][:13]
    rjobs = [(c['id'], float(mpre(c['u'])), float(mpim(c['u'])))
             for c in catalog if c['id'] in rate_ids]
    rate_recs = []
    t_D = time.time()
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for k, row in enumerate(ex.map(rate_job, rjobs)):
            r89 = (row['d_u']['8'] / row['d_u']['9']
                   if '8' in row['d_u'] and '9' in row['d_u'] and row['d_u']['9']
                   else None)
            r910 = (row['d_u']['9'] / row['d_u']['10']
                    if '9' in row['d_u'] and '10' in row['d_u'] and row['d_u']['10']
                    else None)
            row['ratio_89'] = r89
            row['ratio_910'] = r910
            rate_recs.append(row)
            log('  %-6s |u|=%6.3f  d8=%.2e d9=%.2e d10=%.2e  d8/d9=%s d9/d10=%s'
                % (row['id'], row['absu'], row['d_u'].get('8', float('nan')),
                   row['d_u'].get('9', float('nan')), row['d_u'].get('10', float('nan')),
                   ('%.3f' % r89) if r89 else '-', ('%.3f' % r910) if r910 else '-'))
    r89s = sorted(r['ratio_89'] for r in rate_recs if r['ratio_89'])
    r910s = sorted(r['ratio_910'] for r in rate_recs if r['ratio_910'])

    # candidate-ratio table
    inv = lambda x: float(1 / x)
    cands_table = [
        ('1/lambda2', inv(LAM2)),
        ('1/|lambda3|', inv(fabs(CERT[1]) if len(CERT) > 1 else mpf('1'))),
        ('1/|lambda4|', inv(fabs(CERT[2]) if len(CERT) > 2 else mpf('1'))),
        ('1/lambda5', inv(fabs(CERT[3]) if len(CERT) > 3 else mpf('1'))),
        ('alpha6^2', float(ALPHA6 ** 2)),
        ('delta6/alpha6^2', float(DELTA6 / ALPHA6 ** 2)),
        ('delta6', float(DELTA6)),
        ('1/delta6', float(1 / DELTA6)),
    ]
    m89 = r89s[len(r89s) // 2] if r89s else None
    m910 = r910s[len(r910s) // 2] if r910s else None
    log('--- candidate-ratio table (observed medians: d8/d9=%s, d9/d10=%s) ---'
        % (('%.4f' % m89) if m89 else '-', ('%.4f' % m910) if m910 else '-'))
    for name, val in cands_table:
        log('   %-18s %12.6f' % (name, val))
    if r89s:
        log('d_u(8)/d_u(9): median %.4f  [model 1/lambda2 = %.4f]'
            % (m89, float(1 / LAM2)))
    if r910s:
        log('d_u(9)/d_u(10): median %.4f' % m910)

    # ---- verdict data -----------------------------------------------------------------------
    verdict = {
        'catalog_total': len(catalog),
        'interior_total': n_int,
        'interior_accepted': n_acc_int,
        'interior_failed': n_fail_int,
        'edge_total': len(catalog) - n_int,
        'edge_accepted': n_acc - n_acc_int,
        'd_u_median': med,
        'd_u_min': dus_s[0] if dus_s else None,
        'd_u_max': dus_s[-1] if dus_s else None,
        'controls_accepted': sum(1 for r in ctrl_results if r['accepted']),
        'controls_total': N_CONTROLS,
        'specificity_median': spec[len(spec) // 2] if spec else None,
        'specificity_min': spec[0] if spec else None,
        'rate_median_89': m89,
        'rate_median_910': m910,
        'rate_scatter_89': (r89s[-1] - r89s[0]) if r89s else None,
        'model_1_over_lambda2': float(1 / LAM2),
    }
    log('=== VERDICT DATA ===')
    for k, v in verdict.items():
        log('  %s = %s' % (k, v))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'frozen': {'M_MATCH': M_MATCH, 'PERIOD': PERIOD,
                            'R_ACCEPT_U': float(R_ACCEPT_U),
                            'EDGE_U': float(EDGE_U),
                            'EXACT_PERIOD_TOL': float(EXACT_PERIOD_TOL),
                            'N_CONTROLS': N_CONTROLS, 'CONTROL_SEED': CONTROL_SEED,
                            'DENS_GRID': DENS_GRID, 'DENS_SIDE_U': float(DENS_SIDE_U),
                            'family': '1 - a*x^6', 'workers': MAX_WORKERS},
                 'elapsed_sec': time.time() - T0},
        'a_F': MS['a_F'], 'A6': MS['A6'], 'delta6': FP['delta6'],
        'lambda2': FP['checks']['eigen2_deflated'],
        'catalog_results': results,
        'controls': ctrl_results,
        'density': dens,
        'rate_study': rate_recs,
        'candidate_ratios': cands_table,
        'verdict': verdict,
    }
    with open('match_results.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote match_results.json')
    log('DONE phase 4 (match), total elapsed %.1fs' % (time.time() - T0))


def _run_attempts(ex, jobs):
    """Map attempt_point over jobs preserving order; yields (pid, hit, attempts)."""
    from concurrent.futures import as_completed
    futs = {ex.submit(attempt_point, pid, re_, im_): pid for pid, re_, im_ in jobs}
    for fut in as_completed(futs):
        pid = futs[fut]
        hit, attempts = fut.result()
        yield pid, hit, attempts


if __name__ == '__main__':
    main()

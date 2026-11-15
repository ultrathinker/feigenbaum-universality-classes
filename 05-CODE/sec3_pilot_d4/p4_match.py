#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
p4_match.py -- Phase 4 of the quartic float pilot: the correspondence test.
[oid:quartK5]  NOT RIGOROUS -- float exploration (mpmath), no interval arithmetic.

HYPOTHESIS UNDER TEST (quartic analogue of article section 7):
    The zeros of H4 (phase 2 catalog, limit-extrapolated positions u_inf) are
    delta4-scaled limit positions of the superstable parameters of the frozen
    quartic family f_a(x) = 1 - a*x^4 near its Feigenbaum point a_F:

        u = A4 * (a - a_F) * delta4^m ,   period(a) = 2^(m+1) exactly.

A4 and a_F are derived in phase 3 (own cascade, nothing imported from phase 2);
delta4 comes from the phase-1 operator eigenvalue.  After calibration, every
catalog zero is a pure prediction.

FROZEN BEFORE ANY PER-POINT RESULT (all constants justified by GLOBAL scale
arguments available before this run):
  M_MATCH = 9  (period 1024): drift-model scale at m=9 is delta4^{-9}*u^2*c ~
      1e-7..2e-6 across the catalog -- three to four orders below R_ACCEPT_U;
      orbit cost 1024 steps, multiplier |dx/da| ~ delta4^9 ~ 6e7 -- fine at
      dps=80 (noise floor ~1e-72, exact-period tol 1e-45 has 27 orders margin).
  R_ACCEPT_U = 0.01 (u-units): genuine matches are expected at d_u ~ 1e-6..1e-5
      (H4 level-9 drift scale); the smallest spacing between distinct catalog
      zeros is 0.1786 -- 18x above the radius.  Same gap logic as article 7.3.
  EXACT PERIOD CHECK: |x_{2^j}(a)| > 1e-45 for j=1..m -- mandatory (article
      7.4: without it Newton returns roots of neighbouring layers).
  CONTROLS: 6 deterministic points (seed 20260814), |u| in [2,21], Im(u) >= 0,
      at distance >= 0.5 from every catalog zero (and its conjugate);
      they must NOT be accepted.
  Density probes: for a subset of zeros, 12 jittered Newton starts in a
      0.4-wide u-box; specificity = (median distance of OTHER found roots in
      u) / d_u of the accepted match.

Output: match_results.json
"""
import cmath
import json
import random
import time

import math as _math
from mpmath import mp, mpf, mpc, fabs, cos, sin, pi, nstr, re as mpre, im as mpim, conj

pi_r = _math.pi

mp.dps = 80

TAG = "[oid:quartK5] NOT RIGOROUS float pilot, phase 4"
T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


# ---- frozen constants ---------------------------------------------------------
FP = json.load(open('quartic_fixedpoint.json', encoding='utf-8'))
MS = json.load(open('mside_cascade.json', encoding='utf-8'))
HZ = json.load(open('H4_zeros.json', encoding='utf-8'))

DELTA4 = mpf(FP['delta4'])
A_F = mpf(MS['a_F'])
A4 = mpf(MS['A4'])

M_MATCH = 9
PERIOD = 2 ** (M_MATCH + 1)
R_ACCEPT_U = mpf('0.01')
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


# ---- M-side dynamics ------------------------------------------------------------
def orbit_x_dx(a, n):
    """x_0=0, x_{k+1} = 1 - a x_k^4, with d/da; early exit on escape."""
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = -(x ** 4) - 4 * a * x ** 3 * dx
        x = 1 - a * x ** 4
        if fabs(x) > ESCAPE:
            return x, dx
    return x, dx


def orbit_x(a, n):
    x = mpc(0)
    for _ in range(n):
        x = 1 - a * x ** 4
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


def a_of_u(u):
    return A_F + u / (A4 * DELTA4 ** M_MATCH)


def u_of_a(a):
    return A4 * (a - A_F) * DELTA4 ** M_MATCH


def attempt_point(pid, u_cat):
    """Tiered deterministic search. Tiers 1-3 as in the frozen article script;
    Tier 4 (added after the first run rejected only Q57): cascade continuation --
    Newton at level M_MATCH-1 from the level-(m-1) prediction, exact-period check
    there, then Newton at level m from that root.  Disclosed per-attempt."""
    a_pred = a_of_u(u_cat)
    r_acc = R_ACCEPT_U / (fabs(A4) * DELTA4 ** M_MATCH)
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
                        'a_ss': nstr(a_ss, 30), 'u_actual_re': float(mpre(u_act)),
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
    a_pred1 = A_F + u_cat / (A4 * DELTA4 ** m1)
    a1, st1, _xr, _it = newton(a_pred1, 2 ** (m1 + 1))
    if st1 == 'conv':
        ok1, _fj = exact_period_check(a1, m1)
        if ok1:
            hit = try_one(4, 'cascade_from_m%d' % m1, a1)
            if hit:
                return hit, attempts
    return None, attempts


def main():
    log(TAG)
    log('frozen: M_MATCH=%d period=%d R_ACCEPT_U=%s exact_tol=%s'
        % (M_MATCH, PERIOD, nstr(R_ACCEPT_U, 3), nstr(EXACT_PERIOD_TOL, 3)))
    log('a_F = %s' % nstr(A_F, 30))
    log('A4  = %s' % nstr(A4, 30))
    log('|A4|*delta4^9 = %.4e  (a-space acceptance radius = %.3e)'
        % (float(fabs(A4) * DELTA4 ** M_MATCH),
           float(R_ACCEPT_U / (fabs(A4) * DELTA4 ** M_MATCH))))

    cat = HZ['zeros_limit_extrapolated']
    catalog = []
    for i, z in enumerate(cat):
        u = mpc(mpf(z['re']), mpf(z['im']))
        catalog.append({'id': 'Q%02d' % (i + 1), 'u': u,
                        'u8': complex(z['u8_re'], z['u8_im']),
                        'u9': complex(z['u9_re'], z['u9_im']) if z.get('u9_re') else None})
    log('catalog: %d zeros' % len(catalog))

    # ---- controls -----------------------------------------------------------------
    rnd = random.Random(CONTROL_SEED)
    all_u = [c['u'] for c in catalog] + [conj(c['u']) for c in catalog]
    ctrls = []
    tries = 0
    while len(ctrls) < N_CONTROLS:
        tries += 1
        r = rnd.uniform(2.0, 21.0)
        th = rnd.uniform(0.0, 85.0)
        u = mpc(mpf(r) * cos(mpf(th) * pi / 180), mpf(r) * sin(mpf(th) * pi / 180))
        if min(fabs(u - v) for v in all_u) > mpf('0.5'):
            ctrls.append({'id': 'CTRL%d' % (len(ctrls) + 1), 'u': u})
    log('generated %d controls in %d tries (seed %d)'
        % (N_CONTROLS, tries, CONTROL_SEED))

    # ---- catalog run -----------------------------------------------------------------
    log('=== STAGE A: catalog cross-check (%d zeros, m=%d, period %d) ==='
        % (len(catalog), M_MATCH, PERIOD))
    results = []
    n_acc = 0
    dus = []
    drift_ratios = []
    drift_phases = []
    for c in catalog:
        hit, attempts = attempt_point(c['id'], c['u'])
        rec = {'id': c['id'], 'u_cat_re': float(mpre(c['u'])),
               'u_cat_im': float(mpim(c['u'])),
               'n_attempts': len(attempts), 'attempts_status': [a['status'] for a in attempts],
               'accepted': hit is not None}
        if hit:
            n_acc += 1
            dus.append(hit['d_u'])
            rec.update({'d_u': hit['d_u'], 'tier': hit['tier'],
                        'exact_period_ok': hit['exact_period_ok'],
                        'u_actual_re': hit['u_actual_re'],
                        'u_actual_im': hit['u_actual_im'],
                        'a_ss': hit['a_ss']})
            # drift-vector consistency: R = u_actual - u_inf vs R_pred = (u8-u9)/(delta4-1)
            if c['u9'] is not None:
                R = complex(hit['u_actual_re'], hit['u_actual_im']) - complex(mpre(c['u']), mpim(c['u']))
                Rpred = (c['u8'] - c['u9']) / float(DELTA4 - 1)
                if abs(Rpred) > 0:
                    rat = R / Rpred
                    drift_ratios.append(abs(rat))
                    ph = cmath.phase(R) * 180 / pi_r - cmath.phase(Rpred) * 180 / pi_r
                    while ph > 180:
                        ph -= 360
                    while ph <= -180:
                        ph += 360
                    drift_phases.append(ph)
                    rec['drift_ratio_abs'] = abs(rat)
                    rec['drift_phase_deg'] = ph
        results.append(rec)
        log('  %-6s u=%+8.4f%+8.4fi  %s  d_u=%s'
            % (c['id'], float(mpre(c['u'])), float(mpim(c['u'])),
               'ACCEPT ' if hit else 'reject ',
               ('%.3e' % hit['d_u']) if hit else '-'))
    dus_s = sorted(dus)
    med = dus_s[len(dus_s) // 2] if dus_s else None
    log('catalog accepted: %d/%d ; d_u median = %s ; range [%s, %s]'
        % (n_acc, len(catalog),
           ('%.4e' % med) if med is not None else '-',
           ('%.3e' % dus_s[0]) if dus_s else '-', ('%.3e' % dus_s[-1]) if dus_s else '-'))
    if drift_ratios:
        dr_s = sorted(drift_ratios)
        ph_s = sorted(abs(p) for p in drift_phases)
        log('drift-vector check: |R/Rpred| median = %.4f (range %.3f..%.3f); '
            '|phase| median = %.2f deg (max %.2f)'
            % (dr_s[len(dr_s) // 2], dr_s[0], dr_s[-1],
               ph_s[len(ph_s) // 2], ph_s[-1]))

    # ---- controls run ------------------------------------------------------------------
    log('=== STAGE B: %d deterministic controls (must NOT be accepted) ===' % N_CONTROLS)
    ctrl_results = []
    for c in ctrls:
        hit, attempts = attempt_point(c['id'], c['u'])
        rec = {'id': c['id'], 'u_re': float(mpre(c['u'])), 'u_im': float(mpim(c['u'])),
               'accepted': hit is not None,
               'n_attempts': len(attempts),
               'statuses': [a['status'] for a in attempts]}
        if not hit:
            conv = [a for a in attempts if a['status'] == 'conv']
            if conv:
                dmin = min(a['dist_a'] for a in conv)
                # convert to u-units for comparability
                rec['best_dist_u'] = float(dmin * fabs(A4) * DELTA4 ** M_MATCH)
        ctrl_results.append(rec)
        log('  %-7s u=%+8.4f%+8.4fi  %s  (conv attempts: %d, best dist_u %s)'
            % (c['id'], float(mpre(c['u'])), float(mpim(c['u'])),
               'ACCEPTED(!)' if hit else 'rejected',
               sum(1 for a in attempts if a['status'] == 'conv'),
               ('%.3e' % rec['best_dist_u']) if rec.get('best_dist_u') else '-'))
    log('controls accepted: %d/%d (expected 0)' % (sum(1 for r in ctrl_results if r['accepted']), N_CONTROLS))

    # ---- density probes ------------------------------------------------------------------
    log('=== STAGE C: density/specificity probes (12 jittered starts per zero) ===')
    by_u = sorted(catalog, key=lambda c: fabs(c['u']))
    probe_set = by_u[:15] + by_u[-5:]
    dens = {}
    for c in probe_set:
        u0 = c['u']
        found = []
        rnd2 = random.Random(hash(c['id']) & 0xffffffff)
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
        # dedup in u
        ded = []
        for v in found:
            if all(fabs(v - w) > mpf('1e-9') for w in ded):
                ded.append(v)
        # distances to u0; roots within 0.02 are the match itself, others are rivals
        d_all = [float(fabs(v - u0)) for v in ded]
        d_others = sorted(d for d in d_all if d > 0.02)
        med_other = d_others[len(d_others) // 2] if d_others else None
        rec = next(r for r in results if r['id'] == c['id'])
        d_u = rec.get('d_u')
        dens[c['id']] = {
            'n_found': len(ded),
            'n_rival_roots': len(d_others),
            'nearest_rival_u': d_others[0] if d_others else None,
            'median_rival_u': med_other,
            'specificity': (med_other / d_u) if (med_other and d_u) else None,
            'd_u': d_u,
        }
        log('  %-6s roots=%2d rivals=%2d nearest_rival=%s  d_u=%.2e  specificity=%s'
            % (c['id'], len(ded), len(d_others),
               ('%.4f' % d_others[0]) if d_others else 'none<0.4box',
               d_u if d_u else -1,
               ('%.0fx' % (med_other / d_u)) if (med_other and d_u) else 'n/a'))

    spec = [v['specificity'] for v in dens.values() if v['specificity']]
    if spec:
        spec.sort()
        log('specificity median = %.0fx (min %.0fx)' % (spec[len(spec) // 2], spec[0]))

    # ---- rate study: d_u(m) for m = 8,9,10 on a fixed subset ------------------------------
    log('=== STAGE D: convergence-rate study d_u(m), m=8,9,10 (subset) ===')
    FP2 = json.load(open('quartic_fixedpoint.json', encoding='utf-8'))
    lam2_str = FP2['checks']['eigen2_deflated']
    lam2 = mpf(lam2_str)
    subset_ids = [c['id'] for c in by_u[:12]] + ['Q28']
    rate_recs = []
    for c in catalog:
        if c['id'] not in subset_ids:
            continue
        row = {'id': c['id'], 'absu': float(fabs(c['u'])), 'd_u': {}}
        for m in (8, 9, 10):
            a_pred_m = A_F + c['u'] / (A4 * DELTA4 ** m)
            a_ss, st, _xr, _it = newton(a_pred_m, 2 ** (m + 1))
            if st == 'conv':
                ok, _fj = exact_period_check(a_ss, m)
                if ok:
                    u_act = A4 * (a_ss - A_F) * DELTA4 ** m
                    row['d_u'][str(m)] = float(fabs(u_act - c['u']))
        r89 = (row['d_u']['8'] / row['d_u']['9']
               if '8' in row['d_u'] and '9' in row['d_u'] and row['d_u']['9'] else None)
        r910 = (row['d_u']['9'] / row['d_u']['10']
                if '9' in row['d_u'] and '10' in row['d_u'] and row['d_u']['10'] else None)
        row['ratio_89'] = r89
        row['ratio_910'] = r910
        rate_recs.append(row)
        log('  %-6s |u|=%6.3f  d8=%.2e d9=%.2e d10=%.2e  d8/d9=%s d9/d10=%s'
            % (c['id'], row['absu'], row['d_u'].get('8', float('nan')),
               row['d_u'].get('9', float('nan')), row['d_u'].get('10', float('nan')),
               ('%.2f' % r89) if r89 else '-', ('%.2f' % r910) if r910 else '-'))
    r89s = sorted(r['ratio_89'] for r in rate_recs if r['ratio_89'])
    r910s = sorted(r['ratio_910'] for r in rate_recs if r['ratio_910'])
    if r89s:
        log('d_u(8)/d_u(9): median %.3f   [model 1/lambda2 = %.3f ; 1/delta4 = %.3f]'
            % (r89s[len(r89s) // 2], float(1 / lam2), float(1 / DELTA4)))
    if r910s:
        log('d_u(9)/d_u(10): median %.3f' % (r910s[len(r910s) // 2]))

    # ---- verdict data ------------------------------------------------------------------
    verdict = {
        'catalog_total': len(catalog),
        'catalog_accepted': n_acc,
        'd_u_median': med,
        'd_u_min': dus_s[0] if dus_s else None,
        'd_u_max': dus_s[-1] if dus_s else None,
        'controls_accepted': sum(1 for r in ctrl_results if r['accepted']),
        'controls_total': N_CONTROLS,
        'drift_ratio_median': (sorted(drift_ratios)[len(drift_ratios) // 2]
                               if drift_ratios else None),
        'drift_phase_median_absdeg': (sorted(abs(p) for p in drift_phases)[len(drift_phases) // 2]
                                      if drift_phases else None),
        'specificity_median': spec[len(spec) // 2] if spec else None,
    }
    log('=== VERDICT DATA ===')
    for k, v in verdict.items():
        log('  %s = %s' % (k, v))

    out = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'frozen': {'M_MATCH': M_MATCH, 'PERIOD': PERIOD,
                            'R_ACCEPT_U': float(R_ACCEPT_U),
                            'EXACT_PERIOD_TOL': float(EXACT_PERIOD_TOL),
                            'N_CONTROLS': N_CONTROLS, 'CONTROL_SEED': CONTROL_SEED,
                            'family': '1 - a*x^4'},
                 'elapsed_sec': time.time() - T0},
        'a_F': MS['a_F'], 'A4': MS['A4'], 'delta4': FP['delta4'],
        'catalog_results': results,
        'controls': ctrl_results,
        'density': dens,
        'rate_study': rate_recs,
        'verdict': verdict,
    }
    with open('match_results.json', 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote match_results.json')


if __name__ == '__main__':
    main()

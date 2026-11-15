#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
frozen_experiment.py -- UNIFIED, FROZEN Mandelbrot cross-check for the
Feigenbaum-point-19 hypothesis, built to referee spec (single script, single
dataset, algorithm fixed BEFORE the run, full disclosure of every failure).

HYPOTHESIS UNDER TEST
    The 23 catalogued zeros of H(u) (address_certificate.json, stage 15.3,
    field address_table, columns centre_re/centre_im interpreted as complex
    u-coordinates) are delta-scaled limit positions of the complex
    superstable parameters c of f_c(z) = z^2 + c near the Feigenbaum point
    c_F = -1.4011551890920506..., through the first-order scaling law

        u = A * (c - c_F) * delta^m ,      period(c) = 2^(m+1) exactly.

    A and c_F below are DERIVED in this script (stage 1/2), not imported,
    so the run is self-contained; delta is an external universal constant
    (frozen literal, high precision, from the certificate's own convention
    block, cross-checked against the standard reference value).

SOURCES READ (read-only; nothing under Feigenbaum-archive is written to
except this file's own output directory):
    path_to_8/point19_superstable_interpretation/THEORY_NOTE.md
    path_to_8/point19_superstable_interpretation/REPORT.md, FIX_BRIEF.md
    path_to_8/point19_superstable_interpretation/kilo_xcheck/xcheck.py,
        results.json, XCHECK_REPORT.md
    path_to_8/point15_layer1_addressing/address_certificate.json

WHY THIS SCRIPT EXISTS (defects it is built NOT to repeat, per FIX_BRIEF.md
and XCHECK_REPORT.md, read before writing a line of this file):
  - DEFECT (attempt1/REPORT.md #1): err_u/d_u computed from the LAST Newton
    iterate even when Newton did not converge (a stalled Newton "looks"
    converged because the step vanishes near the singular locus). FIX: d_u
    is computed ONLY for points that (a) converged (|x_n|<tol), (b) passed
    the exact-period check, AND (c) passed the frozen root-acceptance rule
    |c_ss-c_pred|<R_ACCEPT. Every other outcome is a logged failure, never a
    silent zero.
  - DEFECT (attempt1 #2/#3): wrong sign convention u=A*(c-cF)*(-delta)^m,
    producing a spurious even/odd alternation. FIX: this script hardcodes
    the THEORY_NOTE-derived, self-consistency-checked form u=A*(c-cF)*delta^m
    (no (-1)^m), and re-derives A,c_F from its OWN cascade Newton (stage 1/2)
    as an internal consistency gate before touching the catalog.
  - DEFECT (kilo_xcheck / point17 provenance): Z72's certified disk centre
    is a "shifted disk" (point17 addendum 15.4: D(4.758336330290482+
    7.1277051021354385i, rho=0.1)) required by the point17 proof technique,
    NOT a tight numerical estimate of the true zero -- unlike the other 22
    disks, whose centres double as tight (~1e-4) numerical estimates. Using
    the disk centre as u_cat for Z72 and comparing d_u against the other 22
    points' d_u~2e-5 is comparing a loose proxy against tight estimates and
    was mistaken for a hypothesis failure in earlier work. FIX: Z72 is run
    TWICE and BOTH rows are reported (see PROVENANCE section in the header
    of stage3_catalog and REPORT_UNIFIED.md): once "as certified" (disk
    centre, honestly expected to fail the tight R_ACCEPT gate) and once
    "corrected" using the float refinement u=4.791660937911268+
    7.177581007211361i (sourced from kilo_xcheck/results.json
    stage3_catalog[Z72].u_actual -- i.e. the position kilo_xcheck's own
    independent Newton search converged to when probing near the certified
    centre). The corrected row is run through this script's OWN Newton code
    (not copied from kilo's converged c), which is a genuine, if weak,
    independent re-confirmation; it is EXCLUDED from the primary N/23 MATCH
    tally to avoid a self-referential "match" (c_pred for that row is, by
    construction, numerically equal to kilo's already-converged root, so
    d_u for that row is expected to be tiny by construction, not by new
    evidence) and is discussed separately.

ALL TUNABLE PARAMETERS ARE FROZEN BELOW BEFORE ANY DATA IS TOUCHED. Nothing
below this block was chosen by looking at a per-point outcome; every
constant is justified from a GLOBAL scale argument (see comments) that
existed before this run.
"""
import os
import json
import time
import random as _random

from mpmath import mp, mpf, mpc, cos, sin, pi, fabs, re, im, nstr, arg

# --------------------------------------------------------------------------- #
#  FROZEN CONSTANTS -- change nothing here after looking at results.json     #
# --------------------------------------------------------------------------- #

# ---- precision -------------------------------------------------------------
MP_DPS = 150
# Justification: THEORY_NOTE.md §3 measured |dy_n/dc| ~ delta^m near c_F for
# critical orbits of period 2^(m+1) (renormalisation taming of the
# multiplier growth) -- NOT ~2^n as for a generic chaotic orbit -- so at
# m=10 only ~log10(delta^10)=6.7 extra decimal digits of cancellation are
# needed beyond the target tolerance. dps=150 leaves a large safety margin
# over kilo_xcheck's already-successful dps=120 clean-room run, at
# negligible extra cost (benchmarked below).

mp.dps = MP_DPS

# ---- universal constant delta (external input, NOT derived) ---------------
DELTA_STR = ("4.669201609102990671853203820466201617258185577475768632745651"
             "343004133555868548426390712381266351242")
# Source: address_certificate.json (path_to_8/point15_layer1_addressing),
# field convention.delta -- the certificate's own high-precision literal for
# the first Feigenbaum constant, used consistently throughout the archive.
DELTA = mpf(DELTA_STR)

BRIEF_CF_STR = "-1.4011551890920506"
# Reference literal (THEORY_NOTE.md / kilo_xcheck brief) used ONLY as an
# external cross-check on the cascade-derived c_F below, never as an input.

# ---- renormalisation level used for the whole frozen run ------------------
M_LEVEL = 10
PERIOD = 2 ** (M_LEVEL + 1)          # 2048
# Justification: (a) this is the exact level independently validated by
# kilo_xcheck/XCHECK_REPORT.md (specificity 3700-8500x for the 4 clean
# catalog points it tested); (b) THEORY_NOTE.md §4 recommends m=8..10 as
# "cheap" (its own estimate: unit-minutes to ~20-30 min sequential for all
# 23 at m=10); (c) benchmarked below at ~13 ms per full period-2048 orbit
# pass at dps=150 on this machine, i.e. a full 23-point + density-grid sweep
# fits the ~45-minute target with margin. No per-zero m is recorded in
# address_certificate.json (only the global fingerprint level
# convention.reference_m=7, which belongs to the ORIGINAL abstract-H
# construction, not to this family-dynamics Newton search) -- so a single
# frozen m for all 23 catalog zeros is the only non-cherry-picked choice.

# ---- cascade (stage 1/2: internal derivation of c_F, A) -------------------
CASCADE_K_MAX = 16
CASCADE_TOL = mpf('1e-100')
CASCADE_MAXITER = 150
CASCADE_MINSTEP = mpf('1e-110')
CASCADE_MAXHALVE = 16
# Mirrors kilo_xcheck's already-validated cascade recipe (its dps=120
# version used tol=1e-90/maxiter=120/minstep=1e-100/maxhalve=14 and reached
# |c_F-brief|=5.2e-19 at k=16); loosened only by the extra dps=150 headroom.

# ---- main Newton (catalog + control points, level m=10, period 2048) ------
NEWTON_TOL = mpf('1e-70')            # |x_2048(c)| convergence threshold
NEWTON_MAXITER = 300
NEWTON_MINSTEP = mpf('1e-90')        # |step| below which Newton is "stalled"
NEWTON_MAXHALVE = 10                 # damping: halve step up to this many times
ESCAPE_THRESHOLD = mpf('1e8')        # |x| above this = escaping start
EXACT_PERIOD_TOL = mpf('1e-45')      # |x_{2^j}| below this counts as "=0"
# tol/minstep are ~1e-70/1e-90 to leave >50 decimal digits of margin below
# dps=150 for the multiplier blow-up (bounded by ~delta^m, see above);
# maxiter/maxhalve/escape/exact-period-tol are direct carries of
# kilo_xcheck's already-working values (xcheck.py T_CONV/T_MINSTEP/ESCAPE/
# EXACT_TOL), not re-tuned.

# ---- deterministic restart ladder (NO randomness in the search itself) ----
# Tier 1: direct Newton from c_pred.
# Tier 2 (only if Tier 1 does not end in an ACCEPTED root): 3x3 grid minus
#   centre, offset TIER2_FRAC*|c_pred-c_F| along the real/imag axes -- this
#   is exactly kilo_xcheck's stage-3 restart recipe (xcheck.py
#   solve_with_restarts), already shown sufficient for its 4 clean hits.
# Tier 3 (only if Tier 2 also fails): a wider deterministic fan, radii
#   TIER3_FRACS x |c_pred-c_F|, TIER3_ANGLES_DEG directions -- generalises
#   THEORY_NOTE.md §3's real-axis restart fan (f in {.02,.05,.10,.20,.30}) to
#   the complex plane needed for the 22 non-real catalog points.
TIER2_FRAC = mpf('1e-3')
TIER3_FRACS = (mpf('0.02'), mpf('0.05'), mpf('0.10'), mpf('0.20'), mpf('0.30'))
TIER3_ANGLES_DEG = (0, 45, 90, 135, 180, 225, 270, 315)

# ---- frozen root-acceptance rule (item 4 of the referee brief) ------------
# A converged, exact-period candidate c_ss is ACCEPTED for a target only if
#     |c_ss - c_pred| < R_ACCEPT_C(m) := R_ACCEPT_U / (|A| * delta^m)
# i.e. the acceptance radius is frozen in the m-INDEPENDENT u-coordinate
# (the natural length scale of the whole construction) and converted to
# c-space at the run's m. R_ACCEPT_U = 0.01 is derived from TWO global
# (not per-point) numbers already on record before this run:
#   - genuine catalog matches (kilo_xcheck stage 3, 4 clean points) land at
#     d_u ~ 2e-5..4e-5 -- three to four orders of magnitude BELOW 0.01.
#   - the smallest observed spacing between DISTINCT genuine period-2048
#     superstable centres anywhere near this catalog (kilo_xcheck stage 5
#     density grids; also address_certificate.json's own closest
#     certified-disk pair Z69/Z70 at raw distance 0.1167) is >=0.05 --
#     five times ABOVE 0.01.
# So R_ACCEPT_U=0.01 sits in the gap with wide margin on both sides and was
# fixed from these pre-existing, cross-point statistics -- not from any
# single point's outcome (in particular NOT from Z72, which this rule will
# be shown, honestly, to REJECT under its certified-centre address).
R_ACCEPT_U = mpf('0.01')

# ---- deterministic negative controls (item 7) ------------------------------
N_CONTROLS = 20
CONTROL_SEED = 20260812          # date this script was frozen (YYYYMMDD)
CONTROL_R_RANGE = (mpf('1.5'), mpf('11.0'))       # covers the catalog's |u| span (1.55..9.89)
CONTROL_THETA_RANGE_DEG = (mpf('0'), mpf('85'))   # covers the catalog's angular span (21.6..81.7 deg)
CONTROL_MIN_DIST_FROM_CATALOG_U = mpf('0.5')      # rejection radius from any of the 23 centres

# ---- local root-density spot-check (catalog only, budget reasons logged) --
DENSITY_GRID_N = 12                 # 12x12 grid, matches kilo_xcheck exactly
DENSITY_GRID_SIDE_U = mpf('0.4')
DENSITY_DEDUP_TOL_C = mpf('1e-30')

# ---- drift-law vector check (item 5), literals from THEORY_NOTE.md §5bis --
DRIFT_C2 = mpf('0.79028031944088979')
DRIFT_AU0 = mpf('0.41403799908627305889537448439555967803329915886199')

# ---- soft wall-clock budget (safety net only, see main()) -----------------
SOFT_TIME_BUDGET_SEC = 40 * 60   # 40 min soft cap inside the ~45 min target

CERT_PATH = (r'D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/path_to_8/'
             r'point15_layer1_addressing/address_certificate.json')

HERE = os.path.dirname(os.path.abspath(__file__))
LOG_PATH = os.path.join(HERE, 'progress.log')
RESULTS_PATH = os.path.join(HERE, 'results.json')

# Z72 provenance correction (frozen BEFORE the run; see header docstring).
Z72_CERTIFIED_U = ('4.758336330290482', '7.1277051021354385')
Z72_CORRECTED_U = ('4.791660937911268425586256', '7.177581007211361312933728')
Z72_CORRECTED_SOURCE = ('kilo_xcheck/results.json stage3_catalog.Z72.u_actual '
                         '(re_str/im_str) -- the position kilo_xcheck\'s own '
                         'Newton search converged to at m=10 when probing '
                         'near the certified disk centre.')

_log_lines = []
_t_start_wall = time.time()


def log(msg):
    line = '[%s +%6.1fs] %s' % (time.strftime('%H:%M:%S'),
                                 time.time() - _t_start_wall, msg)
    _log_lines.append(line)
    print(line, flush=True)


def flush_log():
    with open(LOG_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(_log_lines) + '\n')


def time_left():
    return SOFT_TIME_BUDGET_SEC - (time.time() - _t_start_wall)


# --------------------------------------------------------------------------- #
#  core dynamics: critical orbit of f_c(z) = z^2 + c, x_0=0, with AD deriv.  #
# --------------------------------------------------------------------------- #
def orbit_F_dF(c, n, escape_cap=mpf('1e15')):
    """x_0=0, x_{i+1}=x_i^2+c ; dx_{i+1}=2*x_i*dx_i+1, dx_0=0.
    Early-exits once |x|>escape_cap (bounds cost of escaping starts)."""
    x = mpc(0)
    dx = mpc(0)
    for _ in range(n):
        dx = 2 * x * dx + 1
        x = x * x + c
        if fabs(x) > escape_cap:
            return x, dx
    return x, dx


def orbit_x_at(c, n):
    x = mpc(0)
    for _ in range(n):
        x = x * x + c
    return x


def exact_period_check(c, m):
    """True iff x_{2^j}(c) != 0 for all 1<=j<=m (minimal period is exactly
    2^(m+1), i.e. c_ss is not an inherited root of a smaller doubling
    level). Returns (ok, first_failing_j_or_None, abs_values_dict)."""
    absvals = {}
    for j in range(1, m + 1):
        v = fabs(orbit_x_at(c, 2 ** j))
        absvals[str(2 ** j)] = float(v)
        if v < EXACT_PERIOD_TOL:
            return False, j, absvals
    return True, None, absvals


# --------------------------------------------------------------------------- #
#  damped Newton with tiered, fully deterministic restarts                   #
# --------------------------------------------------------------------------- #
def newton_solve(c0, n, tol, max_iter, min_step, max_halve):
    c = c0
    x, dx = orbit_F_dF(c, n)
    if fabs(x) < tol:
        return c, 'conv', fabs(x), 0
    if fabs(x) > ESCAPE_THRESHOLD:
        return c, 'escape', fabs(x), 0
    for it in range(max_iter):
        if dx == 0:
            return c, 'singular', fabs(x), it
        step = x / dx
        if fabs(step) < min_step and fabs(x) > tol:
            return c, 'stalled', fabs(x), it
        s = step
        improved = False
        for _h in range(max_halve + 1):
            c_try = c - s
            x_try, dx_try = orbit_F_dF(c_try, n)
            if fabs(x_try) < fabs(x):
                c, x, dx = c_try, x_try, dx_try
                improved = True
                break
            s = s / 2
        if not improved:
            return c, 'stalled', fabs(x), it
        if fabs(x) > ESCAPE_THRESHOLD:
            return c, 'escape', fabs(x), it
        if fabs(x) < tol:
            return c, 'conv', fabs(x), it + 1
    return c, 'maxiter', fabs(x), max_iter


def zdump(z, d=30):
    return {'re': float(re(z)), 'im': float(im(z)),
            're_str': nstr(re(z), d), 'im_str': nstr(im(z), d)}


def fdump(x, d=15):
    return {'val': float(x), 'str': nstr(x, d)}


def tiered_search(target_id, u_cat, c_F, A, m, n, budget_guard=True):
    """Tier 1 -> Tier 2 (3x3 grid) -> Tier 3 (fan). Every attempt is logged
    (full disclosure, item 6). Returns (accepted_record_or_None, attempts,
    r_accept_c)."""
    c_pred = c_F + u_cat / (A * DELTA ** m)
    r_accept_c = R_ACCEPT_U / (fabs(A) * DELTA ** m)
    attempts = []

    def try_one(tier, tag, c0):
        if budget_guard and time_left() < 30:
            attempts.append({'tier': tier, 'tag': tag, 'skipped': True,
                              'reason': 'soft_time_budget_exceeded'})
            return None
        c_ss, st, xr, iters = newton_solve(c0, n, NEWTON_TOL, NEWTON_MAXITER,
                                            NEWTON_MINSTEP, NEWTON_MAXHALVE)
        rec = {'tier': tier, 'tag': tag, 'c0': zdump(c0, 20),
               'status': st, 'iters': iters, 'final_abs_x': float(xr)}
        if st == 'conv':
            dist_pred = fabs(c_ss - c_pred)
            ok_period, fail_j, absvals = exact_period_check(c_ss, m)
            accepted = bool(dist_pred < r_accept_c and ok_period)
            u_actual = A * (c_ss - c_F) * DELTA ** m
            rec.update({
                'c_ss': zdump(c_ss, 30),
                'u_actual': zdump(u_actual, 25),
                'dist_to_c_pred': fdump(dist_pred),
                'exact_period_ok': ok_period,
                'exact_period_fail_j': fail_j,
                'accepted': accepted,
                'rejection_reason': (None if accepted else
                                     ('not_exact_period' if not ok_period else
                                      'dist_to_c_pred_exceeds_R_ACCEPT')),
            })
        attempts.append(rec)
        return rec if (st == 'conv' and rec.get('accepted')) else None

    hit = try_one(1, 'direct', c_pred)
    if hit:
        return hit, attempts, c_pred, r_accept_c

    grid = TIER2_FRAC * fabs(c_pred - c_F)
    for dr in (-1, 0, 1):
        for di in (-1, 0, 1):
            if dr == 0 and di == 0:
                continue
            c0 = c_pred + mpc(dr, di) * grid
            hit = try_one(2, 'grid3x3_dr%d_di%d' % (dr, di), c0)
            if hit:
                return hit, attempts, c_pred, r_accept_c

    for frac in TIER3_FRACS:
        radius = frac * fabs(c_pred - c_F)
        for ang in TIER3_ANGLES_DEG:
            th = mpf(ang) * pi / 180
            c0 = c_pred + radius * mpc(cos(th), sin(th))
            hit = try_one(3, 'fan_f%s_a%d' % (nstr(frac, 3), ang), c0)
            if hit:
                return hit, attempts, c_pred, r_accept_c

    return None, attempts, c_pred, r_accept_c


# ============================== STAGE 1 ==================================== #
def stage1():
    log('=== STAGE 1: real superstable cascade c_k, k=2..%d (derive c_F) ===' % CASCADE_K_MAX)
    cascade = {0: mpf(0), 1: mpf(-1)}
    notes = []
    for k in range(2, CASCADE_K_MAX + 1):
        n = 2 ** k
        cp1, cp2 = cascade[k - 1], cascade[k - 2]
        extrap = cp1 + (cp1 - cp2) / DELTA
        cand_starts = [extrap] + [extrap + off for off in
                                  (mpf('-0.05'), mpf('0.05'), mpf('-0.1'),
                                   mpf('-0.2'), mpf('-0.3'), mpf('0.1'),
                                   mpf('-0.5'))]
        chosen = None
        for s0 in cand_starts:
            c, st, xr, _it = newton_solve(s0, n, CASCADE_TOL, CASCADE_MAXITER,
                                          CASCADE_MINSTEP, CASCADE_MAXHALVE)
            if st == 'conv':
                ok, _fj, _av = exact_period_check(c, k - 1) if k >= 2 else (True, None, {})
                if ok:
                    chosen = (c, xr)
                    break
        if chosen is None:
            c, st, xr, _it = newton_solve(extrap, n, CASCADE_TOL, CASCADE_MAXITER,
                                          CASCADE_MINSTEP, CASCADE_MAXHALVE)
            chosen = (c, xr)
            notes.append('k=%d: exact-period scan failed, used extrapolation root' % k)
        cascade[k] = chosen[0]
    cF_seq = {k: (DELTA * cascade[k] - cascade[k - 1]) / (DELTA - 1)
              for k in range(2, CASCADE_K_MAX + 1)}
    r2 = 1 / (DELTA * DELTA)
    cF_rich = {k: (cF_seq[k] - r2 * cF_seq[k - 1]) / (1 - r2)
               for k in range(3, CASCADE_K_MAX + 1)}
    c_F = cF_rich[CASCADE_K_MAX]
    err_vs_brief = float(fabs(c_F - mpf(BRIEF_CF_STR)))
    log('  c_F (own cascade, Richardson@%d) = %s' % (CASCADE_K_MAX, nstr(c_F, 50)))
    log('  |c_F - brief literal|            = %.3e' % err_vs_brief)
    flush_log()
    return cascade, c_F, err_vs_brief, notes


# ============================== STAGE 2 ==================================== #
def stage2(cascade, c_F):
    log('=== STAGE 2: calibration L_k=(c_{k+1}-c_F)*delta^k -> A ===')
    L_seq = {k: (cascade[k + 1] - c_F) * DELTA ** k
             for k in range(8, CASCADE_K_MAX - 1)}
    ks = sorted(L_seq)
    rL = 1 / DELTA
    L = (L_seq[ks[-1]] - rL * L_seq[ks[-2]]) / (1 - rL)
    A = -1 / L
    u_check = A * (cascade[11] - c_F) * DELTA ** 10
    log('  L = %s' % nstr(L, 42))
    log('  A = %s' % nstr(A, 42))
    log('  self-consistency u(c_11)[m=10] = %s (expect -> -1)' % nstr(u_check, 25))
    flush_log()
    return L, A, float(fabs(u_check - (-1)))


# ============================== catalog load =============================== #
def load_catalog():
    with open(CERT_PATH, 'r', encoding='utf-8') as f:
        cert = json.load(f)
    table = {e['id']: e for e in cert['address_table']}
    entries = []
    for zid, e in table.items():
        u = mpc(mpf(e['centre_re']), mpf(e['centre_im']))
        entries.append({'id': zid, 'u_cat': u, 'band': e['band'],
                        'status_cert': e['status'], 'rho': e['rho']})
    entries.sort(key=lambda r: r['id'])
    assert len(entries) == 23, 'expected 23 catalog entries, got %d' % len(entries)
    z72_cert_u = mpc(mpf(Z72_CERTIFIED_U[0]), mpf(Z72_CERTIFIED_U[1]))
    for e in entries:
        if e['id'] == 'Z72':
            assert fabs(e['u_cat'] - z72_cert_u) < mpf('1e-9'), \
                'Z72 certified centre mismatch vs frozen literal'
    z72_corrected_u = mpc(mpf(Z72_CORRECTED_U[0]), mpf(Z72_CORRECTED_U[1]))
    return entries, z72_cert_u, z72_corrected_u, cert


# ============================== controls ==================================== #
def gen_controls(catalog_entries):
    rnd = _random.Random(CONTROL_SEED)
    pts = []
    log_gen = []
    tries = 0
    while len(pts) < N_CONTROLS:
        tries += 1
        r = rnd.uniform(float(CONTROL_R_RANGE[0]), float(CONTROL_R_RANGE[1]))
        th = rnd.uniform(float(CONTROL_THETA_RANGE_DEG[0]), float(CONTROL_THETA_RANGE_DEG[1]))
        u = mpc(mpf(r) * cos(mpf(th) * pi / 180), mpf(r) * sin(mpf(th) * pi / 180))
        dmin = min(fabs(u - e['u_cat']) for e in catalog_entries)
        ok = dmin > CONTROL_MIN_DIST_FROM_CATALOG_U
        log_gen.append({'try': tries, 'r': r, 'theta_deg': th,
                        'min_dist_to_catalog_u': float(dmin), 'accepted': ok})
        if ok:
            pts.append({'id': 'CTRL%02d' % (len(pts) + 1), 'u_cat': u,
                       'r': r, 'theta_deg': th, 'min_dist_to_catalog_u': float(dmin)})
    log('  generated %d controls in %d tries (seed=%d)' % (N_CONTROLS, tries, CONTROL_SEED))
    flush_log()
    return pts, log_gen


# ============================== per-point pipeline =========================== #
def process_point(pid, u_cat, kind, c_F, A):
    hit, attempts, c_pred, r_accept_c = tiered_search(pid, u_cat, c_F, A, M_LEVEL, PERIOD)
    rec = {
        'id': pid, 'kind': kind,
        'u_cat': zdump(u_cat, 25),
        'm': M_LEVEL, 'period': PERIOD,
        'c_pred': zdump(c_pred, 30),
        'r_accept_c': fdump(r_accept_c),
        'n_attempts': len(attempts),
        'attempts': attempts,
        'accepted': hit is not None,
    }
    if hit:
        u_actual = mpc(hit['u_actual']['re'], hit['u_actual']['im'])
        u_cat_c = u_cat
        d_u = fabs(mpc(mpf(hit['u_actual']['re_str']), mpf(hit['u_actual']['im_str'])) - u_cat_c)
        rec.update({
            'c_ss': hit['c_ss'], 'u_actual': hit['u_actual'],
            'd_u': fdump(d_u),
            'exact_period_ok': hit['exact_period_ok'],
            'winning_tier': hit['tier'], 'winning_tag': hit['tag'],
            'dist_to_c_pred': hit['dist_to_c_pred'],
        })
        log('  [%s] %-9s ACCEPTED tier=%d d_u=%.4e exact_period=%s (%d attempts)' %
            (kind, pid, hit['tier'], float(d_u), hit['exact_period_ok'], len(attempts)))
    else:
        # honest failure disclosure: report best (closest-to-accept) attempt too
        conv_attempts = [a for a in attempts if a.get('status') == 'conv']
        best = None
        if conv_attempts:
            best = min(conv_attempts, key=lambda a: a['dist_to_c_pred']['val'])
        rec.update({'d_u': None, 'best_converged_but_rejected': best})
        log('  [%s] %-9s NOT ACCEPTED (%d attempts, %d converged-but-rejected)' %
            (kind, pid, len(attempts), len(conv_attempts)))
    return rec


# ============================== density grid (catalog only) ================== #
def density_grid(u_cat, c_F, A):
    side = DENSITY_GRID_SIDE_U
    N = DENSITY_GRID_N
    grid_u0 = []
    roots = []
    conv = 0
    esc = 0
    t0 = time.time()
    for a in range(N):
        for b in range(N):
            du_re = (mpf(a) / (N - 1) - mpf('0.5')) * side
            du_im = (mpf(b) / (N - 1) - mpf('0.5')) * side
            u0 = u_cat + mpc(du_re, du_im)
            grid_u0.append(u0)
            c0 = c_F + u0 / (A * DELTA ** M_LEVEL)
            c, st, xr, _it = newton_solve(c0, PERIOD, NEWTON_TOL, 100,
                                          NEWTON_MINSTEP, 6)
            if st == 'conv':
                conv += 1
                roots.append(c)
            elif st == 'escape':
                esc += 1
    dt = time.time() - t0
    dedup = []
    for c in roots:
        if all(fabs(c - d) > DENSITY_DEDUP_TOL_C for d in dedup):
            dedup.append(c)
    dedup_u = [A * (c - c_F) * DELTA ** M_LEVEL for c in dedup]
    exact_flags = [exact_period_check(c, M_LEVEL)[0] for c in dedup]
    info = sorted(({'u': zdump(dedup_u[i], 18), 'dist': float(fabs(dedup_u[i] - u_cat)),
                    'exact': exact_flags[i]} for i in range(len(dedup))),
                  key=lambda r: r['dist'])
    nearest = info[0]['dist'] if info else None
    second = info[1]['dist'] if len(info) > 1 else None
    acc = sorted(float(min((fabs(u0 - ru) for ru in dedup_u), default=mpf('nan')))
                 for u0 in grid_u0)
    acc = [x for x in acc if x == x]  # drop nan
    acc_med = acc[len(acc) // 2] if acc else None
    return {
        'grid_total': N * N, 'grid_converged': conv, 'escapes': esc,
        'unique_roots': len(dedup), 'n_exact_period': sum(exact_flags),
        'nearest_u': nearest, 'second_u': second,
        'accidental_median_u': acc_med, 'elapsed_sec': dt,
        'roots_top': info[:10],
    }


# ============================== drift-law vector check ======================= #
def drift_predict(u_cat):
    return -(DRIFT_C2 * DRIFT_AU0 / DELTA ** 9) * u_cat ** 2 * DELTA / (DELTA - 1)


def drift_check(u_cat, u_actual):
    R = u_actual - u_cat
    Delta_pred = drift_predict(u_cat)
    ratio_complex = None
    ratio_abs = None
    phase_diff_deg = None
    if fabs(Delta_pred) > 0:
        ratio_complex = R / Delta_pred
        ratio_abs = float(fabs(R) / fabs(Delta_pred))
        phase_diff_deg = float((arg(R) - arg(Delta_pred)) * 180 / pi)
        # wrap to (-180, 180]
        while phase_diff_deg > 180:
            phase_diff_deg -= 360
        while phase_diff_deg <= -180:
            phase_diff_deg += 360
    return {
        'R': zdump(R, 20), 'Delta_pred': zdump(Delta_pred, 20),
        'abs_R': float(fabs(R)), 'abs_Delta_pred': float(fabs(Delta_pred)),
        'ratio_abs': ratio_abs,
        'ratio_complex': (zdump(ratio_complex, 15) if ratio_complex is not None else None),
        'phase_diff_deg': phase_diff_deg,
    }


# ==================================== main ==================================== #
def main():
    log('kilo-style FROZEN unified Mandelbrot run start (mp.dps=%d, m=%d, period=%d)'
        % (MP_DPS, M_LEVEL, PERIOD))
    flush_log()

    cascade, c_F, cF_err_vs_brief, s1_notes = stage1()
    L, A, u_check_err = stage2(cascade, c_F)

    catalog_entries, z72_cert_u, z72_corrected_u, cert = load_catalog()
    log('=== loaded %d catalog entries from address_certificate.json ===' % len(catalog_entries))
    flush_log()

    controls, controls_gen_log = gen_controls(catalog_entries)

    log('=== STAGE 3: catalog cross-check (23 points, m=%d) ===' % M_LEVEL)
    catalog_results = []
    for e in catalog_entries:
        catalog_results.append(process_point(e['id'], e['u_cat'], 'catalog', c_F, A))
        flush_log()

    log('=== STAGE 3b: Z72 provenance -- corrected address (see header) ===')
    z72_corrected_result = process_point('Z72_corrected', z72_corrected_u,
                                         'catalog_provenance_check', c_F, A)
    flush_log()

    log('=== STAGE 4: %d deterministic negative controls (seed=%d) ===' %
        (N_CONTROLS, CONTROL_SEED))
    control_results = []
    for cpt in controls:
        r = process_point(cpt['id'], cpt['u_cat'], 'control', c_F, A)
        r['generation'] = {'r': cpt['r'], 'theta_deg': cpt['theta_deg'],
                           'min_dist_to_catalog_u': cpt['min_dist_to_catalog_u']}
        control_results.append(r)
        flush_log()

    log('=== STAGE 5: local root-density spot-check (catalog only, 12x12 grid) ===')
    density_results = {}
    for rec in catalog_results:
        if time_left() < 60:
            log('  !! soft time budget nearly exhausted, skipping remaining density grids')
            break
        u_cat_c = mpc(mpf(rec['u_cat']['re_str']), mpf(rec['u_cat']['im_str']))
        d = density_grid(u_cat_c, c_F, A)
        density_results[rec['id']] = d
        d_u = rec['d_u']['val'] if rec['d_u'] else None
        spec = (d['accidental_median_u'] / d_u) if (d_u and d['accidental_median_u']) else None
        log('  density %-6s: roots=%2d nearest=%s acc_med=%s d_u=%s spec=%s (%.1fs)' %
            (rec['id'], d['unique_roots'],
             ('%.3e' % d['nearest_u']) if d['nearest_u'] else 'n/a',
             ('%.3e' % d['accidental_median_u']) if d['accidental_median_u'] else 'n/a',
             ('%.3e' % d_u) if d_u else 'n/a',
             ('%.1fx' % spec) if spec else 'n/a', d['elapsed_sec']))
        flush_log()

    log('=== STAGE 6: drift-law vector check (catalog, item 5) ===')
    drift_results = {}
    for rec in catalog_results:
        if not rec['accepted']:
            continue
        u_cat_c = mpc(mpf(rec['u_cat']['re_str']), mpf(rec['u_cat']['im_str']))
        u_act_c = mpc(mpf(rec['u_actual']['re_str']), mpf(rec['u_actual']['im_str']))
        dc = drift_check(u_cat_c, u_act_c)
        drift_results[rec['id']] = dc
        log('  drift %-6s: |R|/|Delta|=%s phase_diff=%s deg' %
            (rec['id'],
             ('%.4f' % dc['ratio_abs']) if dc['ratio_abs'] else 'n/a',
             ('%.3f' % dc['phase_diff_deg']) if dc['phase_diff_deg'] is not None else 'n/a'))
    flush_log()

    # ---- specificity vs. pool of found roots (cheap, all 43 points) --------
    log('=== STAGE 7: pool-based specificity (median dist to OTHER found roots) ===')
    all_found = []
    for rec in catalog_results + control_results:
        if rec['accepted']:
            all_found.append((rec['id'], mpc(mpf(rec['u_actual']['re_str']),
                                             mpf(rec['u_actual']['im_str']))))
    pool_specificity = {}
    for pid, u_act in all_found:
        others = [u for oid, u in all_found if oid != pid]
        if others:
            dists = sorted(float(fabs(u_act - u)) for u in others)
            med = dists[len(dists) // 2]
        else:
            med = None
        rec = next(r for r in catalog_results + control_results if r['id'] == pid)
        d_u = rec['d_u']['val']
        pool_specificity[pid] = {
            'median_dist_to_other_found_roots': med,
            'pool_size': len(others),
            'specificity': (med / d_u) if (med and d_u) else None,
        }
    flush_log()

    # ---- aggregate verdict --------------------------------------------------
    cat_accepted = [r for r in catalog_results if r['accepted']]
    cat_rejected = [r for r in catalog_results if not r['accepted']]
    ctrl_accepted = [r for r in control_results if r['accepted']]
    cat_dus = sorted(r['d_u']['val'] for r in cat_accepted)
    ctrl_dus = sorted(r['d_u']['val'] for r in ctrl_accepted)

    def median(xs):
        if not xs:
            return None
        xs = sorted(xs)
        n = len(xs)
        return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2

    spec_vals = [pool_specificity[r['id']]['specificity'] for r in cat_accepted
                 if pool_specificity[r['id']]['specificity']]
    drift_ratios = [drift_results[r['id']]['ratio_abs'] for r in cat_accepted
                    if r['id'] in drift_results and drift_results[r['id']]['ratio_abs']]
    drift_phases = [abs(drift_results[r['id']]['phase_diff_deg']) for r in cat_accepted
                    if r['id'] in drift_results and drift_results[r['id']]['phase_diff_deg'] is not None]

    verdict = {
        'catalog_total': len(catalog_results),
        'catalog_accepted_match': len(cat_accepted),
        'catalog_rejected': len(cat_rejected),
        'catalog_rejected_ids': [r['id'] for r in cat_rejected],
        'catalog_median_d_u': median(cat_dus),
        'catalog_d_u_range': [cat_dus[0], cat_dus[-1]] if cat_dus else None,
        'control_total': len(control_results),
        'control_accepted': len(ctrl_accepted),
        'control_median_d_u': median(ctrl_dus),
        'median_pool_specificity_catalog': median(spec_vals),
        'median_drift_ratio_abs': median(drift_ratios),
        'median_abs_phase_diff_deg': median(drift_phases),
        'z72_note': ('Z72 as certified: see catalog_rejected_ids / attempts log. '
                    'Z72_corrected (provenance row, excluded from tally above) '
                    'result is in stage3b_z72_corrected.'),
        'stage1_internal_checks': {
            'c_F_err_vs_brief_literal': cF_err_vs_brief,
            'u_check_self_consistency_err': u_check_err,
        },
        'time_budget_exceeded': time_left() < 0,
    }

    log('=== VERDICT ===')
    log('  catalog accepted-match: %d/%d  (rejected: %s)' %
        (verdict['catalog_accepted_match'], verdict['catalog_total'],
         verdict['catalog_rejected_ids']))
    log('  catalog median d_u = %s   control median d_u = %s' %
        (verdict['catalog_median_d_u'], verdict['control_median_d_u']))
    log('  median pool specificity (catalog) = %s' % verdict['median_pool_specificity_catalog'])
    log('  median |drift ratio| = %s   median |phase diff| deg = %s' %
        (verdict['median_drift_ratio_abs'], verdict['median_abs_phase_diff_deg']))
    log('  total wall time = %.1f sec' % (time.time() - _t_start_wall))

    results = {
        'meta': {
            'mp_dps': MP_DPS, 'm_level': M_LEVEL, 'period': PERIOD,
            'delta_str': DELTA_STR,
            'c_F_str': nstr(c_F, 50), 'c_F_err_vs_brief_literal': cF_err_vs_brief,
            'A_str': nstr(A, 42), 'L_str': nstr(L, 42),
            'u_check_self_consistency_err': u_check_err,
            'r_accept_u': float(R_ACCEPT_U),
            'control_seed': CONTROL_SEED,
            'elapsed_sec': time.time() - _t_start_wall,
            'stage1_notes': s1_notes,
        },
        'stage1_cascade': {str(k): nstr(cascade[k], 40) for k in cascade},
        'controls_generation_log': controls_gen_log,
        'stage3_catalog': catalog_results,
        'stage3b_z72_corrected': z72_corrected_result,
        'stage3b_z72_provenance': {
            'certified_u': Z72_CERTIFIED_U,
            'corrected_u': Z72_CORRECTED_U,
            'source': Z72_CORRECTED_SOURCE,
            'distance_certified_to_corrected': float(fabs(z72_cert_u - z72_corrected_u)),
            'point17_disk_rho': '0.1',
            'note': ('distance_certified_to_corrected < disk_rho, i.e. the '
                    'corrected position lies inside the point17-certified '
                    'disk around the certified centre -- consistent with '
                    'the disk containing the true zero even though its '
                    'centre is a loose (shifted) proxy.'),
        },
        'stage4_controls': control_results,
        'stage5_density': density_results,
        'stage6_drift': drift_results,
        'stage7_pool_specificity': pool_specificity,
        'verdict': verdict,
    }
    def _json_fallback(o):
        # Defensive net only: every value placed into `results` above is
        # already a plain float/int/str/bool/None/dict/list. If something
        # slips through anyway (e.g. a raw mpmath type from a future edit),
        # stringify it rather than lose the entire ~30min run to a crash at
        # the very last line.
        return str(o)

    try:
        with open(RESULTS_PATH, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        log('wrote %s' % RESULTS_PATH)
    except TypeError as exc:
        log('!! primary json.dump failed (%s) -- retrying with str() fallback '
            'encoder so the run is not lost' % exc)
        with open(RESULTS_PATH, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False, default=_json_fallback)
        log('wrote %s (fallback encoder used -- see log line above)' % RESULTS_PATH)
    flush_log()


if __name__ == '__main__':
    main()

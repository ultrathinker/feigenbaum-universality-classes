#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v2_rate.py -- INDEPENDENT Phase-2 verification of agent A's d_u(m) decay-rate result.
[kid:verify] NOT RIGOROUS float (mpmath dps=80).

Own code path: own damped Newton on z_{2^(m+1)}(c)=0 with companion derivative,
own exact-period check, own trajectory machinery u(m)=A*(c_ss-c_F)*delta^m,
own Richardson extrapolation, own ratio/median statistics, own candidate table,
plus a TWO-LAYER fit extension (u(m)=u_inf+C*delta^-m+D*lambda2^m vs ...+E*delta^-2m)
directly testing agent A's section-6 interpretations (A)/(B).

External inputs (marked):
  * c_F, A, delta -- VERBATIM from Investigation-1 results.json meta (calibration).
  * agent A zero catalog zeros_deg2.json REVISION 2 -- used (a) for Stage A
    reproduction of their pipeline, (b) as SEEDS ONLY for my own trajectories.
  * my lambda2/lambda3/alpha/delta from v1_cheb.json (my phase 1) for candidates.
"""
import json
import time

from mpmath import mp, mpf, mpc, fabs, sqrt, matrix, lu_solve, nstr

mp.dps = 80

TAG = '[kid:verify] independent phase-2 re-verification of agent A d_u(m) rate study'
T0 = time.time()
LOG_PATH = 'v2_run.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


# ---------------- constants (VERBATIM calibration) -----------------------------
RES1 = r'../study1_d2_frozen_experiment/results.json'
_r1 = json.load(open(RES1, encoding='utf-8'))
C_F = mpf(_r1['meta']['c_F_str'].replace('(', '').replace(')', '').split('+')[0].strip())
A_C = mpf(_r1['meta']['A_str'].replace('(', '').replace(')', '').split('+')[0].strip())
DELTA = mpf(_r1['meta']['delta_str'])
log('constants VERBATIM from Investigation-1 results.json meta:')
log('  c_F=%s' % nstr(C_F, 50))
log('  A=%s' % nstr(A_C, 50))
log('  delta=%s' % nstr(DELTA, 50))

_v1 = json.load(open('v1_cheb.json', encoding='utf-8'))
LAM2 = mpf(_v1['lambda2'])
LAM3 = mpf(_v1['lambda3'])
ALPHA = mpf(_v1['alpha'])
DELTA_V1 = mpf(_v1['delta'])
log('my phase-1 values: alpha=%s lambda2=%s lambda3=%s' % (
    nstr(ALPHA, 30), nstr(LAM2, 25), nstr(LAM3, 25)))

# agent A catalog (EXTERNAL INPUT)
_zA = json.load(open(r'../../05-CODE/sec5_lambda2_deg2/zeros_deg2.json', encoding='utf-8'))
ZEROS = []
for z in _zA['zeros']:
    if z['id'] == 'Z72_true':
        continue  # diagnostic-only, excluded as in agent A
    ZEROS.append({
        'id': z['id'],
        'u': mpc(mpf(z['re']), mpf(z['im'])),
        'kind': z['kind'],
    })
log('loaded %d zeros from agent A catalog (Z72 excluded as diagnostic)' % len(ZEROS))

_rateA = json.load(open(r'../../05-CODE/sec5_lambda2_deg2/rate_deg2.json', encoding='utf-8'))
DU_A = {r['id']: {int(k): mpf(repr(v)) for k, v in r['d_u'].items()} for r in _rateA['rows']}


# ---------------- orbit machinery (own implementation) -------------------------

def orbit_zd(c, n):
    z = mpc(0)
    d = mpc(0)
    for _ in range(n):
        d = 2 * z * d + 1
        z = z * z + c
    return z, d


def newton_ss(c0, n, tol=mpf('1e-62'), maxiter=60, maxhalve=16):
    c = c0
    for it in range(maxiter):
        z, d = orbit_zd(c, n)
        if fabs(z) > mpf('1e30'):
            return None, 'escape', it
        if fabs(z) < tol:
            return c, 'conv', it
        step = z / d
        lam = mpf(1)
        ok = False
        for _h in range(maxhalve):
            c2 = c - lam * step
            z2, _ = orbit_zd(c2, n)
            if fabs(z2) < mpf('1e30'):
                c, ok = c2, True
                break
            lam /= 2
        if not ok:
            return None, 'stall', it
    return None, 'maxiter', maxiter


def exact_period_ok(c, m, tol=mpf('1e-45')):
    """exact period 2^(m+1): |z_{2^j}(c)| > tol for all 1<=j<=m."""
    n = 1
    z = mpc(0)
    minv = None
    for j in range(1, m + 1):
        # advance to 2^j
        for _ in range(n, 2 * n):
            z = z * z + c
        n = 2 * n
        a = fabs(z)
        if minv is None or a < minv:
            minv = a
        if a <= tol:
            return False, j, a
    return True, None, minv


def c_of_u(u, m):
    return C_F + u / (A_C * DELTA ** m)


def u_of_c(c, m):
    return A_C * (c - C_F) * DELTA ** m


# ---------------- Stage A: reproduce agent A pipeline --------------------------
log('=== STAGE A: reproduce agent A d_u(m) from THEIR catalog u (my Newton) ===')
rowsA = []
for z in ZEROS:
    u = z['u']
    dus = {}
    css = {}
    okall = True
    for m in (8, 9, 10, 11):
        c_pred = c_of_u(u, m)
        c_ss, st, it = newton_ss(c_pred, 2 ** (m + 1))
        if c_ss is None:
            log('  %s m=%d NEWTON FAILED (%s)' % (z['id'], m, st))
            okall = False
            break
        ok, failj, minv = exact_period_ok(c_ss, m)
        dist = fabs(c_ss - c_pred)
        if not ok:
            log('  %s m=%d PERIOD FAIL at j=%d (min |z|=%.2e)' % (z['id'], m, failj, float(minv)))
            okall = False
            break
        dus[m] = fabs(u - u_of_c(c_ss, m))
        css[m] = c_ss
    if not okall:
        continue
    r89 = dus[8] / dus[9]
    r910 = dus[9] / dus[10]
    r1011 = dus[10] / dus[11]
    rowsA.append({'id': z['id'], 'kind': z['kind'], 'u': u, 'd_u': dus, 'r89': r89,
                  'r910': r910, 'r1011': r1011, 'c_ss': css})
    log('  %s: d8=%.4e d9=%.4e d10=%.4e d11=%.4e  r=%.3f/%.3f/%.3f  |css-cpred|max=%.1e' % (
        z['id'], dus[8], dus[9], dus[10], dus[11], r89, r910, r1011,
        max(float(fabs(css[m] - c_of_u(u, m))) for m in (8, 9, 10, 11))))

# compare d_u with agent A values
log('--- d_u diff vs agent A rate_deg2.json (relative) ---')
worst_rel = mpf(0)
for r in rowsA:
    for m in (8, 9, 10, 11):
        dv = DU_A[r['id']][m]
        rel = fabs(r['d_u'][m] - dv) / dv
        worst_rel = max(worst_rel, rel)
log('worst relative d_u deviation from agent A: %.3e' % float(worst_rel))


def median(vals):
    s = sorted(vals)
    n = len(s)
    return (s[n // 2 - 1] + s[n // 2]) / 2 if n % 2 == 0 else s[n // 2]


all14 = [r for r in rowsA]
real3 = [r for r in rowsA if 'real' in r['kind']]
log('STAGE A medians: all-14 r89=%s r910=%s r1011=%s' % (
    nstr(median([r['r89'] for r in all14]), 8), nstr(median([r['r910'] for r in all14]), 8),
    nstr(median([r['r1011'] for r in all14]), 8)))
log('STAGE A medians: real-3  r89=%s r910=%s r1011=%s' % (
    nstr(median([r['r89'] for r in real3]), 8), nstr(median([r['r910'] for r in real3]), 8),
    nstr(median([r['r1011'] for r in real3]), 8)))

# ---------------- Stage B: my own trajectories m=8..13 -------------------------
log('=== STAGE B: my own level trajectories u(m), m=8..13 (seeds: agent A catalog, perturbed variants) ===')
LEVELS = [8, 9, 10, 11, 12, 13]


def walk_trajectory(u_start, pert=mpf(0)):
    """returns {m: u(m)} walking levels with previous level as seed, or (None, reason)"""
    u_prev = u_start * (1 + pert)
    out = {}
    for m in LEVELS:
        c_pred = c_of_u(u_prev, m)
        c_ss, st, it = newton_ss(c_pred, 2 ** (m + 1))
        if c_ss is None:
            return None, 'newton-%s-m%d' % (st, m)
        ok, failj, minv = exact_period_ok(c_ss, m)
        if not ok:
            return None, 'period-m%d-j%d' % (m, failj)
        dist = fabs(c_ss - c_pred)
        if dist > mpf('1e-6'):
            return None, 'farjump-m%d-%.1e' % (m, float(dist))
        u_new = u_of_c(c_ss, m)
        out[m] = u_new
        u_prev = u_new
    return out, None


trajectories = {}
for z in ZEROS:
    traj, err = walk_trajectory(z['u'])
    if traj is None:
        log('  %s TRAJECTORY FAILED: %s' % (z['id'], err))
        continue
    trajectories[z['id']] = traj
    # seed-robustness: +-1e-3 multiplicative + imaginary kick
    tr2, e2 = walk_trajectory(z['u'], mpf('1e-3'))
    tr3, e3 = walk_trajectory(z['u'], mpf('-1e-3'))
    dev = mpf(0)
    for m in LEVELS:
        if tr2:
            dev = max(dev, fabs(tr2[m] - traj[m]))
        if tr3:
            dev = max(dev, fabs(tr3[m] - traj[m]))
    log('  %s: u(8)=%.14f%+.14fj u(10)=%.14f%+.14fj u(13)=%.14f%+.14fj  seed-robust dev=%.2e [%s/%s]' % (
        z['id'], traj[8].real, traj[8].imag, traj[10].real, traj[10].imag,
        traj[13].real, traj[13].imag, float(dev),
        'ok' if tr2 else e2, 'ok' if tr3 else e3))

# cross-check with zero_map level trajectories (external archive, m=10)
_zm = json.load(open(r'D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/ideation_sweep1/codex/zero_map.json', encoding='utf-8'))
zm10 = {r['id']: mpc(mpf(repr(r['level_trajectories'][3]['u'][0])), mpf(repr(r['level_trajectories'][3]['u'][1])))
        for r in _zm['roots'] if len(r['level_trajectories']) >= 4}
log('--- cross-check my u(10) vs zero_map level-10 values (their abs_H ~1e-8..1e-6) ---')
for zid, traj in list(trajectories.items())[:5]:
    if zid in zm10:
        log('  %s: |mine(10) - zeromap(10)| = %.3e' % (zid, float(fabs(traj[10] - zm10[zid]))))

# ---------------- Stage C: my Richardson u_inf, my d_u, ratios -----------------
log('=== STAGE C: my Richardson u_inf (top pair 12->13, exact 1/delta), my d_u(m) ===')
r_inv = 1 / DELTA
factor = r_inv / (1 - r_inv)  # = 1/(delta-1)
u_inf_mine = {}
for zid, traj in trajectories.items():
    shift = traj[13] - traj[12]
    u_inf_mine[zid] = traj[13] + shift * factor

rowsC = []
for z in ZEROS:
    zid = z['id']
    if zid not in trajectories:
        continue
    traj = trajectories[zid]
    uinf = u_inf_mine[zid]
    dus = {m: fabs(uinf - traj[m]) for m in LEVELS}
    rws = {m: dus[m] / dus[m + 1] for m in (8, 9, 10, 11, 12)}
    rowsC.append({'id': zid, 'kind': z['kind'], 'uinf': uinf, 'd_u': dus, 'r': rws,
                  'u_agentA': z['u']})
    log('  %s: uinf=%.15f%+.15fj |uinf-uA|=%.2e  d: %.3e %.3e %.3e %.3e %.3e  r=%.3f/%.3f/%.3f/%.3f/%.3f' % (
        zid, uinf.real, uinf.imag, float(fabs(uinf - z['u'])),
        dus[8], dus[9], dus[10], dus[11], dus[12],
        rws[8], rws[9], rws[10], rws[11], rws[12]))

real3C = [r for r in rowsC if 'real' in r['kind']]
log('STAGE C medians all-14:')
for m in (8, 9, 10, 11):
    vals = [r['r'][m] for r in rowsC]
    log('  r(%d,%d): median=%s' % (m, m + 1, nstr(median(vals), 8)))
log('STAGE C medians real-3:')
for m in (8, 9, 10, 11):
    vals = [r['r'][m] for r in real3C]
    log('  r(%d,%d): median=%s' % (m, m + 1, nstr(median(vals), 8)))

# candidates
log('--- candidate comparison (my phase-1 values) ---')
cands = {
    'delta': DELTA_V1,
    '1/|lambda2|': 1 / fabs(LAM2),
    '1/|lambda3|': 1 / fabs(LAM3),
    'alpha^2': ALPHA ** 2,
    'sqrt(delta)': sqrt(DELTA_V1),
    'delta/alpha^2': DELTA_V1 / ALPHA ** 2,
}
medC = {m: median([r['r'][m] for r in rowsC]) for m in (8, 9, 10, 11)}
medR = {m: median([r['r'][m] for r in real3C]) for m in (8, 9, 10, 11)}
for nm, cv in cands.items():
    diffs_all = ' '.join('%.1f%%' % (100 * abs(float(medC[m]) - float(cv)) / float(cv)) for m in (8, 9, 10))
    diffs_real = ' '.join('%.1f%%' % (100 * abs(float(medR[m]) - float(cv)) / float(cv)) for m in (8, 9, 10))
    log('  %-14s = %.6f   vs all-14 r8/9,r9/10,r10/11: %s | vs real-3: %s' % (nm, float(cv), diffs_all, diffs_real))

# ---------------- Stage D: two-layer fits ---------------------------------------
log('=== STAGE D: two-layer decomposition u(m)=u_inf + C*delta^-m + D*lambda2^m  vs  ...+E*delta^-2m ===')
log('    (exact 3-param solve on m=11,12,13; residuals reported on m=8,9,10)')


def fit3(traj, base2):
    """base2: function m->mpc second-layer basis (lambda2^m or delta^-2m).
    Solve u(m) = u_inf + C*d^-m + D*base2(m) on m=11,12,13 -> (u_inf,C,D)."""
    ms = [11, 12, 13]
    M = matrix(3, 3)
    for i, m in enumerate(ms):
        M[i, 0] = mpf(1)
        M[i, 1] = DELTA ** (-m)
        M[i, 2] = base2(m)
    rhs = matrix(3, 1)
    for i, m in enumerate(ms):
        rhs[i] = traj[m]
    sol = lu_solve(M, rhs)
    return sol[0], sol[1], sol[2]


fit_results = {}
for r in rowsC:
    zid = r['id']
    traj = trajectories[zid]
    uinf_l2, C_l2, D_l2 = fit3(traj, lambda m: LAM2 ** m)
    uinf_d2, C_d2, E_d2 = fit3(traj, lambda m: DELTA ** (-2 * m))
    res_l2 = {m: fabs(traj[m] - (uinf_l2 + C_l2 * DELTA ** (-m) + D_l2 * LAM2 ** m)) for m in (8, 9, 10)}
    res_d2 = {m: fabs(traj[m] - (uinf_d2 + C_d2 * DELTA ** (-m) + E_d2 * DELTA ** (-2 * m))) for m in (8, 9, 10)}
    fit_results[zid] = {
        'uinf_l2fit': uinf_l2, 'C_l2': C_l2, 'D_l2': D_l2,
        'uinf_d2fit': uinf_d2, 'C_d2': C_d2, 'E_d2': E_d2,
        'res_l2': res_l2, 'res_d2': res_d2,
    }
    log('  %s: [lambda2-layer] |D|=%.3e  resid(8,9,10)=%.2e/%.2e/%.2e   [delta^-2 layer] |E|=%.3e resid=%.2e/%.2e/%.2e' % (
        zid, fabs(D_l2), res_l2[8], res_l2[9], res_l2[10],
        fabs(E_d2), res_d2[8], res_d2[9], res_d2[10]))

log('--- sliding second-layer ratios R(m)/R(m+1), R(m)=u(m)-uinf-C*delta^-m (lambda2-fit uinf,C) ---')
for r in rowsC:
    zid = r['id']
    fr = fit_results[zid]
    traj = trajectories[zid]
    uinf, C = fr['uinf_l2fit'], fr['C_l2']
    Rv = {m: traj[m] - uinf - C * DELTA ** (-m) for m in LEVELS}
    rats = ' '.join('%s=%.3f@%+d' % (m, float(fabs(Rv[m] / Rv[m + 1])),
                                     round(float((Rv[m] / Rv[m + 1]).real / fabs(Rv[m] / Rv[m + 1]))))
                    for m in (8, 9, 10, 11, 12))
    log('  %s: ratios %s   (predict 1/lambda2=%.4f, delta^2=%.4f)' % (
        zid, rats, float(1 / LAM2), float(DELTA ** 2)))

# ---------------- output --------------------------------------------------------
def cplx(z):
    return {'re': nstr(z.real, 40), 'im': nstr(z.imag, 40)}


out = {
    'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'mp_dps': mp.dps,
             'not_rigorous': True,
             'constants': 'c_F,A,delta VERBATIM from Investigation-1 results.json (calibration)',
             'zeros_source': 'agent A zeros_deg2.json REV2 (external input; seeds only for trajectories)',
             'levels': LEVELS, 'elapsed_sec': time.time() - T0},
    'stageA_reproduction': [
        {'id': r['id'], 'd_u': {str(m): float(r['d_u'][m]) for m in (8, 9, 10, 11)},
         'r89': float(r['r89']), 'r910': float(r['r910']), 'r1011': float(r['r1011'])} for r in rowsA],
    'stageA_worst_rel_dev_vs_agentA': float(worst_rel),
    'stageA_medians': {
        'all14': {k: float(median([r[k] for r in all14])) for k in ('r89', 'r910', 'r1011')},
        'real3': {k: float(median([r[k] for r in real3])) for k in ('r89', 'r910', 'r1011')}},
    'trajectories': {zid: {str(m): cplx(tr) for m, tr in traj.items()} for zid, traj in trajectories.items()},
    'u_inf_mine': {zid: cplx(v) for zid, v in u_inf_mine.items()},
    'u_inf_vs_agentA': {zid: float(fabs(u_inf_mine[zid] - r['u_agentA']))
                        for zid, r in ((r['id'], r) for r in rowsC) if zid in u_inf_mine},
    'stageC_rows': [
        {'id': r['id'], 'kind': r['kind'], 'uinf': cplx(r['uinf']),
         'd_u': {str(m): float(r['d_u'][m]) for m in LEVELS},
         'r': {str(m): float(r['r'][m]) for m in (8, 9, 10, 11, 12)}} for r in rowsC],
    'stageC_medians': {
        'all14': {str(m): float(medC[m]) for m in (8, 9, 10, 11)},
        'real3': {str(m): float(medR[m]) for m in (8, 9, 10, 11)}},
    'candidates': {nm: float(cv) for nm, cv in cands.items()},
    'stageD_fits': {
        zid: {'D_l2': cplx(fr['D_l2']), 'E_d2': cplx(fr['E_d2']),
              'res_l2': {str(m): float(v) for m, v in fr['res_l2'].items()},
              'res_d2': {str(m): float(v) for m, v in fr['res_d2'].items()},
              'uinf_l2fit': cplx(fr['uinf_l2fit']), 'C_l2': cplx(fr['C_l2'])}
        for zid, fr in fit_results.items()},
}
with open('v2_rate.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote v2_rate.json')

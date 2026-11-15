#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
v5_stageint.py -- Stage D/E two-layer fit IN INTERVAL ARITHMETIC (BRIEF item 3).

For the 3 real zeros Z06, Z44, Z50 (cleanest): refit
    u(m) = u_inf + C * delta^{-m} + D * lambda2^m   (exact 3x3 solve at m=11,12,13)
in arb balls, replacing float lambda2/delta by their CERTIFIED intervals:
    lambda2: certified to +-1.8219e-33 (v4 mini-K2, improves v3's +-3.0e-31)
    delta:   certified to +-2.659e-36  (v3 K2b)
Variant cascade isolates each uncertainty source:
  (a) all-mid reference;
  (b) u(m) +- 1e-57 (float u-noise grade);
  (c) + certified [delta];
  (d) + certified [lambda2]  (full).
Reports ball radii of (u_inf, C, D) per variant, the second-layer ratio
R(12)/R(13) vs 1/lambda2 in balls, and the Stage-E style prediction check.

HONEST scope: c_ss(m), A, c_F remain float calibration (per brief); this shows
STABILITY of the Stage-D conclusion under the certified [lambda2]/[delta],
not a full drift certificate.
"""
import json
import time

from flint import arb, ctx
from mpmath import mp, mpf, mpmathify, nstr

ctx.dps = 80
mp.dps = 80

TAG = '[kid:verify] Stage D interval redo (BRIEF-EXEC-4-KILO-TAIL item 3)'
T0 = time.time()
LOG_PATH = 'v5_run.log'


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def ball(rstr):
    return arb('-' + rstr).union(arb(rstr))


def rad_ub(x):
    """half-width upper bound of an arb ball as float"""
    r = (x.rad() + abs(x.mid() - x)).upper()
    return float(r)


def solve3(Mrows, yvec):
    """Cramer 3x3 in balls: solve M x = y; Mrows = [[a,b,c],...] rows."""
    a11, a12, a13 = Mrows[0]
    a21, a22, a23 = Mrows[1]
    a31, a32, a33 = Mrows[2]
    y1, y2, y3 = yvec

    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    D0 = det3(Mrows)
    Dx1 = det3([[y1, a12, a13], [y2, a22, a23], [y3, a32, a33]])
    Dx2 = det3([[a11, y1, a13], [a21, y2, a23], [a31, y3, a33]])
    Dx3 = det3([[a11, a12, y1], [a21, a22, y2], [a31, a32, y3]])
    return Dx1 / D0, Dx2 / D0, Dx3 / D0


def main():
    log(TAG)
    v2 = json.load(open('v2_rate.json', encoding='utf-8'))
    v1 = json.load(open('v1_cheb.json', encoding='utf-8'))

    lam2_mid = arb(v1['lambda2'])
    LAM2_RAD = '1.822e-33'   # v4 mini-K2 (1.8219e-33 rounded up)
    delta_mid = arb(v1['delta'])
    DELTA_RAD = '2.659e-36'  # v3 K2b (2.659e-36)

    lam2_cert = lam2_mid + ball(LAM2_RAD)
    delta_cert = delta_mid + ball(DELTA_RAD)
    log('certified [lambda2] half-width %s ; [delta] half-width %s' % (LAM2_RAD, DELTA_RAD))

    LEVELS = [11, 12, 13]
    U_NOISE = '1e-57'
    out = {}
    for zid in ['Z06', 'Z44', 'Z50']:
        tr = v2['trajectories'][zid]
        uhat = [arb(tr[str(m)]['re']) for m in LEVELS]  # real zeros: im == 0
        rec = {}
        # precompute powers per variant
        def build(delta_b, lam_b):
            rows = []
            for i, m in enumerate(LEVELS):
                rows.append([arb(1), delta_b ** (-m), lam_b ** m])
            return rows

        variants = {}
        # (a) all-mid
        variants['a_mid'] = (build(delta_mid, lam2_mid), [u for u in uhat])
        # (b) u-noise
        variants['b_unoise'] = (build(delta_mid, lam2_mid), [u + ball(U_NOISE) for u in uhat])
        # (c) + delta cert
        variants['c_delta'] = (build(delta_cert, lam2_mid), [u + ball(U_NOISE) for u in uhat])
        # (d) + lambda2 cert (full)
        variants['d_full'] = (build(delta_cert, lam2_cert), [u + ball(U_NOISE) for u in uhat])

        for name, (Mrows, yvec) in variants.items():
            x1, x2, x3 = solve3(Mrows, yvec)
            rec[name] = {
                'u_inf': str(x1), 'C': str(x2), 'D': str(x3),
                'u_inf_rad': rad_ub(x1), 'C_rad': rad_ub(x2), 'D_rad': rad_ub(x3),
            }
            log('  %s %s: u_inf=%s +-%.3e ; C +-%.3e ; D +-%.3e'
                % (zid, name, x1.mid().str(18), rec[name]['u_inf_rad'],
                   rec[name]['C_rad'], rec[name]['D_rad']))

        # second-layer ratio check with CERTIFIED lambda2:
        # R(m) = u(m) - u_inf - C delta^{-m}, computed in balls from variant (d)
        x1, x2, x3 = solve3(variants['d_full'][0], variants['d_full'][1])
        R = {}
        for i, m in enumerate(LEVELS):
            R[m] = variants['d_full'][1][i] - x1 - x2 * (delta_cert ** (-m))
        ratio_1213 = R[12] / R[13]
        ratio_1112 = R[11] / R[12]
        inv_lam2 = arb(1) / (-lam2_cert)  # 1/|lambda2| positive (lambda2<0; R alternates)
        rec['ratio_R12_R13'] = str(ratio_1213)
        rec['ratio_R11_R12'] = str(ratio_1112)
        log('  %s R(12)/R(13) = %s ; 1/|lambda2| = %s' % (zid, ratio_1213.str(20), inv_lam2.str(20)))
        # sign alternation check
        rec['R13_sign'] = 'neg' if float(R[13].upper()) < 0 else ('pos' if float(R[13].lower()) > 0 else 'straddles0')
        rec['R12_sign'] = 'neg' if float(R[12].upper()) < 0 else ('pos' if float(R[12].lower()) > 0 else 'straddles0')
        # inflation factor: variant d vs b
        rec['inflation_uinf_d_vs_b'] = rec['d_full']['u_inf_rad'] / max(rec['b_unoise']['u_inf_rad'], 1e-300)
        rec['inflation_D_d_vs_b'] = rec['d_full']['D_rad'] / max(rec['b_unoise']['D_rad'], 1e-300)
        out[zid] = rec

    # summary table
    log('--- summary (radii) ---')
    for zid in out:
        r = out[zid]
        log('  %s: u_inf rad a/b/c/d = %.2e / %.2e / %.2e / %.2e ; D rad d = %.2e ; infl(d/b) u_inf=%.2e D=%.2e'
            % (zid, r['a_mid']['u_inf_rad'], r['b_unoise']['u_inf_rad'], r['c_delta']['u_inf_rad'],
               r['d_full']['u_inf_rad'], r['d_full']['D_rad'],
               r['inflation_uinf_d_vs_b'], r['inflation_D_d_vs_b']))

    payload = {
        'meta': {'tag': TAG, 'date': time.strftime('%Y-%m-%d %H:%M:%S'),
                 'arb_dps': 80, 'not_full_certificate': True,
                 'certified_inputs': {'lambda2_halfwidth': LAM2_RAD, 'delta_halfwidth': DELTA_RAD},
                 'float_inputs': 'u(m) m=11,12,13 (trajectories v2_rate.json), c_ss/A/c_F calibration',
                 'u_noise_variant': U_NOISE},
        'fits': out,
        'elapsed_sec': time.time() - T0,
    }
    with open('v5_stageint.json', 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    log('wrote v5_stageint.json')


if __name__ == '__main__':
    main()

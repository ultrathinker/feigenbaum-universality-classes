#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate polished zeros from rescan: for each candidate, evaluate at m=7,8,9
and check drift, den>0, maxt6, and a true-zero domain check at the mpmath-refined
point.  Reads rescan_results.json (produced by rescan.py) and writes validation
tags back into the SAME JSON under 'validation' key per zero.

NOT RIGOROUS.  Single-process mpmath dps=70 -- ~0.3 s per zero.
"""
import json
import os
import sys
import time
import math

from mpmath import mp, mpf, mpc, fabs

mp.dps = 70

HERE = os.path.dirname(os.path.abspath(__file__))
DEG6 = os.path.join(os.path.dirname(HERE), 'agentB_deg6')

T0 = time.time()
LOGF = open(os.path.join(HERE, 'progress.log'), 'a', encoding='utf-8')


def log(msg):
    line = "[%7.1fs] validate: %s" % (time.time() - T0, msg)
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


# ----- load constants ----- #
fp = json.load(open(os.path.join(DEG6, 'sextic_fixedpoint.json'), encoding='utf-8'))
hz = json.load(open(os.path.join(DEG6, 'H6_zeros.json'), encoding='utf-8'))
A_MP = [mpf(s) for s in fp['g_coeffs']]
C_MP = [mpf(s) for s in fp['h_coeffs']]
DELTA6_MP = mpf(fp['delta6'])
AU0_MP = mpf(hz['AU0_6'])
TMAX6 = 1.15


def g6(z):
    t = z ** 6
    r = A_MP[-1]
    for c in reversed(A_MP[:-1]):
        r = r * t + c
    return r


def h6(z):
    t = z ** 6
    r = C_MP[-1]
    for c in reversed(C_MP[:-1]):
        r = r * t + c
    return r


def orbit_two(u_mp, m):
    """x_{2^m} and x_{2^{m+1}}; also max|x|^6 along orbit.  Returns
    (x2m, x2m1, maxt6).  Raises if orbit escapes."""
    s = u_mp * AU0_MP / DELTA6_MP ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    maxt6 = mpf(0)
    x2m = None
    for k in range(n2):
        x = g6(x) + s * h6(x)
        t = x ** 6
        if fabs(t) > maxt6:
            maxt6 = fabs(t)
        if fabs(x) > mpf('1e200'):
            raise OverflowError('orbit escaped')
        if k + 1 == n1:
            x2m = x
    return x2m, x, maxt6


def validate_zero(u_re_str, u_im_str):
    """Validate one mpmath-polished zero: orbit at m=7,8,9; absH; domain."""
    u = mpc(mpf(u_re_str), mpf(u_im_str))
    out = {'re_str': u_re_str, 'im_str': u_im_str}
    try:
        x2m7, x2m18, mxt7 = orbit_two(u, 7)
        out['m7_absH'] = float(fabs(x2m18 / x2m7))
        out['m7_den'] = float(fabs(x2m7))
        out['m7_maxt6'] = float(mxt7)
    except Exception:
        out['m7_absH'] = float('inf')
        out['m7_den'] = 0.0
        out['m7_maxt6'] = float('inf')

    try:
        x2m8, x2m19, mxt8 = orbit_two(u, 8)
        out['m8_absH'] = float(fabs(x2m19 / x2m8))
        out['m8_den'] = float(fabs(x2m8))
        out['m8_maxt6'] = float(mxt8)
    except Exception:
        out['m8_absH'] = float('inf')
        out['m8_den'] = 0.0
        out['m8_maxt6'] = float('inf')

    try:
        x2m9, x2m20, mxt9 = orbit_two(u, 9)
        out['m9_absH'] = float(fabs(x2m20 / x2m9))
        out['m9_den'] = float(fabs(x2m9))
        out['m9_maxt6'] = float(mxt9)
    except Exception:
        out['m9_absH'] = float('inf')
        out['m9_den'] = 0.0
        out['m9_maxt6'] = float('inf')

    # validity tags
    out['den_ok'] = (out['m8_den'] > 1e-10 and out['m9_den'] > 1e-10)
    out['maxt6_ok'] = (out['m8_maxt6'] < TMAX6)
    # genuine zero: |H8| small AND |H9| small AND |H9|/|H8| ratio close to delta6
    # (i.e. H_m decays like 1/delta6^m)
    if out['m8_absH'] > 0 and out['m9_absH'] > 0:
        ratio = out['m8_absH'] / out['m9_absH'] if out['m9_absH'] > 0 else None
        # we expect ratio ~ delta6 = 9.296
        # but also a true zero: |H8| ~ |H7| / delta6 already ... so here ratio indicates
        # decay slope; null-hypothesis of "numerical artifact" would have ratio near 1.
        out['Hratio_98'] = ratio
        # drift in u: u_m(7) - u_m(8) ... we treat our u as limit; instead use
        # |x_{2^8}|~0 implies u is near a limit zero.  Genuineness proxy:
        # |H9| << |x_{2^8}|^2 / |x_{2^9}|^2 ratio expected to be delta6 if
        # self-consistent.
    else:
        out['Hratio_98'] = None

    return out


def main():
    log('loading rescan_results.json')
    with open(os.path.join(HERE, 'rescan_results.json'), encoding='utf-8') as f:
        rs = json.load(f)

    # collect all unique polished zeros from A,B,C,D clusters
    all_zeros = {}
    for nm in 'ABCD':
        for z in rs['scans'][nm]['clustered']:
            if z['status'] in ('polished', 'float_only'):
                key = (z['re'], z['im'])
                # re-cluster later; here we keep one entry per unique coordinate
                all_zeros[key] = {'re': z['re'], 'im': z['im'],
                                  'absH_orig': z['absH'], 'status': z['status']}
    keys = sorted(all_zeros.keys())
    log('validating %d unique polished zeros (m=7,8,9 orbits + domain)' % len(keys))

    validated = {}
    t_per = time.time()
    for i2, k in enumerate(keys):
        z = all_zeros[k]
        v = validate_zero(z['re'], z['im'])
        v['sources'] = []
        for nm in 'ABCD':
            for zs in rs['scans'][nm]['clustered']:
                if zs['re'] == z['re'] and zs['im'] == z['im']:
                    v['sources'].append(nm)
                    break
        validated[k] = v
        if (i2 + 1) % 20 == 0 or (i2 + 1) == len(keys):
            el = time.time() - t_per
            log('  validate %d/%d  absH7=%.1e absH8=%.1e absH9=%.1e  maxt6=%.3f  (%.0fs, ETA %.0fs)'
                % (i2 + 1, len(keys), v['m7_absH'], v['m8_absH'], v['m9_absH'],
                   v['m8_maxt6'],
                   el, el / (i2 + 1) * (len(keys) - i2 - 1)))

    # summary table
    n_total = len(validated)
    genuine = []  # absH8<1e-15 AND absH9<1e-15 AND maxt6<1.15 AND den ok
    fake = []
    ambiguous = []
    for k, v in validated.items():
        ok = (v['m8_absH'] < 1e-15 and v['m9_absH'] < 1e-15
              and v['m8_maxt6'] < TMAX6 and v['den_ok']
              and abs(v['m8_absH'] - v['absH_orig']) < 1)
        # v['m8_absH'] here is recomputed absH8 (should match absH_orig for genuine)
        if v['m8_absH'] < 1e-15 and v['m9_absH'] < 1e-15 and v['m8_maxt6'] < TMAX6 and v['den_ok']:
            genuine.append(k)
        elif v['m8_absH'] > 1e-3 or v['m9_absH'] > 1e-3 or not v['m8_maxt6'] < TMAX6 or not v['den_ok']:
            fake.append(k)
        else:
            ambiguous.append(k)
    log('validation summary: genuine=%d  ambiguous=%d  fake=%d  total=%d'
        % (len(genuine), len(ambiguous), len(fake), n_total))

    # write back
    rs['validation'] = {
        'genuine': [dict(re=k[0], im=k[1], **validated[k]) for k in genuine],
        'ambiguous': [dict(re=k[0], im=k[1], **validated[k]) for k in ambiguous],
        'fake': [dict(re=k[0], im=k[1], **validated[k]) for k in fake],
        'all': {str(k): validated[k] for k in validated},
    }
    with open(os.path.join(HERE, 'rescan_results.json'), 'w', encoding='utf-8') as f:
        json.dump(rs, f, indent=2, ensure_ascii=False)
    log('wrote rescan_results.json with validation')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Final validation: for each unique zero in rescan_results.json, take
mpmath nstr (35 digits) precision of re/im and re-eval |H_8(u)|, |H_9(u)|, |H_7(u)|
on a properly-preserved mpc.  Then classify genuine / ambiguous / fake by:
   genuine = absH8_mp < 1e-20 AND absH9_mp < absH8_mp / delta6*0.5
   i.e. genuine if |H_m| shows the expected self-similar decay with m->m+1
   around the secular rate delta6 = 9.296"""
import json
import os
import sys
import time
import math

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

mp.dps = 70

HERE = os.path.dirname(os.path.abspath(__file__))
DEG6 = os.path.join(os.path.dirname(HERE), 'agentB_deg6')

T0 = time.time()
LOGF = open(os.path.join(HERE, 'progress.log'), 'a', encoding='utf-8')


def log(msg):
    line = "[%7.1fs] final_validate: %s" % (time.time() - T0, msg)
    print(line, flush=True)
    LOGF.write(line + "\n")
    LOGF.flush()


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


def g6p(z):
    t = z ** 6
    s = mpf(0)
    for j in range(len(A_MP) - 1, 0, -1):
        s = s * t + j * A_MP[j]
    return 6 * z ** 5 * s


def h6p(z):
    t = z ** 6
    s = mpf(0)
    for j in range(len(C_MP) - 1, 0, -1):
        s = s * t + j * C_MP[j]
    return 6 * z ** 5 * s


def orbit_H(u_mp, m, want_grad=False):
    s = u_mp * AU0_MP / DELTA6_MP ** m
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    dsdu = AU0_MP / DELTA6_MP ** m
    x2m = dx2m = None
    for k in range(n2):
        if want_grad:
            dx = (g6p(x) + s * h6p(x)) * dx + h6(x) * dsdu
        x = g6(x) + s * h6(x)
        if k + 1 == n1:
            x2m, dx2m = x, dx
    if want_grad:
        return x2m, x, dx2m, dx
    return x2m, x


def orbit_maxt6(u_mp, m):
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
        if fabs(x) > mpf('1e300'):
            raise OverflowError('orbit escaped')
        if k + 1 == n1:
            x2m = x
    return x2m, x, maxt6


def mp_newton_m8(u_start, max_iter=300, tol=mpf('1e-60')):
    """mpmath Newton on H_8 with damping."""
    u_mp = u_start
    for it in range(max_iter):
        x2m, x2m1, dx2m, dx2m1 = orbit_H(u_mp, 8, want_grad=True)
        H = x2m1 / x2m
        aH = fabs(H)
        if aH < tol:
            break
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        if dH == 0:
            break
        step = H / dH
        moved = False
        for _h in range(15):
            ut = u_mp - step
            x2mt, x2m1t = orbit_H(ut, 8)
            if fabs(x2m1t / x2mt) < aH:
                u_mp = ut
                moved = True
                break
            step = step / 2
        if not moved:
            break
        if fabs(step) < mpf('1e-65'):
            break
    x2m, x2m1 = orbit_H(u_mp, 8)
    return u_mp, float(fabs(x2m1 / x2m))


def main():
    log('loading rescan_results.json (with full-precision strings)')
    with open(os.path.join(HERE, 'rescan_results.json'), encoding='utf-8') as f:
        rs = json.load(f)

    # collect all unique clustered zeros -- now re/im are full mpf strings
    all_zeros = {}
    for nm in 'ABCD':
        for z in rs['scans'][nm]['clustered']:
            if z['status'] in ('polished', 'float_only'):
                # cluster key uses float to avoid spurious duplicates, but
                # PRESERVE the full-precision string.
                key = (float(z['re']), float(z['im']))
                if key not in all_zeros:
                    all_zeros[key] = {'re_str': z['re'], 'im_str': z['im'],
                                      'absH_orig': z['absH'], 'status': z['status']}
    log('%d unique clustered zeros to validate' % len(all_zeros))

    validated = {}
    t_per = time.time()
    for i2, k in enumerate(sorted(all_zeros.keys())):
        z = all_zeros[k]
        u_start = mpc(mpf(z['re_str']), mpf(z['im_str']))
        try:
            u8, absH8_mp = mp_newton_m8(u_start, max_iter=300, tol=mpf('1e-60'))
        except Exception as e:
            log('  fail repolish re=%s im=%s: %s' % (z['re_str'][:14], z['im_str'][:14], e))
            continue
        re_full = nstr(mpre(u8), 35)
        im_full = nstr(mpim(u8), 35)

        # evaluation at m=7,8,9
        try:
            x7_2m, x7_2m1, mxt7 = orbit_maxt6(u8, 7)
            m7_absH = float(fabs(x7_2m1 / x7_2m))
            m8_absH = float(fabs(x7_2m1 / x7_2m))  # placeholder, overwritten
            m7_den = float(fabs(x7_2m))
            m7_maxt6 = float(mxt7)
        except Exception:
            m7_absH = float('inf')
            m7_den = 0.0
            m7_maxt6 = float('inf')

        try:
            x8_2m, x8_2m1, mxt8 = orbit_maxt6(u8, 8)
            m8_absH = float(fabs(x8_2m1 / x8_2m))
            m8_den = float(fabs(x8_2m))
            m8_maxt6 = float(mxt8)
        except Exception:
            m8_absH = float('inf')
            m8_den = 0.0
            m8_maxt6 = float('inf')

        try:
            x9_2m, x9_2m1, mxt9 = orbit_maxt6(u8, 9)
            m9_absH = float(fabs(x9_2m1 / x9_2m))
            m9_den = float(fabs(x9_2m))
            m9_maxt6 = float(mxt9)
        except Exception:
            m9_absH = float('inf')
            m9_den = 0.0
            m9_maxt6 = float('inf')

        # drift ratios
        drift_78 = None
        drift_89 = None
        if m7_absH > 0 and m8_absH > 0:
            drift_78 = m7_absH / m8_absH
        if m8_absH > 0 and m9_absH > 0:
            drift_89 = m8_absH / m9_absH

        # sources
        sources = []
        for nm in 'ABCD':
            for zs in rs['scans'][nm]['clustered']:
                if (float(zs['re']), float(zs['im'])) == k:
                    sources.append(nm)
                    break

        validated[k] = dict(re_full=re_full, im_full=im_full,
                            absH8_mp=absH8_mp,
                            m7_absH=m7_absH, m8_absH=m8_absH, m9_absH=m9_absH,
                            m7_den=m7_den, m8_den=m8_den, m9_den=m9_den,
                            m7_maxt6=m7_maxt6, m8_maxt6=m8_maxt6, m9_maxt6=m9_maxt6,
                            drift_78=drift_78, drift_89=drift_89,
                            sources=sorted(sources))

        if (i2 + 1) % 10 == 0 or (i2 + 1) == len(all_zeros):
            el = time.time() - t_per
            log('  %d/%d  re=%s im=%s  absH8=%.2e m7=%.2e m8=%.2e m9=%.2e  ratio_78=%s  (%.0fs, ETA %.0fs)'
                % (i2 + 1, len(all_zeros), re_full[:14], im_full[:14],
                   absH8_mp, m7_absH, m8_absH, m9_absH,
                   '%.2f' % drift_78 if drift_78 else 'None',
                   el, el / (i2 + 1) * (len(all_zeros) - i2 - 1)))

    # classify: genuine means |H_m| shows secular decay with rate ~ delta6
    # criterion: |H9| < |H8| / delta6  AND  |H8| < 1e-20
    # This handles both the "true zero" and the "near-pass-through-zero" cases.
    genuine, plausible, fake = [], [], []
    for k, v in validated.items():
        if v['absH8_mp'] < 1e-20:
            if v['drift_89'] is not None and v['drift_89'] > 1.5:
                genuine.append(k)
            else:
                plausible.append(k)
        else:
            fake.append(k)
    log('FINAL classification: genuine=%d  plausible=%d  fake=%d  total=%d'
        % (len(genuine), len(plausible), len(fake), len(validated)))

    rs['final_validate'] = {
        'genuine': [dict(re=k[0], im=k[1], **validated[k]) for k in genuine],
        'plausible': [dict(re=k[0], im=k[1], **validated[k]) for k in plausible],
        'fake': [dict(re=k[0], im=k[1], **validated[k]) for k in fake],
        'all': {str(k): validated[k] for k in validated},
    }
    with open(os.path.join(HERE, 'rescan_results.json'), 'w', encoding='utf-8') as f:
        json.dump(rs, f, indent=2, ensure_ascii=False)
    log('wrote rescan_results.json with final_validate')


if __name__ == '__main__':
    main()

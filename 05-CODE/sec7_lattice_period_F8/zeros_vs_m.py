#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 Task 4 (rewritten for speed): zeros of H_m (d=6) for late zeros at m=7,8,9,10.

Optimization: dps=40, no full line search (simple Newton with damped halving).
Start from previous level's refined zero for warm-start.
"""
import os
import json
import time
import math
import numpy as np
from mpmath import mp, mpf, mpc, nstr, re as mpre, im as mpim, fabs

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'
T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')
mp.dps = 40


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def dist_to_lattice16(re, im):
    re_m = float(re); im_m = float(im)
    k_re = round(re_m * 16) / 16.0
    k_im = round(im_m * 16) / 16.0
    d = math.hypot(re_m - k_re, im_m - k_im)
    return d, k_re * 16, k_im * 16, k_re, k_im


# Polynomial evaluators using jit-friendly mpmath
def make_g_h(g_coeffs_str, h_coeffs_str, d):
    g_coeffs = [mpf(s) for s in g_coeffs_str]
    h_coeffs = [mpf(s) for s in h_coeffs_str]

    def g(z):
        t = z ** d
        r = g_coeffs[-1]
        for c in reversed(g_coeffs[:-1]):
            r = r * t + c
        return r

    def h(z):
        t = z ** d
        r = h_coeffs[-1]
        for c in reversed(h_coeffs[:-1]):
            r = r * t + c
        return r

    def g1(z):
        t = z ** d
        s = mpf(0)
        for j in range(len(g_coeffs) - 1, 0, -1):
            s = s * t + j * g_coeffs[j]
        return d * z ** (d - 1) * s

    def h1(z):
        t = z ** d
        s = mpf(0)
        for j in range(len(h_coeffs) - 1, 0, -1):
            s = s * t + j * h_coeffs[j]
        return d * z ** (d - 1) * s
    return g, h, g1, h1


def orbit_H_at_m(u_mp, m, AU0_mp, delta_mp, g, h, g1, h1):
    s = u_mp * AU0_mp / delta_mp ** m
    n1 = 2 ** m; n2 = 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    x2m = None; dx2m = None
    dsdu = AU0_mp / delta_mp ** m
    for k in range(n2):
        dx = (g1(x) + s * h1(x)) * dx + h(x) * dsdu
        x = g(x) + s * h(x)
        if k + 1 == n1:
            x2m = x; dx2m = dx
    return x2m, x, dx2m, dx


def newton_zero(u_start_re, u_start_im, m, AU0_mp, delta_mp, g, h, g1, h1,
                max_iter=6, tol=mpf('1e-25')):
    """Plain Newton with damped halving; no line search (each Newton step = 1 orbit eval)."""
    u_mp = mpc(mpf(u_start_re), mpf(u_start_im))
    last_H = None
    history = []
    for it in range(max_iter):
        x2m, x2m1, dx2m, dx2m1 = orbit_H_at_m(u_mp, m, AU0_mp, delta_mp, g, h, g1, h1)
        if x2m == 0:
            return u_mp, float('nan'), history
        H = x2m1 / x2m
        last_H = fabs(H)
        history.append((float(mpre(u_mp)), float(mpim(u_mp)), float(last_H)))
        if last_H < tol:
            return u_mp, float(last_H), history
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        if dH == 0:
            return u_mp, float(last_H), history
        step = H / dH
        if fabs(step) < mpf('1e-42'):
            return u_mp, float(last_H), history
        # Damped Newton: try full step, then half-step if |H| increases
        u_new = u_mp - step
        x2mt, x2m1t, _, _ = orbit_H_at_m(u_new, m, AU0_mp, delta_mp, g, h, g1, h1)
        H_new = None
        if x2mt != 0:
            H_new = fabs(x2m1t / x2mt)
        if H_new is not None and H_new < last_H:
            u_mp = u_new
            continue
        # backtrack
        s_step = step
        moved = False
        for _h in range(8):
            ut = u_mp - s_step
            x2mt, x2m1t, _, _ = orbit_H_at_m(ut, m, AU0_mp, delta_mp, g, h, g1, h1)
            if x2mt == 0:
                s_step = s_step / 2
                continue
            Ht = fabs(x2m1t / x2mt)
            if Ht < last_H:
                u_mp = ut
                last_H = Ht
                moved = True
                break
            s_step = s_step / 2
        if not moved:
            break
    return u_mp, last_H, history


def main():
    log('--- zeros_vs_m (v2): load ---')
    fp = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'sextic_fixedpoint.json'), encoding='utf-8'))
    hz = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'H6_zeros.json'), encoding='utf-8'))
    g, h, g1, h1 = make_g_h(fp['g_coeffs'], fp['h_coeffs'], 6)
    delta6 = mpf(fp['delta6'])
    AU0 = mpf(hz['AU0_6'])
    log('delta6=%.20f  AU0=%.20f  dps=%d' % (float(delta6), float(AU0), mp.dps))

    late_targets = [
        (1.0000, 14.5),
        (20.125, 19.25),
        (41.3125, 10.0),
        (-6.5, 46.9375),
    ]

    cat = hz['zeros_level8']

    out = dict(meta=dict(delta6=float(delta6), AU0=float(AU0),
                          levels=[7, 8, 9, 10], dps=mp.dps))
    zeros_by_target = []

    for tre, tim in late_targets:
        dists = [(z['re'] - tre) ** 2 + (z['im'] - tim) ** 2 for z in cat]
        idx = int(np.argmin(dists))
        z_cat = cat[idx]
        u0_re, u0_im = z_cat['re'], z_cat['im']
        d_lattice, k_re, k_im, re_node, im_node = dist_to_lattice16(u0_re, u0_im)
        log('=== target=(%.4f,%.4f)  cat_m8=(%.10f,%.10f)  dist_lattice=%.3e  (k_re=%d k_im=%d)'
            % (tre, tim, u0_re, u0_im, d_lattice, k_re, k_im))

        target_results = dict(
            target=(tre, tim),
            cat_start=(u0_re, u0_im, z_cat['absH']),
            node=(float(re_node), float(im_node), int(k_re), int(k_im)),
            by_m=[],
        )

        prev_u = (u0_re, u0_im)
        for m in [7, 8, 9, 10]:
            t0 = time.time()
            u_mp, absH, hist = newton_zero(prev_u[0], prev_u[1], m, AU0, delta6,
                                           g, h, g1, h1, max_iter=6, tol=mpf('1e-10'))
            re_f = float(mpre(u_mp)); im_f = float(mpim(u_mp))
            d_l, k_re_l, k_im_l, re_node_l, im_node_l = dist_to_lattice16(re_f, im_f)
            target_results['by_m'].append(dict(
                m=m,
                u=(re_f, im_f),
                u_str=(nstr(mpre(u_mp), 30), nstr(mpim(u_mp), 30)),
                absH=absH,
                dist_to_node=d_l,
                node=(float(re_node_l), float(im_node_l), int(k_re_l), int(k_im_l)),
                n_iters=len(hist),
                elapsed=time.time() - t0,
            ))
            log('  m=%d  u=(%.18f, %.18f)  |H|=%.2e  dist=%.4e  iters=%d (%.1fs)'
                % (m, re_f, im_f, absH, d_l, len(hist), time.time() - t0))
            prev_u = (re_f, im_f)
        zeros_by_target.append(target_results)
    out['targets'] = zeros_by_target

    # Drift analysis
    log('--- drift analysis ---')
    drift_table = []
    for tr in zeros_by_target:
        bms = tr['by_m']
        drifts = []
        for i in range(len(bms) - 1):
            d = math.hypot(bms[i + 1]['u'][0] - bms[i]['u'][0],
                            bms[i + 1]['u'][1] - bms[i]['u'][1])
            drifts.append(d)
        one_over_delta = 1.0 / float(delta6)
        drift_table.append(dict(target=tr['target'], drifts=drifts,
                                drift_times_delta=[d * float(delta6) for d in drifts],
                                expected_one_over_delta=one_over_delta))
        log('  target=%s  drifts=%s  drifts*delta=%s  expected 1/delta=%.4f'
            % (str(tr['target']), drifts,
                ['%.3e' % (d * float(delta6)) for d in drifts], one_over_delta))
    out['drift'] = drift_table

    with open(os.path.join(OUT_DIR, 'zeros_vs_m.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote zeros_vs_m.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('MKL_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
    main()
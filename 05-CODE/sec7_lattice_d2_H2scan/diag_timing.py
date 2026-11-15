# -*- coding: utf-8 -*-
"""Диагностика тайминга: float64-Ньютон (батч, numpy, с градиентом) и
mpmath-полировка на одном известном нуле, чтобы откалибровать бюджет
перед основным прогоном."""
import time
import numpy as np
from mpmath import mp, mpf, mpc, fabs

import h2_common as C

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
known, meta = C.load_known_zeros()

# deriv coeff arrays (numpy)
tc_g_d_np = np.array([(k + 1) * tc_g_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_h_d_np = np.array([(k + 1) * tc_h_np[k + 1] for k in range(N - 1)], dtype=np.float64)


def scan_grad(grid, m):
    s = (C.AU0_F / (C.DELTA_F ** m)) * grid
    dsdu = C.AU0_F / (C.DELTA_F ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    dx = np.zeros(grid.shape, dtype=np.complex128)
    xden = dxden = None
    with np.errstate(over='ignore', invalid='ignore'):
        for k in range(n2):
            xx = x * x - 1.0
            w = xx / C.R_DISC_F
            g = C.horner_np(tc_g_np, w)
            h = C.horner_np(tc_h_np, w)
            dg_w = C.horner_np(tc_g_d_np, w)
            dh_w = C.horner_np(tc_h_d_np, w)
            dw_dx = 2.0 * x / C.R_DISC_F
            dg = dg_w * dw_dx
            dh = dh_w * dw_dx
            dx_new = (dg + s * dh) * dx + h * dsdu
            x = g + s * h
            dx = dx_new
            if k + 1 == n1:
                xden, dxden = x.copy(), dx.copy()
    return x, xden, dx, dxden


# perturb known zeros slightly off-target to simulate "grid seed" and time batch Newton
seeds = np.array([u + complex(0.03, -0.02) for _, u in known], dtype=np.complex128)
log('batch newton (numpy, m=8) on %d seeds, 15 iterations, timing...' % len(seeds))
u = seeds.copy()
t0 = time.time()
for it in range(15):
    xnum, xden, dxnum, dxden = scan_grad(u, 8)
    H = xnum / xden
    dH = (dxnum * xden - xnum * dxden) / (xden * xden)
    with np.errstate(divide='ignore', invalid='ignore'):
        step = H / dH
    step = np.where(np.isfinite(step), step, 0.0)
    mag = np.abs(step)
    cap = 2.0
    scale = np.where(mag > cap, cap / np.maximum(mag, 1e-300), 1.0)
    u = u - step * scale
    if (it + 1) % 5 == 0:
        log('  it=%d max|H|=%.3e median|H|=%.3e' % (it + 1, np.max(np.abs(H)), np.median(np.abs(H))))
dt = time.time() - t0
log('batch newton done in %.2fs for %d seeds x 15 iters (m=8)' % (dt, len(seeds)))
errs = np.abs(u - np.array([uu for _, uu in known]))
log('final distance to true known zeros: max=%.3e median=%.3e' % (errs.max(), np.median(errs)))

# mpmath polish timing on ONE point
log('--- mpmath polish timing test (mp.dps=30, m=8) ---')
mp.dps = 30
u0 = known[0][1]


def orbit_grad_mp(u_mp, m):
    s_mp = u_mp * mpf(C.AU0_STR) / (mpf(C.DELTA_STR) ** m)
    dsdu = mpf(C.AU0_STR) / (mpf(C.DELTA_STR) ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    xden = dxden = None
    tg = [mpf(str(c)) for c in tc_g_np]
    th = [mpf(str(c)) for c in tc_h_np]
    tgd = [mpf(str(c)) for c in tc_g_d_np]
    thd = [mpf(str(c)) for c in tc_h_d_np]
    rdisc = mpf(C.R_DISC_STR)
    for k in range(n2):
        xx = x * x - 1
        w = xx / rdisc
        g = tg[-1]
        for c in tg[-2::-1]:
            g = g * w + c
        h = th[-1]
        for c in th[-2::-1]:
            h = h * w + c
        dg_w = tgd[-1]
        for c in tgd[-2::-1]:
            dg_w = dg_w * w + c
        dh_w = thd[-1]
        for c in thd[-2::-1]:
            dh_w = dh_w * w + c
        dw_dx = 2 * x / rdisc
        dg = dg_w * dw_dx
        dh = dh_w * dw_dx
        dx = (dg + s_mp * dh) * dx + h * dsdu
        x = g + s_mp * h
        if k + 1 == n1:
            xden, dxden = x, dx
    return xden, x, dxden, dx


t0 = time.time()
xden, xnum, dxden, dxnum = orbit_grad_mp(mpc(u0.real, u0.imag), 8)
dt1 = time.time() - t0
log('ONE mpmath orbit+grad eval (m=8, dps=30) took %.3fs' % dt1)
log('estimated for say 6 iters x up to 300 candidates: %.1fs' % (dt1 * 6 * 300))

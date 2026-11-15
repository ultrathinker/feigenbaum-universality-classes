# -*- coding: utf-8 -*-
"""Диагностика-2: (a) сколько итераций float64-Ньютона нужно до предела
точности double; (b) тайминг mpmath-полировки с ОДНОКРАТНЫМ переводом
коэффициентов в mpf (без перепарсинга на каждом вызове)."""
import time
import numpy as np
from mpmath import mp, mpf, mpc, fabs

import h2_common as C

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
known, meta = C.load_known_zeros()

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


seeds = np.array([u + complex(0.03, -0.02) for _, u in known], dtype=np.complex128)
u = seeds.copy()
t0 = time.time()
for it in range(30):
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
errs = np.abs(u - np.array([uu for _, uu in known]))
log('after 30 float64 Newton iters: max dist=%.3e median=%.3e (%.2fs)' % (errs.max(), np.median(errs), time.time() - t0))

# ---------------- mpmath, coefficients precomputed ONCE ----------------
log('--- mpmath polish, coeffs precomputed once, mp.dps=30 ---')
mp.dps = 30
RDISC_MP = mpf(C.R_DISC_STR)
AU0_MP2 = mpf(C.AU0_STR)
DELTA_MP2 = mpf(C.DELTA_STR)
TG = tc_g_mp
TH = tc_h_mp
TGD = [(k + 1) * tc_g_mp[k + 1] for k in range(N - 1)]
THD = [(k + 1) * tc_h_mp[k + 1] for k in range(N - 1)]


def orbit_grad_mp(u_mp, m):
    s_mp = u_mp * AU0_MP2 / (DELTA_MP2 ** m)
    dsdu = AU0_MP2 / (DELTA_MP2 ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    xden = dxden = None
    for k in range(n2):
        xx = x * x - 1
        w = xx / RDISC_MP
        g = TG[-1]
        for c in TG[-2::-1]:
            g = g * w + c
        h = TH[-1]
        for c in TH[-2::-1]:
            h = h * w + c
        dg_w = TGD[-1]
        for c in TGD[-2::-1]:
            dg_w = dg_w * w + c
        dh_w = THD[-1]
        for c in THD[-2::-1]:
            dh_w = dh_w * w + c
        dw_dx = 2 * x / RDISC_MP
        dg = dg_w * dw_dx
        dh = dh_w * dw_dx
        dx = (dg + s_mp * dh) * dx + h * dsdu
        x = g + s_mp * h
        if k + 1 == n1:
            xden, dxden = x, dx
    return xden, x, dxden, dx


u0 = known[0][1]
t0 = time.time()
xden, xnum, dxden, dxnum = orbit_grad_mp(mpc(u0.real, u0.imag), 8)
dt1 = time.time() - t0
log('ONE mpmath orbit+grad eval (m=8, dps=30, precomputed coeffs) took %.3fs' % dt1)
log('estimated for 3 iters x 300 candidates: %.1fs' % (dt1 * 3 * 300))
log('estimated for 3 iters x 150 candidates: %.1fs' % (dt1 * 3 * 150))

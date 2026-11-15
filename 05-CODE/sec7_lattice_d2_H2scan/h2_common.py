# -*- coding: utf-8 -*-
"""Общая библиотека для скана H_2 (класс d=2).

РЕЦЕПТ (проверен на flint 512 бит, воспроизводит H(u_actual)~1e-4..1e-5 при m=7,8
для 23 известных нулей -- см. progress.log сравнение с провалившейся попыткой
MiniMax, которая использовала НЕВЕРНЫЙ рецепт t=x^2 из spectra_grid.json вместо
правильного рецепта из address_certificate.json):

    x_{n+1} = g(w) + s*h(w),   w = (x^2 - 1) / r_disc,   x_0 = 0
    s = u * AU0 / delta^m
    H_m(u) = x_{2^{m+1}}(u) / x_{2^m}(u)

    g(w) = sum_j tc_g[j] w^j ,  h(w) = sum_j tc_h[j] w^j   (degree 99 каждый)

tc_g, tc_h -- ряды Тейлора, полученные из чебышёвских коэффициентов
w33_gh_hiprec.json (M=100, g_str/h_str) через cheb_to_taylor (рецепт
w35/x4_walk_sonnet/x4_common.py): tcoeff -- таблица преобразования Чебышёва
в степенной базис (T_0=1, T_1=w, T_n=2 w T_{n-1} - T_{n-2}); taylor[j] =
sum_{k=j}^{N-1} arb[k] * tcoeff[k][j] * r_disc^j.

Константы (address_certificate.json convention, ТОЛЬКО ЧТЕНИЕ):
    AU0    = 0.41403799908627305889537448439555967803329915886199
    delta  = 4.669201609102990671853203820466201617258185577475768632745651343004133555868548426390712381266351242
    r_disc = 2.5
"""
import json
import os
import time

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')

import numpy as np  # noqa: E402
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim  # noqa: E402

mp.dps = 60

HERE = os.path.dirname(os.path.abspath(__file__))
W33_PATH = r'D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\archive_inputs\w33\t5_cert_sonnet\w33_gh_hiprec.json'
RESULTS_PATH = r'../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json'
CERT_PATH = r'D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\archive_inputs\path_to_8\point15_layer1_addressing\address_certificate.json'

AU0_STR = '0.41403799908627305889537448439555967803329915886199'
DELTA_STR = ('4.669201609102990671853203820466201617258185577475768632745651343004133555868548'
             '426390712381266351242')
R_DISC_STR = '2.5'

AU0_MP = mpf(AU0_STR)
DELTA_MP = mpf(DELTA_STR)
R_DISC_MP = mpf(R_DISC_STR)
AU0_F = float(AU0_STR)
DELTA_F = float(DELTA_STR)
R_DISC_F = 2.5


def _log_stub(msg):
    print(msg, flush=True)


def load_gh_taylor(log=_log_stub):
    """Возвращает (tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N)."""
    t0 = time.time()
    D = json.load(open(W33_PATH, encoding='utf-8'))
    g_arb = [mpf(s) for s in D['g_str']]
    h_arb = [mpf(s) for s in D['h_str']]
    N = len(g_arb)
    assert len(h_arb) == N
    log('  w33_gh_hiprec.json: M(JSON)=%s, N(coeffs)=%d, delta_rayleigh matches cert: %s'
        % (D.get('M'), N, abs(mpf(D['delta_rayleigh']) - DELTA_MP) < mpf('1e-40')))

    tcoeff = [[mpf(0)] * N for _ in range(N)]
    tcoeff[0][0] = mpf(1)
    if N >= 2:
        tcoeff[1][1] = mpf(1)
    for n in range(2, N):
        row_prev = tcoeff[n - 1]
        row_prev2 = tcoeff[n - 2]
        row = tcoeff[n]
        for j in range(n + 1):
            val = mpf(0)
            if j > 0:
                val += 2 * row_prev[j - 1]
            if j <= n - 2:
                val -= row_prev2[j]
            row[j] = val

    def cheb_to_taylor(arr):
        taylor = []
        for j in range(N):
            s = mpf(0)
            rj = R_DISC_MP ** j
            for k in range(j, N):
                s += arr[k] * tcoeff[k][j] * rj
            taylor.append(s)
        return taylor

    tc_g_mp = cheb_to_taylor(g_arb)
    tc_h_mp = cheb_to_taylor(h_arb)
    tc_g_np = np.array([float(c) for c in tc_g_mp], dtype=np.float64)
    tc_h_np = np.array([float(c) for c in tc_h_mp], dtype=np.float64)
    log('  cheb_to_taylor done (%.1fs), N=%d, tc_g[0]=%.6f tc_h[0]=%.6f'
        % (time.time() - t0, N, tc_g_np[0], tc_h_np[0]))
    return tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N


def load_known_zeros():
    """23 известных нуля H_2 (u_actual, id) из results.json (Z72 -> corrected)."""
    R = json.load(open(RESULTS_PATH, encoding='utf-8'))
    catalog = R['stage3_catalog']
    z72 = R['stage3b_z72_corrected']
    known = []
    for z in catalog:
        if z['id'] == 'Z72':
            u = z72['u_actual']
        else:
            u = z['u_actual']
        known.append((z['id'], complex(float(u['re_str']), float(u['im_str']))))
    assert len(known) == 23
    return known, R['meta']


# ----------------------------------------------------------------------- #
#  float64/numpy vectorized scan
# ----------------------------------------------------------------------- #
def horner_np(coeffs, w):
    """coeffs[0] + coeffs[1] w + ... (real float64 array), w complex128 array."""
    r = np.full(w.shape, coeffs[-1], dtype=np.complex128)
    for c in coeffs[-2::-1]:
        r = r * w + c
    return r


def scan_level(grid, m, tc_g_np, tc_h_np, tc_g_np_tr=None, tc_h_np_tr=None):
    """H_m(u) на всей сетке (complex128), плюс диагностика (maxw, maxx, trunc err)."""
    s = (AU0_F / (DELTA_F ** m)) * grid
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    maxw = np.zeros(grid.shape, dtype=np.float64)
    maxx = np.zeros(grid.shape, dtype=np.float64)
    errmax = np.zeros(grid.shape, dtype=np.float64)
    xden = None
    with np.errstate(over='ignore', invalid='ignore'):
        for k in range(n2):
            xx = x * x - 1.0
            w = xx / R_DISC_F
            aw = np.abs(w)
            np.maximum(maxw, np.where(np.isfinite(aw), aw, maxw), out=maxw)
            g = horner_np(tc_g_np, w)
            h = horner_np(tc_h_np, w)
            if tc_g_np_tr is not None:
                gtr = horner_np(tc_g_np_tr, w)
                err = np.abs(g - gtr)
                np.maximum(errmax, np.where(np.isfinite(err), err, errmax), out=errmax)
            x = g + s * h
            ax = np.abs(x)
            np.maximum(maxx, np.where(np.isfinite(ax), ax, maxx), out=maxx)
            if k + 1 == n1:
                xden = x.copy()
        num = np.abs(x)
        den = np.abs(xden)
        with np.errstate(divide='ignore', invalid='ignore'):
            Habs = num / den
    return Habs, xden, x, maxw, maxx, errmax


def detect_local_minima(Habs, valid, xs, ys, thresh):
    Hm = np.where(valid, Habs, np.inf)
    if Hm.size < 9:
        return []
    P = np.full((Hm.shape[0] + 2, Hm.shape[1] + 2), np.inf)
    P[1:-1, 1:-1] = Hm
    core = P[1:-1, 1:-1]
    nb = np.stack([P[:-2, :-2], P[:-2, 1:-1], P[:-2, 2:],
                   P[1:-1, :-2], P[1:-1, 2:],
                   P[2:, :-2], P[2:, 1:-1], P[2:, 2:]])
    ismin = (core <= nb).all(axis=0) & (core < thresh)
    ii, jj = np.nonzero(ismin)
    return [(float(xs[i]), float(ys[j]), float(core[i, j])) for i, j in zip(ii, jj)]


# ----------------------------------------------------------------------- #
#  plain-python complex evaluator (Newton refine, fast-ish)
# ----------------------------------------------------------------------- #
def make_py_evaluators(tc_g_np, tc_h_np):
    G = [float(c) for c in tc_g_np]
    H = [float(c) for c in tc_h_np]
    NG = len(G) - 1
    NH = len(H) - 1

    def gh_and_deriv(x):
        """g(w), h(w), dg/dx, dh/dx at point x (python complex)."""
        if abs(x) > 1e8:
            raise OverflowError('escaped')
        xx = x * x - 1.0
        w = xx / R_DISC_F
        g = G[NG]
        for c in G[NG - 1::-1]:
            g = g * w + c
        h = H[NH]
        for c in H[NH - 1::-1]:
            h = h * w + c
        # dg/dw, dh/dw via derivative Horner
        dg_w = 0.0
        for j in range(NG, 0, -1):
            dg_w = dg_w * w + j * G[j]
        dh_w = 0.0
        for j in range(NH, 0, -1):
            dh_w = dh_w * w + j * H[j]
        dw_dx = 2.0 * x / R_DISC_F
        dg = dg_w * dw_dx
        dh = dh_w * dw_dx
        return g, h, dg, dh

    def H_at(u, m):
        s = u * AU0_F / (DELTA_F ** m)
        x = 0j
        n1 = 2 ** m
        n2 = 2 ** (m + 1)
        xden = None
        try:
            for k in range(n2):
                g, h, _dg, _dh = gh_and_deriv(x)
                x = g + s * h
                if k + 1 == n1:
                    xden = x
        except OverflowError:
            return None
        if xden is None or abs(xden) < 1e-290:
            return None
        return x / xden

    def newton_refine(u0, m, iters=80, rmax=60.0, tol=1e-10):
        s_coef = AU0_F / (DELTA_F ** m)
        n1 = 2 ** m
        n2 = 2 ** (m + 1)
        u = u0
        Hval = None
        for _ in range(iters):
            s = u * s_coef
            x = 0j
            dx = 0j
            xden = None
            dxden = None
            try:
                for k in range(n2):
                    g, h, dg, dh = gh_and_deriv(x)
                    dx = (dg + s * dh) * dx + h * s_coef
                    x = g + s * h
                    if k + 1 == n1:
                        xden, dxden = x, dx
            except OverflowError:
                return None, float('inf')
            if xden is None or abs(xden) < 1e-290 or not np.isfinite(abs(x) + abs(dx)):
                return None, float('inf')
            Hval = x / xden
            dH = (dx * xden - x * dxden) / (xden * xden)
            if abs(Hval) < tol:
                break
            if dH == 0 or not np.isfinite(abs(dH)):
                return None, float('inf')
            step = Hval / dH
            moved = False
            for _h in range(20):
                u_try = u - step
                if abs(u_try) > rmax:
                    step = step / 2
                    continue
                Ht = H_at(u_try, m)
                if Ht is not None and abs(Ht) < abs(Hval):
                    u = u_try
                    moved = True
                    break
                step = step / 2
            if not moved:
                break
            if abs(step) < 1e-15:
                break
        Hf = H_at(u, m)
        if Hf is None:
            return None, float('inf')
        return u, abs(Hf)

    def domain_check(u, m):
        s = u * AU0_F / (DELTA_F ** m)
        x = 0j
        n1 = 2 ** m
        n2 = 2 ** (m + 1)
        maxx = 0.0
        maxw = 0.0
        xden = None
        try:
            for k in range(n2):
                g, h, _dg, _dh = gh_and_deriv(x)
                w = (x * x - 1.0) / R_DISC_F
                aw = abs(w)
                if aw > maxw:
                    maxw = aw
                x = g + s * h
                ax = abs(x)
                if ax > maxx:
                    maxx = ax
                if k + 1 == n1:
                    xden = x
        except OverflowError:
            return float('inf'), float('inf'), 0.0
        return maxx, maxw, (abs(xden) if xden is not None else 0.0)

    return H_at, newton_refine, domain_check


# ----------------------------------------------------------------------- #
#  mpmath orbit + polish
# ----------------------------------------------------------------------- #
def orbit_H_mpmath(u_mp, m, tc_g_mp, tc_h_mp, want_grad=False):
    N = len(tc_g_mp)
    s_mp = u_mp * AU0_MP / (DELTA_MP ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    xden = dxden = None
    for k in range(n2):
        xx = x * x - 1
        w = xx / R_DISC_MP
        g = tc_g_mp[N - 1]
        for c in tc_g_mp[N - 2::-1]:
            g = g * w + c
        h = tc_h_mp[N - 1]
        for c in tc_h_mp[N - 2::-1]:
            h = h * w + c
        if want_grad:
            dg_w = mpf(0)
            for j in range(N - 1, 0, -1):
                dg_w = dg_w * w + j * tc_g_mp[j]
            dh_w = mpf(0)
            for j in range(N - 1, 0, -1):
                dh_w = dh_w * w + j * tc_h_mp[j]
            dw_dx = 2 * x / R_DISC_MP
            dg = dg_w * dw_dx
            dh = dh_w * dw_dx
            dx = (dg + s_mp * dh) * dx + h * (AU0_MP / (DELTA_MP ** m))
        x = g + s_mp * h
        if k + 1 == n1:
            xden, dxden = x, dx
    if want_grad:
        return xden, x, dxden, dx
    return xden, x


def mp_polish(u_start, m, tc_g_mp, tc_h_mp, max_iter=60, tol=mpf('1e-40')):
    u_mp = mpc(mpf(str(u_start.real)), mpf(str(u_start.imag)))
    last_h = None
    for _it in range(max_iter):
        xden, xnum, dxden, dxnum = orbit_H_mpmath(u_mp, m, tc_g_mp, tc_h_mp, want_grad=True)
        Hv = xnum / xden
        last_h = fabs(Hv)
        if last_h < tol:
            break
        dH = (dxnum * xden - xnum * dxden) / (xden * xden)
        if dH == 0:
            break
        step = Hv / dH
        moved = False
        for _h in range(15):
            ut = u_mp - step
            xden_t, xnum_t = orbit_H_mpmath(ut, m, tc_g_mp, tc_h_mp)
            ht = fabs(xnum_t / xden_t)
            if ht < last_h:
                u_mp = ut
                moved = True
                break
            step = step / 2
        if not moved:
            break
        if fabs(step) < mpf('1e-45'):
            break
    xden, xnum = orbit_H_mpmath(u_mp, m, tc_g_mp, tc_h_mp)
    return u_mp, float(fabs(xnum / xden))

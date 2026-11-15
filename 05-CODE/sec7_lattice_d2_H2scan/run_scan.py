# -*- coding: utf-8 -*-
"""S2/H2SCAN -- Задание 1+2: вычислитель H_2, валидация на 23 известных нулях,
расширенный скан Re u in [-16,16], Im u in [0,16], шаг 1/16, спаривание m=7/8,
каталог h2_zeros_scan.json.

Рецепт (см. h2_common.py, docstring): x_{n+1}=g(w)+s*h(w), w=(x^2-1)/r_disc,
g,h -- ряды Тейлора (степень 99) из w33_gh_hiprec.json (Чебышёв M=100) через
cheb_to_taylor (r_disc=2.5), AU0/delta -- из address_certificate.json.
Проверено flint 512 бит (test_h2_recipe.py в этой же папке) и согласовано
с numpy-сканом (diag_calibrate.py): H(u_actual) ~ 1e-4..1e-3 на m=7/8 --
рецепт РАБОЧИЙ (в отличие от ../b1_2_F11_d2, где использовался НЕВЕРНЫЙ
рецепт t=x^2 из spectra_grid.json и восстановился только 1/23).

Домен: max|x_n| вдоль орбиты < DOMAIN_MAXX=10.0 (диагностика: бесэскейпные
орбиты держат max|x|~1.0-1.01, эскейп -- скачок на много порядков, чистое
разделение при пороге 10).
"""
import json
import time

import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

import h2_common as C

T0 = time.time()
LOG = open('progress.log', 'a', encoding='utf-8')


def log(msg):
    s = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(s, flush=True)
    LOG.write(s + '\n')
    LOG.flush()


log('=== run_scan.py START ===')

# ---------------------------------------------------------------------- #
#  0. load coefficients + known zeros
# ---------------------------------------------------------------------- #
tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
tc_g_d_np = np.array([(k + 1) * tc_g_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_h_d_np = np.array([(k + 1) * tc_h_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_g_d_mp = [(k + 1) * tc_g_mp[k + 1] for k in range(N - 1)]
tc_h_d_mp = [(k + 1) * tc_h_mp[k + 1] for k in range(N - 1)]

known, meta_res = C.load_known_zeros()
log('known zeros loaded: n=%d (Z72 -> stage3b_z72_corrected.u_actual)' % len(known))

DOMAIN_MAXX = 10.0
CAND_THRESH = 2.5

# ---------------------------------------------------------------------- #
#  1. Task 1: quick sanity -- H(u_actual) at m=8 directly (no refinement)
# ---------------------------------------------------------------------- #
log('--- Task1 sanity: H_8(u_actual) for all 23 known (direct eval, no refine) ---')
u_known_arr = np.array([u for _, u in known], dtype=np.complex128)
Habs0, xden0, xnum0, maxw0, maxx0, _ = C.scan_level(
    u_known_arr.reshape(1, -1), 8, tc_g_np, tc_h_np)
for i, (cid, u) in enumerate(known):
    log('  %-5s u=(%.6f,%.6f)  |H_8|=%.3e  maxx(orbit)=%.4f' %
        (cid, u.real, u.imag, Habs0[0, i], maxx0[0, i]))

# ---------------------------------------------------------------------- #
#  2. full grid scan, step=1/16, m=8 and m=7
# ---------------------------------------------------------------------- #
RE_LO, RE_HI = -16.0, 16.0
IM_LO, IM_HI = 0.0, 16.0
STEP = 1.0 / 16.0
xs = np.arange(RE_LO, RE_HI + STEP / 2, STEP)
ys = np.arange(IM_LO, IM_HI + STEP / 2, STEP)
log('grid: Re[%.1f,%.1f] x Im[%.1f,%.1f] step=1/%d -> %d x %d = %d nodes'
    % (RE_LO, RE_HI, IM_LO, IM_HI, int(round(1 / STEP)), xs.size, ys.size, xs.size * ys.size))
XX, YY = np.meshgrid(xs, ys, indexing='ij')
grid = XX + 1j * YY

cands_by_level = {}
for m in (8, 7):
    t0 = time.time()
    Habs, xden, xnum, maxw, maxx, _ = C.scan_level(grid, m, tc_g_np, tc_h_np)
    valid = np.isfinite(Habs) & (np.abs(xden) > 1e-12) & (maxx < DOMAIN_MAXX)
    n_valid = int(valid.sum())
    log('scan m=%d done (%.1fs): valid %d/%d (%.1f%%), min|H| on valid = %.3e'
        % (m, time.time() - t0, n_valid, grid.size, 100.0 * n_valid / grid.size,
           float(Habs[valid].min()) if n_valid else float('nan')))
    cand = C.detect_local_minima(Habs, valid, xs, ys, CAND_THRESH)
    log('  candidates (local minima |H|<%.2f): %d' % (CAND_THRESH, len(cand)))
    cands_by_level[m] = cand

seeds = sorted(set((round(re, 6), round(im, 6))
                    for m in (7, 8) for re, im, _h in cands_by_level[m]))
seeds = [complex(re, im) for re, im in seeds]
log('union of m=7/m=8 candidate seeds: %d' % len(seeds))

# ---------------------------------------------------------------------- #
#  3. batch float64 Newton refine at m=8 (vectorized, with analytic grad)
# ---------------------------------------------------------------------- #
def scan_grad(u_arr, m):
    s = (C.AU0_F / (C.DELTA_F ** m)) * u_arr
    dsdu = C.AU0_F / (C.DELTA_F ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = np.zeros(u_arr.shape, dtype=np.complex128)
    dx = np.zeros(u_arr.shape, dtype=np.complex128)
    maxx = np.zeros(u_arr.shape, dtype=np.float64)
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
            dx = (dg + s * dh) * dx + h * dsdu
            x = g + s * h
            ax = np.abs(x)
            np.maximum(maxx, np.where(np.isfinite(ax), ax, maxx), out=maxx)
            if k + 1 == n1:
                xden, dxden = x.copy(), dx.copy()
    return x, xden, dx, dxden, maxx


def batch_newton(u0_arr, m, iters=30, cap=2.0):
    u = u0_arr.copy()
    for _it in range(iters):
        xnum, xden, dxnum, dxden, _maxx = scan_grad(u, m)
        H = xnum / xden
        dH = (dxnum * xden - xnum * dxden) / (xden * xden)
        with np.errstate(divide='ignore', invalid='ignore'):
            step = H / dH
        step = np.where(np.isfinite(step), step, 0.0)
        mag = np.abs(step)
        scale = np.where(mag > cap, cap / np.maximum(mag, 1e-300), 1.0)
        u = u - step * scale
    xnum, xden, dxnum, dxden, maxx = scan_grad(u, m)
    Hf = np.abs(xnum / xden)
    valid = np.isfinite(Hf) & (np.abs(xden) > 1e-12) & (maxx < DOMAIN_MAXX)
    return u, Hf, valid, maxx


log('--- batch float64 Newton refine at m=8, %d seeds, 30 iters ---' % len(seeds))
u0 = np.array(seeds, dtype=np.complex128)
t0 = time.time()
u8, H8, valid8, maxx8 = batch_newton(u0, 8, iters=30, cap=2.0)
log('  done in %.1fs; valid+converged (|H|<1e-6): %d / %d'
    % (time.time() - t0, int((valid8 & (H8 < 1e-6)).sum()), len(seeds)))

ok = valid8 & (H8 < 1e-6)
u8_ok = u8[ok]
H8_ok = H8[ok]

# cluster duplicates (eps=1e-5), keep smallest |H|
log('clustering m=8 refined zeros (eps=1e-5)...')
order = np.argsort(H8_ok)
clustered = []
for idx in order:
    uu, hh = u8_ok[idx], H8_ok[idx]
    placed = False
    for c in clustered:
        if abs(uu - c['u']) < 1e-5:
            placed = True
            break
    if not placed:
        clustered.append({'u': uu, 'absH': float(hh)})
log('  unique m=8 zeros after clustering: %d' % len(clustered))

# restrict to the requested scan box (Newton can wander outside)
in_box = [c for c in clustered
          if RE_LO - 0.5 <= c['u'].real <= RE_HI + 0.5 and IM_LO - 0.5 <= c['u'].imag <= IM_HI + 0.5]
log('  zeros inside scan box (+-0.5 margin): %d' % len(in_box))
clustered = in_box

# ---------------------------------------------------------------------- #
#  4. pairing m=7 <-> m=8 (start from same seeds, converge independently)
# ---------------------------------------------------------------------- #
log('--- pairing: Newton refine same positions at m=7, compare drift ---')
u0_7 = np.array([c['u'] for c in clustered], dtype=np.complex128)
t0 = time.time()
u7, H7, valid7, maxx7 = batch_newton(u0_7, 7, iters=30, cap=2.0)
log('  m=7 refine done in %.1fs' % (time.time() - t0))
for i, c in enumerate(clustered):
    c['u7'] = complex(u7[i]) if (valid7[i] and H7[i] < 1e-4) else None
    c['absH7'] = float(H7[i]) if valid7[i] else None
    c['drift78'] = abs(complex(u7[i]) - c['u']) if c['u7'] is not None else None
n_paired = sum(1 for c in clustered if c['u7'] is not None and c['drift78'] < 0.05)
log('  paired (m=7 partner within 0.05): %d / %d' % (n_paired, len(clustered)))
for c in clustered:
    c['paired'] = bool(c['u7'] is not None and c['drift78'] < 0.05)

# ---------------------------------------------------------------------- #
#  5. light mpmath polish (mp.dps=25, few iterations) for extra digits
# ---------------------------------------------------------------------- #
mp.dps = 25
RDISC_MP = C.R_DISC_MP
AU0_MP = C.AU0_MP
DELTA_MP = C.DELTA_MP


def orbit_grad_mp(u_mp, m):
    s_mp = u_mp * AU0_MP / (DELTA_MP ** m)
    dsdu = AU0_MP / (DELTA_MP ** m)
    n1 = 2 ** m
    n2 = 2 ** (m + 1)
    x = mpc(0)
    dx = mpc(0)
    xden = dxden = None
    for k in range(n2):
        xx = x * x - 1
        w = xx / RDISC_MP
        g = tc_g_mp[-1]
        for c in tc_g_mp[-2::-1]:
            g = g * w + c
        h = tc_h_mp[-1]
        for c in tc_h_mp[-2::-1]:
            h = h * w + c
        dg_w = tc_g_d_mp[-1]
        for c in tc_g_d_mp[-2::-1]:
            dg_w = dg_w * w + c
        dh_w = tc_h_d_mp[-1]
        for c in tc_h_d_mp[-2::-1]:
            dh_w = dh_w * w + c
        dw_dx = 2 * x / RDISC_MP
        dg = dg_w * dw_dx
        dh = dh_w * dw_dx
        dx = (dg + s_mp * dh) * dx + h * dsdu
        x = g + s_mp * h
        if k + 1 == n1:
            xden, dxden = x, dx
    return xden, x, dxden, dx


def mp_polish_light(u_start, m=8, iters=3):
    u_mp = mpc(mpf(str(u_start.real)), mpf(str(u_start.imag)))
    for _ in range(iters):
        xden, xnum, dxden, dxnum = orbit_grad_mp(u_mp, m)
        Hv = xnum / xden
        dH = (dxnum * xden - xnum * dxden) / (xden * xden)
        if dH == 0:
            break
        u_mp = u_mp - Hv / dH
    xden, xnum = orbit_grad_mp(u_mp, m)[0:2]
    absH = float(fabs(xnum / xden))
    return u_mp, absH


log('--- mpmath light polish (dps=25, 3 iters) on %d zeros ---' % len(clustered))
t0 = time.time()
for i, c in enumerate(clustered):
    u_mp, absH_mp = mp_polish_light(c['u'], m=8, iters=3)
    c['u_mp_re'] = nstr(mpre(u_mp), 18)
    c['u_mp_im'] = nstr(mpim(u_mp), 18)
    c['absH_mp'] = absH_mp
    if (i + 1) % 25 == 0 or (i + 1) == len(clustered):
        el = time.time() - t0
        eta = el / (i + 1) * (len(clustered) - i - 1)
        log('  polish %d/%d (%.0fs, ETA %.0fs)' % (i + 1, len(clustered), el, eta))
log('mpmath polish done in %.1fs' % (time.time() - t0))

# ---------------------------------------------------------------------- #
#  6. validation against 23 known zeros
# ---------------------------------------------------------------------- #
log('--- VALIDATION: recover 23 known zeros from the scan catalog (tol=1e-3) ---')
cat_u = np.array([complex(mpf(c['u_mp_re']), mpf(c['u_mp_im'])) for c in clustered])
recovered = 0
val_rows = []
for cid, u in known:
    d = np.abs(cat_u - u)
    j = int(np.argmin(d)) if len(d) else -1
    dist = float(d[j]) if j >= 0 else float('inf')
    ok_v = dist <= 1e-3
    recovered += int(ok_v)
    found_u = clustered[j]['u'] if j >= 0 else None
    val_rows.append(dict(id=cid, u_actual_re=u.real, u_actual_im=u.imag,
                          found_re=(found_u.real if found_u is not None else None),
                          found_im=(found_u.imag if found_u is not None else None),
                          dist=dist, ok=ok_v))
    log('  %-5s u_actual=(%.6f,%.6f) -> nearest found dist=%.3e %s'
        % (cid, u.real, u.imag, dist, 'OK' if ok_v else 'MISS'))
log('VALIDATION RESULT: %d / 23 recovered (tol=1e-3)' % recovered)

# ---------------------------------------------------------------------- #
#  7. build final catalog with conjugates, save
# ---------------------------------------------------------------------- #
final = []
for c in clustered:
    re_s, im_s = c['u_mp_re'], c['u_mp_im']
    final.append(dict(re=re_s, im=im_s, absH_mp=c['absH_mp'], absH_float=c['absH'],
                       paired_m7=c['paired'],
                       drift78=(c['drift78'] if c['drift78'] is not None else None)))
    if abs(float(im_s)) > 1e-9:
        final.append(dict(re=re_s, im=str(-float(im_s)), absH_mp=c['absH_mp'],
                           absH_float=c['absH'], paired_m7=c['paired'],
                           drift78=(c['drift78'] if c['drift78'] is not None else None),
                           conjugate=True))
log('final catalog (with conjugates for Im<0): %d entries' % len(final))

out = dict(
    meta=dict(
        recipe='x_{n+1}=g(w)+s*h(w), w=(x^2-1)/r_disc, address_certificate.json convention',
        AU0=C.AU0_STR, delta=C.DELTA_STR, r_disc=C.R_DISC_STR,
        w33_source=C.W33_PATH,
        N_taylor=N, domain_maxx=DOMAIN_MAXX, cand_thresh=CAND_THRESH,
        grid=dict(re_lo=RE_LO, re_hi=RE_HI, im_lo=IM_LO, im_hi=IM_HI, step=STEP,
                  nodes=int(xs.size * ys.size)),
        levels_scanned=[7, 8], mp_dps_polish=25,
        elapsed_sec=time.time() - T0,
    ),
    candidates=dict(m8=len(cands_by_level[8]), m7=len(cands_by_level[7]), union_seeds=len(seeds)),
    counts=dict(refined_valid=int((valid8 & (H8 < 1e-6)).sum()),
                clustered_unique=len(clustered), paired_m7=n_paired,
                final_with_conjugates=len(final)),
    validation=dict(recovered=recovered, total=23, tolerance=1e-3, rows=val_rows),
    zeros=final,
)
with open('h2_zeros_scan.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote h2_zeros_scan.json (%d zeros incl. conjugates)' % len(final))
log('=== run_scan.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

# -*- coding: utf-8 -*-
"""S2/H2SCAN v2, часть 1/2 -- смешанная сетка (грубая 1/16 по всему боксу +
точная 1/64 в плотной подобласти band III-V, где известные нули стоят
теснее ~0.12-0.26) + Ньютон (float64, векторизовано) + спаривание m=7.

ПРИЧИНА переделки: первый прогон (run_scan.py, шаг 1/16 везде) нашёл 175
уникальных нулей, но восстановил только 11/23 известных -- 12 пропущенных
(Z05,Z38,Z57,Z60,Z65,Z67,Z68,Z69,Z70,Z71,Z72,Z73) все лежат в плотном
кластере Re~[4,8.6] Im~[3.4,8.0] (band III-V адресации), где соседние
известные нули стоят теснее шага 1/16=0.0625 сетки эффективно (при том что
номинальный зазор 0.12-0.26 больше шага, дискретный ландшафт |H| у тесных
соседей сливается в одну котловину -- нужна более тонкая сетка).

Результат этого скрипта: scan_candidates_v2.json (float64 позиции + |H|,
БЕЗ mpmath-полировки -- она в part2, отдельным прогоном, т.к. mpmath дороже
таймбюджета в одном synchронном вызове).
"""
import json
import time

import numpy as np

import h2_common as C

T0 = time.time()
LOG = open('progress.log', 'a', encoding='utf-8')


def log(msg):
    s = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(s, flush=True)
    LOG.write(s + '\n')
    LOG.flush()


log('=== run_scan2_part1.py START (mixed-resolution grid v2) ===')

tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
tc_g_d_np = np.array([(k + 1) * tc_g_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_h_d_np = np.array([(k + 1) * tc_h_np[k + 1] for k in range(N - 1)], dtype=np.float64)

known, meta_res = C.load_known_zeros()
log('known zeros loaded: n=%d' % len(known))

DOMAIN_MAXX = 10.0
CAND_THRESH = 2.5

# ---------------------------------------------------------------------- #
#  grid 1: coarse, full box, step 1/16
# ---------------------------------------------------------------------- #
RE_LO, RE_HI = -16.0, 16.0
IM_LO, IM_HI = 0.0, 16.0
STEP_C = 1.0 / 16.0
xs_c = np.arange(RE_LO, RE_HI + STEP_C / 2, STEP_C)
ys_c = np.arange(IM_LO, IM_HI + STEP_C / 2, STEP_C)
XXc, YYc = np.meshgrid(xs_c, ys_c, indexing='ij')
grid_c = XXc + 1j * YYc
log('coarse grid: Re[%.1f,%.1f] x Im[%.1f,%.1f] step=1/16 -> %d nodes'
    % (RE_LO, RE_HI, IM_LO, IM_HI, grid_c.size))

t0 = time.time()
Habs_c, xden_c, xnum_c, maxw_c, maxx_c, _ = C.scan_level(grid_c, 8, tc_g_np, tc_h_np)
valid_c = np.isfinite(Habs_c) & (np.abs(xden_c) > 1e-12) & (maxx_c < DOMAIN_MAXX)
log('coarse scan m=8 done (%.1fs): valid %d/%d' % (time.time() - t0, int(valid_c.sum()), grid_c.size))
cand_c = C.detect_local_minima(Habs_c, valid_c, xs_c, ys_c, CAND_THRESH)
log('  coarse candidates: %d' % len(cand_c))

# ---------------------------------------------------------------------- #
#  grid 2: fine, dense sub-box (band III-V neighbourhood), step 1/64
# ---------------------------------------------------------------------- #
FRE_LO, FRE_HI = 1.0, 10.0
FIM_LO, FIM_HI = 3.0, 9.0
STEP_F = 1.0 / 64.0
xs_f = np.arange(FRE_LO, FRE_HI + STEP_F / 2, STEP_F)
ys_f = np.arange(FIM_LO, FIM_HI + STEP_F / 2, STEP_F)
XXf, YYf = np.meshgrid(xs_f, ys_f, indexing='ij')
grid_f = XXf + 1j * YYf
log('fine grid: Re[%.1f,%.1f] x Im[%.1f,%.1f] step=1/64 -> %d nodes'
    % (FRE_LO, FRE_HI, FIM_LO, FIM_HI, grid_f.size))

t0 = time.time()
Habs_f, xden_f, xnum_f, maxw_f, maxx_f, _ = C.scan_level(grid_f, 8, tc_g_np, tc_h_np)
valid_f = np.isfinite(Habs_f) & (np.abs(xden_f) > 1e-12) & (maxx_f < DOMAIN_MAXX)
log('fine scan m=8 done (%.1fs): valid %d/%d' % (time.time() - t0, int(valid_f.sum()), grid_f.size))
cand_f = C.detect_local_minima(Habs_f, valid_f, xs_f, ys_f, CAND_THRESH)
log('  fine candidates: %d' % len(cand_f))

seeds = sorted(set((round(re, 6), round(im, 6)) for re, im, _h in cand_c + cand_f))
seeds = [complex(re, im) for re, im in seeds]
log('union coarse+fine candidate seeds: %d' % len(seeds))


# ---------------------------------------------------------------------- #
#  batch float64 Newton refine (m=8) + m=7 pairing
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

in_box = [c for c in clustered
          if RE_LO - 0.5 <= c['u'].real <= RE_HI + 0.5 and IM_LO - 0.5 <= c['u'].imag <= IM_HI + 0.5]
log('  zeros inside scan box (+-0.5 margin): %d' % len(in_box))
clustered = in_box

log('--- pairing: Newton refine same positions at m=7 ---')
u0_7 = np.array([c['u'] for c in clustered], dtype=np.complex128)
t0 = time.time()
u7, H7, valid7, maxx7 = batch_newton(u0_7, 7, iters=30, cap=2.0)
log('  m=7 refine done in %.1fs' % (time.time() - t0))
for i, c in enumerate(clustered):
    c['u7_re'] = float(u7[i].real) if valid7[i] else None
    c['u7_im'] = float(u7[i].imag) if valid7[i] else None
    c['absH7'] = float(H7[i]) if valid7[i] else None
    c['drift78'] = abs(complex(u7[i]) - c['u']) if valid7[i] and H7[i] < 1e-4 else None
n_paired = sum(1 for c in clustered if c['drift78'] is not None and c['drift78'] < 0.05)
log('  paired (m=7 partner within 0.05): %d / %d' % (n_paired, len(clustered)))
for c in clustered:
    c['paired'] = bool(c['drift78'] is not None and c['drift78'] < 0.05)

# quick float64-only validation preview (mpmath polish happens in part2)
log('--- preview validation (float64 only, tol=1e-3) ---')
cat_u = np.array([c['u'] for c in clustered])
recovered = 0
for cid, u in known:
    d = np.abs(cat_u - u)
    dist = float(d.min()) if len(d) else float('inf')
    ok_v = dist <= 1e-3
    recovered += int(ok_v)
    log('  %-5s dist=%.3e %s' % (cid, dist, 'OK' if ok_v else 'MISS'))
log('PREVIEW VALIDATION (float64): %d / 23' % recovered)

# ---------------------------------------------------------------------- #
#  save intermediate (float64) for part2 mpmath polish
# ---------------------------------------------------------------------- #
out = dict(
    meta=dict(
        grids=[dict(kind='coarse', re=[RE_LO, RE_HI], im=[IM_LO, IM_HI], step=STEP_C),
               dict(kind='fine', re=[FRE_LO, FRE_HI], im=[FIM_LO, FIM_HI], step=STEP_F)],
        domain_maxx=DOMAIN_MAXX, cand_thresh=CAND_THRESH,
        candidates=dict(coarse=len(cand_c), fine=len(cand_f), union_seeds=len(seeds)),
        elapsed_sec=time.time() - T0,
    ),
    zeros=[dict(re=float(c['u'].real), im=float(c['u'].imag), absH=c['absH'],
                paired_m7=c['paired'], drift78=c['drift78']) for c in clustered],
)
with open('scan_candidates_v2.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote scan_candidates_v2.json (%d zeros, float64)' % len(clustered))
log('=== run_scan2_part1.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

# -*- coding: utf-8 -*-
"""Диагностика: локальный зум вокруг Z38, Z72, Z73 (не найдены каталогом),
чтобы понять, есть ли там резкая узкая яма |H_8|, которую пропускает даже
сетка 1/64, или локальный минимум существует но выше порога (=> проблема
в самой посадке Ньютона, а не в сетке)."""
import time
import numpy as np
import h2_common as C

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
known, meta_res = C.load_known_zeros()
kd = dict(known)

targets = {'Z38': kd['Z38'], 'Z72': kd['Z72'], 'Z73': kd['Z73']}

for name, u0 in targets.items():
    log('--- zoom around %s, u_actual=(%.6f,%.6f) ---' % (name, u0.real, u0.imag))
    half = 0.6
    step = 1.0 / 256.0
    xs = np.arange(u0.real - half, u0.real + half + step / 2, step)
    ys = np.arange(u0.imag - half, u0.imag + half + step / 2, step)
    XX, YY = np.meshgrid(xs, ys, indexing='ij')
    grid = XX + 1j * YY
    Habs, xden, xnum, maxw, maxx, _ = C.scan_level(grid, 8, tc_g_np, tc_h_np)
    valid = np.isfinite(Habs) & (np.abs(xden) > 1e-12) & (maxx < 10.0)
    hv = np.where(valid, Habs, np.inf)
    imin = np.unravel_index(np.argmin(hv), hv.shape)
    log('  grid %dx%d step=1/256, min|H8|=%.4e at u=(%.6f,%.6f)  dist_to_u_actual=%.4e'
        % (grid.shape[0], grid.shape[1], hv[imin], xs[imin[0]], ys[imin[1]],
           abs(complex(xs[imin[0]], ys[imin[1]]) - u0)))
    # count how many grid points have |H|<0.5, <2.5
    for thr in (0.1, 0.5, 2.5):
        log('    fraction with |H8|<%.1f: %.4f (%d points)' % (thr, float((hv < thr).mean()), int((hv < thr).sum())))
    # also report H8 exactly AT u_actual (nearest grid point) and slightly around
    Habs_at, xden_at, xnum_at, maxw_at, maxx_at, _ = C.scan_level(
        np.array([[u0]], dtype=np.complex128), 8, tc_g_np, tc_h_np)
    log('   H8(u_actual) = %.4e (maxx=%.4f)' % (Habs_at[0, 0], maxx_at[0, 0]))
log('DONE diag_zoom')

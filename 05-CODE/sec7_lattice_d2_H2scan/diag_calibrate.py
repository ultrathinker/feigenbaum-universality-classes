# -*- coding: utf-8 -*-
"""Диагностика: проверить, что numpy-скан согласуется с flint-проверкой
(test_h2_recipe.py), и оценить диапазон max|w|, max|x| по региону скана,
чтобы откалибровать порог 'орбита не улетает' перед основным скриптом.
Быстрый прогон (не более пары минут)."""
import time
import numpy as np
import h2_common as C

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


log('loading tc_g, tc_h ...')
tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)

log('known zeros ...')
known, meta = C.load_known_zeros()
log('n known = %d' % len(known))

# numpy scan sanity check against flint reference for first 5 known zeros
log('--- numpy H_m(u_actual) for first 5 known, m=7,8 (compare with flint test) ---')
for cid, u in known[:5]:
    grid = np.array([[u]], dtype=np.complex128)
    for m in (7, 8):
        Habs, xden, xnum, maxw, maxx, errmax = C.scan_level(grid, m, tc_g_np, tc_h_np)
        log('%s m=%d |H|=%.6e maxw=%.4f maxx=%.4f' % (cid, m, Habs[0, 0], maxw[0, 0], maxx[0, 0]))

# scan a coarse grid over the FULL target region to see maxw/maxx distribution and
# where things blow up (if anywhere), at m=8 (cheaper: use step=1.0 coarse probe)
log('--- coarse probe over Re[-16,16] x Im[0,16] step=1.0, m=8 ---')
xs = np.arange(-16.0, 16.0 + 0.5, 1.0)
ys = np.arange(0.0, 16.0 + 0.5, 1.0)
XX, YY = np.meshgrid(xs, ys, indexing='ij')
grid = XX + 1j * YY
t0 = time.time()
Habs, xden, xnum, maxw, maxx, errmax = C.scan_level(grid, 8, tc_g_np, tc_h_np)
log('coarse probe m=8 done in %.1fs, grid=%dx%d' % (time.time() - t0, grid.shape[0], grid.shape[1]))
finite_mask = np.isfinite(Habs)
log('finite H fraction: %.3f (%d/%d)' % (finite_mask.mean(), int(finite_mask.sum()), finite_mask.size))
log('maxw stats: min=%.4f median=%.4f p90=%.4f p99=%.4f max=%.4f'
    % (np.nanmin(maxw), np.nanmedian(maxw), np.nanpercentile(maxw, 90),
       np.nanpercentile(maxw, 99), np.nanmax(maxw)))
log('maxx stats: min=%.4f median=%.4f p90=%.4f p99=%.4f max=%.4e'
    % (np.nanmin(maxx), np.nanmedian(maxx), np.nanpercentile(maxx, 90),
       np.nanpercentile(maxx, 99), np.nanmax(maxx)))
finite_maxx = maxx[np.isfinite(maxx)]
log('finite maxx count=%d, count with maxx<10: %d, <100: %d, <1e6: %d'
    % (finite_maxx.size, int((finite_maxx < 10).sum()), int((finite_maxx < 100).sum()),
       int((finite_maxx < 1e6).sum())))
log('non-finite H count: %d' % int((~finite_mask).sum()))

# timing estimate for full fine scan at step=1/16 (513 x 257 nodes), m=8
log('--- timing estimate: fine grid step=1/16 subgrid 40x40 sample, m=8 ---')
xs2 = np.linspace(-16.0, 16.0, 40)
ys2 = np.linspace(0.0, 16.0, 40)
XX2, YY2 = np.meshgrid(xs2, ys2, indexing='ij')
grid2 = XX2 + 1j * YY2
t0 = time.time()
Habs2, *_ = C.scan_level(grid2, 8, tc_g_np, tc_h_np)
dt = time.time() - t0
n_full = 513 * 257
n_sample = grid2.size
log('sample 40x40=%d nodes took %.3fs -> extrapolated full 513x257=%d nodes ~ %.1fs (m=8)'
    % (n_sample, dt, n_full, dt * n_full / n_sample))

log('DONE diag_calibrate')

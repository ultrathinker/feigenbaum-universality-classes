# -*- coding: utf-8 -*-
"""Диагностика: один прогон скана m=8 (шаг 1/16) + подсчёт числа локальных
минимумов |H| при разных порогах, чтобы выбрать CAND_THRESH разумно (не
пересчитывая скан по многу раз)."""
import time
import numpy as np
import h2_common as C

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)

RE_LO, RE_HI = -16.0, 16.0
IM_LO, IM_HI = 0.0, 16.0
STEP = 1.0 / 16.0
xs = np.arange(RE_LO, RE_HI + STEP / 2, STEP)
ys = np.arange(IM_LO, IM_HI + STEP / 2, STEP)
XX, YY = np.meshgrid(xs, ys, indexing='ij')
grid = XX + 1j * YY
log('grid %d x %d = %d' % (xs.size, ys.size, xs.size * ys.size))

t0 = time.time()
Habs, xden, xnum, maxw, maxx, _ = C.scan_level(grid, 8, tc_g_np, tc_h_np)
log('scan m=8 done in %.1fs' % (time.time() - t0))
valid = np.isfinite(Habs) & (np.abs(xden) > 1e-12) & (maxx < 10.0)
log('valid fraction: %.3f' % valid.mean())

for thr in (0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
    cand = C.detect_local_minima(Habs, valid, xs, ys, thr)
    log('  thresh=%.2f -> %d local minima candidates' % (thr, len(cand)))

# distribution of Habs on valid points
hv = Habs[valid]
for p in (0.001, 0.01, 0.1, 1, 5, 10, 25, 50):
    log('  percentile %5.3f%% of |H| (valid) = %.4e' % (p, np.percentile(hv, p)))
log('DONE diag_threshold')

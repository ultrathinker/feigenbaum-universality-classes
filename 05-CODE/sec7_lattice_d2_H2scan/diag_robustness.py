# -*- coding: utf-8 -*-
"""Диагностика робастности: не является ли сигнал q=16 (p~1e-22 на полном
каталоге) артефактом того, что сама сетка скана имеет шаг 1/16 (и 1/64,
кратное 1/16)? Разбиваем каталог по географии (внутри плотного under-box
Re[1,10]xIm[3,9] -- где работала ТОНКАЯ сетка 1/64 -- против снаружи, где
только ГРУБАЯ сетка 1/16) и повторяем тест независимо на каждой половине.
Если сигнал держится в ОБЕИХ половинах (включая ту, что видела только
1/64-сетку, отличную от 1/16), это довод против чистого grid-артефакта."""
import json
import math
import time

import numpy as np
from scipy import stats as sstats

T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


scan = json.load(open('h2_zeros_scan.json', encoding='utf-8'))
pts = [complex(float(z['re']), float(z['im'])) for z in scan['zeros']]
log('total catalog (with conj): %d' % len(pts))

FRE_LO, FRE_HI = 1.0, 10.0
FIM_LO, FIM_HI = -9.0, 9.0  # conjugates included, mirror the fine box too

inside = [p for p in pts if FRE_LO <= p.real <= FRE_HI and FIM_LO <= p.imag <= FIM_HI]
outside = [p for p in pts if not (FRE_LO <= p.real <= FRE_HI and FIM_LO <= p.imag <= FIM_HI)]
log('inside fine-box region (could be found by 1/64 grid): %d' % len(inside))
log('outside fine-box region (ONLY the 1/16 coarse grid could find these): %d' % len(outside))


def dist_1d(x, q):
    f = q * x
    g = f - math.floor(f)
    return min(g, 1.0 - g) / q


def analyze(pts, q, label):
    n = len(pts)
    if n < 5:
        log('  %s: n=%d too small, skip' % (label, n))
        return
    eps = 0.125 / (2.0 * q)
    re = np.array([p.real for p in pts])
    im = np.array([p.imag for p in pts])
    re_d = np.array([dist_1d(x, q) for x in re])
    im_d = np.array([dist_1d(x, q) for x in im])
    re_near = re_d < eps
    im_near = im_d < eps
    both = re_near & im_near
    k_re, k_im, k_both = int(re_near.sum()), int(im_near.sum()), int(both.sum())
    p_re = float(sstats.binom.sf(k_re - 1, n, 0.125))
    p_im = float(sstats.binom.sf(k_im - 1, n, 0.125))
    p_both = float(sstats.binom.sf(k_both - 1, n, 0.015625))
    # Rayleigh at exactly q
    s = np.sum(np.exp(2j * np.pi * q * re))
    Rre = abs(s) / n
    s2 = np.sum(np.exp(2j * np.pi * q * im))
    Rim = abs(s2) / n
    log('  %-28s q=%2d n=%3d Re_near=%d/%.1f%%(p=%.2e) Im_near=%d/%.1f%%(p=%.2e) both=%d/%.1f%%(p=%.2e) R_re=%.3f R_im=%.3f'
        % (label, q, n, k_re, 100 * k_re / n, p_re, k_im, 100 * k_im / n, p_im,
           k_both, 100 * k_both / n, p_both, Rre, Rim))


for q in (16, 32):
    analyze(inside, q, 'inside fine-box (1/64-eligible)')
    analyze(outside, q, 'outside fine-box (1/16-only)')

log('DONE diag_robustness')

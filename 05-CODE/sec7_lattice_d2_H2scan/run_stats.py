# -*- coding: utf-8 -*-
"""S2/H2SCAN -- Задание 3+4: статистика решётки 1/16 (и 1/32) на:
  (A) 23 известных нуля (u_actual),
  (B) расширенном каталоге скана (h2_zeros_scan.json, с сопряжениями),
  (C) объединении (A union B, с дедупликацией совпадающих точек).

Протокол (как для d=4/6, НО с исправленным симметричным определением
"у узла", в отличие от ../b1_2_F11_d2/partA_lattice23.py, где тест
re_frac16<1/16 был асимметричным и давал expected=6.25% вместо заявленных
в отчёте 12.5% -- здесь используется правильная симметричная функция
dist_to_lattice_1d с порогом eps(q)=0.125/(2q), что ДАЁТ РОВНО 12.5%
ожидания по каждой координате и 1.5625% совместно, как того требует
бриф S2 п.3):
  - доля точек на расстоянии < eps(q) от узла k/q по Re, Im, оба (joint)
  - 2D евклидово расстояние до ближайшего узла (k/q, l/q): медиана/мин/макс
  - хи-квадрат по q бинам на frac(q*x) (Re, Im отдельно)
  - KS тест на равномерность frac(q*x)
  - точные биномиальные p (одностороннее P[X>=k_obs])
  - Рэлей R_re(q)=|mean(exp(2pi i q Re(u)))|, R_im(q) аналогично, на сетке
    q in [0.5,40] шаг 0.05 + тонкая сетка вокруг 16 и 32 (шаг 0.01)
"""
import json
import math
import time

import numpy as np
from scipy import stats as sstats

import h2_common as C

T0 = time.time()
LOG = open('progress.log', 'a', encoding='utf-8')


def log(msg):
    s = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(s, flush=True)
    LOG.write(s + '\n')
    LOG.flush()


log('=== run_stats.py START ===')

known, meta_res = C.load_known_zeros()
scan = json.load(open('h2_zeros_scan.json', encoding='utf-8'))
scan_pts = [complex(float(z['re']), float(z['im'])) for z in scan['zeros']]
log('known=%d, scan(with conj)=%d' % (len(known), len(scan_pts)))

known_pts = [u for _, u in known]

# dedup combined: known + scan, drop scan points within 1e-4 of a known point
combined = list(known_pts)
n_dupe = 0
for z in scan_pts:
    if any(abs(z - k) < 1e-4 for k in known_pts):
        n_dupe += 1
        continue
    combined.append(z)
log('combined (known + scan, dedup dist<1e-4): %d  (dupes dropped: %d)' % (len(combined), n_dupe))

DATASETS = {'known23': known_pts, 'scan_full': scan_pts, 'combined': combined}


# ----------------------------------------------------------------------- #
#  lattice protocol
# ----------------------------------------------------------------------- #
def dist_1d(x, q):
    f = q * x
    g = f - math.floor(f)
    return min(g, 1.0 - g) / q


def dist_2d(u, q):
    return math.hypot(dist_1d(u.real, q), dist_1d(u.imag, q))


def frac(x, q):
    f = q * x
    return f - math.floor(f)


def analyze(pts, q, label):
    n = len(pts)
    eps = 0.125 / (2.0 * q)   # symmetric window s.t. expected fraction = 0.125 exactly
    re = np.array([p.real for p in pts])
    im = np.array([p.imag for p in pts])
    re_d = np.array([dist_1d(x, q) for x in re])
    im_d = np.array([dist_1d(x, q) for x in im])
    d2 = np.array([dist_2d(p, q) for p in pts])
    re_near = re_d < eps
    im_near = im_d < eps
    both_near = re_near & im_near
    k_re, k_im, k_both = int(re_near.sum()), int(im_near.sum()), int(both_near.sum())
    p_re = 0.125
    p_both = 0.125 * 0.125
    out = dict(
        q=q, n=n, eps=eps, label=label,
        re_near=k_re, re_near_frac=k_re / n, im_near=k_im, im_near_frac=k_im / n,
        both_near=k_both, both_near_frac=k_both / n,
        p_binom_re=float(sstats.binom.sf(k_re - 1, n, p_re)) if n else None,
        p_binom_im=float(sstats.binom.sf(k_im - 1, n, p_re)) if n else None,
        p_binom_both=float(sstats.binom.sf(k_both - 1, n, p_both)) if n else None,
        d2_median=float(np.median(d2)), d2_min=float(np.min(d2)), d2_max=float(np.max(d2)),
    )
    # chi2 on q bins of frac(q*x)
    re_f = np.array([frac(x, q) for x in re])
    im_f = np.array([frac(x, q) for x in im])
    nb = min(q, max(4, n // 2)) if n >= 8 else 0
    if nb >= 4:
        re_bins = np.clip((re_f * nb).astype(int), 0, nb - 1)
        im_bins = np.clip((im_f * nb).astype(int), 0, nb - 1)
        e = n / nb
        chi2_re = float(((np.bincount(re_bins, minlength=nb) - e) ** 2 / e).sum())
        chi2_im = float(((np.bincount(im_bins, minlength=nb) - e) ** 2 / e).sum())
        out['chi2_nbins'] = nb
        out['chi2_re'] = chi2_re
        out['chi2_p_re'] = float(sstats.chi2.sf(chi2_re, nb - 1))
        out['chi2_im'] = chi2_im
        out['chi2_p_im'] = float(sstats.chi2.sf(chi2_im, nb - 1))
        ks_re = sstats.kstest(re_f, 'uniform')
        ks_im = sstats.kstest(im_f, 'uniform')
        out['ks_re_stat'] = float(ks_re.statistic)
        out['ks_re_p'] = float(ks_re.pvalue)
        out['ks_im_stat'] = float(ks_im.statistic)
        out['ks_im_p'] = float(ks_im.pvalue)
    else:
        out['chi2_nbins'] = None
    return out


def rayleigh_1d(vals, q):
    n = len(vals)
    s = np.sum(np.exp(2j * np.pi * q * np.asarray(vals)))
    R = abs(s) / n
    chi2 = 2 * n * R * R
    p = float(sstats.chi2.sf(chi2, 2))
    return float(R), p


def rayleigh_scan(pts, qs):
    re = np.array([p.real for p in pts])
    im = np.array([p.imag for p in pts])
    rows = []
    for q in qs:
        Rre, pre = rayleigh_1d(re, q)
        Rim, pim = rayleigh_1d(im, q)
        rows.append((float(q), Rre, pre, Rim, pim))
    return rows


qs_coarse = np.arange(0.5, 40.0001, 0.05)
qs_fine16 = np.arange(15.5, 16.5001, 0.01)
qs_fine32 = np.arange(31.5, 32.5001, 0.01)
qs_all = np.unique(np.concatenate([qs_coarse, qs_fine16, qs_fine32]))
log('Rayleigh grid: %d q-values (coarse 0.5..40 step .05 + fine near 16,32 step .01)' % len(qs_all))

results = {}
for name, pts in DATASETS.items():
    results[name] = {}
    for q in (16, 32):
        stat = analyze(pts, q, '%s@q=%d' % (name, q))
        results[name][q] = stat
        log('%-10s q=%2d n=%3d  Re_near=%d/%.1f%% Im_near=%d/%.1f%% both=%d/%.1f%%  '
            'p_re=%.3g p_im=%.3g p_both=%.3g  2Dmed=%.4f'
            % (name, q, stat['n'], stat['re_near'], 100 * stat['re_near_frac'],
               stat['im_near'], 100 * stat['im_near_frac'],
               stat['both_near'], 100 * stat['both_near_frac'],
               stat['p_binom_re'], stat['p_binom_im'], stat['p_binom_both'],
               stat['d2_median']))
        if stat.get('chi2_nbins'):
            log('           chi2_re=%.2f(p=%.3f) chi2_im=%.2f(p=%.3f) KS_re p=%.3f KS_im p=%.3f'
                % (stat['chi2_re'], stat['chi2_p_re'], stat['chi2_im'], stat['chi2_p_im'],
                   stat['ks_re_p'], stat['ks_im_p']))

log('--- Rayleigh spectra (this can take a few seconds) ---')
rayleigh_results = {}
for name, pts in DATASETS.items():
    rows = rayleigh_scan(pts, qs_all)
    rayleigh_results[name] = rows
    # top peaks
    top_re = sorted(rows, key=lambda r: -r[1])[:8]
    top_im = sorted(rows, key=lambda r: -r[3])[:8]
    log('%s: top Re-peaks: %s'
        % (name, ', '.join('q=%.2f R=%.3f p=%.1e' % (q, R, p) for q, R, p, _, _ in top_re)))
    log('%s: top Im-peaks: %s'
        % (name, ', '.join('q=%.2f R=%.3f p=%.1e' % (q, Rim, pim) for q, _, _, Rim, pim in top_im)))
    # value exactly at q=16 and q=32
    for qtag in (16.0, 32.0):
        row = min(rows, key=lambda r: abs(r[0] - qtag))
        log('   at q=%.2f: R_re=%.4f (p=%.2e)  R_im=%.4f (p=%.2e)'
            % (row[0], row[1], row[2], row[3], row[4]))

# ----------------------------------------------------------------------- #
#  save
# ----------------------------------------------------------------------- #
out = dict(
    meta=dict(n_known=len(known_pts), n_scan_full=len(scan_pts), n_combined=len(combined),
              n_scan_dupe_of_known=n_dupe,
              note='near-node test uses SYMMETRIC dist_1d with eps(q)=0.125/(2q) '
                   '(expected 12.5% per-coord, 1.5625% joint) -- corrected vs '
                   '../b1_2_F11_d2 asymmetric frac<1/16 bug',
              elapsed_sec=time.time() - T0),
    lattice=results,
    rayleigh_top10=dict(
        (name, dict(
            top_re=[list(r) for r in sorted(rayleigh_results[name], key=lambda r: -r[1])[:10]],
            top_im=[list(r) for r in sorted(rayleigh_results[name], key=lambda r: -r[3])[:10]],
        )) for name in DATASETS
    ),
    rayleigh_at_16_32=dict(
        (name, dict(
            q16=[list(min(rayleigh_results[name], key=lambda r: abs(r[0] - 16.0)))],
            q32=[list(min(rayleigh_results[name], key=lambda r: abs(r[0] - 32.0)))],
        )) for name in DATASETS
    ),
)
with open('stats_lattice.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote stats_lattice.json')
log('=== run_stats.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

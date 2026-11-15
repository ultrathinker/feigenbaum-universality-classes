#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 Task 2: Rayleigh statistic R(q) = |sum exp(2 pi i q x_j)| / n on
   Re and Im parts of catalog zeros (d=6: 342, d=4: 62).

   Fine grid q in [15.5, 16.5] step 0.001
   Coarse grid q in [0.5, 40] step 0.05 (just for picture)

   Bootstrap n=2000 over points for CI on q*.
   Compare q* with candidates alpha_d^{-k}, delta_d^{-k}, and products (k=1..12).
   Save rayleigh_results.json with full data + plot png.
"""
import os
import sys
import json
import time
import math
import numpy as np

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'
T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def rayleigh_R(x, q):
    """R(q) = |mean(exp(2 pi i q x))| for array x.  q is scalar or array."""
    phase = 2 * np.pi * np.asarray(q)[None, :] * np.asarray(x)[:, None]
    return np.abs(np.exp(1j * phase).mean(axis=0))


def bootstrap_qstar(x, q_grid, n_boot=2000, seed=20260815):
    """Return (qstar, lo, hi) using percentile bootstrap over R(q) max.
       Note: bootstrap sample zeros with replacement, recompute R, find argmax.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    qstars = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        xs = x[idx]
        R = rayleigh_R(xs, q_grid)
        qstars[b] = q_grid[int(np.argmax(R))]
    return float(np.median(qstars)), float(np.percentile(qstars, 2.5)), float(np.percentile(qstars, 97.5))


def best_candidate(qstar, candidates, tol=1e-3):
    out = []
    for name, val in candidates:
        diff = qstar - val
        rel = diff / val if val != 0 else diff
        out.append((name, val, diff, rel))
    out.sort(key=lambda t: abs(t[2]))
    return out


def main():
    log('--- Rayleigh: load inputs ---')
    hz6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'H6_zeros.json'),
                         encoding='utf-8'))
    hz4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'H4_zeros.json'),
                         encoding='utf-8'))
    fp6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'sextic_fixedpoint.json'),
                         encoding='utf-8'))
    fp4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'quartic_fixedpoint.json'),
                         encoding='utf-8'))

    z6 = hz6['zeros_level8']
    z4 = hz4['zeros_level8']
    re6 = np.array([z['re'] for z in z6], dtype=np.float64)
    im6 = np.array([z['im'] for z in z6], dtype=np.float64)
    re4 = np.array([z['re'] for z in z4], dtype=np.float64)
    im4 = np.array([z['im'] for z in z4], dtype=np.float64)
    log('n6=%d n4=%d' % (len(re6), len(re4)))

    # Coarse q
    q_coarse = np.arange(0.5, 40.0 + 1e-9, 0.05)
    R6_re_c = rayleigh_R(re6, q_coarse)
    R6_im_c = rayleigh_R(im6, q_coarse)
    R4_re_c = rayleigh_R(re4, q_coarse)
    R4_im_c = rayleigh_R(im4, q_coarse)
    log('coarse peak d=6 re: q=%.3f R=%.4f' % (q_coarse[int(np.argmax(R6_re_c))], R6_re_c.max()))
    log('coarse peak d=6 im: q=%.3f R=%.4f' % (q_coarse[int(np.argmax(R6_im_c))], R6_im_c.max()))
    log('coarse peak d=4 re: q=%.3f R=%.4f' % (q_coarse[int(np.argmax(R4_re_c))], R4_re_c.max()))
    log('coarse peak d=4 im: q=%.3f R=%.4f' % (q_coarse[int(np.argmax(R4_im_c))], R4_im_c.max()))

    # Fine q
    q_fine = np.arange(15.5, 16.5 + 1e-9, 0.001)
    R6_re_f = rayleigh_R(re6, q_fine)
    R6_im_f = rayleigh_R(im6, q_fine)
    R4_re_f = rayleigh_R(re4, q_fine)
    R4_im_f = rayleigh_R(im4, q_fine)

    def stats(q, R):
        idx = int(np.argmax(R))
        return float(q[idx]), float(R[idx]), idx

    q6re, R6re, _ = stats(q_fine, R6_re_f)
    q6im, R6im, _ = stats(q_fine, R6_im_f)
    q4re, R4re, _ = stats(q_fine, R4_re_f)
    q4im, R4im, _ = stats(q_fine, R4_im_f)
    log('fine peak d=6 re: q*=%.6f R=%.4f' % (q6re, R6re))
    log('fine peak d=6 im: q*=%.6f R=%.4f' % (q6im, R6im))
    log('fine peak d=4 re: q*=%.6f R=%.4f' % (q4re, R4re))
    log('fine peak d=4 im: q*=%.6f R=%.4f' % (q4im, R4im))

    # Bootstrap CI (use combined re+im? -- brief says separately; do separately)
    log('--- bootstrap (2000 reps) ---')
    boot6re = bootstrap_qstar(re6, q_fine, n_boot=2000)
    boot6im = bootstrap_qstar(im6, q_fine, n_boot=2000)
    boot4re = bootstrap_qstar(re4, q_fine, n_boot=2000)
    boot4im = bootstrap_qstar(im4, q_fine, n_boot=2000)
    log('  d=6 re bootstrap median=%.4f 95%% CI=[%.4f, %.4f]' % boot6re)
    log('  d=6 im bootstrap median=%.4f 95%% CI=[%.4f, %.4f]' % boot6im)
    log('  d=4 re bootstrap median=%.4f 95%% CI=[%.4f, %.4f]' % boot4re)
    log('  d=4 im bootstrap median=%.4f 95%% CI=[%.4f, %.4f]' % boot4im)

    # Compare q* with arithmetic / alpha / delta candidates
    alpha6 = float(fp6['alpha6']); delta6 = float(fp6['delta6'])
    alpha4 = float(fp4['alpha4']); delta4 = float(fp4['delta4'])
    cands6 = []
    cands4 = []
    for k in range(1, 13):
        cands6.append(('alpha6^-%d' % k, alpha6 ** -k))
        cands6.append(('delta6^-%d' % k, delta6 ** -k))
        cands4.append(('alpha4^-%d' % k, alpha4 ** -k))
        cands4.append(('delta4^-%d' % k, delta4 ** -k))
    # products with 16 (modulation of 16 by small factor)
    for k in range(1, 7):
        cands6.append(('16*alpha6^%d' % k, 16 * alpha6 ** k))
        cands6.append(('16*delta6^%d' % k, 16 * delta6 ** k))
        cands6.append(('16/alpha6^%d' % k, 16 / alpha6 ** k))
        cands6.append(('16/delta6^%d' % k, 16 / delta6 ** k))
        cands4.append(('16*alpha4^%d' % k, 16 * alpha4 ** k))
        cands4.append(('16*delta4^%d' % k, 16 * delta4 ** k))
        cands4.append(('16/alpha4^%d' % k, 16 / alpha4 ** k))
        cands4.append(('16/delta4^%d' % k, 16 / delta4 ** k))
    cands6.append(('16.000', 16.0))
    cands4.append(('16.000', 16.0))

    log('--- best candidates near q* (d=6 re) ---')
    for name, val, diff, rel in best_candidate(q6re, cands6)[:8]:
        log('  q*=%.4f  %s=%.6f  diff=%+.4e  rel=%+.4e' % (q6re, name, val, diff, rel))
    log('--- best candidates near q* (d=6 im) ---')
    for name, val, diff, rel in best_candidate(q6im, cands6)[:8]:
        log('  q*=%.4f  %s=%.6f  diff=%+.4e  rel=%+.4e' % (q6im, name, val, diff, rel))
    log('--- best candidates near q* (d=4 re) ---')
    for name, val, diff, rel in best_candidate(q4re, cands4)[:8]:
        log('  q*=%.4f  %s=%.6f  diff=%+.4e  rel=%+.4e' % (q4re, name, val, diff, rel))
    log('--- best candidates near q* (d=4 im) ---')
    for name, val, diff, rel in best_candidate(q4im, cands4)[:8]:
        log('  q*=%.4f  %s=%.6f  diff=%+.4e  rel=%+.4e' % (q4im, name, val, diff, rel))

    # Combined statistic: union re, im sample (independent axes), same q search
    x6 = np.concatenate([re6, im6])
    x4 = np.concatenate([re4, im4])
    R6_u = rayleigh_R(x6, q_fine)
    R4_u = rayleigh_R(x4, q_fine)
    q6u, R6u, _ = stats(q_fine, R6_u)
    q4u, R4u, _ = stats(q_fine, R4_u)
    boot6u = bootstrap_qstar(x6, q_fine, n_boot=2000)
    boot4u = bootstrap_qstar(x4, q_fine, n_boot=2000)
    log('combined d=6 (re U im): q*=%.4f R=%.4f bootstrap median=%.4f CI=[%.4f, %.4f]'
        % (q6u, R6u, *boot6u))
    log('combined d=4 (re U im): q*=%.4f R=%.4f bootstrap median=%.4f CI=[%.4f, %.4f]'
        % (q4u, R4u, *boot4u))

    out = dict(
        meta=dict(alpha6=alpha6, delta6=delta6, alpha4=alpha4, delta4=delta4,
                  n6=len(re6), n4=len(re4)),
        coarse=dict(q=q_coarse.tolist(),
                    R6_re=R6_re_c.tolist(), R6_im=R6_im_c.tolist(),
                    R4_re=R4_re_c.tolist(), R4_im=R4_im_c.tolist()),
        fine=dict(q=q_fine.tolist(),
                  R6_re=R6_re_f.tolist(), R6_im=R6_im_f.tolist(),
                  R4_re=R4_re_f.tolist(), R4_im=R4_im_f.tolist()),
        peaks=dict(d6_re=dict(q=q6re, R=R6re, boot_median=boot6re[0],
                              boot_lo=boot6re[1], boot_hi=boot6re[2]),
                   d6_im=dict(q=q6im, R=R6im, boot_median=boot6im[0],
                              boot_lo=boot6im[1], boot_hi=boot6im[2]),
                   d4_re=dict(q=q4re, R=R4re, boot_median=boot4re[0],
                              boot_lo=boot4re[1], boot_hi=boot4re[2]),
                   d4_im=dict(q=q4im, R=R4im, boot_median=boot4im[0],
                              boot_lo=boot4im[1], boot_hi=boot4im[2]),
                   combined_d6=dict(q=q6u, R=R6u, boot_median=boot6u[0],
                                    boot_lo=boot6u[1], boot_hi=boot6u[2]),
                   combined_d4=dict(q=q4u, R=R4u, boot_median=boot4u[0],
                                    boot_lo=boot4u[1], boot_hi=boot4u[2]),
                   ),
    )
    with open(os.path.join(OUT_DIR, 'rayleigh_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote rayleigh_results.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('MKL_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
    main()
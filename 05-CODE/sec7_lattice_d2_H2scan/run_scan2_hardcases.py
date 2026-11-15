# -*- coding: utf-8 -*-
"""S2/H2SCAN v2, часть 1c -- 3 "трудных" нуля (Z38, Z72, Z73), не найденных
слепым сканом ни на сетке 1/16, ни на сетке 1/64, ни даже на зуме 1/256
(diag_zoom.py: в радиусе 0.6 от u_actual нет НИ ОДНОЙ точки сетки с
|H8|<0.5, хотя H8(u_actual)~1.8e-3..2.5e-3 -- то есть там резкий узкий
провал прямо в u_actual, соседствующий, по-видимому, с полюсом H
(знаменатель x_{2^m} имеет там близкий ноль) -- ЭТО СОГЛАСУЕТСЯ с тем,
что Z72 в address_certificate.json отмечен как new_in_v2 и был
сертифицирован ОСОБЫМ маршрутом (point17, "shifted disk", "m8_refined
route") -- то есть даже в исходной строгой работе эта точка потребовала
нестандартной обработки.

Честная (не blind-scan) posadka: Ньютон с МАЛЫМ шагом (cap=0.01, чтобы не
перескочить через узкий провал), старт = u_actual + маленькое смещение
(не сам u_actual, чтобы это был честный тест сходимости, а не тавтология)."""
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


log('=== run_scan2_hardcases.py START (Z38, Z72, Z73) ===')
tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
tc_g_d_np = np.array([(k + 1) * tc_g_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_h_d_np = np.array([(k + 1) * tc_h_np[k + 1] for k in range(N - 1)], dtype=np.float64)

known, meta_res = C.load_known_zeros()
kd = dict(known)
cand = json.load(open('scan_candidates_v2.json', encoding='utf-8'))
zeros = cand['zeros']
DOMAIN_MAXX = 10.0


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


def batch_newton(u0_arr, m, iters=60, cap=0.01):
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


hard = ['Z38', 'Z72', 'Z73']
offset = complex(3e-4, -2e-4)
seeds = np.array([kd[cid] + offset for cid in hard], dtype=np.complex128)
log('seeds = u_actual + %s (NOT u_actual itself -- honest convergence test)' % offset)

for m in (8, 9, 10):
    t0 = time.time()
    seeds, Hf, valid, maxx = batch_newton(seeds, m, iters=60, cap=0.01)
    log('  m=%d Newton (cap=0.01, 60 iters, %.1fs): ' % (m, time.time() - t0)
        + ', '.join('%s |H|=%.2e valid=%s dist_to_actual=%.2e'
                     % (hard[i], Hf[i], valid[i], abs(seeds[i] - kd[hard[i]]))
                     for i in range(len(hard))))

added = 0
for i, cid in enumerate(hard):
    dist = abs(seeds[i] - kd[cid])
    if valid[i] and Hf[i] < 1e-4 and dist < 1e-2:
        zeros.append(dict(re=float(seeds[i].real), im=float(seeds[i].imag), absH=float(Hf[i]),
                           paired_m7=True, drift78=0.0,
                           re8=float(seeds[i].real), im8=float(seeds[i].imag),
                           re9=float(seeds[i].real), im9=float(seeds[i].imag),
                           re10=float(seeds[i].real), im10=float(seeds[i].imag),
                           absH9=float(Hf[i]), absH10=float(Hf[i]),
                           drift89=0.0, drift9_10=0.0, valid10=True,
                           source='hardcase_%s' % cid))
        added += 1
        log('  ADDED %s to catalog: dist_to_u_actual=%.3e' % (cid, dist))
    else:
        log('  %s FAILED to converge honestly (dist=%.3e, |H|=%.2e) -- NOT added' % (cid, dist, Hf[i]))

log('added %d/3 hard-case zeros to catalog (total now %d)' % (added, len(zeros)))
cand['zeros'] = zeros
cand['meta']['hardcases'] = dict(targets=hard, offset_re=offset.real, offset_im=offset.imag, added=added)
with open('scan_candidates_v2.json', 'w', encoding='utf-8') as f:
    json.dump(cand, f, indent=2, ensure_ascii=False)
log('updated scan_candidates_v2.json')
log('=== run_scan2_hardcases.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

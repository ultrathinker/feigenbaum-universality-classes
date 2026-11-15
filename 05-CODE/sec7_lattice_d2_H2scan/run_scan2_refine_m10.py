# -*- coding: utf-8 -*-
"""S2/H2SCAN v2, часть 1b -- довести float64-Ньютон с m=8 до m=9 и m=10
(m_level=10 -- ИМЕННО тот уровень, на котором построены u_actual в
results.json, meta.m_level=10). Для 3 нулей (Z38, Z72, Z73) m=8-приближение
давало расстояние 0.3-0.9 от u_actual при том что m=7<->m=8 дрейф был мал
(<0.05) -- то есть конечно-m сходимость H_m -> H для этих точек медленнее,
и апелляция к тому же m=10, на котором сертифицирован u_actual, должна
снять расхождение (либо честно показать, что это разные нули)."""
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


log('=== run_scan2_refine_m10.py START ===')
tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
tc_g_d_np = np.array([(k + 1) * tc_g_np[k + 1] for k in range(N - 1)], dtype=np.float64)
tc_h_d_np = np.array([(k + 1) * tc_h_np[k + 1] for k in range(N - 1)], dtype=np.float64)

known, meta_res = C.load_known_zeros()
cand = json.load(open('scan_candidates_v2.json', encoding='utf-8'))
zeros = cand['zeros']
log('loaded %d candidate zeros (m=8) from scan_candidates_v2.json' % len(zeros))

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


def batch_newton(u0_arr, m, iters=25, cap=0.5):
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


u8 = np.array([complex(z['re'], z['im']) for z in zeros], dtype=np.complex128)

log('--- Newton refine m=9 (seed = m=8 result), cap=0.5 (small steps: track same root) ---')
t0 = time.time()
u9, H9, valid9, maxx9 = batch_newton(u8, 9, iters=25, cap=0.5)
log('  m=9 done in %.1fs; valid %d/%d, median|H9|=%.2e' % (time.time() - t0, int(valid9.sum()), len(zeros), float(np.median(H9))))

log('--- Newton refine m=10 (seed = m=9 result), cap=0.5 ---')
t0 = time.time()
u10, H10, valid10, maxx10 = batch_newton(u9, 10, iters=25, cap=0.5)
log('  m=10 done in %.1fs; valid %d/%d, median|H10|=%.2e' % (time.time() - t0, int(valid10.sum()), len(zeros), float(np.median(H10))))

drift89 = np.abs(u9 - u8)
drift9_10 = np.abs(u10 - u9)
log('drift m8->m9: median=%.3e max=%.3e ; drift m9->m10: median=%.3e max=%.3e'
    % (np.median(drift89), np.max(drift89), np.median(drift9_10), np.max(drift9_10)))

for i, z in enumerate(zeros):
    z['re8'], z['im8'] = z['re'], z['im']
    z['re9'], z['im9'] = float(u9[i].real), float(u9[i].imag)
    z['re10'], z['im10'] = float(u10[i].real), float(u10[i].imag)
    z['absH9'] = float(H9[i])
    z['absH10'] = float(H10[i])
    z['drift89'] = float(drift89[i])
    z['drift9_10'] = float(drift9_10[i])
    z['valid10'] = bool(valid10[i] and H10[i] < 1e-4)
    # use m=10 as the authoritative position (matches u_actual's own m_level=10)
    z['re'] = z['re10']
    z['im'] = z['im10']

log('--- validation (float64 m=10, tol=1e-3) ---')
cat_u = np.array([complex(z['re10'], z['im10']) for z in zeros])
recovered = 0
for cid, u in known:
    d = np.abs(cat_u - u)
    j = int(np.argmin(d))
    dist = float(d[j])
    ok_v = dist <= 1e-3
    recovered += int(ok_v)
    log('  %-5s dist=%.3e %s  (H8=%.2e H9=%.2e H10=%.2e)'
        % (cid, dist, 'OK' if ok_v else 'MISS', zeros[j]['absH'], zeros[j]['absH9'], zeros[j]['absH10']))
log('VALIDATION (m=10 float64): %d / 23' % recovered)

cand['zeros'] = zeros
cand['meta']['m910_refine'] = dict(drift89_median=float(np.median(drift89)),
                                    drift89_max=float(np.max(drift89)),
                                    drift9_10_median=float(np.median(drift9_10)),
                                    drift9_10_max=float(np.max(drift9_10)),
                                    validation_m10=recovered)
with open('scan_candidates_v2.json', 'w', encoding='utf-8') as f:
    json.dump(cand, f, indent=2, ensure_ascii=False)
log('updated scan_candidates_v2.json with m=9,m=10 positions (authoritative u = m=10)')
log('=== run_scan2_refine_m10.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

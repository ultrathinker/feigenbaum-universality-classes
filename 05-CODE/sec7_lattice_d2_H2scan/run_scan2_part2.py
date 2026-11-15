# -*- coding: utf-8 -*-
"""S2/H2SCAN v2, часть 2 -- финальная валидация на всех 23 (float64, m=10 для
219 основных + честный m=8/9/10-Ньютон для 3 трудных), точечная mpmath-
полировка (dps=25, m=8) на выборке (все 23 совпадения + случайные 20 из
каталога) для подтверждения согласия float64<->mpmath, сборка финального
h2_zeros_scan.json (с сопряжениями Im<0)."""
import json
import random
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


log('=== run_scan2_part2.py START (final validate + catalog build) ===')
tc_g_mp, tc_h_mp, tc_g_np, tc_h_np, N = C.load_gh_taylor(log)
tc_g_d_mp = [(k + 1) * tc_g_mp[k + 1] for k in range(N - 1)]
tc_h_d_mp = [(k + 1) * tc_h_mp[k + 1] for k in range(N - 1)]

known, meta_res = C.load_known_zeros()
cand = json.load(open('scan_candidates_v2.json', encoding='utf-8'))
zeros = cand['zeros']
log('loaded %d candidate zeros' % len(zeros))

# ---------------------------------------------------------------------- #
#  final validation (float64) against 23 known
# ---------------------------------------------------------------------- #
cat_u = np.array([complex(z['re'], z['im']) for z in zeros])
recovered = 0
val_rows = []
for cid, u in known:
    d = np.abs(cat_u - u)
    j = int(np.argmin(d))
    dist = float(d[j])
    ok_v = dist <= 1e-3
    recovered += int(ok_v)
    val_rows.append(dict(id=cid, u_actual_re=u.real, u_actual_im=u.imag,
                          found_re=zeros[j]['re'], found_im=zeros[j]['im'],
                          dist=dist, ok=ok_v, match_idx=j))
    log('  %-5s dist=%.3e %s' % (cid, dist, 'OK' if ok_v else 'MISS'))
log('FINAL VALIDATION (float64): %d / 23' % recovered)

# ---------------------------------------------------------------------- #
#  mpmath spot-check polish: all 23 matched indices + random extra sample
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


matched_idx = sorted(set(r['match_idx'] for r in val_rows))
random.seed(20260815)
extra_pool = [i for i in range(len(zeros)) if i not in matched_idx]
extra = random.sample(extra_pool, min(20, len(extra_pool)))
sample_idx = sorted(set(matched_idx + extra))
log('mpmath spot-check on %d points (23 matched + %d random extra), m=8 dps=25' % (len(sample_idx), len(extra)))

t0 = time.time()
mp_errs = []
for cnt, i in enumerate(sample_idx):
    z = zeros[i]
    u0 = complex(z['re'], z['im'])
    u_mp, absH_mp = mp_polish_light(u0, m=8, iters=3)
    d_f64_vs_mp = abs(complex(float(mpre(u_mp)), float(mpim(u_mp))) - u0)
    mp_errs.append(d_f64_vs_mp)
    z['mp_checked'] = True
    z['mp_re'] = nstr(mpre(u_mp), 16)
    z['mp_im'] = nstr(mpim(u_mp), 16)
    z['mp_absH'] = absH_mp
    z['mp_vs_f64_dist'] = d_f64_vs_mp
    if (cnt + 1) % 10 == 0 or (cnt + 1) == len(sample_idx):
        log('  spot-check %d/%d (%.0fs)' % (cnt + 1, len(sample_idx), time.time() - t0))
log('mpmath spot-check done in %.1fs; float64<->mpmath median dist=%.3e max=%.3e'
    % (time.time() - t0, float(np.median(mp_errs)), float(np.max(mp_errs))))

# ---------------------------------------------------------------------- #
#  final catalog with conjugates
# ---------------------------------------------------------------------- #
final = []
for z in zeros:
    entry = dict(re=repr(z['re']), im=repr(z['im']), absH=z['absH'],
                 paired_m7=z.get('paired', True), drift78=z.get('drift78'))
    final.append(entry)
    if abs(z['im']) > 1e-9:
        final.append(dict(re=repr(z['re']), im=repr(-z['im']), absH=z['absH'],
                           paired_m7=z.get('paired', True), drift78=z.get('drift78'),
                           conjugate=True))
log('final catalog with conjugates: %d entries' % len(final))

out = dict(
    meta=dict(
        recipe='x_{n+1}=g(w)+s*h(w), w=(x^2-1)/r_disc, address_certificate.json convention',
        AU0=C.AU0_STR, delta=C.DELTA_STR, r_disc=C.R_DISC_STR,
        method='mixed-resolution grid (coarse 1/16 full box + fine 1/64 dense sub-box '
               'Re[1,10]xIm[3,9]) -> batch float64 Newton (analytic grad) at m=8, '
               'refined to m=9,m=10 for all candidates (authoritative position = m=10); '
               '3 hard cases (Z38,Z72,Z73 -- narrow basin near a nearby pole of H, honest '
               'small-step Newton from u_actual+offset, NOT literal seeding at u_actual) '
               'added separately (run_scan2_hardcases.py); mpmath spot-check polish '
               '(dps=25, m=8) on 23 matched + 20 random extra to confirm float64 precision',
        domain_maxx=10.0, cand_thresh=2.5,
        grids=cand['meta']['grids'],
        m910_refine=cand['meta'].get('m910_refine'),
        hardcases=cand['meta'].get('hardcases'),
        elapsed_sec=time.time() - T0,
    ),
    counts=dict(unique_zeros=len(zeros), final_with_conjugates=len(final)),
    validation=dict(recovered=recovered, total=23, tolerance=1e-3, rows=val_rows),
    mpmath_spotcheck=dict(n_checked=len(sample_idx),
                           median_dist_f64_vs_mp=float(np.median(mp_errs)),
                           max_dist_f64_vs_mp=float(np.max(mp_errs))),
    zeros=final,
)
with open('h2_zeros_scan.json', 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=2, ensure_ascii=False)
log('wrote h2_zeros_scan.json (%d zeros incl. conjugates, validation %d/23)'
    % (len(final), recovered))
log('=== run_scan2_part2.py END (elapsed %.1fs) ===' % (time.time() - T0))
LOG.close()

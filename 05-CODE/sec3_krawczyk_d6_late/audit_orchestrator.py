# -*- coding: utf-8 -*-
"""audit_orchestrator.py -- third code path (orchestrator) for the 9 late boxes.
Own orbit code (mpmath, complex a, f(x)=1-a*x^6), independent of Arb pipeline.
Checks per box: |x_P(c)| tiny; min_j |x_{2^j}(c)| matches table order; a
finite-difference derivative magnitude agrees in order with df_abs_lb; the
center is far from every one of the 287 wave-2 boxes and from other 8.
"""
import json, sys
from mpmath import mp, mpc, mpf, fabs

mp.dps = 90
D = 6

def orbit_stats(a, P, m):
    x = mpf(0)
    mins = {}
    for k in range(1, P + 1):
        x = 1 - a * x**D
        if k & (k - 1) == 0 and k <= 2**m:  # power of two up to 2^m
            mins[k] = fabs(x)
    return x, mins

cert = json.load(open('krawczyk6_late_certificates.json', encoding='utf-8'))
w2 = json.load(open('../sec3_krawczyk_d6/krawczyk6_certificates.json', encoding='utf-8'))
w2roots = w2['roots'] if isinstance(w2.get('roots'), list) else list(w2['roots'].values())
import re as _re
def _num(s):
    """parse plain decimal or Arb-style '[mid +/- rad]' -> upper bound as mpf"""
    s = str(s).strip()
    if s.startswith('['):
        mm = _re.match(r'\[\s*([-+0-9.eE]+)\s*\+/-\s*([-+0-9.eE]+)\s*\]', s)
        if mm:
            return fabs(mpf(mm.group(1))) + mpf(mm.group(2))
        s = s.strip('[]').split()[0]
    return mpf(s)
w2c = [(r['label'], mpc(_num(r['box_center_re']), _num(r['box_center_im'])), _num(r['box_radius'])) for r in w2roots]

ok_all = True
rows = []
centers = []
for r in cert['roots']:
    lab, m, P = r['label'], int(r['m']), int(r['n'])
    assert P == 2**(m+1), (lab, m, P)
    a = mpc(mpf(r['box_center_re']), mpf(r['box_center_im']))
    rad = _num(r['box_radius'])
    xP, mins = orbit_stats(a, P, m)
    res = fabs(xP)
    minlb = min(mins.values())
    # finite-difference derivative (own path)
    h = mpf('1e-30')
    xPh, _ = orbit_stats(a + h, P, m)
    dF = fabs((xPh - xP) / h)
    df_lb = _num(r['df_abs_lb'])
    newton = res / dF
    # isolation vs wave2 + others
    dmin = min(fabs(a - c2) - rad - r2 for _, c2, r2 in w2c)
    centers.append((lab, a, rad))
    good = (res < mpf('1e-40')) and (minlb > mpf('1e-3')) and (dF > df_lb*mpf('0.5')) and (newton < rad*mpf('1e-6')) and (dmin > 0)
    ok_all &= good
    rows.append((lab, m, P, res, minlb, _num(r['period_min_lb']), dF, df_lb, newton, dmin, good))

# pairwise among 9
pmin = None
for i in range(len(centers)):
    for j in range(i+1, len(centers)):
        g = fabs(centers[i][1]-centers[j][1]) - centers[i][2] - centers[j][2]
        pmin = g if pmin is None or g < pmin else pmin

print('label m P | |x_P(c)| | min|x_2^j| (mine / table lb) | |F\'| fd / table lb | newton_dist | gap to nearest w2 box | OK')
for lab, m, P, res, minlb, tlb, dF, df_lb, newton, dmin, good in rows:
    print('%s %2d %5d | %.2e | %.4f / %.4f | %.2e / %.2e | %.2e | %.2e | %s' % (
        lab, m, P, res, minlb, tlb, dF, df_lb, newton, dmin, 'OK' if good else 'FAIL'))
print('pairwise min gap among 9 (minus radii): %.3e' % pmin)
print('ALL_OK =', ok_all)
sys.exit(0 if ok_all and pmin > 0 else 1)

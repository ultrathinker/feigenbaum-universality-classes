# -*- coding: utf-8 -*-
"""
f4_predcheck.py -- проверка предсказаний P1-P5 (записаны в progress.log
19:36:30 ДО этого прогона): наличие/отсутствие членов башни alpha^{-dk}
по закону "dk чётное целое". Обратная итерация на N=40 и N=48, два усечения.
[oid:f4GLM] NOT RIGOROUS.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C
from f4_test3_dictfit import build_DT, inv_iter

mp.dps = 45
T0 = time.time()

SG = json.load(open(os.path.join(C.ROOT, 'wave1_spectra_d', 'spectra_grid.json'),
                   encoding='utf-8'))
RECS = {r['d']: r for r in SG['points']}


def spow(y, p):
    if y == 0:
        return mpf(0)
    return (mpf(1) if y > 0 else -1) * abs(y) ** p


CASES = [  # (d_str, k, ожидание)
    ('3.0', 2, 'P2-like: есть (нечётное d, чётное k)'),
    ('3.0', 4, 'P1: есть'),
    ('5.0', 2, 'P2: есть'),
    ('7.0', 2, 'P3: есть'),
    ('2.5', 4, 'P4: есть (рациональное d=5/2, k0=4)'),
    ('2.5', 1, 'P4: НЕТ'),
    ('2.5', 2, 'P5: НЕТ (dk=5 нечётное)'),
    ('2.5', 3, 'P4: НЕТ'),
]

out = {}
for dstr, k, expect in CASES:
    rec = RECS[dstr]
    a = [mpf(s) for s in rec['_a48']]
    alpha = mpf(rec['alpha48'])
    d = mpf(dstr)
    mu = alpha ** (-d * k)
    row = {'expect': expect, 'alpha': mp.nstr(alpha, 24), 'mu': mp.nstr(mu, 20),
           'dk': float(d * k)}
    # кандидат (для косинуса)
    def Phi(t):
        G_ = C.horner(a, t)
        return (d * t ** (k + 1) * C.Gp(a, t) - spow(G_, d * k + 1) + G_
                - d * t * C.Gp(a, t))
    for N in (40, 48):
        A = build_DT(a, alpha, d, N)
        lam, v, mres = inv_iter(A, mu, iters=70)
        mx = max(abs(x) for x in v)
        vn = [x / mx for x in v]

        def Vf(t, vn=vn):
            return C.horner([mpf(0)] + vn, t)
        nodes = C.grid_t(50, mpf('1.0'))
        num = sum(Phi(t) * Vf(t) for t in nodes)
        n1 = mp.sqrt(sum(Phi(t) ** 2 for t in nodes))
        n2 = mp.sqrt(sum(Vf(t) ** 2 for t in nodes))
        row['N%d' % N] = {'lam': mp.nstr(lam, 20),
                          'rel_to_mu': float(abs(lam - mu) / mu),
                          'mat_resid': float(mres),
                          'cos_phi': float(abs(num / (n1 * n2)))}
    row['N_drift'] = float(abs(mpf(row['N40']['lam']) - mpf(row['N48']['lam']))
                           / max(abs(mpf(row['N48']['lam'])), mpf('1e-30')))
    out['d=%s_k=%d' % (dstr, k)] = row
    C.plog('PRED d=%s k=%d (dk=%s): %s | mu=%s | N40 rel=%.3e N48 rel=%.3e drift=%.1e cos=%.6f'
           % (dstr, k, row['dk'], expect, mp.nstr(mu, 10),
              row['N40']['rel_to_mu'], row['N48']['rel_to_mu'], row['N_drift'],
              row['N48']['cos_phi']))

with open(os.path.join(C.HERE, 'f4_predcheck.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
C.plog('PRED DONE (%.0fs) -> f4_predcheck.json' % (time.time() - T0))

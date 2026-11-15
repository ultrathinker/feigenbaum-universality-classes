# -*- coding: utf-8 -*-
"""
f4_extras.py -- финальные числа для отчёта:
 (1) остатки тождества на [0,1] (уровень точности G) vs [0,1.2];
 (2) d=6: собственный вектор при N=48, косинус в коэф. пространстве с Phi-рядом(48);
 (3) контроль: DT_bare(K)/K = 1 и полный DT(K) = alpha^d K.
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C
from f4_test3_dictfit import build_DT, inv_iter

mp.dps = 50
out = {}

C.plog('EXTRAS START')

# (1) остатки на [0,1]
out['resid_0_1'] = {}
for tag in ('2.0', '4.0', '6.0'):
    a, alpha, d, meta = C.load_fp(tag)
    out['resid_0_1'][tag] = []
    for k in (1, 2):
        Phi = C.cand_C2x if False else (lambda t, k=k:
                                        d * t ** (k + 1) * C.Gp(a, t)
                                        - C.horner(a, t) ** (d * k + 1)
                                        + C.horner(a, t) - d * t * C.Gp(a, t))
        mu = alpha ** (-d * k)
        rr1, _, _ = C.resid_table(a, alpha, d, Phi, mu, npts=200, tmax=mpf(1))
        rr2, _, _ = C.resid_table(a, alpha, d, Phi, mu, npts=200, tmax=mpf('1.2'))
        out['resid_0_1'][tag].append({'k': k, 'rel_0_1': float(rr1),
                                      'rel_0_12': float(rr2)})
        C.plog('  d=%s k=%d: |r|/max|Phi| на [0,1] = %.3e ; на [0,1.2] = %.3e'
               % (tag, k, float(rr1), float(rr2)))

# (2) d=6, N=48 eigenvector + косинус в коэф. пространстве
a6, alpha6, d6, _ = C.load_fp('6.0')
A48 = build_DT(a6, alpha6, d6, 48)
lam6, v6, mres = inv_iter(A48, alpha6 ** (-6), iters=90)
from f4_main import phi_series  # noqa: E402
Phis = phi_series(a6, d6, 1, 48)
mx = max(abs(x) for x in v6)
v6n = [x / mx for x in v6]
num = sum(Phis[i + 1] * v6n[i] for i in range(48))
n1 = mp.sqrt(sum(x * x for x in Phis[1:]))
n2 = mp.sqrt(sum(x * x for x in v6n))
ccoef = abs(num / (n1 * n2))
out['d6_N48'] = {'lam': mp.nstr(lam6, 18),
                 'rel_to_mu': float(abs(lam6 - alpha6 ** (-6)) / alpha6 ** (-6)),
                 'cos_coef_48': float(ccoef)}
C.plog('  d=6 N=48: lam=%s (rel %.2e), cos_коэф(48)=%.8f'
       % (mp.nstr(lam6, 14), out['d6_N48']['rel_to_mu'], float(ccoef)))

# (3) контроль координатной моды
for tag in ('2.0', '6.0'):
    a, alpha, d, _ = C.load_fp(tag)
    K = C.K_fun(a, d)
    nodes = C.grid_t(60, mpf(1))
    bare = [C.DT_unproj_t(a, alpha, d, K, t) / K(t) for t in nodes]
    full = [C.DT_t(a, alpha, d, K, t) / K(t) for t in nodes]
    out['coord_%s' % tag] = {
        'bare_ratio_minmax': [float(min(bare)), float(max(bare))],
        'full_ratio_minmax': [float(min(full)), float(max(full))],
        'alpha_pow_d': mp.nstr(alpha ** d, 16)}
    C.plog('  d=%s: DT_bare(K)/K = %.17f..%.17f ; DT(K)/K = %.6f..%.6f (alpha^d=%s)'
           % (tag, min(bare), max(bare), min(full), max(full), mp.nstr(alpha ** d, 10)))

import json
with open(os.path.join(C.HERE, 'f4_extras.json'), 'w', encoding='utf-8') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
C.plog('EXTRAS DONE')

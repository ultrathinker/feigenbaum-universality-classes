# -*- coding: utf-8 -*-
"""
f4_test6_layers.py -- многослойная комбинация:
  Phi~_k = d t^{k+1} G' - G^{dk+1} + G - d t G'   (слои k и 0, спаривание b=-a)
Вопрос: остаток DT[Phi~] - alpha^{-dk} Phi~ точно пропорционален K?
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C

mp.dps = 50

C.plog('TEST6 START: слоистая комбинация, d=2/4/6, k=1,2')

for tag in ('2.0', '4.0', '6.0'):
    a, alpha, d, meta = C.load_fp(tag)
    G = lambda t: C.horner(a, t)
    Gp = lambda t: C.Gp(a, t)
    K = lambda t: G(t) - d * t * Gp(t)
    for k in (1, 2):
        mu = alpha ** (-d * k)

        def PhiT(t):
            return (d * t ** (k + 1) * Gp(t) - G(t) ** (d * k + 1) + G(t)
                    - d * t * Gp(t))

        # остаток
        res = [C.DT_t(a, alpha, d, PhiT, t) - mu * PhiT(t)
               for t in C.grid_t(80, mpf('1.05'))]
        rmax = max(abs(x) for x in res)
        pmax = max(abs(PhiT(t)) for t in C.grid_t(80, mpf('1.05')))
        # регрессия res = c*K
        nodes = C.grid_t(80, mpf('1.05'))
        Kn = [K(t) for t in nodes]
        num = sum(res[i] * Kn[i] for i in range(80))
        den = sum(Kn[i] * Kn[i] for i in range(80))
        cc = num / den
        fitmax = max(abs(res[i] - cc * Kn[i]) for i in range(80))
        C.plog('d=%s k=%d: max|Phi|=%s  max|res|=%.3e  коэф.при K=%s  непрокрытый остаток=%.3e'
               % (tag, k, mp.nstr(pmax, 6), float(rmax), mp.nstr(cc, 10), float(fitmax)))
        # прогноз C_k из алгебры: (1-mu) + g'(1)(1-alpha) + mu/alpha - 1/alpha
        gp1 = d * Gp(mpf(1))
        Ck = (1 - mu) + gp1 * (1 - alpha) + mu / alpha - 1 / alpha
        C.plog('         алгебраический C_k=%s (g\'(1)=%s=-alpha^{d-1}? %s)'
               % (mp.nstr(Ck, 10), mp.nstr(gp1, 10),
                  mp.nstr(-alpha ** (d - 1), 10)))

C.plog('TEST6 DONE')

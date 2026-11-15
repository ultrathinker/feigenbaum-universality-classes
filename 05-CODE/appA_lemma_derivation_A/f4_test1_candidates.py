# -*- coding: utf-8 -*-
"""
f4_test1_candidates.py -- разведка: какой явный кандидат даёт тождество
DT Phi = alpha^{-dk} Phi при чётных d. Кандидаты:
  C1: d t^{k+1} G'                  (касание g(psi(x)), psi=x(1+eps x^{dk}))
  C2: G^{dk+1} - d t^{k+1} G' - G   (сопряжение с нормировкой f(0)=1)
  C3: t^k K,  C4: t^k G             (полиномиальные умножения)
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C

mp.dps = 60

C.plog('TEST1 START: кандидаты башенной моды, dps=%d' % mp.dps)

for tag in ('2.0', '4.0', '6.0'):
    a, alpha, d, meta = C.load_fp(tag)
    C.plog('d=%s: %s, len(a)=%d, alpha=%s' % (tag, meta['src'], len(a), mp.nstr(alpha, 20)))
    for k in (1, 2):
        lam = alpha ** (-d * k)
        row = []
        for name, Phi in (('C1', C.cand_C1(a, d, k)),
                          ('C2', C.cand_C2(a, d, k)),
                          ('C3', lambda t, k=k: t ** k * C.K_fun(a, d)(t)),
                          ('C4', lambda t, k=k: t ** k * C.horner(a, t))):
            rr, wt, pm = C.resid_table(a, alpha, d, Phi, lam, npts=200)
            row.append('%s: %.3e (worst t=%s)' % (name, float(rr), mp.nstr(wt, 4)))
        C.plog('  d=%s k=%d lam=alpha^-%d=%s | %s' %
               (tag, k, d * k, mp.nstr(lam, 16), ' | '.join(row)))

C.plog('TEST1 DONE')

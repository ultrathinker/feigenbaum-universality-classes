# -*- coding: utf-8 -*-
"""
f4_test2_structure.py -- реверс-инжиниринг: структура численного собственного
вектора башенной моды (d=2, k=1). Отношения к G', K, G; проверка Phi(1)=0;
подгонка Phi/G' полиномом.
[oid:f4GLM] NOT RIGOROUS.
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C

mp.dps = 40

C.plog('TEST2 START: структура собственного вектора башенной моды d=2')

a, alpha, d, meta = C.load_fp('2.0')
j = json.load(open(os.path.join(C.ROOT, 'wave2_vector_homotopy',
                                'spectrum_window.json'), encoding='utf-8'))
p = j['points']['2.0']
lams = [mpf(s) for s in p['lams']]
tow = alpha ** (-2)
i = min(range(len(lams)), key=lambda ii: abs(lams[ii] - tow))
C.plog('мода idx %d: lam=%s (alpha^-2=%s)' % (i, mp.nstr(lams[i], 12), mp.nstr(tow, 12)))
v = [mpf(s) for s in p['rvecs'][i]]


def Phi(t):
    return C.horner([mpf(0)] + v, t)


G1 = C.horner(a, mpf(1))
C.plog('Phi(1) = %s   (G(1)=%s, K(1)=%s)' %
       (mp.nstr(Phi(mpf(1)), 8), mp.nstr(G1, 8),
        mp.nstr(C.K_fun(a, d)(mpf(1)), 8)))
C.plog('c_j (j=1..24): %s' % [mp.nstr(x, 4) for x in v[:24]])

print()
C.plog('t      Phi/G''             Phi/(t*G'')        Phi/K            Phi/G')
for t in [mpf(x) / 10 for x in range(1, 13)]:
    Gp_ = C.Gp(a, t)
    K_ = C.K_fun(a, d)(t)
    G_ = C.horner(a, t)
    C.plog('%5s  %s  %s  %s  %s' % (mp.nstr(t, 4), mp.nstr(Phi(t) / Gp_, 8),
                                     mp.nstr(Phi(t) / (t * Gp_), 8),
                                     mp.nstr(Phi(t) / K_, 8), mp.nstr(Phi(t) / G_, 8)))

# fit Phi/G' по полиномам степени 1..6 (МНК по 40 узлам, через реш. Вандермонда)
import functools

nodes = [mpf('0.05') + mpf('1.1') * i / 39 for i in range(40)]


def lstsq_poly(f, deg):
    A = [[tt ** q for q in range(deg + 1)] for tt in nodes]
    b = [f(tt) for tt in nodes]
    n = deg + 1
    # normal equations (mpf)
    AtA = [[sum(A[r][i] * A[r][j2] for r in range(len(nodes))) for j2 in range(n)] for i in range(n)]
    Atb = [sum(A[r][i] * b[r] for r in range(len(nodes))) for i in range(n)]
    M = [[AtA[i][j2] + (Atb[i] if j2 == n else 0)] for i in range(n) for j2 in range(n + 1)]
    # gaussian elim
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col:
                fac = M[r][col] / M[col][col]
                for cc in range(col, n + 1):
                    M[r][cc] -= fac * M[col][cc]
    return [M[i][n] / M[i][i] for i in range(n)]


for deg in (2, 3, 4, 5):
    co = lstsq_poly(lambda tt: Phi(tt) / C.Gp(a, tt), deg)
    resid = max(abs(functools.reduce(lambda acc, cq: acc * tt + cq, co[::-1])
                    - Phi(tt) / C.Gp(a, tt)) for tt in nodes)
    C.plog('deg %d: Phi/G'' ~ poly, coefs=%s, maxresid=%.3e'
           % (deg, [mp.nstr(x, 5) for x in co], float(resid)))

C.plog('TEST2 DONE')

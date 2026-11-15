# -*- coding: utf-8 -*-
"""
f4_test4_vg_form.py -- тест гипотезы: башенный собственный вектор имеет вид
Phi = v*g' (сопряжённый тип), где v решает (E1):
    v(g(y)) + g'(y) v(y) = alpha^{-dk-1} g'(y) v(alpha y),   v нечётна, v(1)=0.
Проверяем (E1) численно на v_num = Phi_num/g' (свой eigenvector, d=2, k=1),
и заодно k=2 (alpha^-4).
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C
from f4_test3_dictfit import build_DT, inv_iter

mp.dps = 55

C.plog('TEST4 START: гипотеза Phi = v*g'' (E1), d=2, k=1,2')

a, alpha, d, meta = C.load_fp('2.0')
A40 = build_DT(a, alpha, d, 40)


def g_x(x):
    # g(x) для x>=0: G(x^d)
    return C.horner(a, x ** d)


def gp_x(x):
    # g'(x), x>0: d x^{d-1} G'(x^d)
    return d * x ** (d - 1) * C.Gp(a, x ** d)


for k in (1, 2):
    mu = alpha ** (-d * k)
    lam, v40, mres = inv_iter(A40, mu)
    # нормировка по максимуму
    mx = max(abs(x) for x in v40)
    coef = [mpf(0)] + [x / mx for x in v40]

    def Phi(t):
        return C.horner(coef, t)

    # порядок обращения в нуль: c_1, c_2
    nz = next(j for j in range(1, 12) if abs(coef[j]) > mpf('1e-8'))
    C.plog('k=%d: lam=%s  первый значимый коэффициент c_%d=%s'
           % (k, mp.nstr(lam, 16), nz, mp.nstr(coef[nz], 6)))

    # v_num(x) = Phi(x^d)/g'(x)
    def v_num(x):
        return Phi(x ** d) / gp_x(x)

    # проверка (E1): v(g(y)) + g'(y) v(y) - alpha^{-dk-1} g'(y) v(alpha y) = 0
    lamE1 = alpha ** (-d * k - 1)
    worst = mpf(0)
    worst_y = None
    vn_max = mpf(0)
    for i in range(40):
        # зона надёжности: (alpha*y)^d <= 1.2  =>  y <= (1.2)^(1/d)/alpha ~ 0.44
        y = mpf('0.05') + mpf('0.35') * i / 39
        gy = g_x(y)
        if gy < 0:
            gy_abs = abs(gy)
            # v нечётная: v(gy) = -v(|gy|); g'(gy) = d gy^{d-1} G'(|gy|^d) (знак через степень)
            vg = -v_num(gy_abs)
            gpg = d * gy ** (d - 1) * C.Gp(a, gy_abs ** d)
        else:
            vg = v_num(gy)
            gpg = gp_x(gy)
        lhs = vg + gp_x(y) * v_num(y) - lamE1 * gp_x(y) * v_num(alpha * y)
        if abs(lhs) > worst:
            worst, worst_y = abs(lhs), y
        vn_max = max(vn_max, abs(v_num(y)))
    C.plog('  (E1) невязка: max|lhs|/max|v| = %.3e (worst y=%s)'
           % (float(worst / vn_max), mp.nstr(worst_y, 4)))
    C.plog('  v(1) численно = %s (max|v|=%s)' % (mp.nstr(v_num(1), 5), mp.nstr(vn_max, 5)))

C.plog('TEST4 DONE')

# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (ворота 2, задача F): третий кодовый путь для знака λ₃(4).
# Независимость: неподвижная точка — из ПИЛОТА (quartic_fixedpoint.json, не kilo-2);
# оператор — ненормированная форма (LΦ)(x) = -α[Φ(g(x/α)) + g'(g(x/α))Φ(x/α)]
# (у kilo-2 — нормированный DT с проекцией); базис x^{4j}, generalized eig mpmath.
import json, sys
from mpmath import mp, mpf, matrix, eig, fabs

mp.dps = 40

FP = json.load(open('../sec3_pilot_d4/quartic_fixedpoint.json', encoding='utf-8'))
al = mpf(FP['alpha4'])
gc = [mpf(repr(c)) if isinstance(c, float) else mpf(c) for c in FP['g_coeffs']]
print(f"[конст] alpha4={mp.nstr(al,20)}  N_coeffs={len(gc)}", flush=True)

def G(t):
    s = mpf(0)
    for c in reversed(gc): s = s*t + c
    return s
def Gp(t):
    s = mpf(0)
    for j in range(len(gc)-1, 0, -1): s = s*t + j*gc[j]
    return s
def g(x):  return G(x**4)
def gp(x): return 4*x**3*Gp(x**4)

# санити неподвижной точки: g(g(x/al))*(-al) == g(x)
for xt in (mpf('0.3'), mpf('0.7'), mpf('0.95')):
    r = -al*g(g(xt/al)) - g(xt)
    assert fabs(r) < mpf('1e-30'), f"fixed point residual {r}"
print("[санити] уравнение неподвижной точки: остатки < 1e-30", flush=True)

M = 26
# узлы: t_i чебышёвские на (0,1), x_i = t^(1/4)
from mpmath import cos, pi
ts = [(1+cos(pi*(2*i+1)/(2*M)))/2 for i in range(M)]
xs = [t**mpf('0.25') for t in ts]

A = matrix(M, M); B = matrix(M, M)
for i in range(M):
    x = xs[i]
    u = g(x/al)        # аргумент внешней Φ
    w = gp(u)          # множитель при Φ(x/α)
    for j in range(M):
        A[i, j] = -al*(u**(4*j) + w*(x/al)**(4*j))
        B[i, j] = ts[i]**j
E, _ = eig(B**-1 * A)
E = sorted(E, key=lambda z: -abs(z))
print("[спектр] топ-8 по модулю:")
for z in E[:8]:
    print(f"   {mp.nstr(z.real, 12)} {'+' if z.imag>=0 else '-'} {mp.nstr(abs(z.imag), 3)}i   |λ|={mp.nstr(abs(z), 12)}", flush=True)

# ищем моду с |λ| ~ 0.2557
target = mpf('0.255664558474')
cand = [z for z in E if fabs(abs(z) - target) < mpf('0.01')]
ok = False
if cand:
    z = cand[0]
    is_real_neg = fabs(z.imag) < 1e-8 and z.real < 0
    print(f"\n[λ₃(4)] найдено {mp.nstr(z.real, 14)} (im {mp.nstr(z.imag, 2)}); "
          f"заявка kilo-2: -0.255664558474; пилот был: +0.255664558474")
    ok = is_real_neg and fabs(fabs(z.real) - target) < mpf('1e-6')
else:
    print("\n[λ₃(4)] мода около 0.2557 не найдена!")
# заодно δ и λ₂
d4 = [z for z in E if fabs(z - mpf(FP['delta4'])) < mpf('1e-6')]
l2 = [z for z in E if fabs(z.imag) < 1e-8 and fabs(z.real - mpf('0.291838408984')) < mpf('1e-6')]
print(f"[якоря] delta4 найден: {bool(d4)}; lambda2=+0.29184 найден: {bool(l2)}")
print('\nВЕРДИКТ: ' + ('ЗНАК МИНУС ПОДТВЕРЖДЁН, якоря на месте' if (ok and d4 and l2) else 'ПРОВЕРИТЬ ВРУЧНУЮ'))
sys.exit(0 if (ok and d4 and l2) else 1)

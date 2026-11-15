# -*- coding: utf-8 -*-
# АУДИТ ОРКЕСТРАТОРА (задача M): третий кодовый путь.
# Моя ненормированная линеаризация (как в аудите d=4) + mpmath.eig:
# (1) d=2.0 (fp = _a48 сетки kilo-2): сколько мод в окне |λ-α⁻²|<0.02? (M: ровно одна)
# (2) d=2.05 (fp = НЕЗАВИСИМАЯ точка агента K): присутствуют ли ОБЕ моды +0.1645 и -0.1281.
import json, sys
from mpmath import mp, mpf, matrix, eig, fabs, cos, pi

mp.dps = 40

def spectrum(alpha, gc, d_int_pow, M=30):
    # g(x) = G(x^p), p = d (целочисленная степень как t=x^p; для 2.05 берём p как mpf)
    al = alpha
    def G(t):
        s = mpf(0)
        for c in reversed(gc): s = s*t + c
        return s
    def Gp(t):
        s = mpf(0)
        for j in range(len(gc)-1, 0, -1): s = s*t + j*gc[j]
        return s
    p = d_int_pow
    def g(x):  return G(x**p)
    def gp(x): return p*x**(p-1)*Gp(x**p)
    for xt in (mpf('0.3'), mpf('0.7')):
        r = -al*g(g(xt/al)) - g(xt)
        assert fabs(r) < mpf('1e-25'), f"fp residual {r}"
    ts = [(1+cos(pi*(2*i+1)/(2*M)))/2 for i in range(M)]
    xs = [t**(1/p) for t in ts]
    A = matrix(M, M); B = matrix(M, M)
    for i in range(M):
        x = xs[i]; u = g(x/al); w = gp(u)
        for j in range(M):
            A[i, j] = -al*(u**(p*j) + w*(x/al)**(p*j))
            B[i, j] = ts[i]**j
    E, _ = eig(B**-1 * A)
    return sorted([z for z in E if fabs(z.imag) < 1e-6], key=lambda z: -abs(z))

# (1) d=2.0 из сетки kilo-2
gpt = None
for pt in json.load(open('../sec2_spectra_grid_d/spectra_grid.json', encoding='utf-8'))['points']:
    if abs(float(pt['d']) - 2.0) < 1e-12: gpt = pt; break
al2 = mpf(pt['alpha48'] if isinstance(gpt['alpha48'], str) else repr(gpt['alpha48']))
gc2 = [mpf(c) if isinstance(c, str) else mpf(repr(c)) for c in gpt['_a48']]
E2 = spectrum(al2, gc2, mpf(2))
tower = al2**-2
win = [z for z in E2 if fabs(abs(z) - tower) < mpf('0.02')]
print(f"[d=2.0] α⁻²={mp.nstr(tower,10)}; мод в окне ±0.02: {len(win)} -> {[mp.nstr(z.real,10) for z in win]}")
near = sorted(E2, key=lambda z: fabs(abs(z)-tower))[:4]
print(f"        ближайшие: {[mp.nstr(z.real,8) for z in near]}")
has_l2 = any(fabs(z.real + mpf('0.123652712553')) < 1e-8 for z in E2)
has_l3 = any(fabs(z.real + mpf('0.057307021067')) < 1e-8 for z in E2)
print(f"        λ₂=-0.1237 присутствует: {has_l2}; λ₃=-0.0573: {has_l3}")

# (2) d=2.05 (fp из фланга сетки; моя сборка с |x|^d и моя eig)
gpt5 = None
for pt in json.load(open('../sec2_spectra_grid_d/spectra_grid.json', encoding='utf-8'))['points']:
    if abs(float(pt['d']) - 2.05) < 1e-12: gpt5 = pt; break
assert gpt5 is not None, "нет точки 2.05 в сетке"
al5 = mpf(gpt5['alpha48'] if isinstance(gpt5['alpha48'], str) else repr(gpt5['alpha48']))
gc5 = [mpf(c) if isinstance(c, str) else mpf(repr(c)) for c in gpt5['_a48']]
p = mpf('2.05')
def G5(t):
    s = mpf(0)
    for c in reversed(gc5): s = s*t + c
    return s
def G5p(t):
    s = mpf(0)
    for j in range(len(gc5)-1, 0, -1): s = s*t + j*gc5[j]
    return s
def g5(x):  return G5(fabs(x)**p)
def g5p(x):
    sgn = 1 if x >= 0 else -1
    return p*sgn*fabs(x)**(p-1)*G5p(fabs(x)**p)
for xt in (mpf('0.3'), mpf('0.7')):
    r = -al5*g5(g5(xt/al5)) - g5(xt)
    assert fabs(r) < mpf('1e-25'), f"fp5 residual {r}"
M5 = 30
ts = [(1+cos(pi*(2*i+1)/(2*M5)))/2 for i in range(M5)]
xs = [t**(1/p) for t in ts]
A = matrix(M5, M5); B = matrix(M5, M5)
for i in range(M5):
    x = xs[i]; u = g5(x/al5); w = g5p(u); ua = fabs(u)
    for j in range(M5):
        A[i, j] = -al5*(ua**(p*j) + w*(x/al5)**(p*j))
        B[i, j] = ts[i]**j
E5, _ = eig(B**-1 * A)
E5r = sorted([z for z in E5 if fabs(z.imag) < 1e-6], key=lambda z: -abs(z))
print(f"[d=2.05] топ-7 вещественных: {[mp.nstr(z.real, 9) for z in E5r[:7]]}")
t5 = al5**(-p)
pos_ok = any(fabs(z.real - mpf('0.164454844')) < 1e-4 for z in E5r)
neg_ok = any(fabs(z.real + mpf('0.128069980')) < 1e-4 for z in E5r)
near_t = [z for z in E5r if fabs(abs(z) - t5) < t5*mpf('1e-4')]
print(f"[d=2.05] +0.16445 присутствует: {pos_ok}; -0.12807 присутствует: {neg_ok}; "
      f"мод точно на α^-d={mp.nstr(t5,8)}: {len(near_t)}")
print('\nВЕРДИКТ:', 'ЗЕЛЁНЫЙ (кратность d=2 простая; обе спорные моды при 2.05 — честные простые)'
      if (len(win) == 1 and pos_ok and neg_ok) else 'ПРОВЕРИТЬ ВРУЧНУЮ')
sys.stdout.flush()

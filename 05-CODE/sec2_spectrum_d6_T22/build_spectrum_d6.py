# -*- coding: utf-8 -*-
# T2.2 d=6: полный вещественный спектр линеаризации DT при d=6, поиск калибровочной
# ветви башни alpha^{-6k} у alpha6^{-6}=0.10002322865..., два усечения N=32/48,
# функциональный (не только матричный) остаток, тест координатной моды.
# NOT RIGOROUS. Один процесс, mpmath dps=50.
import json, sys, time, os
from mpmath import mp, mpf, matrix, eig, fabs, cos, pi, nstr

mp.dps = 50
os.chdir(os.path.dirname(os.path.abspath(__file__)))
LOG = open('progress.log', 'a', encoding='utf-8')

def log(msg):
    t = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"{t} {msg}"
    print(line)
    LOG.write(line + '\n')
    LOG.flush()

log("build_spectrum_d6.py: старт")

GRID_PATH = r'../sec2_spectra_grid_d/spectra_grid.json'
pts = json.load(open(GRID_PATH, encoding='utf-8'))['points']
pt6 = None
for pt in pts:
    if abs(float(pt['d']) - 6.0) < 1e-9:
        pt6 = pt
        break
assert pt6 is not None, "d=6.0 не найдено в сетке"

alpha = mpf(pt6['alpha48'])
gc = [mpf(c) for c in pt6['_a48']]
p = 6  # d=6, целая чётная степень
log(f"alpha6 = {nstr(alpha, 25)}  (бриф даёт 1.4677424503199009)")
log(f"len(_a48) coeffs = {len(gc)}")

def G(t):
    s = mpf(0)
    for c in reversed(gc):
        s = s * t + c
    return s

def Gp(t):
    s = mpf(0)
    for j in range(len(gc) - 1, 0, -1):
        s = s * t + j * gc[j]
    return s

def g(x):
    return G(x**p)

def gp(x):
    return p * x**(p - 1) * Gp(x**p)

# --- Санити: неподвижная точка g(x) = -alpha*g(g(x/alpha)) ---
max_fp_resid = mpf(0)
for xt in (mpf('0.1'), mpf('0.3'), mpf('0.5'), mpf('0.7'), mpf('0.9'), mpf('0.99')):
    r = -alpha * g(g(xt / alpha)) - g(xt)
    max_fp_resid = max(max_fp_resid, fabs(r))
log(f"fixed-point residual (6 точек x=0.1..0.99): max={nstr(max_fp_resid, 6)}")
assert max_fp_resid < mpf('1e-25'), "неподвижная точка не воспроизводится"

# --- Функциональная (не матричная) форма DT, независимая проверка ---
alpha_p = alpha**p

def Phi_eval(v, s):
    r = mpf(0)
    for c in reversed(v):
        r = r * s + c
    return r

def DTPhi_at_t(v, t):
    # (DTPhi)(t) = -alpha*( Phi(W) + g'(U)*Phi(t/alpha^p) ),  U=G(t/alpha^p), W=U^p
    U = G(t / alpha_p)
    W = U**p
    gpU = p * U**(p - 1) * Gp(U**p)
    return -alpha * (Phi_eval(v, W) + gpU * Phi_eval(v, t / alpha_p))

def build_matrices(M):
    ts = [(1 + cos(pi * (2 * i + 1) / (2 * M))) / 2 for i in range(M)]
    xs = [t**(mpf(1) / p) for t in ts]
    A = matrix(M, M)
    B = matrix(M, M)
    for i in range(M):
        x = xs[i]
        u = g(x / alpha)
        w = gp(u)
        for j in range(M):
            A[i, j] = -alpha * (u**(p * j) + w * (x / alpha)**(p * j))
            B[i, j] = ts[i]**j
    return A, B, ts

# --- Проверка формулы DT (матрица) против независимой функциональной формы на j=0..3 ---
Mchk = 12
Achk, Bchk, tschk = build_matrices(Mchk)
max_form_diff = mpf(0)
for j in range(4):
    v = [mpf(0)] * Mchk
    v[j] = mpf(1)
    for i in range(Mchk):
        via_matrix = Achk[i, j]  # (DT t^j)(t_i) по матричной формуле
        via_func = DTPhi_at_t(v, tschk[i])
        max_form_diff = max(max_form_diff, fabs(via_matrix - via_func))
log(f"проверка формулы DT (матрица vs функционал) на j=0..3, M={Mchk}: max diff={nstr(max_form_diff, 6)}")
assert max_form_diff < mpf('1e-30'), "формулы DT расходятся"

# --- Координатная мода K(t) = G(t) - p*t*G'(t): DT(K)=K без проекции ---
Kvec = [gc[j] * (1 - p * j) for j in range(len(gc))]
# нормировка для устойчивого сравнения
Knorm = max(fabs(c) for c in Kvec)
Kvec_n = [c / Knorm for c in Kvec]
dense_ts_K = [mpf(k) / 200 for k in range(1, 200)]
max_K_resid = mpf(0)
for t in dense_ts_K:
    lhs = DTPhi_at_t(Kvec_n, t)
    rhs = Phi_eval(Kvec_n, t)
    denom = max(fabs(rhs), mpf('1e-30'))
    max_K_resid = max(max_K_resid, fabs(lhs - rhs) / denom)
log(f"координатная мода K=G-p*t*G': DT(K)/K max relative residual (199 точек, без проекции) = {nstr(max_K_resid, 6)}")

results = {'d': 6.0, 'alpha': nstr(alpha, 30), 'checks': {
    'fp_residual_max': nstr(max_fp_resid, 10),
    'DT_formula_diff_matrix_vs_functional': nstr(max_form_diff, 10),
    'coord_mode_K_DT_over_K_residual': nstr(max_K_resid, 10),
}}
json.dump(results, open('spectrum_d6_partial.json', 'w', encoding='utf-8'), indent=2)
log("частичный JSON (проверки) сохранён -> spectrum_d6_partial.json")
log("этап 1 (проверки формул) завершён; переходим к eig() на M=32")
LOG.close()

# -*- coding: utf-8 -*-
# Этап 3: сдвинутая обратная итерация (уточнение) для всех кандидатов |λ|>=0.05,
# сравнение M=32/48 (дрейф), сравнение tower1 с alpha^{-6} точно, тест координатной
# моды (cos-подобие собственного вектора tower1 с K=G-6tG'), итоговый spectrum_d6.json.
import json, time, os
from mpmath import mp, mpf, matrix, eig, fabs, nstr, cos, pi, lu_solve

mp.dps = 60
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def log(msg):
    t = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"{t} [refine] {msg}"
    print(line)
    with open('progress.log', 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()

GRID_PATH = r'../sec2_spectra_grid_d/spectra_grid.json'
pts = json.load(open(GRID_PATH, encoding='utf-8'))['points']
pt6 = next(pt for pt in pts if abs(float(pt['d']) - 6.0) < 1e-9)
alpha = mpf(pt6['alpha48'])
gc = [mpf(c) for c in pt6['_a48']]
p = 6
alpha_p = alpha**p

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

def Phi_eval(v, s):
    r = mpf(0)
    for c in reversed(v):
        r = r * s + c
    return r

def DTPhi_at_t(v, t):
    U = G(t / alpha_p)
    W = U**p
    gpU = p * U**(p - 1) * Gp(U**p)
    return -alpha * (Phi_eval(v, W) + gpU * Phi_eval(v, t / alpha_p))

def dense_residual(v, lam, npts=161):
    vnorm = max(fabs(c) for c in v)
    vv = [c / vnorm for c in v]
    maxr = mpf(0)
    maxphi = mpf(0)
    for k in range(1, npts):
        t = mpf(k) / npts
        lhs = DTPhi_at_t(vv, t)
        rhs = lam * Phi_eval(vv, t)
        maxr = max(maxr, fabs(lhs - rhs))
        maxphi = max(maxphi, fabs(Phi_eval(vv, t)))
    return maxr / max(maxphi, mpf('1e-30'))

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
    return A, B

def shifted_inverse_iteration(C, lam0, M, iters=8):
    # C = B^-1 A (уже собранная плотная матрица); (C - lam0*I)^-1 степенная итерация
    Cs = matrix(M, M)
    for i in range(M):
        for j in range(M):
            Cs[i, j] = C[i, j] - (lam0 if i == j else 0)
    v = matrix(M, 1)
    for i in range(M):
        v[i] = mpf(1) / (i + 1)
    lam = lam0
    for it in range(iters):
        w = lu_solve(Cs, v)
        nrm = max(fabs(w[i]) for i in range(M))
        w = w / nrm
        # Rayleigh quotient через C напрямую: lam_new = (w^T C w)/(w^T w)
        Cw = C * w
        num = sum(w[i] * Cw[i] for i in range(M))
        den = sum(w[i] * w[i] for i in range(M))
        lam_new = num / den
        v = w
        lam = lam_new
    return lam, [v[i] for i in range(M)]

M = 48
log(f"строим C=B^-1*A при M={M} для обратной итерации...")
t0 = time.time()
A, B = build_matrices(M)
C = B**-1 * A
log(f"C готова, {time.time()-t0:.1f}s")

seeds = {
    'delta': mpf('9.29624683277137'),
    'coord_K': mpf('1.0'),
    'lambda2_pos': mpf('0.37107146768193'),
    'lambda3_neg': mpf('-0.34196832254780'),
    'lambda4_m3': mpf('-0.21051852081427'),
    'lambda5_new': mpf('0.20784669102300'),
    'lambda6_new': mpf('-0.14488609906251'),
    'tower1_k1': mpf('0.10002322865338'),
    'lambda7_new': mpf('-0.06922351219536'),
}

refined = {}
for name, lam0 in seeds.items():
    lam_ref, vec = shifted_inverse_iteration(C, lam0, M, iters=10)
    dres = dense_residual(vec, lam_ref)
    refined[name] = {'lam': lam_ref, 'vec': vec, 'dense_resid': dres}
    log(f"{name}: refined lam={nstr(lam_ref,25)}  dense_resid={nstr(dres,6)}  (seed={nstr(lam0,10)})")

# --- Сравнение M=32 vs M=48 (дрейф) ---
log("повторяем на M=32 для дрейфа...")
A32, B32 = build_matrices(32)
C32 = B32**-1 * A32
drift = {}
for name, lam0 in seeds.items():
    lam32, vec32 = shifted_inverse_iteration(C32, lam0, 32, iters=10)
    d = fabs(lam32 - refined[name]['lam'])
    drift[name] = {'lam_M32': lam32, 'drift_abs': d}
    log(f"{name}: M32={nstr(lam32,20)}  drift|M48-M32|={nstr(d,6)}")

# --- alpha^-6 точное сравнение с tower1 ---
alpha_m6 = alpha**(-6)
tower_val = refined['tower1_k1']['lam']
match_rel = fabs(tower_val - alpha_m6) / alpha_m6
log(f"alpha^-6 = {nstr(alpha_m6,25)}")
log(f"tower1 refined = {nstr(tower_val,25)}")
log(f"|tower1 - alpha^-6|/alpha^-6 = {nstr(match_rel,6)}")

# --- тест координатной моды: cos-подобие tower1 vs K=G-6tG' (в монономном базисе) ---
Kvec = [gc[j] * (1 - p * j) for j in range(min(len(gc), M))]
while len(Kvec) < M:
    Kvec.append(mpf(0))
def cosine(v1, v2):
    n1 = sum(c * c for c in v1) ** mpf('0.5')
    n2 = sum(c * c for c in v2) ** mpf('0.5')
    dot = sum(v1[i] * v2[i] for i in range(len(v1)))
    return fabs(dot / (n1 * n2))

cos_tower_K = cosine(refined['tower1_k1']['vec'], Kvec)
log(f"тест координатной моды: cos(tower1_vec, K_vec) = {nstr(cos_tower_K,10)} (ожидание: << 1, мода отличается от координатной)")

# для контроля -- cos координатной моды с самой собой (через её же refined-вектор, если он сошёлся к K)
cos_K_self = cosine(refined['coord_K']['vec'], Kvec)
log(f"контроль: cos(refined coord_K vector, K_vec analytic) = {nstr(cos_K_self,10)} (ожидание: ~1, это ТА ЖЕ мода)")

# --- также cos прочих новых веток с K, для полноты ---
for name in ('lambda5_new', 'lambda6_new', 'lambda7_new', 'lambda4_m3'):
    c = cosine(refined[name]['vec'], Kvec)
    log(f"cos({name}_vec, K_vec) = {nstr(c,10)}")

# --- сохранение итогового JSON ---
out = {
    'd': 6.0,
    'alpha6': nstr(alpha, 30),
    'alpha6_pow_minus6': nstr(alpha_m6, 30),
    'M_truncations': [32, 48],
    'branches': []
}
order = ['delta', 'coord_K', 'lambda2_pos', 'lambda3_neg', 'lambda4_m3', 'lambda5_new',
         'lambda6_new', 'tower1_k1', 'lambda7_new']
for name in order:
    r = refined[name]
    dr = drift[name]
    out['branches'].append({
        'name': name,
        'lam_M48_refined': nstr(r['lam'], 25),
        'lam_M48_f': float(r['lam']),
        'dense_resid_M48': nstr(r['dense_resid'], 6),
        'lam_M32_refined': nstr(dr['lam_M32'], 25),
        'drift_M32_M48_abs': nstr(dr['drift_abs'], 6),
    })
out['tower1_vs_alpha_pow_minus6'] = {
    'alpha_pow_minus6': nstr(alpha_m6, 30),
    'tower1_refined': nstr(tower_val, 30),
    'match_rel': nstr(match_rel, 6),
}
out['coordinate_mode_test'] = {
    'cos_tower1_vs_K': nstr(cos_tower_K, 10),
    'cos_coordK_selfcheck_vs_K': nstr(cos_K_self, 10),
    'note': 'K(t)=G(t)-6*t*G_t(t), монономный базис, векторы из обратной итерации на M=48'
}
json.dump(out, open('spectrum_d6.json', 'w', encoding='utf-8'), indent=2)
log("итоговый JSON сохранён -> spectrum_d6.json")
log("refine_and_report.py завершён")

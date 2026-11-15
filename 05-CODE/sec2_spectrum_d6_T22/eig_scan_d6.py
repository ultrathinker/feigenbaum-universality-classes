# -*- coding: utf-8 -*-
# Этап 2: полный eig() на M=32 и M=48, вещественные собственные значения |λ|>=0.05,
# функциональный (дальний от узлов) остаток, дрейф по N, сопоставление с тем, что
# нашёл kilo-2 (spectra_grid.json: delta, pos/lam2/lam3, tower1).
import json, time, os
from mpmath import mp, mpf, matrix, eig, fabs, cos, pi, nstr

mp.dps = 50
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def log(msg):
    t = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"{t} {msg}"
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

def dense_residual(v, lam, npts=161):
    # относительный функциональный остаток на дальней (не узловой) сетке
    vnorm = max(fabs(c) for c in v)
    if vnorm == 0:
        return mpf('nan')
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

def run_scan(M, tag):
    t0 = time.time()
    log(f"[{tag}] строим A,B (M={M})...")
    A, B = build_matrices(M)
    log(f"[{tag}] матрицы готовы ({time.time()-t0:.1f}s), считаем eig(B^-1 A)...")
    E, EV = eig(B**-1 * A)
    log(f"[{tag}] eig() готов ({time.time()-t0:.1f}s), всего {len(E)} собственных значений")
    reals = []
    for k in range(len(E)):
        z = E[k]
        if fabs(z.imag) < mpf('1e-8') * max(fabs(z.real), mpf(1)):
            reals.append((z.real, k))
    reals.sort(key=lambda pr: -abs(pr[0]))
    out = []
    for lam, k in reals:
        if abs(lam) < mpf('0.05'):
            continue
        vec = [EV[i, k] for i in range(M)]
        # приводим к вещественному (мнимая часть подавлена)
        vec_r = [c.real if hasattr(c, 'real') else c for c in vec]
        dres = dense_residual(vec_r, lam)
        out.append({'lam': nstr(lam, 20), 'lam_f': float(lam), 'dense_resid': nstr(dres, 6), 'dense_resid_f': float(dres)})
        log(f"[{tag}] lambda={nstr(lam,15)}  dense_resid={nstr(dres,4)}")
    log(f"[{tag}] сканирование |λ|>=0.05 завершено: {len(out)} вещественных кандидатов, {time.time()-t0:.1f}s всего")
    return out

log("=== старт eig-скана M=32 ===")
out32 = run_scan(32, 'M32')
json.dump(out32, open('eig_scan_M32.json', 'w', encoding='utf-8'), indent=2)
log("M32 сохранён -> eig_scan_M32.json")

log("=== старт eig-скана M=48 ===")
out48 = run_scan(48, 'M48')
json.dump(out48, open('eig_scan_M48.json', 'w', encoding='utf-8'), indent=2)
log("M48 сохранён -> eig_scan_M48.json")

log("этап 2 (сырой eig-скан) завершён")

# -*- coding: utf-8 -*-
# Кросс-проверка через ВТОРОЙ независимый движок (агент M, wave2_vector_homotopy/engine.py):
# другая матричная формулировка (V^{-1}F, с явным gauge-проекционным членом
# -(Phi(1)/G(1))*K(t)), собственные лево+право векторы, dense_eigen_residual.
# Импортируется как read-only модуль, ничего не пишем в чужую папку.
import sys, os, json, time
sys.path.insert(0, r'../sec2_vector_homotopy')
import engine as M  # noqa: E402
from mpmath import mp, mpf, fabs, nstr  # noqa: E402

mp.dps = 50
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def log(msg):
    t = time.strftime('%Y-%m-%d %H:%M:%S')
    line = f"{t} [cross-M] {msg}"
    print(line)
    with open('progress.log', 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()

pts = M.load_wave1_points()
rec = pts[6.0] if 6.0 in pts else next(v for k, v in pts.items() if abs(float(k) - 6.0) < 1e-9)
a48, alpha48 = M.fixedpoint_from_wave1(rec)
log(f"alpha6 (via engine.py loader) = {nstr(alpha48, 20)}")

d = mpf(6)
candidates_named = {
    'delta': mpf('9.29624683277137'),
    'coord_K': mpf('1.0'),
    'lambda2(pos)': mpf('0.37107146768193'),
    'lambda3(neg)': mpf('-0.34196832254780'),
    'm3(lam4?)': mpf('-0.21051852081427'),
    'NEW_0.2078': mpf('0.20784669102300'),
    'NEW_-0.1449': mpf('-0.14488609906251'),
    'tower1(0.10002)': mpf('0.10002322865338'),
    'NEW_-0.0692': mpf('-0.06922351219536'),
}

for N in (32, 48):
    log(f"engine.py full_eig_window(N={N}, window=1.5) ...")
    t0 = time.time()
    out = M.full_eig_window(a48, alpha48, d, N, window=mpf('1.5'), imag_tol=1e-10)
    log(f"  -> {len(out)} вещественных мод в |lam|<=1.5, {time.time()-t0:.1f}s")
    results = []
    for rowd in out:
        lam = rowd['lam']
        if abs(lam) < mpf('0.05'):
            continue
        hcoef = [mpf(0)] + rowd['rvec']  # engine's rvec = coeffs of t^1..t^N; horner0 expects c[0]=const (=0)
        worst, phmax = M.dense_eigen_residual(a48, alpha48, d, hcoef, lam, npts=121)
        rel = worst / max(phmax, mpf('1e-30'))
        # найти ближайший именованный кандидат
        best_name, best_diff = None, None
        for name, val in candidates_named.items():
            diff = abs(lam - val)
            if best_diff is None or diff < best_diff:
                best_name, best_diff = name, diff
        results.append({'lam': nstr(lam, 20), 'dense_resid_rel': nstr(rel, 6), 'closest': best_name, 'diff_to_closest': nstr(best_diff, 6)})
        log(f"  N={N} lam={nstr(lam,15)}  dense_resid_rel={nstr(rel,4)}  closest={best_name} (diff={nstr(best_diff,4)})")
    json.dump(results, open(f'cross_check_engine_N{N}.json', 'w', encoding='utf-8'), indent=2)
    log(f"сохранено -> cross_check_engine_N{N}.json")

log("cross_check_engine.py завершён")

# -*- coding: utf-8 -*-
"""
f4_main.py -- главный счёт F4:
 (A) таблица тождеств d in {2,4,6} x k in {1,2,3,4}: max|DT Phi_k - alpha^{-dk} Phi_k| / max|Phi_k|;
 (B) косинусы Phi_k с численными собственными векторами (коэфф. пространство и сетка t);
 (C) формальный кандидат при d=2.01, 2.05, 3.0: остаток, локализация, косинус с B+;
 (D) ПРОВЕРКА ПРЕДСКАЗАНИЯ (записано в progress.log 19:31): при d=3,5 нет alpha^{-dk} в спектре.
[oid:f4GLM] NOT RIGOROUS -- mpmath, no interval certificates.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf

import f4_common as C
from f4_test3_dictfit import build_DT, inv_iter

mp.dps = 50
RES = {}
T0 = time.time()


def spow(y, p):
    """Вещественная ветка y^p: sign(y)|y|^p (конвенция кодов для |U|^d)."""
    if y == 0:
        return mpf(0)
    return (mpf(1) if y > 0 else -1) * abs(y) ** p


def phi_factory(a, d, k):
    """Phi_k = d t^{k+1} G' - G^{dk+1} + G - d t G'."""
    def Phi(t):
        G_ = C.horner(a, t)
        Gp_ = C.Gp(a, t)
        return (d * t ** (k + 1) * Gp_ - spow(G_, d * k + 1) + G_
                - d * t * Gp_)
    return Phi


def phi_series(a, d, k, nmax):
    """Коэффициенты t-ряда Phi_k (целые d): через произведение рядов."""
    Gs = a[:nmax + 1] + [mpf(0)] * max(0, nmax + 1 - len(a))
    # G^{dk+1} свёртками
    q = int(d * k + 1)
    P = [mpf(1)] + [mpf(0)] * nmax
    for _ in range(q):
        Pn = [mpf(0)] * (nmax + 1)
        for i in range(nmax + 1):
            if P[i] == 0:
                continue
            for j in range(0, nmax + 1 - i):
                if Gs[j] != 0:
                    Pn[i + j] += P[i] * Gs[j]
        P = Pn
    # d t^{k+1} G' - d t G' + G - G^{dk+1}
    out = [mpf(0)] * (nmax + 1)
    for j in range(1, min(nmax, len(a) - 1) + 1):
        if j + k + 1 <= nmax:
            out[j + k + 1] += d * j * a[j]
        if j + 1 <= nmax:
            out[j + 1] -= d * j * a[j]
    for j in range(0, min(nmax, len(Gs) - 1) + 1):
        out[j] += Gs[j] - P[j]
    return out


def cos_grid(f, g, nodes):
    num = sum(f(t) * g(t) for t in nodes)
    nf = mp.sqrt(sum(f(t) ** 2 for t in nodes))
    ng = mp.sqrt(sum(g(t) ** 2 for t in nodes))
    return abs(num / (nf * ng))


def cos_coef(u, v):
    num = sum(u[i] * v[i] for i in range(min(len(u), len(v))))
    nu = mp.sqrt(sum(x * x for x in u))
    nv = mp.sqrt(sum(x * x for x in v))
    return abs(num / (nu * nv))


# ============================== (A) тождество ==============================
C.plog('MAIN(A) START: таблица тождеств d=2,4,6 x k=1..4 (200 точек, dps=%d)' % mp.dps)
RES['identity'] = {}
for tag in ('2.0', '4.0', '6.0'):
    a, alpha, d, meta = C.load_fp(tag)
    RES['identity'][tag] = {'src': meta['src'], 'alpha': mp.nstr(alpha, 30), 'rows': []}
    for k in (1, 2, 3, 4):
        Phi = phi_factory(a, d, k)
        mu = alpha ** (-d * k)
        rr, wt, pm = C.resid_table(a, alpha, d, Phi, mu, npts=200)
        RES['identity'][tag]['rows'].append(
            {'k': k, 'mu': mp.nstr(mu, 25), 'rel_resid': float(rr),
             'worst_t': mp.nstr(wt, 6)})
        C.plog('  d=%s k=%d: alpha^-%d=%s  max|r|/max|Phi| = %.3e'
               % (tag, k, d * k, mp.nstr(mu, 16), float(rr)))

# ============================== (B) косинусы ===============================
C.plog('MAIN(B) START: косинусы Phi_k с численными eigenvector')
WIN = json.load(open(os.path.join(C.ROOT, 'wave2_vector_homotopy',
                                  'spectrum_window.json'), encoding='utf-8'))
RES['cosines'] = {}
nodesB = C.grid_t(60, mpf('1.0'))

for tag, k, idx in (('2.0', 1, 0), ('2.0', 2, 3), ('4.0', 1, 3), ('4.0', 2, 7)):
    a, alpha, d, meta = C.load_fp(tag)
    p = WIN['points'][tag]
    lam = mpf(p['lams'][idx])
    v = [mpf(s) for s in p['rvecs'][idx]]
    mu = alpha ** (-d * k)
    Phi = phi_factory(a, d, k)
    Phis = phi_series(a, d, k, len(v))
    cvals = cos_grid(Phi, lambda t: C.horner([mpf(0)] + v, t), nodesB)
    ccoef = cos_coef(Phis[1:], v)
    # косинус с K для контекста
    K = C.K_fun(a, d)
    cK = cos_grid(Phi, K, nodesB)
    RES['cosines']['%s_k%d' % (tag, k)] = {
        'lam_window': mp.nstr(lam, 20), 'mu_alpha': mp.nstr(mu, 20),
        'rel_lam_mu': float(abs(lam - mu) / mu),
        'cos_values_grid': float(cvals), 'cos_coef_space': float(ccoef),
        'cos_with_K': float(cK)}
    C.plog('  d=%s k=%d: lam=%s vs mu=%s (rel %.1e): cos(сетка)=%.12f  cos(коэф.)=%.12f  cos(K)=%.4f'
           % (tag, k, mp.nstr(lam, 12), mp.nstr(mu, 12),
              float(abs(lam - mu) / mu), float(cvals), float(ccoef), float(cK)))

# d=6: свой eigenvector (inverse iteration), window не имеет d=6
a6, alpha6, d6, meta6 = C.load_fp('6.0')
A40_6 = build_DT(a6, alpha6, d6, 40)
lam6, v6, mres6 = inv_iter(A40_6, alpha6 ** (-6))
mx6 = max(abs(x) for x in v6)
v6n = [x / mx6 for x in v6]
Phis6 = phi_series(a6, d6, 1, 40)
Phi6 = phi_factory(a6, d6, 1)
nodesB6 = C.grid_t(60, mpf('1.0'))
cvals6 = cos_grid(Phi6, lambda t: C.horner([mpf(0)] + v6n, t), nodesB6)
ccoef6 = cos_coef(Phis6[1:], v6n)
K6 = C.K_fun(a6, d6)
RES['cosines']['6.0_k1'] = {'lam_inviter': mp.nstr(lam6, 20),
                            'mu_alpha': mp.nstr(alpha6 ** (-6), 20),
                            'rel_lam_mu': float(abs(lam6 - alpha6 ** (-6)) / alpha6 ** (-6)),
                            'cos_values_grid': float(cvals6),
                            'cos_coef_space': float(ccoef6),
                            'cos_with_K': float(cos_grid(Phi6, K6, nodesB6))}
C.plog('  d=6.0 k=1 (свой eig): lam=%s vs mu (rel %.1e): cos(сетка)=%.12f cos(коэф.)=%.12f cos(K)=%.4f'
       % (mp.nstr(lam6, 12), float(abs(lam6 - alpha6 ** (-6)) / alpha6 ** (-6)),
          float(cvals6), float(ccoef6), float(cos_grid(Phi6, K6, nodesB6))))

# ============== (C) дробные/нечётные d: формальный кандидат ================
C.plog('MAIN(C) START: формальный кандидат при d=2.01, 2.05, 3.0')
SG = json.load(open(os.path.join(C.ROOT, 'wave1_spectra_d', 'spectra_grid.json'),
                    encoding='utf-8'))
RES['non_even'] = {}
for tagd in ('2.01', '2.05', '3.0'):
    rec = next(r for r in SG['points'] if r['d'] == tagd)
    a = [mpf(s) for s in rec['_a48']]
    alpha = mpf(rec['alpha48'])
    d = mpf(tagd)
    row = {'alpha': mp.nstr(alpha, 25), 'rows': []}
    for k in (1, 2):
        Phi = phi_factory(a, d, k)
        mu = alpha ** (-d * k)
        # остаток и его локализация
        rmax = mpf(0)
        pmax = mpf(0)
        prof = []
        for t in C.grid_t(200, mpf('1.2')):
            rv = abs(C.DT_t(a, alpha, d, Phi, t) - mu * Phi(t))
            prof.append((float(t), float(rv)))
            rmax = max(rmax, rv)
            pmax = max(pmax, abs(Phi(t)))
        half = [tt for tt, rv in prof if rv > 0.5 * rmax]
        tzeroG = None
        # где G(tau)=0: граница ветки |U|^d
        lo, hi = mpf(0), mpf('1.2')
        for _ in range(60):
            mid = (lo + hi) / 2
            if C.horner(a, mid / alpha ** d) > 0:
                lo = mid
            else:
                hi = mid
        tzeroG = mp.nstr(hi, 8)
        row['rows'].append({
            'k': k, 'mu_formal': mp.nstr(mu, 16),
            'rel_resid': float(rmax / pmax),
            'worst_t': '%.4f' % max(prof, key=lambda z: z[1])[0],
            't_half_band': ['%.4f' % min(half), '%.4f' % max(half)],
            't_where_U_zero': tzeroG})
        C.plog('  d=%s k=%d: max|r|/max|Phi|=%.3e (худший t=%s, полуполоса %s, U=0 при t=%s)'
               % (tagd, k, float(rmax / pmax), '%.4f' % max(prof, key=lambda z: z[1])[0],
                  ['%.4f' % min(half), '%.4f' % max(half)], tzeroG))
    # косинус формального кандидата с ветвью B+ (rvecs[0]) при 2.01/2.05
    if tagd in WIN['points']:
        p = WIN['points'][tagd]
        vB = [mpf(s) for s in p['rvecs'][0]]
        lamB = mpf(p['lams'][0])
        Phi1 = phi_factory(a, d, 1)
        cB = cos_grid(Phi1, lambda t: C.horner([mpf(0)] + vB, t), nodesB)
        row['Bplus'] = {'lam': mp.nstr(lamB, 16),
                        'alpha_pow': mp.nstr(alpha ** (-d), 16),
                        'rel_diff': float(abs(lamB - alpha ** (-d)) / alpha ** (-d)),
                        'cos_formalPhi1_Bplus': float(cB)}
        C.plog('  d=%s: B+=%s vs alpha^{-d}=%s (rel %.2e), cos(формальный Phi_1, B+)=%.6f'
               % (tagd, mp.nstr(lamB, 10), mp.nstr(alpha ** (-d), 10),
                  float(abs(lamB - alpha ** (-d)) / alpha ** (-d)), float(cB)))
    RES['non_even'][tagd] = row

# ============== (D) проверка предсказания: d=3, 5 ==========================
C.plog('MAIN(D) START: предсказание -- при d=3.0 и 5.0 нет alpha^{-dk} в спектре DT')
from mpmath import eig as mp_eig, matrix

RES['prediction'] = {'stated_before_check':
                     '2026-08-15 19:31:00 progress.log: при нечётных d=3.0,5.0 '
                     'нет собственных значений alpha^{-dk} (k=1,2,3) с rel<1e-6',
                     'checks': {}}
for tagd in ('3.0', '5.0'):
    rec = next(r for r in SG['points'] if r['d'] == tagd)
    a = [mpf(s) for s in rec['_a48']]
    alpha = mpf(rec['alpha48'])
    d = mpf(tagd)
    towers = {k: alpha ** (-d * k) for k in (1, 2, 3)}
    chk = {'alpha': mp.nstr(alpha, 25), 'towers': {str(k): mp.nstr(v, 16)
                                                   for k, v in towers.items()}}
    for N in (20, 28):
        A = build_DT(a, alpha, d, N)
        E = mp_eig(matrix(A), right=False)
        E = [z.real if abs(z.imag) < abs(z) * mpf('1e-12') else z for z in E]
        E = [z for z in E if abs(z) > mpf('0.002')]
        best = {}
        for k, tw in towers.items():
            cand = min(E, key=lambda z: abs(z - tw))
            best['k%d' % k] = {'nearest': mp.nstr(cand, 14),
                               'rel_miss': float(abs(cand - tw) / tw)}
        chk['N%d' % N] = {'n_eigs': len(E),
                          'top': [mp.nstr(z, 10) for z in sorted(E, key=lambda z: -abs(z))[:10]],
                          'best_tower_match': best}
        C.plog('  d=%s N=%d: %d соб.знач.; ближайшие к башне: %s'
               % (tagd, N, len(E),
                  {kk: ('%.3e' % vv['rel_miss']) for kk, vv in best.items()}))
    RES['prediction']['checks'][tagd] = chk

RES['meta'] = {'tag': '[oid:f4GLM] NOT RIGOROUS', 'dps': mp.dps,
               'elapsed_s': round(time.time() - T0, 1)}
with open(os.path.join(C.HERE, 'f4_results.json'), 'w', encoding='utf-8') as f:
    json.dump(RES, f, indent=1, ensure_ascii=False)
C.plog('MAIN DONE (%.0fs) -> f4_results.json' % (time.time() - T0))

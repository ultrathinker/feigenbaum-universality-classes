# -*- coding: utf-8 -*-
"""
f4_test3_dictfit.py -- свой eigenvector башенной моды (d=2, dps 60, обратная
итерация на коллокационной матрице N=40) + словарный фит явными формулами.
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, cos, pi

import f4_common as C

mp.dps = 60

# ---- LU (Doolittle, partial pivoting) -- из spectra_grid.py ---------------
def lu_decomp(A):
    n = len(A)
    LU = [row[:] for row in A]
    piv = list(range(n))
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(LU[i][k]))
        if p != k:
            LU[k], LU[p] = LU[p], LU[k]
            piv[k], piv[p] = piv[p], piv[k]
        pv = LU[k][k]
        rk = LU[k]
        for i in range(k + 1, n):
            m = LU[i][k] / pv
            LU[i][k] = m
            if m != 0:
                ri = LU[i]
                for j in range(k + 1, n):
                    ri[j] -= m * rk[j]
    return LU, piv


def lu_solve(LU, piv, b):
    n = len(b)
    x = [b[piv[i]] for i in range(n)]
    for i in range(1, n):
        s = x[i]
        ri = LU[i]
        for j in range(i):
            s -= ri[j] * x[j]
        x[i] = s
    for i in range(n - 1, -1, -1):
        s = x[i]
        ri = LU[i]
        for j in range(i + 1, n):
            s -= ri[j] * x[j]
        x[i] = s / ri[i]
    return x


def colloc(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def build_DT(a, alpha, d, N):
    pts = colloc(N)
    ad = alpha ** d
    G1 = C.horner(a, mpf(1))
    V, F = [], []
    for t in pts:
        ta = t / ad
        U = C.horner(a, ta)
        aU = abs(U)
        sU = 1 if U >= 0 else -1
        W = aU ** d
        Kt = C.horner(a, t) - d * t * C.Gp(a, t)
        gauge = -(1 / G1) * Kt
        c3 = -alpha * d * sU * aU ** (d - 1) * C.Gp(a, W)
        vr = [mpf(0)] * N
        fr = [mpf(0)] * N
        pt = pw = pta = mpf(1)
        for j in range(1, N + 1):
            pt *= t
            pw *= W
            pta *= ta
            vr[j - 1] = pt
            fr[j - 1] = gauge - alpha * pw + c3 * pta
        V.append(vr)
        F.append(fr)
    LU, piv = lu_decomp(V)
    cols = [lu_solve(LU, piv, [F[i][j] for i in range(N)]) for j in range(N)]
    return [[cols[j][i] for j in range(N)] for i in range(N)]


def inv_iter(A, sigma, iters=80):
    n = len(A)
    M = [[A[i][j] - (sigma if i == j else 0) for j in range(n)] for i in range(n)]
    LU, piv = lu_decomp(M)
    v = [mpf(1) / (i + 2) for i in range(n)]
    for it in range(iters):
        w = lu_solve(LU, piv, v)
        nrm = max(abs(x) for x in w)
        v = [x / nrm for x in w]
    Av = [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]
    kk = max(range(n), key=lambda i: abs(v[i]))
    lam = Av[kk] / v[kk]
    resid = max(abs(Av[i] - lam * v[i]) for i in range(n))
    return lam, v, resid


def lstsq(basis_f, target_f, nodes):
    """МНК: target ~ sum c_i basis_i. Возвращает (coefs, maxresid)."""
    m = len(nodes)
    nb = len(basis_f)
    B = [[f(t) for f in basis_f] for t in nodes]
    y = [target_f(t) for t in nodes]
    AtA = [[sum(B[r][i] * B[r][j] for r in range(m)) for j in range(nb)] for i in range(nb)]
    Atb = [sum(B[r][i] * y[r] for r in range(m)) for i in range(nb)]
    n = nb
    M = [AtA[i][:] + [Atb[i]] for i in range(n)]
    for col in range(n):
        p = max(range(col, n), key=lambda r: abs(M[r][col]))
        M[col], M[p] = M[p], M[col]
        for r in range(n):
            if r != col and M[col][col] != 0:
                fac = M[r][col] / M[col][col]
                for cc in range(col, n + 1):
                    M[r][cc] -= fac * M[col][cc]
    co = [M[i][n] / M[i][i] for i in range(n)]
    resid = max(abs(sum(co[i] * B[r][i] for i in range(nb)) - y[r]) for r in range(m))
    return co, resid


C.plog('TEST3 START: свой eigenvector d=2 k=1 + словарный фит (dps=%d)' % mp.dps)

a, alpha, d, meta = C.load_fp('2.0')
N = 40
A40 = build_DT(a, alpha, d, N)
C.plog('матрица DT N=40 построена')

mu = alpha ** (-2)
lam, v, mres = inv_iter(A40, mu)
C.plog('inverse iter: lam=%s (mu=%s, rel.diff %.2e), mat_resid %.2e'
       % (mp.nstr(lam, 20), mp.nstr(mu, 20), float(abs(lam - mu) / mu), float(mres)))

# собственный вектор как функция (нормировка: c1 = 1)
ss = sum(v)
coef = [mpf(0)] + [x / ss for x in v]      # сумма коэффициентов = 1 -> Phi(1)=1


def Phi(t):
    return C.horner(coef, t)


G = lambda t: C.horner(a, t)
Gp_ = lambda t: C.Gp(a, t)
K = lambda t: C.horner(a, t) - d * t * C.Gp(a, t)

# плотная проверка собственного тождества
rr, wt, pm = C.resid_table(a, alpha, d, Phi, lam, npts=200)
C.plog('плотный остаток своего eigenvector: %.3e (pmax=%s)' % (float(rr), mp.nstr(pm, 6)))
C.plog('коэффициенты c_j (норм. sum=1): %s' % [mp.nstr(x, 6) for x in coef[1:16]])

nodes = [mpf('1.2') * (1 - cos(pi * (i + mpf('0.5')) / 60)) / 2 for i in range(60)]

# словарь (k=1)
Gt = G
dicts = {
    'G,Gp,K,tG,tGp,tK,t2K,t2G,t2Gp,G2,G3,K2,GK': None,
}

basis1 = [G, Gp_, K,
          lambda t: t * G(t), lambda t: t * Gp_(t), lambda t: t * K(t),
          lambda t: t ** 2 * K(t), lambda t: t ** 2 * G(t), lambda t: t ** 2 * Gp_(t),
          lambda t: G(t) ** 2, lambda t: G(t) ** 3,
          lambda t: K(t) ** 2, lambda t: G(t) * K(t)]
co, resid = lstsq(basis1, Phi, nodes)
C.plog('фит deg-словарь: maxresid=%.3e  коэф=%s'
       % (float(resid), [mp.nstr(x, 5) for x in co]))

# малые словари по 3-4 функции (гипотеза смешанного семейства)
small = {
    'A*G+(B+t)*tGp': [G, lambda t: t * Gp_(t), lambda t: t ** 2 * Gp_(t)],
    'A*G+B*tGp+C*tK': [G, lambda t: t * Gp_(t), lambda t: t * K(t)],
    'A*G+B*K+C*t2Gp': [G, K, lambda t: t ** 2 * Gp_(t)],
    'A*G3+B*G+C*t2Gp': [lambda t: G(t) ** 3, G, lambda t: t ** 2 * Gp_(t)],
}
for name, bs in small.items():
    co, resid = lstsq(bs, Phi, nodes)
    C.plog('  %s: maxresid=%.3e co=%s' % (name, float(resid), [mp.nstr(x, 6) for x in co]))

C.plog('TEST3 DONE')

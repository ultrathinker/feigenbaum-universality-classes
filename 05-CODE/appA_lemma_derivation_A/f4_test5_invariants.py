# -*- coding: utf-8 -*-
"""
f4_test5_invariants.py -- ищем инвариантное подперство DT, натянутое на явные
формулы (t-пространство, d=2, k=1,2). Базисы строятся из S- и J-собственных
форм: b1 = d t^{k+1} G' (S-диагональ), b2 = G^{dk+1}-1 (J-диагональ),
b3 = K-1, b4 = t^k(K-1) и т.п. Считаем матрицу сужения МНК и её спектр.
[oid:f4GLM] NOT RIGOROUS.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp, mpf, cos, pi

import f4_common as C

mp.dps = 50


def lstsq_proj(basis, F_nodes, nodes):
    """коэффициенты представления F_nodes в span(basis) + макс.невязка."""
    m = len(nodes)
    nb = len(basis)
    B = [[f(t) for f in basis] for t in nodes]
    y = F_nodes
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


def mp_eig_2x2(M22):
    (a, b), (c, dd) = M22
    tr = a + dd
    det = a * dd - b * c
    disc = tr * tr / 4 - det
    if disc < 0:
        return None
    sq = mp.sqrt(disc)
    return (tr / 2 + sq, tr / 2 - sq)


C.plog('TEST5 START: инвариантные подпространства DT (d=2), dps=%d' % mp.dps)

a, alpha, d, meta = C.load_fp('2.0')
G = lambda t: C.horner(a, t)
Gp = lambda t: C.Gp(a, t)
K = lambda t: G(t) - d * t * Gp(t)
tau = lambda t: t / alpha ** d

nodes = [mpf('1.05') * (1 - cos(pi * (i + mpf('0.5')) / 40)) / 2 for i in range(40)]

for k in (1, 2):
    mu = alpha ** (-d * k)
    C.plog('--- k=%d: цель mu=alpha^-%d=%s' % (k, d * k, mp.nstr(mu, 16)))
    bases = {
        'A=(dt^{k+1}Gp, G^{dk+1}-1)': [
            (lambda t: d * t ** (k + 1) * Gp(t)), (lambda t: G(t) ** (d * k + 1) - 1)],
        'A+(K-1)': [
            (lambda t: d * t ** (k + 1) * Gp(t)), (lambda t: G(t) ** (d * k + 1) - 1),
            (lambda t: K(t) - 1)],
        'A+(K-1,t^{k}(K-1))': [
            (lambda t: d * t ** (k + 1) * Gp(t)), (lambda t: G(t) ** (d * k + 1) - 1),
            (lambda t: K(t) - 1), (lambda t: t ** k * (K(t) - 1))],
        'A+(K-1,t^k(K-1),tGp)': [
            (lambda t: d * t ** (k + 1) * Gp(t)), (lambda t: G(t) ** (d * k + 1) - 1),
            (lambda t: K(t) - 1), (lambda t: t ** k * (K(t) - 1)),
            (lambda t: t * Gp(t))],
    }
    for bname, basis in bases.items():
        nb = len(basis)
        cols = []
        worst_resid = mpf(0)
        Bmat = [[mpf(0)] * nb for _ in range(nb)]
        for jcol in range(nb):
            Fnodes = [C.DT_t(a, alpha, d, basis[jcol], t) for t in nodes]
            co, resid = lstsq_proj(basis, Fnodes, nodes)
            worst_resid = max(worst_resid, resid / max(abs(x) for x in Fnodes))
            for i in range(nb):
                Bmat[i][jcol] = co[i]
        C.plog('  базис %s: относ.невязка инвариантности %.2e'
               % (bname, float(worst_resid)))
        if nb == 2:
            ev = mp_eig_2x2(Bmat)
            if ev:
                C.plog('    соб.значения сужения: %s , %s (цель %s)'
                       % (mp.nstr(ev[0], 14), mp.nstr(ev[1], 14), mp.nstr(mu, 14)))
        else:
            # след и детерминант малой матрицы для ориентира
            tr = sum(Bmat[i][i] for i in range(nb))
            C.plog('    след матрицы сужения: %s' % mp.nstr(tr, 14))

C.plog('TEST5 DONE')

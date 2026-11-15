#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
search_alpha_candidates.py -- [oid:F4-M3]

Try candidates that explicitly involve alpha:
- G(t/alpha^d) - G(t)
- d·t^k·G'(t/alpha^d)
- d·t^k·(G(t) - G(t/alpha^d))
- t^k · (G(t) - alpha · G(t/alpha^d))   -- vanishes by FP eq if t replaced by W
"""
import json
import os
from mpmath import mp, mpf, pi, cos, fsum

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))


def G_eval(a, t):
    r = a[-1]
    for c in reversed(a[:-1]):
        r = r * t + c
    return r


def Gp_eval(a, t):
    s = mpf(0)
    for j in range(1, len(a)):
        s += j * a[j] * t ** (j - 1)
    return s


def colloc_points(npts):
    return [(1 - cos(pi * (i + mpf('0.5')) / npts)) / 2 for i in range(npts)]


def DT_scalar_t(a, alpha, d, c, t):
    ad = alpha ** d
    G1 = G_eval(a, mpf(1))
    ta = t / ad
    U = G_eval(a, ta)
    aU = abs(U)
    sU = 1 if U >= 0 else -1
    W = aU ** d
    Kt = G_eval(a, t) - d * t * Gp_eval(a, t)

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    Phi1 = sum(c)
    return (-(Phi1 / G1) * Kt - alpha * horner0(c, W)
            - alpha * d * sU * aU ** (d - 1) * Gp_eval(a, W) * horner0(c, ta))


def load_fp_d2():
    with open(os.path.join(HERE, '..', '..', 'wave1_spectra_d', 'spectra_grid.json'), encoding='utf-8') as f:
        sg = json.load(f)
    for r in sg['points']:
        if r['d'] == '2.0':
            return [mpf(s) for s in r['_a48']], mpf(r['alpha48'])


def main():
    a, alpha = load_fp_d2()
    d = mpf(2)
    N = 48

    def horner0(cc, tt):
        r = cc[-1]
        for x in reversed(cc[:-1]):
            r = r * tt + x
        return r

    # Get numerical rvec for tower1
    with open(os.path.join(HERE, '..', '..', 'wave2_vector_homotopy', 'spectrum_window.json'), encoding='utf-8') as f:
        sw = json.load(f)
    rec = sw['points']['2.0']
    tv = [mpf(s) for s in rec['rvecs'][0]]
    print(f'Numerical tower1 rvec (lam=alpha^-2={float(rec["lams"][0])}), first 8:')
    print(f'  {[float(tv[i]) for i in range(8)]}')
    print(f'normalized to tv[1]=1:')
    rn = [tv[i]/tv[1] for i in range(len(tv))]
    print(f'  {[float(rn[i]) for i in range(8)]}')

    def test(candidate, name, target):
        c = [mpf(0)] * N
        for j in range(1, N):
            c[j] = candidate(j)
        tpts = colloc_points(201)
        maxr = mpf(0); maxp = mpf(0)
        for t in tpts:
            lhs = DT_scalar_t(a, alpha, d, c, t)
            Phi = horner0(c, t)
            rhs = target * Phi
            ar = abs(lhs - rhs)
            if ar > maxr: maxr = ar
            if abs(Phi) > maxp: maxp = abs(Phi)
        rel = maxr / maxp if maxp > 0 else mpf(0)
        # cosine
        L = min(len(c), len(tv))
        n1 = fsum(x*x for x in c[:L]) ** mpf('0.5')
        n2 = fsum(x*x for x in tv[:L]) ** mpf('0.5')
        cos = abs(fsum(c[i]*tv[i] for i in range(L)) / (n1*n2)) if n1 > 0 and n2 > 0 else mpf(0)
        print(f'  {name:50s} | rel={float(rel):.4e} | cos={float(cos):.4e}')

    # G(t/alpha^d) - G(t) candidates
    target = alpha ** -2
    print(f'\nTarget = alpha^-2 = {float(target):.6f}')
    print('\nCandidates involving alpha:')

    # G(t/alpha^d) - G(t)
    test(lambda j: a[j] * (alpha ** (-d*j) - mpf(1)), 'G(t/alpha^d) - G(t)', target)

    # t * (G(t/alpha^d) - G(t))
    test(lambda j: a[j - 1] * (alpha ** (-d*(j-1)) - mpf(1)) if j - 1 >= 0 else mpf(0), 't * (G(t/alpha^d) - G(t))', target)

    # G(t/alpha^d) - alpha^d * G(t) -- using FP equation: G(t) + alpha*G(W) = 0 with W = G(t/alpha^d)^d
    # If t = W, then G(W) + alpha*G(G(W)/alpha^d) = 0, etc. Not obvious.
    test(lambda j: a[j] * alpha ** (-d*j) - a[j] * alpha**(-d*j + 1), 'G(t/alpha^d) - alpha*G(t/alpha^d)', target)

    # d * t^k * G'(t/alpha^d) -- G' at the rescaled argument
    test(lambda j: d * a[j - 1] * (j - 1) * alpha ** (-d*(j - 1)) if j - 1 >= 1 else mpf(0),
         'd * t * Gprime(t/alpha^d)', target)

    # d * t^k * (G(t) - G(t/alpha^d))
    test(lambda j: d * (a[j] - a[j] * alpha**(-d*j)) if j - 1 >= 0 else d * a[j],
         'd * t * (G(t) - G(t/alpha^d))', target)

    # alpha^{-d} * G(t) - G(t/alpha^d) -- this might relate to FP eq
    test(lambda j: a[j] * alpha**(-d*j) * alpha**(-d) - a[j] * alpha**(-d*j),
         'alpha^{-d} * G(t/alpha^d) - G(t/alpha^d)', target)

    # Combination: G(t) - alpha^d * G(t/alpha^d)
    test(lambda j: a[j] - a[j] * alpha**(-d*j) * alpha**d,
         'G(t) - alpha^d * G(t/alpha^d)', target)

    # Just the gauge mode K(t) = G - d t G'
    test(lambda j: a[j] - d * j * a[j], 'K(t) [gauge mode]', target)

    # t^k * (1 - alpha^{-d})^k -- trivial
    test(lambda j: mpf(1) if j == 1 else mpf(0), 't (monomial)', target)

    # Try: alpha^d * G(t) - G(alpha^d * t)
    test(lambda j: alpha**d * a[j] - a[j] * alpha**(d*j),
         'alpha^d * G(t) - G(alpha^d * t)', target)

    # Try: d*t*G'(alpha^d * t)
    test(lambda j: d * j * a[j] * alpha**(d*(j - 1)) if j >= 1 else mpf(0),
         'd*t*Gprime(alpha^d * t)', target)


if __name__ == '__main__':
    main()

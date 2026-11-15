import time
from mpmath import mp, mpf, nstr
mp.dps = 50
import core

d = mpf(2)
for N in [12, 20, 30]:
    a, al, it, res = core.solve_fixedpoint(N, d)
    A = core.build_DT_matrix(a, al, d, N)
    lam1, v1, _ = core.power_iter(A, N)
    D1 = core.deflate(A, lam1, v1, N)
    lam2, v2, _ = core.power_iter(D1, N)
    D2 = core.deflate(D1, lam2, v2, N)
    lam3, v3, _ = core.power_iter(D2, N)
    D3 = core.deflate(D2, lam3, v3, N)
    lam4, v4, _ = core.power_iter(D3, N)

    alpha_neg_d = al ** (-d)
    print("N=%d alpha=%s" % (N, nstr(al, 15)))
    print("  delta      =", nstr(lam1, 20))
    print("  eig#2      =", nstr(lam2, 20), " alpha^-d=", nstr(alpha_neg_d, 20), " diff=", nstr(lam2-alpha_neg_d,5))
    print("  eig#3      =", nstr(lam3, 20))
    print("  eig#4      =", nstr(lam4, 20))

    # coordinate mode K under PROJECTED operator: check if it's an eigenvector
    Kcoef = [a[0]] + [a[j]*(1 - d*j) for j in range(1, N+1)]
    pts = core.colloc_points(7)
    ratios = []
    for t in pts:
        num = core.DT_scalar_t(a, al, d, Kcoef, t)
        den = sum(Kcoef[j]*t**j for j in range(0, N+1))
        ratios.append(num/den)
    print("  K-mode under PROJECTED DT, ratio sample:", [nstr(r,10) for r in ratios])
    print()

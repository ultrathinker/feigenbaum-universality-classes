import time
from mpmath import mp, mpf, nstr
mp.dps = 40
import core

t0 = time.time()
d = mpf(2)
N = 20
a, al, it, res = core.solve_fixedpoint(N, d)
print("alpha=%s res=%.2e" % (nstr(al, 20), float(res)))

A = core.build_DT_matrix(a, al, d, N)
lam1, v1, it1 = core.power_iter(A, N)
print("delta (raw power iter) =", nstr(lam1, 20), "iters", it1)

D1 = core.deflate(A, lam1, v1, N)
lam2, v2, it2 = core.power_iter(D1, N)
print("lambda2 =", nstr(lam2, 20), "iters", it2)

D2 = core.deflate(D1, lam2, v2, N)
lam3, v3, it3 = core.power_iter(D2, N)
print("lambda3 =", nstr(lam3, 20), "iters", it3)

print("elapsed", time.time() - t0)
print("known anchors: delta=4.66920160910299 lambda2=-0.123652712552687 lambda3=-0.057307021067")

# coordinate mode test
Kcoef = [a[0]] + [a[j]*(1 - d*j) for j in range(1, N+1)]
pts = core.colloc_points(11)
ratios = []
for t in pts:
    num = core.DT_unprojected_t(a, al, d, Kcoef, t)
    den = sum(Kcoef[j]*t**j for j in range(0, N+1))
    ratios.append(num/den)
print("coord mode ratio sample:", [nstr(r, 12) for r in ratios[:5]])

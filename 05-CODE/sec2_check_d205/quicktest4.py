import time
from mpmath import mp, mpf, nstr
mp.dps = 50
import core

d = mpf(2)
N = 20
a, al, it, res = core.solve_fixedpoint(N, d)

# M(t) = d * t * G'(t)  (generator of x -> (1+eps) x rescaling), c_0 = 0 automatically
Mcoef = [mpf(0)] + [d * j * a[j] for j in range(1, N + 1)]
pts = core.colloc_points(9)
ratios = []
for t in pts:
    num = core.DT_scalar_t(a, al, d, Mcoef, t)
    den = sum(Mcoef[j] * t ** j for j in range(1, N + 1))
    ratios.append(num / den)
print("M-mode (x-rescaling generator) under PROJECTED DT, ratios:")
for r in ratios:
    print("  ", nstr(r, 15))

print("alpha^-d =", nstr(al ** (-d), 15))
print("alpha^+d =", nstr(al ** d, 15))

# also check the eigenvector v2 (associated to alpha^-d) found via deflation -- does it
# align with Mcoef?
A = core.build_DT_matrix(a, al, d, N)
lam1, v1, _ = core.power_iter(A, N)
D1 = core.deflate(A, lam1, v1, N)
lam2, v2, _ = core.power_iter(D1, N)
print("lam2 =", nstr(lam2, 20))
# normalize both to compare direction
def normalize(v, n):
    s = max(abs(v[i]) for i in range(n))
    return [v[i]/s for i in range(n)]
v2n = normalize(v2, N)
Mn = normalize(Mcoef[1:], N)
diffs = [abs(v2n[i]-Mn[i]) for i in range(N)]
print("max |v2_normalized - M_normalized| (component-wise, up to N) =", nstr(max(diffs), 8))
print("v2 first 5:", [nstr(x,8) for x in v2n[:5]])
print("M  first 5:", [nstr(x,8) for x in Mn[:5]])

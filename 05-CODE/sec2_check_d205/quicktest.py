import time
from mpmath import mp, mpf, nstr
mp.dps = 40
import core

t0 = time.time()
ladder = [4, 8, 12, 16, 20]
a_cur, al_cur = None, None
d = mpf(2)
for N in ladder:
    a, al, it, res = core.solve_fixedpoint(N, d, a_cur, al_cur)
    print("N=%2d alpha=%s iters=%d res=%.2e" % (N, nstr(al, 20), it, float(res)))
    a_cur, al_cur = a, al

print("elapsed", time.time() - t0)
print("known anchor alpha(2.0) = 2.50290787509589")

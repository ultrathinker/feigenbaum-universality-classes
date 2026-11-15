import time
from mpmath import mp, mpf, nstr
mp.dps = 40
import core

N = 24
d0 = mpf(2)
a, al, it, res = core.solve_fixedpoint(N, d0)
d_list = [mpf('2.00'), mpf('2.01'), mpf('2.02'), mpf('2.03'), mpf('2.04'), mpf('2.05')]
a_cur, al_cur = a, al
for d in d_list:
    if d != d0:
        a_cur, al_cur, it, res = core.solve_fixedpoint(N, d, a_cur, al_cur)
    A = core.build_DT_matrix(a_cur, al_cur, d, N)
    eigs = []
    M = A
    for k in range(5):
        lam, v, itc = core.power_iter(M, N)
        eigs.append(lam)
        M = core.deflate(M, lam, v, N)
    andneg = al_cur ** (-d)
    print("d=%s alpha=%s" % (nstr(d, 6), nstr(al_cur, 15)))
    print("  eig1..5:", [nstr(e, 12) for e in eigs])
    print("  alpha^-d=%s  eig2/alpha^-d - 1 = %s" % (nstr(andneg, 12), nstr(eigs[1]/andneg - 1, 8)))

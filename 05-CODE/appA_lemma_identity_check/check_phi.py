# Orchestrator step 0: test candidate Phi_k(t) = G^{dk+1} - d t^{k+1} G' - (G - d t G')  vs DT with eigenvalue alpha^{-dk}
import json, sys
from mpmath import mp, mpf, fabs, sign
mp.dps = 40
D = json.load(open('../sec2_spectra_grid_d/spectra_grid.json', encoding='utf-8'))
pts = {p['d']: p for p in D['points']}
def mkG(a):
    a = [mpf(x) for x in a]
    def G(t):
        r = a[-1]
        for c in reversed(a[:-1]): r = r*t + c
        return r
    def Gp(t):
        s = mpf(0)
        for j in range(len(a)-1, 0, -1): s = s*t + j*a[j]
        return s
    return G, Gp
def run(dstr, k):
    p = pts[dstr]; d = mpf(dstr); al = mpf(p['alpha48']); G, Gp = mkG(p['_a48'])
    def K(t): return G(t) - d*t*Gp(t)
    def Phi(t): return G(t)**(int(d*k)+1) - d*t**(k+1)*Gp(t) - K(t)   # requires d*k integer
    def DT(F, t):
        U = G(t/al**d); W = fabs(U)**d
        return -(F(mpf(1))/G(mpf(1)))*K(t) - al*F(W) - d*al*sign(U)*fabs(U)**(d-1)*Gp(W)*F(t/al**d)
    lam = al**(-d*k)
    mx = mxr = mpf(0)
    for i in range(1, 241):
        t = mpf(i)/200          # (0,1.2]
        ph = Phi(t); r = DT(Phi, t) - lam*ph
        mx = max(mx, fabs(ph)); mxr = max(mxr, fabs(r))
    print("d=%s k=%d lam=alpha^-dk=%s  max|Phi|=%.3e  max|DTPhi-lam Phi|=%.3e  rel=%.3e" % (dstr, k, mp.nstr(lam, 18), float(mx), float(mxr), float(mxr/mx)))
for dstr in ['2.0', '4.0', '6.0']:
    for k in (1, 2):
        run(dstr, k)
for dstr in ['3.0', '5.0']:   # odd d: k=2 -> dk even -> predicted exact; k=1 -> dk odd -> Phi has G^{d+1} with even power... test anyway
    for k in (1, 2):
        try: run(dstr, k)
        except Exception as e: print(dstr, k, 'ERR', e)

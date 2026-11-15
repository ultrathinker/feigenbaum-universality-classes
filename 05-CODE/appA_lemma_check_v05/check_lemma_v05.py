# check_lemma_v05.py — расширенная проверка тождества DT Phi_k = alpha^{-dk} Phi_k
# (i) целые чётные dk; (ii) нечётная вещественная ветвь [G]^{dk+1} для дробных dk;
# (iii) x-представление: Phi_k(t) с t=|x|^d совпадает с w_k o g - g' w_k;
# (iv) g'(1) = -alpha^{d-1} численно.
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

def Kfn(d, G, Gp):
    def K(t): return G(t) - d*t*Gp(t)
    return K

def dt_phi(d, G, Gp, K, F, al):
    def op(t):
        U = G(t/al**d); W = fabs(U)**d
        return -(F(mpf(1))/G(mpf(1)))*K(t) - al*F(W) - d*al*sign(U)*fabs(U)**(d-1)*Gp(W)*F(t/al**d)
    return op

def max_resid(dstr, k, branch='integer'):
    p = pts[dstr]; d = mpf(dstr); al = mpf(p['alpha48'])
    G, Gp = mkG(p['_a48'])
    K = Kfn(d, G, Gp)
    if branch == 'integer':
        n = int(d*k) + 1
        def F(t): return G(t)**n - d*t**(k+1)*Gp(t) - K(t)
    elif branch == 'odd':
        n = d*k + 1
        def F(t):
            v = G(t)
            return sign(v)*fabs(v)**n - d*t**(k+1)*Gp(t) - K(t)
    else:
        raise ValueError(branch)
    op = dt_phi(d, G, Gp, K, F, al)
    lam = al**(-d*k)
    mx = mxr = mpf(0)
    for i in range(1, 241):
        t = mpf(i)/200
        ph = F(t); r = op(t) - lam*ph
        mx = max(mx, fabs(ph)); mxr = max(mxr, fabs(r))
    return float(mx), float(mxr), float(mxr/mx if mx != 0 else 0)

def xrep_resid(dstr, k):
    p = pts[dstr]; d = mpf(dstr); al = mpf(p['alpha48'])
    G, Gp = mkG(p['_a48'])
    K = Kfn(d, G, Gp)
    n = int(d*k) + 1
    def Phi_t(t): return G(t)**n - d*t**(k+1)*Gp(t) - K(t)
    def g(x):
        t = fabs(x)**d
        return G(t)
    def gp(x):
        # g'(x) = d*|x|^{d-2}*x * G'(|x|^d)
        if x == 0: return mpf(0)
        t = fabs(x)**d
        return d*fabs(x)**(d-2)*x*Gp(t)
    def w(x): return x*fabs(x)**(d*k) - x
    def phi_x(x):
        return w(g(x)) - gp(x)*w(x)
    mxr = mpf(0); mx = mpf(0)
    for i in range(1, 101):
        x = mpf(i)/100
        t = fabs(x)**d
        ph_t = Phi_t(t)
        ph_x = phi_x(x)
        mx = max(mx, fabs(ph_x))
        mxr = max(mxr, fabs(ph_x - ph_t))
    return float(mx), float(mxr), float(mxr/mx if mx != 0 else 0)

def gp_at_one(dstr):
    p = pts[dstr]; d = mpf(dstr); al = mpf(p['alpha48'])
    G, Gp = mkG(p['_a48'])
    # g(x)=G(|x|^d), so g'(1)=d*G'(1) (chain rule with t=|x|^d)
    gp1 = float(d*Gp(mpf(1)))
    rhs = float(-al**(d-1))
    return gp1, rhs

def run_all():
    out = []
    out.append("=== (i) INTEGER EVEN dk ===")
    for dstr in ['2.0', '4.0', '6.0']:
        for k in (1, 2, 3):
            mx, mxr, rel = max_resid(dstr, k, 'integer')
            out.append(f"d={dstr} k={k} (dk={float(mpf(dstr)*k):.4g})  max|Phi|={mx:.3e}  max|DTPhi-lam Phi|={mxr:.3e}  rel={rel:.3e}")
    out.append("")
    out.append("=== (ii) ODD BRANCH [G]^{dk+1} for FRACTIONAL dk ===")
    for dstr in ['2.5', '2.05', '3.0', '5.0']:
        for k in (1,):
            mx, mxr, rel = max_resid(dstr, k, 'odd')
            out.append(f"d={dstr} k={k} (dk={float(mpf(dstr)*k):.4g})  odd-branch  max|Phi|={mx:.3e}  max|DTPhi-lam Phi|={mxr:.3e}  rel={rel:.3e}")
    out.append("")
    out.append("=== (iii) x-REPRESENTATION: Phi_k(t=|x|^d) vs w_k o g - g' w_k ===")
    for dstr in ['2.0', '4.0', '6.0']:
        for k in (1, 2):
            mx, mxr, rel = xrep_resid(dstr, k)
            out.append(f"d={dstr} k={k}  max|phi_x|={mx:.3e}  max|phi_x-Phi_t|={mxr:.3e}  rel={rel:.3e}")
    out.append("")
    out.append("=== (iv) g'(1) = -alpha^{d-1} ===")
    for dstr in ['2.0', '4.0', '6.0']:
        gp1, rhs = gp_at_one(dstr)
        out.append(f"d={dstr}  g'(1)={gp1:.12f}  -alpha^(d-1)={rhs:.12f}  diff={gp1-rhs:.3e}")
    return "\n".join(out)

if __name__ == '__main__':
    s = run_all()
    print(s)
    with open('check_lemma_v05.out', 'w', encoding='utf-8') as f:
        f.write(s + "\n")

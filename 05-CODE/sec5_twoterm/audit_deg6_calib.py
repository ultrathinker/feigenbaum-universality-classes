# [oid:twotermC-AUDIT] NOT RIGOROUS - orchestrator: independent deg-6 calibration test of the min-rule.
# Own Newton cascade a_1..a_14 for f_a(x)=1-a*x^6 (NOT reusing agent B's cascade), own Richardson a_F,
# then linear LS fit e(m) = -D*lam2^m - C*delta^-m + G*delta^+m (G absorbs a_F error - lesson of tonight).
# Operator constants delta6, lam2(6) taken from agent B's certified phase 1 (cross-checked vs my cascade).
# min-rule prediction for deg 6: dominant ratio = min(delta6, 1/lam2) = 1/0.371071 = 2.6949.
import json, time, sys
from mpmath import mp, mpf, fabs, matrix, lu_solve

mp.dps = 100
T0 = time.time()
def log(msg):
    print(f"[+{time.time()-T0:6.1f}s] {msg}"); sys.stdout.flush()

d6   = mpf("9.2962468327713700828")        # agent B operator value (certified, res 4.8e-38)
lam2 = mpf("0.37107146768193013247")       # agent B certified

def superstable_k(k, a0):
    P = 2 ** k
    a = mpf(a0)
    for it in range(80):
        x = mpf(0); dx = mpf(0)
        for i in range(P):
            x2 = x * x; x5 = x2 * x2 * x; x6 = x5 * x
            dx = -x6 - 6 * a * x5 * dx
            x = 1 - a * x6
        step = x / dx
        a = a - step
        if fabs(step) < mpf(10) ** (-85):
            break
    x = mpf(0); mins = []; nxt = 2
    for i in range(1, P + 1):
        x = 1 - a * x ** 6
        if i == nxt:
            mins.append(fabs(x)); nxt *= 2
    return a, (min(mins[:-1]) if len(mins) > 1 else mpf(1))

a = {1: superstable_k(1, mpf(1))[0]}
log(f"a_1 = {mp.nstr(a[1], 20)}")
prev_gap = None
for k in range(2, 15):
    if k == 2:
        seeds = [mpf(s) for s in ("1.55", "1.5", "1.6", "1.45", "1.65")]
    else:
        base = a[k - 1] + (a[k - 1] - a[k - 2]) / d6
        h = (a[k - 1] - a[k - 2]) / d6
        seeds = [base, base + h / 4, base - h / 4, base + h / 2]
    ak = None
    for s in seeds:
        cand, minsep = superstable_k(k, s)
        # reject lower-period impostors (exact-period trap, §7.4) and out-of-order values
        if minsep > mpf(10) ** (-8) and cand > a[k - 1]:
            ak = cand; break
    if ak is None:
        log(f"a_{k}: ALL SEEDS FAILED exact-period check - abort"); sys.exit(1)
    a[k] = ak
    gap = a[k] - a[k - 1]
    dexp = (prev_gap / gap) if (prev_gap is not None and gap != 0) else None
    log(f"a_{k} = {mp.nstr(a[k], 25)}  minsep={mp.nstr(minsep,3)}" + (f"  delta_est={mp.nstr(dexp,10)}" if dexp else ""))
    prev_gap = gap

# own a_F via Richardson from top levels (r = 1/d6, two-step)
r = 1 / d6
aF_1 = a[14] + (a[14] - a[13]) * r / (1 - r)
aF_2 = a[13] + (a[13] - a[12]) * r / (1 - r)
a_F = aF_1 + (aF_1 - aF_2) * (r / (1 - r)) * r   # second-order touch-up
log(f"a_F(own) = {mp.nstr(aF_1, 30)}  (2nd-order {mp.nstr(a_F, 30)}; agent B: 1.68326019820464300051197)")
a_F = aF_1
L6 = (a[14] - a_F) * d6 ** 13   # pilot convention: L = lim (a_k - a_F)*delta^(k-1)
A6 = -1 / L6
log(f"L6(own) = {mp.nstr(L6, 20)}  A6(own) = {mp.nstr(A6, 20)}  (agent B: L6=-0.7568191, A6=1.3213197)")

# m=12,13 excluded: identically ~0 by circularity (a13,a14 build a_F and L6)
ms = list(range(6, 12))
e = {m: A6 * (a[m + 1] - a_F) * d6 ** m + 1 for m in ms}
log("signed e6(m) and ratios (m>=12 shown for reference, excluded from fits as circular):")
shown = ms + [12, 13]
for m in shown:
    if m not in e:
        e[m] = A6 * (a[m + 1] - a_F) * d6 ** m + 1
for m in shown[:-1]:
    nxt = e[m + 1]
    log(f"  m={m}: e={mp.nstr(e[m], 10)}  e(m)/e(m+1)={mp.nstr(e[m]/nxt, 8) if nxt != 0 else 'inf'}")
log(f"  m={shown[-1]}: e={mp.nstr(e[shown[-1]], 10)}")

def fit(basis):
    n = len(basis)
    Amat = [[sum(b1(m) * b2(m) / e[m] ** 2 for m in ms) for b2 in basis] for b1 in basis]
    b = [sum(e[m] * b1(m) / e[m] ** 2 for m in ms) for b1 in basis]
    x = lu_solve(matrix(Amat), matrix(b))
    resid = mp.sqrt(sum(((e[m] - sum(x[i] * basis[i](m) for i in range(n))) / e[m]) ** 2 for m in ms) / len(ms))
    return list(x), resid

lam3 = mpf("-0.34196832254780210510")   # agent B certified; |lam3| close to lam2 - near-degenerate slow pair
b_lam  = lambda m: -(lam2 ** m)
b_lam3 = lambda m: -(lam3 ** m)
b_dec  = lambda m: -(d6 ** (-m))
b_grow = lambda m: d6 ** m

x1, r1 = fit([b_lam]);            log(f"fit lam2 only:          D={mp.nstr(x1[0],8)}  relRMS={mp.nstr(r1,4)}")
x2, r2 = fit([b_dec]);            log(f"fit delta^-m only:      C={mp.nstr(x2[0],8)}  relRMS={mp.nstr(r2,4)}")
x3, r3 = fit([b_lam, b_grow]);    log(f"fit lam2 + G*d^m:       D={mp.nstr(x3[0],8)} G={mp.nstr(x3[1],4)}  relRMS={mp.nstr(r3,4)}  eps_aF={mp.nstr(x3[1]/A6,4)}")
x4, r4 = fit([b_lam, b_lam3]);    log(f"fit lam2 + lam3:        D={mp.nstr(x4[0],8)} D3={mp.nstr(x4[1],8)}  relRMS={mp.nstr(r4,4)}")
x5, r5 = fit([b_lam, b_lam3, b_grow]); log(f"fit lam2 + lam3 + G:    D={mp.nstr(x5[0],8)} D3={mp.nstr(x5[1],8)} G={mp.nstr(x5[2],4)}  relRMS={mp.nstr(r5,4)}  eps_aF={mp.nstr(x5[2]/A6,4)}")
log(f"min-rule prediction: dominant ratio = 1/lam2 = {mp.nstr(1/lam2, 8)} (delta6 = {mp.nstr(d6, 8)} is the fast one here)")
log("DONE")

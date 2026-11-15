# [oid:twotermC-AUDIT] NOT RIGOROUS - orchestrator: does an a_F offset explain the deg-4 tail drift?
# Model: e4(m) = -D*lam2^m + G*delta4^m, where G = A4*eps and eps = error in a_F (Richardson tail).
# Growing delta4^m signature = the candidate agent C did NOT test (they tested decaying delta4^-m).
# Uses pilot cascade values directly (verified exact to ~5e-46 by audit_a14_check.py) - no Newton needed.
import json, sys
from mpmath import mp, mpf, fabs, mp

mp.dps = 120
MS = json.load(open(r"../sec3_pilot_d4/mside_cascade.json", encoding="utf-8"))
a_F = mpf(MS["a_F"]); A4 = mpf(MS["A4"]); d4 = mpf(MS["delta4"])
lam2 = mpf("0.29183840898354510453")

ms = list(range(7, 14))  # m=7..13, non-circular band (cascade k=8..14 do not enter a_F's Richardson chain)
e = {m: A4 * (mpf(MS["cascade"][str(m + 1)]) - a_F) * d4 ** m + 1 for m in ms}
print("signed e4(m):")
for m in ms:
    print(f"  m={m}: {mp.nstr(e[m], 12)}")

def fit(basis):  # linear LS with relative weights, basis = list of functions of m
    n = len(basis)
    Amat = [[sum(b1(m) * b2(m) / e[m] ** 2 for m in ms) for b2 in basis] for b1 in basis]
    b = [sum(e[m] * b1(m) / e[m] ** 2 for m in ms) for b1 in basis]
    from mpmath import matrix, lu_solve
    x = lu_solve(matrix(Amat), matrix(b))
    resid = mp.sqrt(sum(((e[m] - sum(x[i] * basis[i](m) for i in range(n))) / e[m]) ** 2 for m in ms) / len(ms))
    return list(x), resid

b_lam = lambda m: -(lam2 ** m)
b_grow = lambda m: d4 ** m
b_decay = lambda m: -(d4 ** (-m))

x1, r1 = fit([b_lam])
print(f"\nM0L pure-lambda2:            D={mp.nstr(x1[0],8)}  relRMS={mp.nstr(r1,4)}")
x2, r2 = fit([b_lam, b_decay])
print(f"M2  lambda2 + C*delta^-m:    D={mp.nstr(x2[0],8)} C={mp.nstr(x2[1],8)}  relRMS={mp.nstr(r2,4)}")
x3, r3 = fit([b_lam, b_grow])
print(f"NEW lambda2 + G*delta^+m:    D={mp.nstr(x3[0],8)} G={mp.nstr(x3[1],4)}  relRMS={mp.nstr(r3,4)}")
if r3 < r1:
    eps = x3[1] / A4
    print(f"  -> implied a_F error eps = G/A4 = {mp.nstr(eps, 6)}")
    print(f"     (pilot's own a_F Richardson last-vs-prev was ~4e-23; plausible if |eps| ~ 1e-22..1e-20)")
    print(f"  -> corrected a_F* = a_F + eps = {mp.nstr(a_F + eps, 30)}")
    # refit ratios with corrected a_F: do ratios flatten to 1/lam2=3.4266?
    e_corr = {m: A4 * (mpf(MS["cascade"][str(m + 1)]) - (a_F + eps)) * d4 ** m + 1 for m in ms}
    print("  corrected ratios e4*(m)/e4*(m+1):")
    for i in range(len(ms) - 1):
        print(f"    {ms[i]}/{ms[i+1]}: {mp.nstr(e_corr[ms[i]] / e_corr[ms[i+1]], 8)}")
print(f"\ncandidate 1/lam2 = {mp.nstr(1/lam2, 8)}")

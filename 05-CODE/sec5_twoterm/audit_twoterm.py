# [oid:twotermC-AUDIT] NOT RIGOROUS - orchestrator's independent audit of agent C phase 1.
# Own Newton (not C's code), signed e(m) for u=-1, m=7..13, M0-residual sign pattern,
# and joint (C,D) fit quality vs M0-only.
import json, time, sys
from mpmath import mp, mpf, fabs, mpc

mp.dps = 120
T0 = time.time()
def log(msg):
    print(f"[+{time.time()-T0:6.1f}s] {msg}"); sys.stdout.flush()

meta = json.load(open(r"../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json", encoding="utf-8"))["meta"]
c_F   = mpf(meta["c_F_str"].strip("()").split()[0].rstrip(","))
A     = mpf(meta["A_str"].strip("()").split()[0].rstrip(","))
delta = mpf(meta["delta_str"])
lam2  = mpf("-0.1236527125526870243353937")

def superstable(m):
    P = 2 ** (m + 1)
    c = c_F + mpf(-1) / (A * delta ** m)
    for it in range(80):
        z = mpf(0); dz = mpf(0)
        for k in range(P):
            dz = 2 * z * dz + 1
            z = z * z + c
        step = z / dz
        c = c - step
        if fabs(step) < mpf(10) ** (-100):
            break
    z = mpf(0); mins = []; nxt = 2
    for k in range(1, P + 1):
        z = z * z + c
        if k == nxt:
            mins.append(fabs(z)); nxt *= 2
    return c, min(mins[:-1])

ms = list(range(7, 14))
e = {}
for m in ms:
    c, minsep = superstable(m)
    e[m] = A * (c - c_F) * delta ** m + 1
    log(f"m={m}: e(m)={mp.nstr(e[m], 12)}  minsep={mp.nstr(minsep,3)}")

# M0-only weighted LS (relative weights, as C did): minimize sum ((e+C*d^-m)/e)^2
num = sum((-e[m]) * (delta ** -m) / e[m] ** 2 for m in ms)
den = sum((delta ** -m) ** 2 / e[m] ** 2 for m in ms)
C_M0 = num / den
log(f"C_M0 (rel-weighted) = {mp.nstr(C_M0, 10)}   (agent C: -0.0700219)")

log("residual r(m) = e(m) + C_M0*delta^-m :")
r = {m: e[m] + C_M0 * delta ** -m for m in ms}
alt = 0
for i, m in enumerate(ms):
    ratio = fabs(r[m] / r[ms[i+1]]) if i + 1 < len(ms) else None
    sgn = "+" if r[m] > 0 else "-"
    log(f"  m={m}: sign={sgn}  |r|={mp.nstr(fabs(r[m]),4)}" + (f"  |r(m)/r(m+1)|={mp.nstr(ratio,5)}" if ratio else ""))
    if i + 1 < len(ms) and (r[m] > 0) != (r[ms[i+1]] > 0):
        alt += 1
ratios = [fabs(r[ms[i]] / r[ms[i+1]]) for i in range(len(ms) - 1)]
gm = mp.exp(sum(mp.log(x) for x in ratios) / len(ratios))
log(f"alternations: {alt}/{len(ms)-1}; geometric-mean ratio = {mp.nstr(gm,6)}  (predict 1/|lam2|=8.0872)")

# joint (C,D) LS fit, relative weights
def fit2():
    a11 = sum((delta ** -m) ** 2 / e[m] ** 2 for m in ms)
    a12 = sum((delta ** -m) * (lam2 ** m) / e[m] ** 2 for m in ms)
    a22 = sum((lam2 ** m) ** 2 / e[m] ** 2 for m in ms)
    b1  = sum((-e[m]) * (delta ** -m) / e[m] ** 2 for m in ms)
    b2  = sum((-e[m]) * (lam2 ** m) / e[m] ** 2 for m in ms)
    det = a11 * a22 - a12 * a12
    return (b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det

C2, D2 = fit2()
log(f"joint fit: C={mp.nstr(C2,10)} D={mp.nstr(D2,10)}  (agent C: C=-0.070018, D=+0.0019849)")
res0 = mp.sqrt(sum(((e[m] + C_M0 * delta ** -m) / e[m]) ** 2 for m in ms) / len(ms))
res2 = mp.sqrt(sum(((e[m] + C2 * delta ** -m + D2 * lam2 ** m) / e[m]) ** 2 for m in ms) / len(ms))
log(f"rel RMS residual: M0-only={mp.nstr(res0,4)}  M2={mp.nstr(res2,4)}  improvement={mp.nstr(res0/res2,4)}x")
log("DONE")

# [oid:twotermC-AUDIT] NOT RIGOROUS - orchestrator check of the m=13 outlier flag:
# recompute quartic cascade a_12..a_14 (family 1 - a*x^4) at dps 120 with own Newton,
# compare vs quartic_pilot/mside_cascade.json, recompute e4(11..13) ratios.
import json, time, sys
from mpmath import mp, mpf, fabs

mp.dps = 120
T0 = time.time()
def log(msg):
    print(f"[+{time.time()-T0:6.1f}s] {msg}"); sys.stdout.flush()

MS = json.load(open(r"../sec3_pilot_d4/mside_cascade.json", encoding="utf-8"))
a_F = mpf(MS["a_F"]); L4 = mpf(MS["L4"]); A4 = mpf(MS["A4"]); d4 = mpf(MS["delta4"])

def superstable_k(k, a0):
    P = 2 ** k
    a = mpf(a0)
    for it in range(80):
        x = mpf(0); dx = mpf(0)
        for i in range(P):
            dx = -x ** 4 - 4 * a * x ** 3 * dx
            x = 1 - a * x ** 4
        step = x / dx
        a = a - step
        if fabs(step) < mpf(10) ** (-100):
            break
    x = mpf(0); mins = []; nxt = 2
    for i in range(1, P + 1):
        x = 1 - a * x ** 4
        if i == nxt:
            mins.append(fabs(x)); nxt *= 2
    return a, (min(mins[:-1]) if len(mins) > 1 else mpf(1)), it + 1

e = {}
for k in (12, 13, 14):
    a_pilot = mpf(MS["cascade"][str(k)])
    a_mine, minsep, nit = superstable_k(k, a_pilot)
    diff = fabs(a_mine - a_pilot)
    m = k - 1
    e[m] = A4 * (a_mine - a_F) * d4 ** m + 1
    e_pilot = A4 * (a_pilot - a_F) * d4 ** m + 1
    log(f"k={k} (m={m}): |a_mine - a_pilot| = {mp.nstr(diff, 4)}  minsep={mp.nstr(minsep,3)} steps={nit}")
    log(f"   e4({m}) mine = {mp.nstr(e[m], 12)}   from pilot a = {mp.nstr(e_pilot, 12)}")

log(f"ratio e4(11)/e4(12) mine = {mp.nstr(e[11]/e[12], 8)}")
log(f"ratio e4(12)/e4(13) mine = {mp.nstr(e[12]/e[13], 8)}  (agent C from pilot artifact: 3.585; band 3.36-3.46; 1/lam2=3.4266)")
log("DONE")

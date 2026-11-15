# [oid:deg6B-AUDIT] NOT RIGOROUS - orchestrator's gate-1 audit of agent B (kilo) deg-6 pilot.
# Independent of kilo's code AND of kilo's audit_self.py:
#  P1: data-level recompute of every verdict statistic straight from match_results.json
#  P2: own complex Newton spot-checks (accepted zeros across |u|, real chaotic-side zero,
#      anomalous-subfamily failure, 2 controls) with exact-period check
# (P3, done earlier: independent calibration-node rate test in agentC_twoterm\audit_deg6_calib.py)
import json, time, sys
from mpmath import mp, mpf, mpc, fabs

mp.dps = 60
T0 = time.time()
def log(msg):
    print(f"[+{time.time()-T0:6.1f}s] {msg}"); sys.stdout.flush()

MR = json.load(open(r"match_results.json", encoding="utf-8"))
MC = json.load(open(r"mside_cascade.json", encoding="utf-8"))

def getconst(name, fallback_key):
    v = MR.get(name)
    if isinstance(v, str):
        return mpf(v)
    fv = MC.get(fallback_key)
    return mpf(fv) if isinstance(fv, str) else mpf(repr(v))

a_F = getconst("a_F", "a_F"); A6 = getconst("A6", "A6"); d6 = getconst("delta6", "delta6")
log(f"constants: a_F={mp.nstr(a_F,25)} A6={mp.nstr(A6,20)} delta6={mp.nstr(d6,20)}")

# ---------- P1: data-level verdict recompute ----------
cat = MR["catalog_results"]; V = MR["verdict"]
interior = [c for c in cat if not c["edge"]]
acc_int = [c for c in interior if c["accepted"]]
edge = [c for c in cat if c["edge"]]
acc_edge = [c for c in edge if c["accepted"]]
dus = sorted(c["d_u"] for c in acc_int)
med = dus[len(dus)//2] if len(dus) % 2 else 0.5*(dus[len(dus)//2-1]+dus[len(dus)//2])
ctrl_acc = sum(1 for c in MR["controls"] if c["accepted"])
r89 = sorted(r["ratio_89"] for r in MR["rate_study"]); r910 = sorted(r["ratio_910"] for r in MR["rate_study"])
med89 = r89[len(r89)//2]; med910 = r910[len(r910)//2]
gm = (med89*med910)**0.5
checks = [
    ("interior_total", len(interior), V["interior_total"]),
    ("interior_accepted", len(acc_int), V["interior_accepted"]),
    ("edge_accepted", len(acc_edge), V["edge_accepted"]),
    ("d_u_median", round(med, 12), round(V["d_u_median"], 12)),
    ("d_u_min", min(dus), V["d_u_min"]),
    ("d_u_max", max(dus), V["d_u_max"]),
    ("controls_accepted", ctrl_acc, V["controls_accepted"]),
    ("rate_median_89", round(med89, 10), round(V["rate_median_89"], 10)),
    ("rate_median_910", round(med910, 10), round(V["rate_median_910"], 10)),
]
ok = True
for name, mine, theirs in checks:
    good = (mine == theirs)
    ok &= good
    log(f"P1 {name}: mine={mine} report={theirs} {'OK' if good else 'MISMATCH'}")
log(f"P1 GM of rate medians = {gm:.6f} (report 2.6893; 1/lambda2=2.694899)")
tier1 = sum(1 for c in cat if c["accepted"] and c.get("tier") == 1)
log(f"P1 tier-1 accepted: {tier1} (report 332); exact_period_ok on all accepted: {all(c['exact_period_ok'] for c in cat if c['accepted'])}")
log(f"P1 verdict-stats: {'ALL OK' if ok else 'HAS MISMATCHES'}")

# ---------- P2: own Newton spot checks ----------
def newton_a(u, m=9):
    P = 2 ** (m + 1)
    seed = a_F + u / (A6 * d6 ** m)
    a = mpc(seed) if isinstance(u, mpc) or (hasattr(u, 'imag') and u.imag != 0) else mpf(u.real) * 0 + seed.real if isinstance(seed, mpc) else seed
    a = mpc(seed)
    for it in range(60):
        x = mpc(0); dx = mpc(0)
        for i in range(P):
            x2 = x * x; x5 = x2 * x2 * x; x6 = x5 * x
            dx = -x6 - 6 * a * x5 * dx
            x = 1 - a * x6
        step = x / dx
        a = a - step
        if abs(step) < mpf(10) ** (-50):
            break
    x = mpc(0); mins = []; nxt = 2
    for i in range(1, P + 1):
        x = 1 - a * x ** 6
        if i == nxt:
            mins.append(abs(x)); nxt *= 2
    res = abs(x)   # |x_P(a)| - must be ~0 for a genuine superstable point
    u_act = A6 * (a - a_F) * d6 ** m
    return a, min(mins[:-1]), abs(u_act - u), it + 1, res

acc_sorted = sorted(acc_int, key=lambda c: abs(complex(c["u_cat_re"], c["u_cat_im"])))
picks = [acc_sorted[1], acc_sorted[len(acc_sorted)//2], acc_sorted[-2]]
reals = [c for c in acc_int if abs(c["u_cat_im"]) < 1e-12 and c["u_cat_re"] > 0]
if reals: picks.append(reals[0])
fails = [c for c in interior if not c["accepted"]]
if fails: picks.append(fails[0])

for c in picks:
    u = mpc(c["u_cat_re"], c["u_cat_im"])
    a, minsep, du, nit, res = newton_a(u)
    kdu = c.get("d_u")
    cmp = f"(kilo {kdu:.4e}) rel-diff={abs(float(du)-kdu)/max(kdu,1e-300):.1e}" if kdu is not None else "(kilo d_u: none - FAILED entry)"
    log(f"P2 {c['id']} |u|={abs(complex(c['u_cat_re'],c['u_cat_im'])):.3f} acc={c['accepted']}: "
        f"my d_u={float(du):.4e} {cmp} minsep={float(minsep):.2e} steps={nit}")
    if c.get("a_ss_re") is not None:
        da = abs(mpc(c["a_ss_re"], c["a_ss_im"]) - a)
        log(f"     |a_mine - a_kilo| = {float(da):.2e}")

for c in MR["controls"]:
    u = mpc(c["u_re"], c["u_im"])
    a, minsep, du, nit, res = newton_a(u)
    converged = res < mpf(10) ** (-20)
    verdict = "FAILS (good)" if (not converged or du > 0.01) else "ACCEPTED (BAD!)"
    log(f"P2 CONTROL {c['id']}: newton_converged={converged} |x_P|={float(res):.2e} steps={nit} "
        f"my d_u={float(du):.4e} (kilo best {c['best_dist_u']:.4e}) accepted_by_kilo={c['accepted']} -> {verdict}")

log("DONE")

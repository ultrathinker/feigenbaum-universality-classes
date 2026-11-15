# [oid:lam2A-AUDIT] NOT RIGOROUS - orchestrator's independent audit (third code path).
# Checks agent A's REPORT_LAMBDA2_DEG2.md without reusing its code:
#  P1: catalog identity claim (u_cat == representative_u == m=7 trajectory) across all zeros
#  P2: calibration test u=-1 (catalog-free!): rate of |A(c_m-c_F)*delta^m + 1| for m=8..11
#  P3: independent recompute of real zeros Z06, Z44, Z50 (own Richardson + own Newton)
#  P4: spot check of one complex zero Z19
import json, time, sys
from mpmath import mp, mpf, mpc, fabs

mp.dps = 80
T0 = time.time()
LOG = open(r"audit_progress.log", "a", buffering=1)
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{time.time()-T0:6.1f}s] {msg}"
    print(line); sys.stdout.flush(); LOG.write(line + "\n")

RES = json.load(open(r"../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json", encoding="utf-8"))
ZM  = json.load(open(r"D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json", encoding="utf-8"))
meta = RES["meta"]
c_F   = mpf(meta["c_F_str"].strip("()").split()[0].rstrip(","))
A     = mpf(meta["A_str"].strip("()").split()[0].rstrip(","))
delta = mpf(meta["delta_str"])
log(f"constants: c_F={mp.nstr(c_F,25)} A={mp.nstr(A,25)} delta={mp.nstr(delta,25)}")

OUT = {"tag": "[oid:lam2A-AUDIT] NOT RIGOROUS", "parts": {}}
def dump():
    json.dump(OUT, open(r"audit_results.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)

# ---------- P1: catalog identity ----------
roots = {r["id"]: r for r in ZM["roots"]}
p1 = {"n_checked": 0, "max_diff_ucat_vs_repr": 0.0, "max_diff_repr_vs_m7": 0.0, "mismatches": []}
for e in RES["stage3_catalog"]:
    zid = e.get("id")
    if zid not in roots or "u_cat" not in e:
        continue
    r = roots[zid]
    uc = complex(e["u_cat"]["re"], e["u_cat"]["im"])
    ru = complex(r["representative_u"][0], r["representative_u"][1])
    m7 = next((t for t in r["level_trajectories"] if t["m"] == 7), None)
    d1 = abs(uc - ru)
    p1["max_diff_ucat_vs_repr"] = max(p1["max_diff_ucat_vs_repr"], d1)
    if m7:
        d2 = abs(ru - complex(m7["u"][0], m7["u"][1]))
        p1["max_diff_repr_vs_m7"] = max(p1["max_diff_repr_vs_m7"], d2)
    if d1 > 1e-12:
        p1["mismatches"].append(zid)
    p1["n_checked"] += 1
OUT["parts"]["P1_catalog_identity"] = p1
log(f"P1: checked={p1['n_checked']} max|u_cat-repr|={p1['max_diff_ucat_vs_repr']:.3e} "
    f"max|repr-m7traj|={p1['max_diff_repr_vs_m7']:.3e} mismatches={p1['mismatches']}")
dump()

# ---------- Newton machinery (own implementation) ----------
def superstable(c0, m, is_complex):
    """Newton for z_{2^(m+1)}(c)=0, z0=0, companion derivative. Returns (c*, min_j|z_{2^j}|, steps)."""
    P = 2 ** (m + 1)
    c = mpc(c0) if is_complex else mpf(c0)
    for it in range(60):
        z = c * 0; dz = c * 0
        for k in range(P):
            dz = 2 * z * dz + 1
            z = z * z + c
        step = z / dz
        c = c - step
        if fabs(step) < mpf(10) ** (-60):
            break
    # exact-period check: min over j<=m of |z_{2^j}(c*)|
    z = c * 0; mins = []
    nxt = 2; k = 0
    for k in range(1, P + 1):
        z = z * z + c
        if k == nxt:
            mins.append(fabs(z)); nxt *= 2
    minsep = min(mins[:-1]) if len(mins) > 1 else mpf(1)  # exclude j=m+1 (that's the zero itself)
    return c, minsep, it + 1

def rate_row(u_inf, is_complex, label):
    row = {"u_inf": str(u_inf), "d": {}, "minsep": {}, "ratios": {}}
    prev = None
    for m in (8, 9, 10, 11):
        seed = c_F + u_inf / (A * delta ** m)
        cst, minsep, nit = superstable(seed, m, is_complex)
        d = fabs(u_inf - A * (cst - c_F) * delta ** m)
        row["d"][m] = float(d); row["minsep"][m] = float(minsep)
        log(f"  {label} m={m}: d_u={float(d):.4e} minsep={float(minsep):.2e} newton_steps={nit}")
        if prev is not None:
            row["ratios"][f"{m-1}/{m}"] = float(prev / d)
        prev = d
    return row

# ---------- P2: calibration u=-1 (catalog-free) ----------
log("P2: calibration zero u=-1 (no catalog involved; limit is exactly -1 by A=-1/L)")
p2 = rate_row(mpf(-1), False, "u=-1")
OUT["parts"]["P2_calibration_u_minus1"] = p2
log(f"P2 ratios: {p2['ratios']}  candidates: delta={float(delta):.4f} 1/|lam2|=8.0872")
dump()

# ---------- P3: real zeros, own Richardson ----------
def richardson(zid):
    r = roots[zid]
    tr = {t["m"]: complex(t["u"][0], t["u"][1]) for t in r["level_trajectories"]}
    u9, u10 = tr[9], tr[10]
    rr = 1.0 / float(delta)
    u_inf = u10 + (u10 - u9) * rr / (1.0 - rr)
    return u_inf

p3 = {}
for zid in ("Z06", "Z44", "Z50"):
    u_inf_c = richardson(zid)
    u_inf = mpf(u_inf_c.real)  # real zeros
    log(f"P3 {zid}: u_inf(own Richardson)={u_inf_c.real:.10f} (im part {u_inf_c.imag:.1e})")
    p3[zid] = rate_row(u_inf, False, zid)
    dump()
OUT["parts"]["P3_real_zeros"] = p3

# ---------- P4: complex spot check Z19 ----------
u_inf_c = richardson("Z19")
log(f"P4 Z19: u_inf(own Richardson)={u_inf_c:.10f}")
p4 = rate_row(mpc(u_inf_c.real, u_inf_c.imag), True, "Z19")
OUT["parts"]["P4_complex_Z19"] = p4
dump()

# ---------- summary ----------
import statistics
all_r = []
for part in (p3["Z06"], p3["Z44"], p3["Z50"]):
    all_r.append(part["ratios"])
log("SUMMARY:")
log(f"  P2 u=-1 ratios: {p2['ratios']}")
for zid in p3: log(f"  P3 {zid} ratios: {p3[zid]['ratios']}")
log(f"  P4 Z19 ratios: {p4['ratios']}")
log(f"  candidates: delta={float(delta):.4f}, 1/|lambda2|=8.0872, alpha^2=6.2646")
log("DONE")
dump()

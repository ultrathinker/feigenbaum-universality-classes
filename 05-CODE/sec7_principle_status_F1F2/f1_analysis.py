#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F1 analysis, stage 1: assemble residual data and test the perturbation form.
NOT RIGOROUS (float). OMP threads capped at 2.
"""
import json, os, math, cmath
import numpy as np

BASE = r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип"
OUT = os.path.dirname(os.path.abspath(__file__))

def c(s):
    return complex(s)

def load(p):
    with open(os.path.join(BASE, p), encoding="utf-8") as f:
        return json.load(f)

# ---------- spectral anchors ----------
D6 = load("agentB_deg6/match_results.json")
D4 = load("quartic_pilot/match_results.json")
H6 = load("agentB_deg6/H6_zeros.json")
H4 = load("quartic_pilot/H4_zeros.json")

SPEC = {"d6": dict(delta=9.2962468327713700828344765663674575503066887803593,
                   lam2=0.3710714676819301324709834,
                   lam3=-0.3419683225478021050998591,
                   lam4=-0.2105185208142716405642349,
                   alpha=1.4677424503199009444538343151089737463687971293968,
                   A=1.3213197054428516594468381447179159106238901404509,
                   AU0=-0.6679695396324661394893857007410712354653),
        "d4": dict(delta=7.2846862170733433643089305679955530694780466189573,
                   lam2=0.291838408983545104530996961457,
                   lam3=-0.255664558474297406325347479489,
                   lam4=-0.155562761843395028402820942804,
                   alpha=1.6903029714052448533437801503241613482282780596328,
                   A=1.4904325463540699872206657264594941418825094844754,
                   AU0=-0.6087812101890336104440950663154768942781)}

# ---------- catalog data ----------
def build_catalog(mr, hz):
    """Return list of dicts per accepted interior zero."""
    lim = hz["zeros_limit_extrapolated"]
    rows = []
    for r in mr["catalog_results"]:
        if not r.get("accepted"):
            continue
        edge = r.get("edge", False)
        if edge:
            continue
        u_cat = complex(r["u_cat_re"], r["u_cat_im"])
        u_act = complex(r.get("u_actual_re", 0.0), r.get("u_actual_im", 0.0))
        e9 = u_act - u_cat
        d9 = abs(e9)
        # level-m zero positions from the extrapolation record (d=6 has u7,u8,u9)
        rec = {"id": r["id"], "u_cat": u_cat, "u_act9": u_act, "e9": e9, "d9": d9,
               "tier": r.get("tier"), "absu": abs(u_cat)}
        z = _best_zero(hz["zeros_limit_extrapolated"], u_cat)
        if z is not None:
            rec["u7"] = complex(z.get("u7_re", z["re"]), z.get("u7_im", z["im"]))
            rec["u8"] = complex(z.get("u8_re", z["re"]), z.get("u8_im", z["im"]))
            rec["u9"] = complex(z.get("u9_re", z["re"]), z.get("u9_im", z["im"]))
        rows.append(rec)
    return rows

def _best_zero(lim, u_cat):
    best, bd = None, 1e9
    for z in lim:
        zz = complex(z["re"], z["im"])
        d = abs(zz - u_cat)
        if d < bd:
            bd, best = d, z
    return best

cat6 = build_catalog(D6, H6)
cat4 = build_catalog(D4, H4)
print("catalog interior accepted: d6=%d d4=%d" % (len(cat6), len(cat4)))

# ---------- cross-section diagnostics at m=9 ----------
def diag(rows, tag, sp):
    d9 = np.array([abs(r["e9"]) for r in rows])
    absu = np.array([abs(r["u_cat"]) for r in rows])
    e9 = np.array([r["e9"] for r in rows])
    u = np.array([r["u_cat"] for r in rows])
    print("\n== %s  n=%d" % (tag, len(rows)))
    print("d9: median=%.3e min=%.3e max=%.3e" % (np.median(d9), d9.min(), d9.max()))
    # correlations
    print("corr(log d9, log|u|)      = %.3f" % np.corrcoef(np.log(d9), np.log(absu))[0, 1])
    print("corr(d9, |u|)             = %.3f" % np.corrcoef(d9, absu)[0, 1])
    print("corr(d9, |u|^2)           = %.3f" % np.corrcoef(d9, absu**2)[0, 1])
    # model: e9 = c0 + c1*u   (the dominant lambda2-term, u-linear + constant)
    #   and add u^2*delta^{-m} term explicitly.
    m9 = 9.0
    lam2m = sp["lam2"]**m9
    deltam = sp["delta"]**(-m9)
    X = np.column_stack([np.ones(len(u)), u, u**2 * deltam])
    coef, res, *_ = np.linalg.lstsq(X, e9, rcond=None)
    pred = X @ coef
    rel = np.abs(pred - e9) / np.maximum(np.abs(e9), 1e-12)
    print("OLS e9 = c0 + c1*u + c2*u^2*delta^{-m}:")
    print("   c0 = %.4e%+.4ei   c1 = %.4e%+.4ei   c2 = %.4e%+.4ei" %
          (coef[0].real, coef[0].imag, coef[1].real, coef[1].imag, coef[2].real, coef[2].imag))
    print("   |c0|/mean|e9| = %.3f   |c1*u|/|e9| median = %.3f   |c2*u^2 d^{-m}|/|e9| median = %.3e"
          % (abs(coef[0])/np.mean(np.abs(e9)),
             np.median(np.abs(coef[1]*u)/np.maximum(np.abs(e9),1e-12)),
             np.median(np.abs(coef[2]*u**2*deltam)/np.maximum(np.abs(e9),1e-12))))
    # pure linear-only fit (no delta term)
    X1 = np.column_stack([np.ones(len(u)), u])
    coef1, *_ = np.linalg.lstsq(X1, e9, rcond=None)
    pred1 = X1 @ coef1
    rel1 = np.abs(pred1 - e9) / np.maximum(np.abs(e9), 1e-12)
    print("OLS e9 = c0 + c1*u (lambda2-term only): median|rel err|=%.3f  pct>10%%=%.2f  pct>100%%=%.2f"
          % (np.median(rel1), 100*np.mean(rel1 > 0.1), 100*np.mean(rel1 > 1.0)))
    # scaling of |c1*u| term vs |u|
    s = np.abs(coef1[1]*u)/np.maximum(np.abs(e9),1e-12)
    print("   median |(c1*u)/e9| = %.3f  (c1*u explains fraction of residual)" % np.median(s))
    return rows, dict(coef=coef, coef1=coef1, d9=d9, absu=absu, e9=e9, u=u)

r6, diag6 = diag(cat6, "d=6 m=9 cross-section", SPEC["d6"])
r4, diag4 = diag(cat4, "d=4 m=9 cross-section", SPEC["d4"])

# ---------- out-of-sample fit (fit on 30% training, predict 70% test) ----------
def oos(rows, sp, frac=0.30, seed=20260815, with_delta=True):
    rng = np.random.RandomState(seed)
    n = len(rows)
    idx = np.arange(n)
    rng.shuffle(idx)
    ntr = int(frac * n)
    tr, te = idx[:ntr], idx[ntr:]
    e9 = np.array([r["e9"] for r in rows]); u = np.array([r["u_cat"] for r in rows])
    m9 = 9.0
    if with_delta:
        X = np.column_stack([np.ones(n), u, u**2*sp["delta"]**(-m9)])
    else:
        X = np.column_stack([np.ones(n), u])
    coef, *_ = np.linalg.lstsq(X[tr], e9[tr], rcond=None)
    pred = X[te] @ coef
    err = np.abs(pred - e9[te]); truth = np.abs(e9[te])
    ratio = np.maximum(err, 1e-12) / np.maximum(truth, 1e-12)
    within2 = np.mean(ratio < 2.0)
    within10 = np.mean(ratio < 10.0)
    print("OOS (%s, with_delta=%s): train=%d test=%d  median ratio(pred/true)=%.2f  "
          "within2x=%.1f%%  within10x=%.1f%%  mean|e|_test=%.2e"
          % (rows[0]["id"][:2], with_delta, ntr, n-ntr, np.median(ratio), 100*within2, 100*within10, np.mean(truth)))
    return dict(median_ratio=np.median(ratio), within2=within2, within10=within10)

print("\n-- out-of-sample prediction of d9 --")
for tag, rows, sp in [("d6", cat6, SPEC["d6"]), ("d4", cat4, SPEC["d4"])]:
    for wd in (True, False):
        oos(rows, sp, with_delta=wd)

# ---------- rate data: two-term model e(m) ~ C*lam2^m + D*delta^{-m} ----------
def fit_rate(d_u, ms, sp):
    """d_u: dict level->value. fit log-linear model |e| = |C*lam2^m + D*delta^{-m}| by brute force phase scan."""
    lam2, delta = sp["lam2"], sp["delta"]
    ms_arr = np.array([float(m) for m in ms], dtype=float)
    ys = np.array([d_u[m] for m in ms], dtype=float)
    # model e(m) = C*lam2^m + D*delta^{-m}  (real scalars, signed)
    # brute-force fit: for each zero, the sign/amplitude structure; we fit by least squares
    # over candidate (C,D) grid, then refine with scipy.
    best = None
    for sgn in (+1.0, -1.0):
        A = np.column_stack([sgn*lam2**ms_arr, delta**(-ms_arr)])
        coef, res, rank, sv = np.linalg.lstsq(A, ys, rcond=None)
        pred = A @ coef
        r = float(np.sum((pred-ys)**2))
        if best is None or r < best[0]:
            best = (r, coef)
    _, coef = best
    pred = np.column_stack([lam2**ms_arr, delta**(-ms_arr)]) @ coef
    return coef, pred

print("\n-- rate study (d6, 13 zeros: m=8,9,10) --")
rs = D6["rate_study"]
rows_rate = []
for row in rs:
    du = row["d_u"]
    ms = sorted(du.keys())
    if len(ms) < 3:
        continue
    coef, pred = fit_rate(du, ms, SPEC["d6"])
    vals = [du[m] for m in ms]
    lam2 = SPEC["d6"]["lam2"]; delta = SPEC["d6"]["delta"]
    # per-level ratio
    r89 = vals[0]/vals[1]; r910 = vals[1]/vals[2]
    rows_rate.append((row["id"], row["absu"], ms, vals, coef, r89, r910))
    print("%-6s |u|=%6.3f  d(m)=%s  r8/9=%.3f r9/10=%.3f  [1/l2=%.3f]  fit C=%.2e D=%.2e"
          % (row["id"], row["absu"], ["%.2e" % v for v in vals], r89, r910,
             1/abs(lam2), coef[0], coef[1]))
med_r89 = np.median([r[5] for r in rows_rate]); med_r910 = np.median([r[6] for r in rows_rate])
print("medians r8/9=%.3f r9/10=%.3f   (1/|lam2|=%.3f, 1/delta=%.3f)" %
      (med_r89, med_r910, 1/abs(SPEC["d6"]["lam2"]), 1/SPEC["d6"]["delta"]))

print("\n-- rate study (d4, 13 zeros) --")
rs4 = D4["rate_study"]
for row in rs4:
    du = row["d_u"]; ms = sorted(du.keys())
    if len(ms) < 3: continue
    vals = [du[m] for m in ms]
    r89 = vals[0]/vals[1]; r910 = vals[1]/vals[2]
    print("%-6s |u|=%6.3f  d(m)=%s  r8/9=%.3f r9/10=%.3f  [1/l2=%.3f]"
          % (row["id"], row["absu"], ["%.2e" % v for v in vals], r89, r910, 1/abs(SPEC["d4"]["lam2"])))
r89s4 = [du[sorted(du)[0]]/du[sorted(du)[1]] for du in [r["d_u"] for r in rs4] if len(du)>=3]
r910s4 = [du[sorted(du)[1]]/du[sorted(du)[2]] for du in [r["d_u"] for r in rs4] if len(du)>=3]
print("medians d4: r8/9=%.3f r9/10=%.3f  (1/|lam2|=%.3f)" %
      (np.median(r89s4), np.median(r910s4), 1/abs(SPEC["d4"]["lam2"])))

print("\n-- rate_full (d6, 9 late zeros: m=7..12) --")
rf = load("agentB_deg6/subfamily/rate_full.json")
for k, v in rf.items():
    du = v["d_u"]; ms = sorted(du.keys())
    vals = [du[m] for m in ms]
    print("%-6s d(m)=%s" % (k, ["%.2e" % vals[i] for i in range(len(ms))]))

print("\ndone stage 1")

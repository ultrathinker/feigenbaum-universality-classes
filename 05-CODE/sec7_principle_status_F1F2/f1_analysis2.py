#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F1 analysis, stage 2: baselines, rate-model fit, phase, u-dependence.
NOT RIGOROUS (float).
"""
import json, os, math, cmath
import numpy as np

BASE = r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип"
def load(p):
    with open(os.path.join(BASE, p), encoding="utf-8") as f:
        return json.load(f)

D6 = load("agentB_deg6/match_results.json"); H6 = load("agentB_deg6/H6_zeros.json")
D4 = load("quartic_pilot/match_results.json"); H4 = load("quartic_pilot/H4_zeros.json")
SPEC = {"d6": dict(delta=9.2962468327713700828344765663674575503066887803593,
                   lam2=0.3710714676819301324709834, lam3=-0.3419683225478021050998591,
                   A=1.3213197054428516594468381447179159106238901404509),
        "d4": dict(delta=7.2846862170733433643089305679955530694780466189573,
                   lam2=0.291838408983545104530996961457, lam3=-0.255664558474297406325347479489,
                   A=1.4904325463540699872206657264594941418825094844754)}

def build_cat(mr):
    rows=[]
    for r in mr["catalog_results"]:
        if not r.get("accepted"): continue
        if r.get("edge"): continue
        u_cat=complex(r["u_cat_re"], r["u_cat_im"])
        e=complex(r.get("u_actual_re",0.0), r.get("u_actual_im",0.0))-u_cat
        rows.append((r["id"], u_cat, e))
    return rows

cat6=build_cat(D6); cat4=build_cat(D4)

def oos_model(rows, sp, model):
    """model in {'const','u','u_del','u2_del','full'} ; 30/70 split, report within-2x."""
    rng=np.random.RandomState(20260815); n=len(rows)
    idx=np.arange(n); rng.shuffle(idx); ntr=int(0.30*n)
    tr,te=idx[:ntr],idx[ntr:]
    u=np.array([r[1] for r in rows]); e=np.array([r[2] for r in rows])
    delta=sp["delta"]
    if model=="const": X=np.ones((n,1))
    elif model=="u": X=np.column_stack([np.ones(n), u])
    elif model=="u_del": X=np.column_stack([np.ones(n), u, u*sp["lam2"]**9])
    elif model=="u2_del": X=np.column_stack([np.ones(n), u, u**2*delta**(-9)])
    else: X=np.column_stack([np.ones(n), u, u**2*delta**(-9)])
    coef,*_=np.linalg.lstsq(X[tr], e[tr], rcond=None)
    pred=X[te]@coef
    ratio=np.maximum(np.abs(pred-e[te]),1e-12)/np.maximum(np.abs(e[te]),1e-12)
    return np.median(ratio), np.mean(ratio<2.0), np.mean(ratio<10.0)

print("-- OOS baselines (within-2x / within-10x / median ratio) --")
for tag, rows, sp in [("d6",cat6,SPEC["d6"]),("d4",cat4,SPEC["d4"])]:
    for mdl in ["const","u","full"]:
        med, w2, w10 = oos_model(rows, sp, mdl)
        print("  %s %-6s: within2x=%.3f within10x=%.3f medratio=%.2f" % (tag, mdl, w2, w10, med))

print("\n-- phase of residual e9 vs u (d6, d4) --")
for tag, rows, sp in [("d6",cat6,SPEC["d6"]),("d4",cat4,SPEC["d4"])]:
    args=[]; rel_ang=[]
    for _id,u,e in rows:
        if abs(u)>1e-9 and abs(e)>1e-9:
            args.append(np.angle(u)); rel_ang.append((np.angle(e)-np.angle(u)+np.pi)%(2*np.pi)-np.pi)
    args=np.array(args); rel_ang=np.array(rel_ang)
    print("  %s: mean(rel angle e vs u)=%.1f deg  std=%.1f deg" %
          (tag, np.degrees(np.mean(rel_ang)), np.degrees(np.std(rel_ang))))
    # sign pattern: for lambda2>0 model e ~ (c0+c1*u), non-alternating; check Re(e/u)
    re_eu = np.array([(e/u).real for _id,u,e in rows if abs(e)>1e-9])
    print("   Re(e9/u): median=%.4f  frac>0=%.2f" % (np.median(re_eu), np.mean(re_eu>0)))

print("\n-- two-term rate fit per zero (rate_study d6 & d4) --")
def fit2(du, sp, m8=True):
    """e(m)=C*lam2^m + D*delta^{-m}; least squares on |e| with sign scan."""
    ms=[float(m) for m in du]; ys=np.array([du[str(int(m))] for m in ms])
    l2,de=sp["lam2"],sp["delta"]
    best=None
    for sgn in (1.0,-1.0):
        A=np.column_stack([sgn*l2**np.array(ms), de**(-np.array(ms))])
        coef,*_=np.linalg.lstsq(A,ys,rcond=None); pred=A@coef
        r=float(np.sum((pred-ys)**2))
        if best is None or r<best[0]: best=(r,np.abs(coef))
    return best[1]

for tag, rs, sp in [("d6",D6["rate_study"],SPEC["d6"]),("d4",D4["rate_study"],SPEC["d4"])]:
    print("  %s: C (lambda2-coef) and D (delta-coef) per zero" % tag)
    C,D,absu=[],[],[]
    for row in rs:
        du=row["d_u"]
        if len(du)<3: continue
        c2,de=fit2(du,sp)
        C.append(c2);D.append(de);absu.append(row["absu"])
        print("    %-5s |u|=%.3f  C=%.2e  D=%.2e  |C*lam2^9|=%.2e  |D*delta^-9|=%.2e" %
              (row["id"],row["absu"],c2,de,c2*abs(sp["lam2"])**9, de*sp["delta"]**(-9)))
    C=np.array(C);D=np.array(D);absu=np.array(absu)
    print("    corr(C, |u|)=%.3f  corr(D,|u|)=%.3f  corr(logC,log|u|)=%.3f" %
          (np.corrcoef(C,absu)[0,1], np.corrcoef(D,absu)[0,1],
           np.corrcoef(np.log(C+1e-12),np.log(absu))[0,1]))

print("\n-- d9 vs |u| log-log slope (d6, d4) --")
for tag, rows in [("d6",cat6),("d4",cat4)]:
    d9=np.array([abs(e) for _,_,e in rows]); u=np.array([abs(u) for _,u,_ in rows])
    sl=np.polyfit(np.log(u), np.log(d9), 1)[0]
    print("  %s: slope(log d9 vs log|u|)=%.3f" % (tag, sl))
    # quantiles of d9 by |u| band
    for lo,hi in [(0,1),(1,3),(3,10),(10,44)]:
        m=(u>=lo)&(u<hi)
        if m.sum()>0:
            print("     |u| in [%2d,%2d): n=%3d  d9 median=%.2e  min=%.2e  max=%.2e" %
                  (lo,hi,m.sum(),np.median(d9[m]),d9[m].min(),d9[m].max()))
print("\ndone stage 2")
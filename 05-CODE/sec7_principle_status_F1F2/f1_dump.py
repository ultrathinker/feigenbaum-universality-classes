#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F1: dump final result tables (predicted vs measured at m=9; rate fits) to JSON."""
import json, os
import numpy as np

BASE = r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип"
OUT = os.path.dirname(os.path.abspath(__file__))
def load(p):
    with open(os.path.join(BASE, p), encoding="utf-8") as f:
        return json.load(f)

D6 = load("agentB_deg6/match_results.json"); H6 = load("agentB_deg6/H6_zeros.json")
D4 = load("quartic_pilot/match_results.json"); H4 = load("quartic_pilot/H4_zeros.json")
RF = load("agentB_deg6/subfamily/rate_full.json")

SPEC = {"d6": dict(delta=9.2962468327713700828344765663674575503066887803593,
                   lam2=0.3710714676819301324709834, lam3=-0.3419683225478021050998591,
                   lam4=-0.2105185208142716405642349, A=1.3213197054428516594468381447179159106238901404509,
                   AU0=-0.6679695396324661394893857007410712354653, aF="1.6832601982046430005119730057818436560674199035301"),
        "d4": dict(delta=7.2846862170733433643089305679955530694780466189573,
                   lam2=0.291838408983545104530996961457, lam3=-0.255664558474297406325347479489,
                   lam4=-0.155562761843395028402820942804, A=1.4904325463540699872206657264594941418825094844754,
                   AU0=-0.6087812101890336104440950663154768942781, aF="1.594901356228820564497830482717353475129697571561")}

def build_cat(mr):
    rows=[]
    for r in mr["catalog_results"]:
        if not r.get("accepted") or r.get("edge"): continue
        u_cat=complex(r["u_cat_re"], r["u_cat_im"])
        u_act=complex(r.get("u_actual_re",0.0), r.get("u_actual_im",0.0))
        e=u_act-u_cat
        rows.append({"id":r["id"],"u_re":u_cat.real,"u_im":u_cat.imag,"absu":abs(u_cat),
                     "u_act_re":u_act.real,"u_act_im":u_act.imag,"d9":abs(e),
                     "e_re":e.real,"e_im":e.imag,"tier":r.get("tier")})
    return rows

def oos_full(rows, sp, frac=0.30, seed=20260815):
    rng=np.random.RandomState(seed); n=len(rows)
    idx=np.arange(n); rng.shuffle(idx); ntr=int(frac*n); tr,te=idx[:ntr],idx[ntr:]
    u=np.array([r["u_re"]+1j*r["u_im"] for r in rows]); e=np.array([r["e_re"]+1j*r["e_im"] for r in rows])
    X=np.column_stack([np.ones(n), u, u**2*sp["delta"]**(-9)])
    coef,*_=np.linalg.lstsq(X[tr],e[tr],rcond=None)
    pred=X[te]@coef; ratio=np.maximum(np.abs(pred-e[te]),1e-12)/np.maximum(np.abs(e[te]),1e-12)
    return dict(train=int(ntr),test=int(n-ntr),median_ratio=float(np.median(ratio)),
                within2=float(np.mean(ratio<2)), within10=float(np.mean(ratio<10)),
                coef=[complex(c).__str__() for c in coef])

out={"SPEC":{k:{kk:(float(vv) if not isinstance(vv,str) else vv) for kk,vv in v.items()} for k,v in SPEC.items()}}
out["catalog_d6"]=build_cat(D6)
out["catalog_d4"]=build_cat(D4)
out["rate_study_d6"]=D6["rate_study"]
out["rate_study_d4"]=D4["rate_study"]
out["rate_full_d6"]={k:{kk:vv for kk,vv in v.items()} for k,v in RF.items()}
out["oos"]={"d6":oos_full(out["catalog_d6"],SPEC["d6"]),"d4":oos_full(out["catalog_d4"],SPEC["d4"])}
with open(os.path.join(OUT,"f1f2_results.json"),"w",encoding="utf-8") as f:
    json.dump(out,f,indent=1,ensure_ascii=False)
print("wrote f1f2_results.json: d6=%d d4=%d rows" % (len(out["catalog_d6"]),len(out["catalog_d4"])))
print("OOS:", json.dumps(out["oos"]))
# quick table excerpt for the report
import cmath
print("\n-- sample: 12 rows of d6 predicted vs measured --")
for r in out["catalog_d6"][:12]:
    print("%s u=%+6.3f%+6.3fi d9=%.3e" % (r["id"], r["u_re"], r["u_im"], r["d9"]))
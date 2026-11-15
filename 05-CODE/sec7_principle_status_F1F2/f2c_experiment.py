#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F2(c): linearized (1st-order-in-s) catalog vs true nonlinear catalog.
Insight: to 1st order in s, x_n = x^{(0)}_n + s x^{(1)}_n with the SAME g-orbit
for all u. So H^{lin}_m(u) = (x0n + s*x1n)/(x0d + s*x1d), s = u*AU0/delta^m, is a
ratio of linear functions of u -> at most ONE zero. The true H_m has MANY.
NOT RIGOROUS (mpmath float).
"""
import json, os, cmath
import numpy as np
from mpmath import mp, mpf, mpc, fabs, nstr
mp.dps = 34

BASE = r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип"
OUT = os.path.dirname(os.path.abspath(__file__))
def load(p):
    with open(os.path.join(BASE, p), encoding="utf-8") as f:
        return json.load(f)

def horner(c, t):
    s = mpf(0)
    for k in range(len(c)-1, -1, -1):
        s = s*t + c[k]
    return s

class UVec:
    def __init__(self, fp, d):
        self.g=[mpf(x) for x in fp["g_coeffs"]]; self.h=[mpf(x) for x in fp["h_coeffs"]]
        self.d=d
    def gval(self,x):
        return horner(self.g, x**self.d)
    def gpx(self,x):
        gp=mpf(0); t=x**self.d
        for k in range(len(self.g)-1,0,-1):
            gp=gp*t + k*self.g[k]
        return gp*self.d*x**(self.d-1)
    def hval(self,x):
        return horner(self.h, x**self.d)

def build_orbits(uv, nmax):
    """x0n, x1n along the shared orbit; return arrays as python complex."""
    x0=mpf(0); x1=mpf(0)
    X0=[0.0]; X1=[0.0]
    for _ in range(nmax):
        x0n=uv.gval(x0); x1n=uv.gpx(x0)*x1 + uv.hval(x0)
        x0,x1=x0n,x1n
        X0.append(float(x0)); X1.append(float(x1))
    return np.array(X0), np.array(X1)

def Hlin_realaxis(X0n,X1n, X0d,X1d, delta, m, AU0):
    """Return u->H^{lin}_m(u) sampled on a dense real grid; u real."""
    s = np.linspace(0, 2, 40001)  # s>=0 side; scan s in [0,2] -> u in [0, 2*delta^m/|AU0|]
    num = X0n + s*X1n
    den = X0d + s*X1d
    H = num/den
    return s, H

def count_zero_crossings(H):
    cnt=0
    h=H
    for i in range(len(h)-1):
        if h[i]*h[i+1] <= 0:
            cnt+=1
    return cnt

for tag, fpname, hzname, AU0, delta, d in [
    ("d6","agentB_deg6/sextic_fixedpoint.json","agentB_deg6/H6_zeros.json",
     -0.6679695396324661394893857007410712354653, 9.2962468327713700828344765663674575503066887803593, 6),
    ("d4","quartic_pilot/quartic_fixedpoint.json","quartic_pilot/H4_zeros.json",
     -0.6087812101890336104440950663154768942781, 7.2846862170733433643089305679955530694780466189573, 4),
]:
    fp=load(fpname); uv=UVec(fp,d)
    hz=load(hzname)
    lim=hz["zeros_limit_extrapolated"]
    for m in (8,9):
        nmax=2**(m+1)
        X0,X1=build_orbits(uv,nmax)
        x0n, x1n = X0[nmax], X1[nmax]
        x0d, x1d = X0[nmax//2], X1[nmax//2]
        s, H = Hlin_realaxis(x0n,x1n, x0d,x1d, delta,m, mpf(AU0))
        ncross=count_zero_crossings(H)
        # true catalog count on real axis with |u|<=uMAX (s-scale)
        uMAX = 2*delta**m/abs(mpf(AU0))
        ntrue = sum(1 for z in lim if abs(z["im"])<1e-9 and abs(z["re"])<=uMAX)
        # single linearized zero (numerator=0)
        if abs(x1n)>0:
            slin = -x0n/x1n
            ulin = slin*delta**m/mpf(AU0)
        else:
            ulin=None
        print("%s m=%d: H^{lin}_m real-axis zeros=%d (true catalog real-axis zeros<=u=%.0f: %d); "
              "single linearized zero u^{lin}=%s" % (tag, m, ncross, uMAX, ntrue,
              (nstr(ulin,6) if ulin is not None else "n/a")))
        # compare linearized single zero to nearest true level-m zero
        if ulin is not None:
            ucat=complex(float(ulin))
            best=None
            for z in lim:
                for lvl in (m,):
                    k="u%d"%lvl
                    zz=complex(float(z[k+"_re"]),float(z[k+"_im"]))
                    dd=abs(zz-ucat)
                    if best is None or dd<best[0]:
                        best=(dd, (z.get("id",""), zz))
            print("   u^{lin}(%d)=%s ; dist to nearest true level%d zero = %.3e" %
                  (m, nstr(ulin,6), m, best[0]))
print("\ndone")
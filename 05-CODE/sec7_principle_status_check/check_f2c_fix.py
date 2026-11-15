#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Follow-up to check_f2c.py: the crossing-counter there counted sign flips of
H_lin = num/den, which conflates true zeros (num=0) with POLES (den=0,
where the ratio also flips sign). Since num and den are each exactly
linear in u (s is linear in u), each has exactly ONE root algebraically.
This script separately tracks num and den sign changes to see whether the
second "zero" found at u~-8.46 (d=6) / u~-6.11 (d=4) in check_f2c.py is a
real zero (num=0) or a pole (den=0) of H^lin.
NOT RIGOROUS (float64). <2 cores.
"""
import os, json
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
import numpy as np

BASE = r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип"
def load(p):
    with open(os.path.join(BASE, p), encoding="utf-8") as f:
        return json.load(f)

def horner_vec(coeffs, t):
    s = np.zeros_like(t)
    for k in range(len(coeffs)-1, -1, -1):
        s = s*t + coeffs[k]
    return s

def run(tag, fp_path, AU0, delta, d, m, u_lo, u_hi, n_u):
    fp = load(fp_path)
    g = [float(x) for x in fp["g_coeffs"]]
    h = [float(x) for x in fp["h_coeffs"]]
    u = np.linspace(u_lo, u_hi, n_u)
    s = u * AU0 / (delta**m)
    x0 = np.zeros_like(u); x1 = np.zeros_like(u)
    nmax = 2**(m+1); half = 2**m
    x_half_x0 = x_half_x1 = None
    for step in range(1, nmax+1):
        t0 = x0**d
        gval0 = horner_vec(g, t0); hval0 = horner_vec(h, t0)
        gprime_t = np.zeros_like(x0)
        for k in range(len(g)-1, 0, -1):
            gprime_t = gprime_t*t0 + k*g[k]
        dgdx = gprime_t * d * np.where(x0 != 0, x0**(d-1), 0.0)
        x1_new = dgdx*x1 + hval0
        x0, x1 = gval0, x1_new
        if step == half:
            x_half_x0, x_half_x1 = x0.copy(), x1.copy()
    num = x0 + s*x1              # numerator: x^(0)_{2^(m+1)} + s*x^(1)_{2^(m+1)}
    den = x_half_x0 + s*x_half_x1  # denominator

    def sign_changes(y):
        pts = []
        for i in range(len(y)-1):
            if y[i] == 0 or y[i]*y[i+1] < 0:
                pts.append(0.5*(u[i]+u[i+1]))
        return pts

    num_roots = sign_changes(num)
    den_roots = sign_changes(den)
    print("%s m=%d: NUM sign changes (true zeros of H^lin) at u=%s" % (tag, m, ["%.6f"%p for p in num_roots]))
    print("%s m=%d: DEN sign changes (POLES of H^lin)         at u=%s" % (tag, m, ["%.6f"%p for p in den_roots]))

for m in (5,6,7,8):
    run("d6", "agentB_deg6/sextic_fixedpoint.json",
        -0.6679695396324661394893857007410712354653,
        9.2962468327713700828344765663674575503066887803593, 6, m, -15, 15, 200001)
print()
for m in (5,6,7,8):
    run("d4", "quartic_pilot/quartic_fixedpoint.json",
        -0.6087812101890336104440950663154768942781,
        7.2846862170733433643089305679955530694780466189573, 4, m, -15, 15, 200001)

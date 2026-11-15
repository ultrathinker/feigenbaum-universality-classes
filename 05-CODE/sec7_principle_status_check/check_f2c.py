#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent reviewer (Sonnet) re-check of F2(c): does the linearized-in-s
catalog really collapse to ~1 zero on the real axis, versus many for the
true nonlinear H_m(u)?  Independent re-implementation (numpy float64,
vectorized over a real-u grid), NOT copied from f2c_experiment.py (which
uses mpmath + a single Newton-style algebraic solve for the linear zero).
This script directly counts sign changes of H_m(u) = x_{2^(m+1)}/x_{2^m}
on a dense real grid, for BOTH the true nonlinear map and its 1st-order
(linear-in-s) truncation, at a modest depth m (float64 stability limits
us to smaller m than the original NOT RIGOROUS-mpmath run; this is a
qualitative cross-check of the MULTIPLICITY claim, not a precision match).

NOT RIGOROUS (float64). <2 cores, read-only inputs, writes only here.
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
    """coeffs: python floats, low-to-high; t: numpy array. Returns sum coeffs[k]*t^k."""
    s = np.zeros_like(t)
    for k in range(len(coeffs)-1, -1, -1):
        s = s*t + coeffs[k]
    return s

def run_class(tag, fp_path, AU0, delta, d, m, u_lo, u_hi, n_u):
    fp = load(fp_path)
    g = [float(x) for x in fp["g_coeffs"]]
    h = [float(x) for x in fp["h_coeffs"]]

    u = np.linspace(u_lo, u_hi, n_u)
    s = u * AU0 / (delta**m)

    # nonlinear orbit: x_{k+1} = g(x_k) + s*h(x_k), vectorized over u (each u own s)
    x = np.zeros_like(u)
    xlin0 = np.zeros_like(u)  # 0th order (s=0) orbit, shared, but keep per-u for uniform code
    x1 = np.zeros_like(u)     # 1st-order-in-s response

    nmax = 2**(m+1)
    half = 2**m
    x_half_nl = None
    x_half_x0 = None
    x_half_x1 = None
    for step in range(1, nmax+1):
        t = x**d
        gval = horner_vec(g, t)
        hval = horner_vec(h, t)
        x_new = gval + s*hval
        # linear-order bookkeeping needs g'(x0)*x1 + h(x0); track separately using s=0 orbit
        t0 = xlin0**d
        gval0 = horner_vec(g, t0)
        hval0 = horner_vec(h, t0)
        # g'(x0): derivative of horner poly in t wrt x, chain rule via t=x^d
        gprime_t = np.zeros_like(xlin0)
        for k in range(len(g)-1, 0, -1):
            gprime_t = gprime_t*t0 + k*g[k]
        dgdx = gprime_t * d * np.where(xlin0 != 0, xlin0**(d-1), 0.0)
        x1_new = dgdx*x1 + hval0
        x, xlin0, x1 = x_new, gval0, x1_new
        if step == half:
            x_half_nl = x.copy(); x_half_x0 = xlin0.copy(); x_half_x1 = x1.copy()
    x_full_nl, x_full_x0, x_full_x1 = x, xlin0, x1

    Hnl = x_full_nl / x_half_nl
    num_lin = x_full_x0 + s*x_full_x1
    den_lin = x_half_x0 + s*x_half_x1
    Hlin = num_lin / den_lin

    def count_crossings(H):
        finite = np.isfinite(H)
        cnt = 0
        pts = []
        idxs = np.where(finite)[0]
        for i in range(len(idxs)-1):
            a, b = idxs[i], idxs[i+1]
            if b != a+1:
                continue
            if H[a] == 0 or (H[a]*H[b] < 0):
                cnt += 1
                pts.append(0.5*(u[a]+u[b]))
        return cnt, pts

    n_nl, pts_nl = count_crossings(Hnl)
    n_lin, pts_lin = count_crossings(Hlin)
    return dict(tag=tag, m=m, n_nonlinear_zeros=n_nl, n_linear_zeros=n_lin,
                nonlinear_zero_positions=pts_nl, linear_zero_positions=pts_lin)

out = []
def log(s):
    print(s)
    out.append(s)

log("=== Independent F2(c) re-check: multiplicity, linear vs nonlinear H_m(u) on real axis ===")
log("(own numpy float64 implementation, vectorized sign-crossing count; modest m for float64 stability)\n")

AU0_6 = -0.6679695396324661394893857007410712354653
delta6 = 9.2962468327713700828344765663674575503066887803593
AU0_4 = -0.6087812101890336104440950663154768942781
delta4 = 7.2846862170733433643089305679955530694780466189573

for m in (5, 6, 7):
    r = run_class("d6", "agentB_deg6/sextic_fixedpoint.json", AU0_6, delta6, 6, m,
                  u_lo=-15, u_hi=15, n_u=200001)
    log("d=6  m=%d  u in [-15,15]:  nonlinear real zeros=%d   linear(1st-order-in-s) real zeros=%d"
        % (m, r["n_nonlinear_zeros"], r["n_linear_zeros"]))
    if r["linear_zero_positions"]:
        log("     linear zero(s) at u ~ %s" % ["%.6f" % p for p in r["linear_zero_positions"]])
    if r["nonlinear_zero_positions"]:
        log("     nonlinear zero(s) at u ~ %s" % ["%.4f" % p for p in r["nonlinear_zero_positions"][:15]])

log("")
for m in (5, 6, 7):
    r = run_class("d4", "quartic_pilot/quartic_fixedpoint.json", AU0_4, delta4, 4, m,
                  u_lo=-15, u_hi=15, n_u=200001)
    log("d=4  m=%d  u in [-15,15]:  nonlinear real zeros=%d   linear(1st-order-in-s) real zeros=%d"
        % (m, r["n_nonlinear_zeros"], r["n_linear_zeros"]))
    if r["linear_zero_positions"]:
        log("     linear zero(s) at u ~ %s" % ["%.6f" % p for p in r["linear_zero_positions"]])
    if r["nonlinear_zero_positions"]:
        log("     nonlinear zero(s) at u ~ %s" % ["%.4f" % p for p in r["nonlinear_zero_positions"][:15]])

log("\ndone check_f2c.py")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_f2c_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent reviewer check (Sonnet), decisive test of the F1 dispute:
does the subleading residual scale as lambda2^m (DeepSeek/draft) or
(lambda2/delta)^m (brief's naive guess)?

Method: use the RATIO of residual d_u between consecutive levels m, m+1
for the SAME zero. This ratio isolates the m-scaling law independent of
any u-dependence assumption (no regression, no model form needed):
  if d_u(m) ~ C * r^m  then d_u(m)/d_u(m+1) -> 1/r
Compare observed ratio to 1/|lambda2| (their claim) vs delta/|lambda2|
(brief's naive (lambda2/delta)^m form).

Data (read-only): agentB_deg6/match_results.json (rate_study, 13 zeros,
levels 8/9/10), agentB_deg6/subfamily/rate_full.json (9 late zeros,
levels 7..12), quartic_pilot/match_results.json (rate_study, 13 zeros).

Writes results to this folder only. NOT RIGOROUS (float). <2s runtime.
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

SPEC = {
    "d6": dict(delta=9.2962468327713700828344765663674575503066887803593,
               lam2=0.3710714676819301324709834),
    "d4": dict(delta=7.2846862170733433643089305679955530694780466189573,
               lam2=0.291838408983545104530996961457),
}

out = []
def log(s):
    print(s)
    out.append(s)

log("=== Decisive ratio test: lambda2^m vs (lambda2/delta)^m ===\n")

for tag, path in [("d6", "agentB_deg6/match_results.json"), ("d4", "quartic_pilot/match_results.json")]:
    sp = SPEC[tag]
    inv_lam2 = 1.0 / abs(sp["lam2"])
    ratio_naive = sp["delta"] / abs(sp["lam2"])   # if d_u(m) ~ (lam2/delta)^m, ratio(m/m+1) -> delta/lam2
    D = load(path)
    rs = D["rate_study"]
    r89, r910 = [], []
    for row in rs:
        du = row["d_u"]
        ms = sorted(du.keys(), key=int)
        if len(ms) < 3:
            continue
        vals = [du[m] for m in ms]
        # ms are strings '8','9','10' typically
        r89.append(vals[0] / vals[1])
        r910.append(vals[1] / vals[2])
    r89 = np.array(r89); r910 = np.array(r910)
    log("[%s] rate_study n=%d" % (tag, len(r89)))
    log("  median ratio d(8)/d(9)  = %.4f" % np.median(r89))
    log("  median ratio d(9)/d(10) = %.4f" % np.median(r910))
    log("  claim 1/|lam2|          = %.4f  (lam2^m hypothesis)" % inv_lam2)
    log("  claim delta/|lam2|      = %.4f  ((lam2/delta)^m hypothesis, brief naive)" % ratio_naive)
    log("  |r910 - 1/lam2| = %.4f   |r910 - delta/lam2| = %.4f" %
        (abs(np.median(r910) - inv_lam2), abs(np.median(r910) - ratio_naive)))
    verdict = "lambda2^m WINS" if abs(np.median(r910) - inv_lam2) < abs(np.median(r910) - ratio_naive) else "(lambda2/delta)^m WINS"
    log("  -> %s (by ~%.1fx margin)\n" % (verdict, abs(np.median(r910)-ratio_naive)/max(abs(np.median(r910)-inv_lam2),1e-9)))

# rate_full (d=6, 9 late zeros, levels 7..12) - deeper levels, less m=8/9 interference
log("=== rate_full (d=6, 9 late zeros, levels up to 12) ===")
rf = load("agentB_deg6/subfamily/rate_full.json")
sp = SPEC["d6"]
inv_lam2 = 1.0/abs(sp["lam2"])
all_ratios = []
for k, v in rf.items():
    du = v["d_u"]
    ms = sorted(du.keys(), key=int)
    vals = [du[m] for m in ms]
    ratios = [vals[i]/vals[i+1] for i in range(len(vals)-1)]
    all_ratios.extend(ratios)
    log("  %-6s levels=%s ratios=%s" % (k, ms, ["%.3f"%r for r in ratios]))
all_ratios = np.array(all_ratios)
log("  median all consecutive-level ratios = %.4f  (1/|lam2|=%.4f, delta/|lam2|=%.4f)" %
    (np.median(all_ratios), inv_lam2, sp["delta"]/abs(sp["lam2"])))

log("\ndone check_rate.py")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_rate_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))

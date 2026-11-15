#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Independent reviewer (Sonnet) re-check of the F1 out-of-sample claim:
  "prediction of d_u(9) within factor 2, out-of-sample: 95% (d=6), 67% (d=4)"

Independently re-implemented (not copied) from raw match_results.json.
Adds a robustness sweep: many random 30/70 splits (not just seed 20260815)
plus k-fold cross-validation, to see whether the reported single-split
numbers are representative or a favorable draw.

NOT RIGOROUS (float). <2 cores, <1 min runtime, read-only on inputs,
writes only to this folder.
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
    "d6": dict(delta=9.2962468327713700828344765663674575503066887803593),
    "d4": dict(delta=7.2846862170733433643089305679955530694780466189573),
}

def build_catalog(mr):
    """Independent re-implementation: accepted, non-edge zeros, u_cat and residual e9=u_actual-u_cat."""
    rows = []
    for r in mr["catalog_results"]:
        if not r.get("accepted"):
            continue
        if r.get("edge"):
            continue
        u_cat = complex(r["u_cat_re"], r["u_cat_im"])
        u_act = complex(r.get("u_actual_re", 0.0), r.get("u_actual_im", 0.0))
        rows.append(dict(id=r["id"], u=u_cat, e=u_act - u_cat))
    return rows

def fit_predict(u_tr, e_tr, u_te, delta, with_delta=True):
    n_tr = len(u_tr)
    if with_delta:
        X = np.column_stack([np.ones(n_tr), u_tr, u_tr**2 * delta**(-9.0)])
    else:
        X = np.column_stack([np.ones(n_tr), u_tr])
    coef, *_ = np.linalg.lstsq(X, e_tr, rcond=None)
    if with_delta:
        Xte = np.column_stack([np.ones(len(u_te)), u_te, u_te**2 * delta**(-9.0)])
    else:
        Xte = np.column_stack([np.ones(len(u_te)), u_te])
    return Xte @ coef

def within2_frac(pred, truth):
    err = np.abs(pred - truth)
    ratio = np.maximum(err, 1e-14) / np.maximum(np.abs(truth), 1e-14)
    return np.mean(ratio < 2.0), np.median(ratio)

def repeated_split_test(rows, delta, frac=0.30, n_reps=500, with_delta=True, base_seed=1):
    u = np.array([r["u"] for r in rows])
    e = np.array([r["e"] for r in rows])
    n = len(rows)
    ntr = int(frac * n)
    w2s = []
    rng_master = np.random.RandomState(base_seed)
    for rep in range(n_reps):
        idx = np.arange(n)
        rs = np.random.RandomState(rng_master.randint(0, 2**31 - 1))
        rs.shuffle(idx)
        tr, te = idx[:ntr], idx[ntr:]
        pred = fit_predict(u[tr], e[tr], u[te], delta, with_delta=with_delta)
        w2, medr = within2_frac(pred, e[te])
        w2s.append(w2)
    w2s = np.array(w2s)
    return w2s

def kfold_test(rows, delta, k=5, with_delta=True, seed=1):
    u = np.array([r["u"] for r in rows])
    e = np.array([r["e"] for r in rows])
    n = len(rows)
    rng = np.random.RandomState(seed)
    idx = np.arange(n)
    rng.shuffle(idx)
    folds = np.array_split(idx, k)
    preds = np.zeros(n, dtype=complex)
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        preds[te] = fit_predict(u[tr], e[tr], u[te], delta, with_delta=with_delta)
    w2, medr = within2_frac(preds, e)
    return w2, medr

out = []
def log(s):
    print(s)
    out.append(s)

log("=== Independent OOS re-check of F1 within-2x claim (95%% d6 / 67%% d4) ===\n")

for tag, path in [("d6", "agentB_deg6/match_results.json"), ("d4", "quartic_pilot/match_results.json")]:
    D = load(path)
    rows = build_catalog(D)
    delta = SPEC[tag]["delta"]
    log("[%s] n=%d accepted interior zeros" % (tag, len(rows)))

    # 1) reproduce their exact split (seed=20260815, frac=0.30) for a sanity check
    u = np.array([r["u"] for r in rows]); e = np.array([r["e"] for r in rows])
    n = len(rows); ntr = int(0.30 * n)
    idx = np.arange(n); rng = np.random.RandomState(20260815); rng.shuffle(idx)
    tr, te = idx[:ntr], idx[ntr:]
    pred = fit_predict(u[tr], e[tr], u[te], delta, with_delta=True)
    w2, medr = within2_frac(pred, e[te])
    log("  reproduction of their seed=20260815 split: within2x=%.1f%%  median_ratio=%.2f  (train=%d test=%d)"
        % (100*w2, medr, ntr, n-ntr))

    # 2) robustness: 500 independent random 30/70 splits, different seeds
    w2s = repeated_split_test(rows, delta, frac=0.30, n_reps=500, with_delta=True, base_seed=777)
    log("  500 random 30/70 splits (fresh seeds): within2x median=%.1f%%  mean=%.1f%%  [p5=%.1f%%, p95=%.1f%%]  min=%.1f%% max=%.1f%%"
        % (100*np.median(w2s), 100*np.mean(w2s), 100*np.percentile(w2s,5), 100*np.percentile(w2s,95),
           100*np.min(w2s), 100*np.max(w2s)))

    # 3) 5-fold and 10-fold cross-validation (every point predicted exactly once, out-of-sample)
    w2_5, medr_5 = kfold_test(rows, delta, k=5, with_delta=True, seed=42)
    w2_10, medr_10 = kfold_test(rows, delta, k=10, with_delta=True, seed=42)
    log("  5-fold CV (all points OOS once): within2x=%.1f%%  median_ratio=%.2f" % (100*w2_5, medr_5))
    log("  10-fold CV (all points OOS once): within2x=%.1f%%  median_ratio=%.2f" % (100*w2_10, medr_10))

    # 4) same but WITHOUT the delta term (pure linear-in-u model) for comparison
    w2s_lin = repeated_split_test(rows, delta, frac=0.30, n_reps=500, with_delta=False, base_seed=778)
    log("  500 splits, LINEAR-ONLY model (c0+c1*u, no u^2*delta^-9 term): within2x median=%.1f%%  [p5=%.1f%%,p95=%.1f%%]"
        % (100*np.median(w2s_lin), 100*np.percentile(w2s_lin,5), 100*np.percentile(w2s_lin,95)))

    # 5) constant-only baseline (does the u-structure add anything beyond a global magnitude rescale?)
    def fit_predict_const(u_tr, e_tr, u_te):
        c0 = np.mean(e_tr)
        return np.full(len(u_te), c0, dtype=complex)
    w2c_list = []
    rng2 = np.random.RandomState(999)
    for rep in range(500):
        idx2 = np.arange(n); rs = np.random.RandomState(rng2.randint(0,2**31-1)); rs.shuffle(idx2)
        tr2, te2 = idx2[:ntr], idx2[ntr:]
        predc = fit_predict_const(u[tr2], e[tr2], u[te2])
        w2c, _ = within2_frac(predc, e[te2])
        w2c_list.append(w2c)
    w2c_arr = np.array(w2c_list)
    log("  500 splits, CONSTANT-ONLY baseline: within2x median=%.1f%%  [p5=%.1f%%,p95=%.1f%%]"
        % (100*np.median(w2c_arr), 100*np.percentile(w2c_arr,5), 100*np.percentile(w2c_arr,95)))
    log("")

log("done check_oos.py")
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_oos_output.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(out))

# -*- coding: utf-8 -*-
"""F10 — honest 1/16 lattice statistics."""
import os, sys, json, time, math, pathlib
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
from scipy.stats import gaussian_kde, kstest, mannwhitneyu

ROOT = pathlib.Path(r"..\..").resolve()
OUT  = pathlib.Path(r".")
LOG  = OUT / "progress.log"
def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%S')}  f10_stats  {msg}"
    with LOG.open("a", encoding="utf-8") as f:
        f.write(line + "\n"); f.flush()
    print(line)

log("imports OK")

# ---------- inputs ----------
H6_Z = json.loads((ROOT/"agentB_deg6/H6_zeros.json").read_text(encoding="utf-8"))
H6_M = json.loads((ROOT/"agentB_deg6/match_results.json").read_text(encoding="utf-8"))
H4_Z = json.loads((ROOT/"quartic_pilot/H4_zeros.json").read_text(encoding="utf-8"))
H4_M = json.loads((ROOT/"quartic_pilot/match_results.json").read_text(encoding="utf-8"))
LATE_C = json.loads((ROOT/"wave3_lateboxes/targets_late.json").read_text(encoding="utf-8"))
CERT   = json.loads((ROOT/"wave3_lateboxes/krawczyk6_late_certificates.json").read_text(encoding="utf-8"))
log("inputs loaded")

# ---------- helpers ----------
def frac16(x):
    f = (16.0 * x) % 1.0
    return np.minimum(f, 1.0 - f)

def lattice_distance(u_re, u_im):
    """Physical half-width distance to nearest 1/16 2D node."""
    p_re = frac16(u_re) / 16.0
    p_im = frac16(u_im) / 16.0
    return np.hypot(p_re, p_im)

def lattice_stats(u_re, u_im):
    d = lattice_distance(u_re, u_im)
    fr = frac16(u_re)
    fi = frac16(u_im)
    one_d = 1.0/16.0
    # per-coordinate fractions
    re_near = float(np.mean(fr < one_d))
    im_near = float(np.mean(fi < one_d))
    joint   = float(np.mean((fr < one_d) & (fi < one_d)))
    # 1/8 and 1/4
    one_8 = 1.0/8.0
    one_4 = 1.0/4.0
    re_8 = float(np.mean(((8.0*u_re) % 1.0 <= 1.0/8.0) | (((8.0*u_re) % 1.0) >= 1.0 - 1.0/8.0)))
    # careful: distance to 1/8 node = min(frac(8x), 1-frac(8x))/8 ; tolerance 1/8 ⇒ frac in [0, 1/8]∪[7/8, 1]
    def frac8(x):
        g = (8.0*x) % 1.0
        return np.minimum(g, 1.0 - g)
    def frac4(x):
        g = (4.0*x) % 1.0
        return np.minimum(g, 1.0 - g)
    re_8 = float(np.mean(frac8(u_re) < 1.0/8.0))
    im_8 = float(np.mean(frac8(u_im) < 1.0/8.0))
    re_4 = float(np.mean(frac4(u_re) < 1.0/4.0))
    im_4 = float(np.mean(frac4(u_im) < 1.0/4.0))
    # line-vs-point: at least one coordinate near 1/16 node
    line_or = float(np.mean((fr < one_d) | (fi < one_d)))
    # Rayleigh per coordinate at q=16
    def rayleigh(x, q):
        ang = 2.0*np.pi*q*x
        R = abs(np.exp(1j*ang).sum()) / len(x)
        z = (R*len(x))**2 / len(x)
        return float(R), float(z), float(np.exp(-z))
    Rr, zr, pr = rayleigh(u_re, 16.0)
    Ri, zi, pi = rayleigh(u_im, 16.0)
    return {
        "n": int(len(u_re)),
        "re_near": re_near, "im_near": im_near, "joint": joint,
        "re_near_8": re_8, "im_near_8": im_8,
        "re_near_4": re_4, "im_near_4": im_4,
        "line_or": line_or,
        "dists": {"min": float(d.min()),
                  "q01": float(np.quantile(d, 0.01)),
                  "q05": float(np.quantile(d, 0.05)),
                  "q10": float(np.quantile(d, 0.10)),
                  "q25": float(np.quantile(d, 0.25)),
                  "median": float(np.median(d)),
                  "mean": float(d.mean()),
                  "q75": float(np.quantile(d, 0.75)),
                  "q90": float(np.quantile(d, 0.90)),
                  "q95": float(np.quantile(d, 0.95)),
                  "q99": float(np.quantile(d, 0.99)),
                  "max": float(d.max())},
        "n_below_5e-4":   int((d < 5.0e-4).sum()),
        "n_below_1/256":  int((d < 1.0/256.0).sum()),
        "n_below_4.4e-4": int((d < 4.4e-4).sum()),
        "rayleigh_re": {"R": Rr, "z": zr, "p": pr},
        "rayleigh_im": {"R": Ri, "z": zi, "p": pi},
    }

# ---------- extract sample sets ----------
def to_xy(rows, re_key, im_key):
    return np.array([r[re_key] for r in rows], float), np.array([r[im_key] for r in rows], float)

# H6 zeros
def zeros_xy(zlist):
    return np.array([z["re"] for z in zlist], float), np.array([z["im"] for z in zlist], float)

H6_lvl8 = zeros_xy(H6_Z["zeros_level8"])
H6_lvl7 = zeros_xy(H6_Z["zeros_level7"])
H6_extrap = zeros_xy(H6_Z["zeros_limit_extrapolated"])

cat_rows = H6_M["catalog_results"]
def filter_cat(rows, accepted=None, edge=None):
    out = []
    for r in rows:
        if accepted is not None and r.get("accepted") != accepted: continue
        if edge is not None and bool(r.get("edge", False)) != edge: continue
        out.append(r)
    return out

good_h6 = filter_cat(cat_rows, accepted=True, edge=False)
good_h6_re, good_h6_im = to_xy(good_h6, "u_actual_re", "u_actual_im")
good_h6_cat_re, good_h6_cat_im = to_xy(good_h6, "u_cat_re", "u_cat_im")

late_ids = {"S031","S036","S070","S104","S187","S235","S279","S319","S339"}
late_rows = [r for r in cat_rows if r.get("id") in late_ids]
late_cat_re, late_cat_im = to_xy(late_rows, "u_cat_re", "u_cat_im")

# Certificate-derived u_actual
A6 = float(H6_M["A6"]); aF = float(H6_M["a_F"]); dlt = float(H6_M["delta6"])
cert_u = []
for box in CERT["roots"]:
    m = int(box["m"])
    c_re = float(box["box_center_re"]); c_im = float(box["box_center_im"])
    u = A6 * ((c_re - aF) + 1j*c_im) * dlt**m
    cert_u.append((box["label"], u.real, u.imag, m))
cert_u.sort()
cert_re = np.array([u[1] for u in cert_u], float)
cert_im = np.array([u[2] for u in cert_u], float)
# late u_actual from match_results.json is not present (unaccepted); use certificate
late_act_re, late_act_im = cert_re, cert_im

# H4
H4_lvl8 = zeros_xy(H4_Z["zeros_level8"])
H4_lvl7 = zeros_xy(H4_Z["zeros_level7"])
H4_extrap = zeros_xy(H4_Z["zeros_limit_extrapolated"])
h4_rows = H4_M["catalog_results"]
h4_good = filter_cat(h4_rows, accepted=True, edge=False)
h4_good_re, h4_good_im = to_xy(h4_good, "u_actual_re", "u_actual_im")
h4_good_cat_re, h4_good_cat_im = to_xy(h4_good, "u_cat_re", "u_cat_im")

log(f"H6 samples: lvl8={len(H6_lvl8[0])}, extrap={len(H6_extrap[0])}, good_actual={len(good_h6_re)}, late={len(late_rows)}, cert={len(cert_re)}")
log(f"H4 samples: lvl8={len(H4_lvl8[0])}, extrap={len(H4_extrap[0])}, good_actual={len(h4_good_re)}")

samples = {
    "H6_lvl8": H6_lvl8,
    "H6_extrap": H6_extrap,
    "H6_287_actual": (good_h6_re, good_h6_im),
    "H6_287_cat": (good_h6_cat_re, good_h6_cat_im),
    "H6_late_cat": (late_cat_re, late_cat_im),
    "H6_late_actual": (late_act_re, late_act_im),
    "H6_late_cert": (cert_re, cert_im),
    "H4_lvl8": H4_lvl8,
    "H4_extrap": H4_extrap,
    "H4_61_actual": (h4_good_re, h4_good_im),
    "H4_61_cat": (h4_good_cat_re, h4_good_cat_im),
}

results = {name: lattice_stats(re, im) for name, (re, im) in samples.items()}
log("per-set lattice stats computed")

# ---------- save cumulative ----------
with (OUT/"f10_results.json").open("w", encoding="utf-8") as f:
    json.dump({"per_set": results}, f, indent=2)

# ---------- radial KDE null ----------
K = 6000
def sample_radial_kde(radii, n, bw=2.0, rng=None):
    if rng is None:
        rng = np.random.default_rng(0xC0FFEE)
    kde = gaussian_kde(radii[:, None].T, bw_method=bw / float(np.std(radii)))
    seed_int = int(rng.integers(0, 2**31 - 1))
    np.random.seed(seed_int)
    new_r = kde.resample(n).ravel()
    new_r = np.clip(new_r, 0.0, None)
    theta = rng.uniform(0, 2*np.pi, n)
    return new_r * np.cos(theta), new_r * np.sin(theta)

def simulate_null(re, im, seed=11):
    r = np.hypot(re, im)
    n = len(re)
    s_re = np.empty((K, n)); s_im = np.empty((K, n))
    base = 0
    chunk = 1000
    rng = np.random.default_rng(seed)
    while base < K:
        c = min(chunk, K-base)
        sr, si = sample_radial_kde(r, n*c, 2.0, rng)
        sr = sr.reshape(c, n); si = si.reshape(c, n)
        s_re[base:base+c] = sr; s_im[base:base+c] = si
        base += c
    return s_re, s_im

def mc_pvalue(observed, sims, direction="excess"):
    sims = np.asarray(sims)
    if direction == "excess":
        return float((np.sum(sims >= observed) + 1) / (sims.size + 1))
    if direction == "deficit":
        return float((np.sum(sims <= observed) + 1) / (sims.size + 1))
    return float((np.sum(np.abs(sims - 0.5) >= np.abs(observed - 0.5)) + 1) / (sims.size + 1))

# Helper: derive per-set null for each statistic
def kde_null_pvals(re, im, observed, key_name):
    s_re, s_im = simulate_null(re, im)
    seen = []
    # rates per coordinate / joint
    fr_obs = frac16(re); fi_obs = frac16(im); one_d = 1.0/16.0
    re_n = np.mean((s_re*16.0)%1.0 <= 1.0/16.0, axis=1)
    # use min(frac, 1-frac) < 1/16
    def per_coord(s):
        f = (16.0*s) % 1.0
        return np.mean(np.minimum(f, 1.0-f) < 1.0/16.0, axis=1)
    re_n = per_coord(s_re); im_n = per_coord(s_im)
    # joint: per-point indicator both < 1/16, then average
    def per_indicator(s):
        f = (16.0*s) % 1.0
        return np.minimum(f, 1.0-f) < 1.0/16.0
    joint_n = np.mean(per_indicator(s_re) & per_indicator(s_im), axis=1)
    d = lattice_distance(re, im)
    n_below = (lattice_distance(s_re, s_im) < d.max()*0 + 5.0e-4).sum(axis=1)
    n_below_1256 = (lattice_distance(s_re, s_im) < 1.0/256.0).sum(axis=1)
    med = np.median(lattice_distance(s_re, s_im), axis=1)
    # quantile-based
    q = observed["dists"]
    q01_n = np.quantile(lattice_distance(s_re, s_im), 0.01, axis=1)
    q05_n = np.quantile(lattice_distance(s_re, s_im), 0.05, axis=1)
    q10_n = np.quantile(lattice_distance(s_re, s_im), 0.10, axis=1)
    q25_n = np.quantile(lattice_distance(s_re, s_im), 0.25, axis=1)
    p_re = mc_pvalue(observed["re_near"], re_n)
    p_im = mc_pvalue(observed["im_near"], im_n)
    p_joint = mc_pvalue(observed["joint"], joint_n)
    p_n5e4 = mc_pvalue(observed["n_below_5e-4"], n_below)
    p_n1256 = mc_pvalue(observed["n_below_1/256"], n_below_1256)
    p_med = mc_pvalue(observed["dists"]["median"], med, "deficit")
    p_q05 = mc_pvalue(observed["dists"]["q05"], q05_n, "deficit")
    p_q01 = mc_pvalue(observed["dists"]["q01"], q01_n, "deficit")
    p_q10 = mc_pvalue(observed["dists"]["q10"], q10_n, "deficit")
    p_q25 = mc_pvalue(observed["dists"]["q25"], q25_n, "deficit")
    return {
        "p_re_near": p_re, "p_im_near": p_im, "p_joint": p_joint,
        "p_n_below_5e-4": p_n5e4, "p_n_below_1/256": p_n1256,
        "p_median": p_med, "p_q05": p_q05, "p_q01": p_q01, "p_q10": p_q10, "p_q25": p_q25,
        "null_re_near_mean": float(re_n.mean()), "null_re_near_std": float(re_n.std()),
        "null_im_near_mean": float(im_n.mean()), "null_im_near_mean": float(im_n.mean()),
        "null_joint_mean": float(joint_n.mean()),
        "null_n_below_5e-4_mean": float(n_below.mean()),
        "null_n_below_1/256_mean": float(n_below_1256.mean()),
        "null_median_mean": float(med.mean()),
    }

needs_null = ["H6_extrap", "H6_287_actual", "H6_late_cat", "H6_late_actual", "H6_late_cert",
              "H4_extrap", "H4_61_actual"]
null_results = {}
for name in needs_null:
    re, im = samples[name]
    obs = results[name]
    log(f"KDE null: {name} n={len(re)}")
    null_results[name] = kde_null_pvals(re, im, obs, name)

results_with_null = {n: {**results[n], "kde_null": null_results[n]} if n in null_results else results[n]
                     for n in results}

# ---------- shift controls ----------
def shift_test(re, im, shift, n_trials=2000, seed=99):
    rng = np.random.default_rng(seed)
    base = frac16(re); base_im = frac16(im)
    joint = float(np.mean((base < 1.0/16.0) & (base_im < 1.0/16.0)))
    sims = np.empty(n_trials)
    for i in range(n_trials):
        dx = rng.uniform(-shift, shift)
        dy = rng.uniform(-shift, shift)
        sr = re + dx; si = im + dy
        fr = frac16(sr); fi = frac16(si)
        sims[i] = np.mean((fr < 1.0/16.0) & (fi < 1.0/16.0))
    p = mc_pvalue(joint, sims, "excess")
    return joint, float(sims.mean()), float(sims.std()), p

re, im = samples["H6_287_actual"]
sub_joint, sub_mean, sub_std, sub_p = shift_test(re, im, 1.0/32.0, n_trials=2000)
log(f"random sub-cell shift (H6 287): joint obs={sub_joint:.4f} null={sub_mean:.4f}±{sub_std:.4f} p={sub_p:.4f}")

# Integer-lattice shift: add (k/16, l/16) — should leave lattice stats unchanged
def int_shift_test(re, im, n_trials=2000, seed=7):
    rng = np.random.default_rng(seed)
    base = frac16(re); base_im = frac16(im)
    joint = float(np.mean((base < 1.0/16.0) & (base_im < 1.0/16.0)))
    sims = np.empty(n_trials)
    for i in range(n_trials):
        k = int(rng.integers(-3, 4)); l = int(rng.integers(-3, 4))
        sr = re + k/16.0; si = im + l/16.0
        fr = frac16(sr); fi = frac16(si)
        sims[i] = np.mean((fr < 1.0/16.0) & (fi < 1.0/16.0))
    return joint, float(sims.mean()), float(sims.std())

joint_int, int_mean, int_std = int_shift_test(re, im)
log(f"integer-lattice shift (H6 287): joint obs={joint_int:.4f} null={int_mean:.4f}±{int_std:.4f}")

results_with_null["shifts_control"] = {
    "sub_cell_shift_1/32": {"joint": sub_joint, "null_mean": sub_mean, "null_std": sub_std, "p": sub_p},
    "integer_shift_1/16": {"joint": joint_int, "null_mean": int_mean, "null_std": int_std},
}

# ---------- Rayleigh curves ----------
def rayleigh_curve(re, im, qs):
    re = np.asarray(re, float).ravel()
    im = np.asarray(im, float).ravel()
    n = len(re)
    out = []
    for q in qs:
        ang_r = 2.0*np.pi*q*re
        ang_i = 2.0*np.pi*q*im
        esum_r = np.exp(1j*ang_r).sum()
        esum_i = np.exp(1j*ang_i).sum()
        Rr = abs(esum_r) / n
        Ri = abs(esum_i) / n
        C = np.cos(2*np.pi*q*re) + np.cos(2*np.pi*q*im)
        S = np.sin(2*np.pi*q*re) + np.sin(2*np.pi*q*im)
        R2 = float(np.hypot(C.sum(), S.sum()) / (2*n))
        Rr = float(Rr); Ri = float(Ri)
        z_r = (Rr*n)**2 / n; z_i = (Ri*n)**2 / n
        z_2 = (R2*2*n)**2 / (2*n)
        out.append({
            "q": float(q),
            "R_re": Rr, "R_im": Ri, "R_combined": R2,
            "p_re": float(np.exp(-z_r)),
            "p_im": float(np.exp(-z_i)),
            "p_combined": float(np.exp(-z_2)),
        })
    return out

qs = np.arange(0.5, 40.0+1e-9, 0.05)
qs_special = [4, 8, 12, 16, 20, 24, 32]
rayleigh = {}
for name in ["H6_extrap", "H6_287_actual", "H6_late_cat", "H4_extrap", "H4_61_actual"]:
    re, im = samples[name]
    full = rayleigh_curve(re, im, qs)
    spec = rayleigh_curve(re, im, qs_special)
    rayleigh[name] = {"full": full, "special": spec}
    log(f"rayleigh: {name} done ({len(full)} qs)")

# Find top-10 peaks in combined R for H6_extrap
def top_peaks(curve, k=10):
    rs = np.array([c["R_combined"] for c in curve])
    qs = np.array([c["q"] for c in curve])
    peaks = []
    for i in range(1, len(rs)-1):
        if rs[i] > rs[i-1] and rs[i] >= rs[i+1]:
            peaks.append((i, rs[i], qs[i]))
    peaks.sort(key=lambda t: -t[1])
    return [{"q": float(q), "R_combined": float(r), "p_combined": next(c["p_combined"] for c in curve if c["q"]==q)} for _, r, q in peaks[:k]]

top_peaks_table = {name: top_peaks(rayleigh[name]["full"], 10) for name in rayleigh}

# ---------- ring decomposition ----------
def ring_stats(re, im, rings):
    out = []
    r = np.hypot(re, im)
    for (lo, hi) in rings:
        mask = (r >= lo) & (r < hi)
        if mask.sum() == 0:
            out.append({"radius": [lo, hi], "n": 0})
            continue
        sub = lattice_stats(re[mask], im[mask])
        out.append({"radius": [lo, hi], **sub})
    return out

rings = [(0,8),(8,16),(16,24),(24,32),(32,48),(48,1e9)]
ring_data = {}
for name in ["H6_extrap", "H6_287_actual", "H4_extrap"]:
    ring_data[name] = ring_stats(samples[name][0], samples[name][1], rings)

# ---------- line vs point: among joint and per-coord ----------
def line_vs_point(re, im):
    fr = frac16(re); fi = frac16(im)
    tol = 1.0/16.0
    only_re = (fr < tol) & (fi >= tol)
    only_im = (fr >= tol) & (fi < tol)
    both = (fr < tol) & (fi < tol)
    n = len(re)
    return {
        "re_only": float(only_re.mean()),
        "im_only": float(only_im.mean()),
        "both": float(both.mean()),
        "either": float((only_re | only_im | both).mean()),
        "expected_per_coord": 0.125,
        "expected_joint": 0.015625,
        "expected_either": 2*0.125 - 0.015625,
    }

lvp = {name: line_vs_point(samples[name][0], samples[name][1]) for name in samples}

# ---------- late uniform-cell probability ----------
r_late = 4.4e-4
one_p = (np.pi * r_late**2) / (1.0/16.0)**2  # circle cell radius r vs square of side 1/16
# paper uses square radius; use that
one_p_sq = (2*r_late)**2 / (1.0/16.0)**2
p_nine = one_p_sq ** 9
log(f"late uniform cell p (square radius {r_late}): single={one_p_sq:.6e}, 9-points={p_nine:.6e}")

# KDE null probability for 9 late points: fraction of 9-point K-simulations with all 9 < 4.4e-4
def kde_late_9pt_p(re, im, r_thresh=4.4e-4, n_trials=200000, seed=23):
    rng = np.random.default_rng(seed)
    r = np.hypot(re, im)
    n = len(re)
    kde = gaussian_kde(r[:, None].T, bw_method=2.0/np.std(r))
    # 9-point block
    sims = n_trials // 9
    blocks = max(1, sims)
    np.random.seed(int(rng.integers(0, 2**31 - 1)))
    radii = kde.resample(9*blocks).ravel()
    radii = np.clip(radii, 0, None)
    theta = rng.uniform(0, 2*np.pi, 9*blocks)
    xs = radii * np.cos(theta); ys = radii * np.sin(theta)
    xs = xs.reshape(blocks, 9); ys = ys.reshape(blocks, 9)
    d = lattice_distance(xs, ys)
    max_per_block = d.max(axis=1)
    return float(np.mean(max_per_block < r_thresh))

log("computing late 9-point KDE null probability ...")
late_9pt_p = kde_late_9pt_p(cert_re, cert_im, r_thresh=4.4e-4, n_trials=200000)
log(f"late 9-point KDE null p (all < 4.4e-4): {late_9pt_p:.6e}")

# ---------- good tail comparison: do good points have a near-node excess? ----------
def good_tail_p(re, im, d_thresh=5.0e-4, n_trials=200000, seed=31):
    rng = np.random.default_rng(seed)
    r = np.hypot(re, im)
    n = len(re)
    kde = gaussian_kde(r[:, None].T, bw_method=2.0/np.std(r))
    sims = n_trials // n
    blocks = max(1, sims)
    np.random.seed(int(rng.integers(0, 2**31 - 1)))
    radii = kde.resample(n*blocks).ravel()
    radii = np.clip(radii, 0, None)
    theta = rng.uniform(0, 2*np.pi, n*blocks)
    xs = radii * np.cos(theta); ys = radii * np.sin(theta)
    xs = xs.reshape(blocks, n); ys = ys.reshape(blocks, n)
    d = lattice_distance(xs, ys)
    counts = (d < d_thresh).sum(axis=1)
    obs = int((lattice_distance(re, im) < d_thresh).sum())
    return obs, float(counts.mean()), float((counts >= obs).mean() + 1.0/(counts.size+1))

log("computing good-tail KDE null for H6 287 ...")
good_tail = good_tail_p(good_h6_re, good_h6_im, d_thresh=5.0e-4, n_trials=200000)
log(f"good 287 tail: obs={good_tail[0]} mean={good_tail[1]:.3f} p={good_tail[2]:.4e}")

# ---------- bimodality test on log-distance for good sample ----------
def bimodality(re, im):
    d = lattice_distance(re, im)
    d = d[d > 0]
    z = np.log10(d)
    # 1- and 2-component Gaussian mixture on log10(d) via EM-style moment fit
    mu1 = z.mean()
    sd1 = z.std()
    # 2-component: split at 0.5 percentile vs rest
    thr = np.quantile(z, 0.05)
    near = z[z < thr]; far = z[z >= thr]
    if len(near) < 5 or len(far) < 5:
        return {"mu1": float(mu1), "sd1": float(sd1), "two_modes": None}
    mu_n = near.mean(); sd_n = near.std()
    mu_f = far.mean(); sd_f = far.std()
    # BIC: 2-component mixture log-likelihood under diagonal Gaussian
    def loglik(z, params):
        m1, s1, m2, s2, w = params
        p1 = w * np.exp(-(z-m1)**2/(2*s1**2))/np.sqrt(2*np.pi)/s1
        p2 = (1-w) * np.exp(-(z-m2)**2/(2*s2**2))/np.sqrt(2*np.pi)/s2
        return float(np.log(p1+p2+1e-300).sum())
    # 1-comp LL
    ll1 = loglik(z, (mu1, sd1, mu1, sd1, 1.0))
    # 2-comp LL (grid)
    from itertools import product
    best = (-1e18, None)
    for w in [0.02, 0.05, 0.10, 0.15, 0.20]:
        for sn in [0.3, 0.5, 0.7, 1.0]:
            best = max(best, (loglik(z, (mu_n, sn, mu_f, sd_f, w)),
                              (mu_n, sn, mu_f, sd_f, w)))
    ll2 = best[0]
    n = len(z); k2 = 5; k1 = 2
    bic2 = -2*ll2 + k2*np.log(n); bic1 = -2*ll1 + k1*np.log(n)
    return {"mu1": float(mu1), "sd1": float(sd1),
            "two_modes": {"near_n": int(len(near)), "far_n": int(len(far)),
                          "mu_near": float(mu_n), "sd_near": float(sd_n),
                          "mu_far": float(mu_f), "sd_far": float(sd_f),
                          "loglik_1": float(ll1), "loglik_2": float(ll2),
                          "BIC_1": float(bic1), "BIC_2": float(bic2),
                          "delta_BIC": float(bic1 - bic2)}}

bim = {}
for name in ["H6_287_actual", "H4_61_actual"]:
    bim[name] = bimodality(samples[name][0], samples[name][1])
log("bimodality tests done")

# ---------- H4 late / residue ---
# No explicit late unaccepted in H4; only Q57 rejected. Q57 row:
q57 = [r for r in h4_rows if r.get("id") == "Q57" or r.get("accepted") is False]
log(f"H4 rejected: {len(q57)}")

# ---------- finalize JSON ----------
final = {
    "per_set": {k: {kk: vv for kk, vv in v.items() if kk != "dists" or True} for k, v in results.items()},
    "nulls": null_results,
    "shifts": results_with_null.pop("shifts_control"),
    "rayleigh": {k: {"special": v["special"],
                     "top10_peaks": top_peaks_table[k]} for k, v in rayleigh.items()},
    "rings": ring_data,
    "line_vs_point": lvp,
    "late_uniform_cell": {"single_pt_p": one_p_sq, "nine_pt_p": p_nine, "r_thresh": r_late},
    "late_kde_p_all_below_4.4e-4": late_9pt_p,
    "good_tail_test": {"obs_n_below_5e-4": good_tail[0], "null_mean": good_tail[1], "p": good_tail[2]},
    "bimodality": bim,
    "ring_counts_H6_extrap": {
        "[0,8)": int(((np.hypot(*H6_extrap)>=0)&(np.hypot(*H6_extrap)<8)).sum()),
        "[8,16)": int(((np.hypot(*H6_extrap)>=8)&(np.hypot(*H6_extrap)<16)).sum()),
        "[16,24)": int(((np.hypot(*H6_extrap)>=16)&(np.hypot(*H6_extrap)<24)).sum()),
        "[24,32)": int(((np.hypot(*H6_extrap)>=24)&(np.hypot(*H6_extrap)<32)).sum()),
        "[32,48)": int(((np.hypot(*H6_extrap)>=32)&(np.hypot(*H6_extrap)<48)).sum()),
        ">=48": int((np.hypot(*H6_extrap)>=48).sum())
    }
}

with (OUT/"f10_results.json").open("w", encoding="utf-8") as f:
    json.dump(final, f, indent=2, ensure_ascii=False)
log("f10_results.json written")

# ---------- plot Rayleigh ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# validated palette (light)
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
plt.rcParams.update({
    "font.family": "system-ui, Segoe UI, sans-serif",
    "axes.axisbelow": True,
    "axes.edgecolor": "#c3c2b7",
    "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.grid": True, "grid.color": "#e1e0d9", "grid.linewidth": 1.0,
    "axes.titlecolor": "#0b0b0b",
})

fig, ax = plt.subplots(figsize=(8.8, 4.6), dpi=140)
fig.patch.set_facecolor("#fcfcfb"); ax.set_facecolor("#fcfcfb")
names = ["H6_extrap", "H6_287_actual", "H4_extrap", "H4_61_actual"]
labels = {"H6_extrap":"H6 342 экстрап.","H6_287_actual":"H6 287 actual",
          "H4_extrap":"H4 62 экстрап.","H4_61_actual":"H4 61 actual"}
for k, name in enumerate(names):
    curve = rayleigh[name]["full"]
    qs = np.array([c["q"] for c in curve])
    R  = np.array([c["R_combined"] for c in curve])
    ax.plot(qs, R, color=PAL[k], linewidth=2.0, label=labels[name])
# mark q=16
ax.axvline(16.0, color="#52514e", linewidth=1.0, linestyle="--")
ax.axvline(8.0, color="#898781", linewidth=1.0, linestyle=":")
ax.set_xlim(0.5, 40.0)
ax.set_xlabel("q")
ax.set_ylabel("R(q) = |Σ exp(2πi·q·u)| / (2n)")
ax.set_title("Совмещённый Рэлей (Re+Im) по 1/16-решётке")
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(OUT/"fig_rayleigh_q.png", dpi=160)
plt.close(fig)
log("fig_rayleigh_q.png written")
log("done")

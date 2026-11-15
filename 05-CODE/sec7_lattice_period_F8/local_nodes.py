#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 Task 5: Local structure near k/16 lattice nodes.

For 4 late targets: compute |H(node)| and |H'(node)| (central diff).
For 20 random k/16 nodes WITHOUT nearby catalog zero: distribution of |H(node)|.
Compare to |H| at random non-node points.  Is the lattice a set of near-zeros?
"""
import os
import json
import time
import math
import numpy as np

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'
T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def horner_d(z, coeffs, d):
    z = np.asarray(z, dtype=np.complex128)
    t = z ** d
    r = np.full(z.shape, float(coeffs[-1]), dtype=np.float64).astype(np.complex128)
    for c in reversed(coeffs[:-1]):
        r = r * t + complex(c)
    return r


def H_at(u_arr, d, g_np, h_np, AU0, delta, m=8, TMAX=1.15, DEN_MIN=1e-12):
    u_arr = np.asarray(u_arr, dtype=np.complex128)
    s_coef = float(AU0) / (float(delta) ** m)
    s = s_coef * u_arr
    n1 = 2 ** m; n2 = 2 ** (m + 1)
    x = np.zeros(u_arr.shape, dtype=np.complex128)
    maxt = np.zeros(u_arr.shape, dtype=np.float64)
    x2m = np.zeros(u_arr.shape, dtype=np.complex128)
    for k in range(n2):
        t_d = x ** d
        x = horner_d(x, g_np, d) + s * horner_d(x, h_np, d)
        np.maximum(maxt, np.abs(t_d), out=maxt)
        if k + 1 == n1:
            x2m = x.copy()
    Habs = np.abs(x) / np.abs(x2m)
    valid = (np.isfinite(Habs) & (np.abs(x2m) > DEN_MIN)
             & (maxt < TMAX) & (np.abs(x) < 1e100))
    return Habs, valid, x


def nearest_cat_dist(re, im, cat_re, cat_im):
    d2 = (cat_re - re) ** 2 + (cat_im - im) ** 2
    return float(np.sqrt(d2.min()))


def main():
    log('--- local_nodes: load ---')
    fp6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'sextic_fixedpoint.json'), encoding='utf-8'))
    hz6 = json.load(open(os.path.join(ROOT_IN, 'agentB_deg6', 'H6_zeros.json'), encoding='utf-8'))
    fp4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'quartic_fixedpoint.json'), encoding='utf-8'))
    hz4 = json.load(open(os.path.join(ROOT_IN, 'quartic_pilot', 'H4_zeros.json'), encoding='utf-8'))

    g6 = np.array(fp6['g_coeffs'], dtype=np.float64)
    h6 = np.array(fp6['h_coeffs'], dtype=np.float64)
    g4 = np.array(fp4['g_coeffs'], dtype=np.float64)
    h4 = np.array(fp4['h_coeffs'], dtype=np.float64)

    cat6_re = np.array([z['re'] for z in hz6['zeros_level8']], dtype=np.float64)
    cat6_im = np.array([z['im'] for z in hz6['zeros_level8']], dtype=np.float64)
    cat4_re = np.array([z['re'] for z in hz4['zeros_level8']], dtype=np.float64)
    cat4_im = np.array([z['im'] for z in hz4['zeros_level8']], dtype=np.float64)

    out = dict(meta=dict(delta6=fp6['delta6'], AU0_6=hz6['AU0_6'],
                         delta4=fp4['delta4'], AU0_4=hz4['AU0_4']))

    # ---- 4 late targets: |H(node)| and |H'(node)| via central diff ----
    late_targets = [(1.0, 14.5), (20.125, 19.25), (41.3125, 10.0), (-6.5, 46.9375)]
    eps = 1e-6  # float64 step
    late_results = []
    for tre, tim in late_targets:
        node_re = round(tre * 16) / 16.0
        node_im = round(tim * 16) / 16.0
        u_node = np.array([node_re + 1j * node_im], dtype=np.complex128)
        u_pp = np.array([node_re + eps + 1j * node_im], dtype=np.complex128)
        u_mm = np.array([node_re - eps + 1j * node_im], dtype=np.complex128)
        u_ppi = np.array([node_re + 1j * (node_im + eps)], dtype=np.complex128)
        u_mmi = np.array([node_re + 1j * (node_im - eps)], dtype=np.complex128)
        H_node, valid_n, _ = H_at(u_node, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        H_pp, valid_pp, _ = H_at(u_pp, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        H_mm, valid_mm, _ = H_at(u_mm, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        H_ppi, valid_ppi, _ = H_at(u_ppi, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        H_mmi, valid_mmi, _ = H_at(u_mmi, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        dH_dre = (H_pp[0] - H_mm[0]) / (2 * eps)
        dH_dim = (H_ppi[0] - H_mmi[0]) / (2 * eps)
        dH_abs = math.sqrt(dH_dre ** 2 + dH_dim ** 2)
        late_results.append(dict(target=(tre, tim), node=(node_re, node_im),
                                 absH_at_node=float(H_node[0]),
                                 dH_dre=float(dH_dre),
                                 dH_dim=float(dH_dim),
                                 abs_grad=float(dH_abs),
                                 valid=bool(valid_n[0])))
        log('target=(%.4f,%.4f)  node=(%.6f,%.6f)  |H(node)|=%.3e  |H\'(node)|=%.3e'
            % (tre, tim, node_re, node_im, H_node[0], dH_abs))

    out['late_local'] = late_results

    # ---- 20 random k/16 nodes WITHOUT nearby catalog zero (d=6) ----
    rng = np.random.default_rng(20260815)
    n_random = 200  # sample more for stat reliability
    while True:
        # Sample k_re in [-50*16, 50*16), k_im in [-50*16, 50*16), so u in [-50, 50]^2
        k_re_int = rng.integers(-40 * 16, 40 * 16, n_random)
        k_im_int = rng.integers(-40 * 16, 40 * 16, n_random)
        nodes_re = k_re_int / 16.0
        nodes_im = k_im_int / 16.0
        # Exclude nodes within 0.05 of any catalog zero
        u = (nodes_re + 1j * nodes_im).astype(np.complex128)
        Habs, valid, _ = H_at(u, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
        # Check distance to nearest catalog zero
        # for each node, nearest cat dist
        # Vectorize: compute all pairwise distances?
        # Simpler: just for each node, search nearest catalog
        logH = np.log10(np.where(valid, Habs, np.nan))
        # Filter out invalid (out of domain)
        keep = np.isfinite(logH) & valid
        # distance to nearest catalog
        d_to_cat = np.array([nearest_cat_dist(r, i, cat6_re, cat6_im) for r, i in zip(nodes_re, nodes_im)])
        # Take those far from any catalog zero (>0.05) AND valid
        mask = keep & (d_to_cat > 0.05)
        if mask.sum() >= 20:
            break
    log('random k/16 nodes test d=6: n=%d valid=%d  far_from_cat=%d'
        % (n_random, keep.sum(), mask.sum()))
    logH_nodes = logH[mask]
    absH_nodes = Habs[mask]
    log('  |H| at k/16 nodes (no cat nearby): min=%.3e p10=%.3e p25=%.3e median=%.3e p75=%.3e p90=%.3e max=%.3e'
        % (absH_nodes.min(), np.percentile(absH_nodes, 10), np.percentile(absH_nodes, 25),
            np.median(absH_nodes), np.percentile(absH_nodes, 75),
            np.percentile(absH_nodes, 90), absH_nodes.max()))

    # Random non-node u points (uniform in [-30, 30]^2, exclude nodes and exclude near-cat)
    n_rand2 = 500
    re_rand = rng.uniform(-30, 30, n_rand2)
    im_rand = rng.uniform(-30, 30, n_rand2)
    # Make sure NOT near k/16 node: |u - round(16u)/16| > some threshold
    fr = 16 * re_rand - np.floor(16 * re_rand)
    fi = 16 * im_rand - np.floor(16 * im_rand)
    fr = np.minimum(fr, 1 - fr); fi = np.minimum(fi, 1 - fi)
    dist_to_node = np.sqrt((fr / 16) ** 2 + (fi / 16) ** 2)
    mask_off = dist_to_node > 0.01  # not within 1/100 of a node
    re_rand = re_rand[mask_off]; im_rand = im_rand[mask_off]
    u = (re_rand + 1j * im_rand).astype(np.complex128)
    Habs, valid, _ = H_at(u, 6, g6, h6, hz6['AU0_6'], fp6['delta6'])
    keep = valid
    absH_off = Habs[keep]
    log('random u points (off-node, no cat nearby): n=%d'
        % absH_off.size)
    log('  |H| off-node: min=%.3e p10=%.3e p25=%.3e median=%.3e p75=%.3e p90=%.3e max=%.3e'
        % (absH_off.min(), np.percentile(absH_off, 10), np.percentile(absH_off, 25),
            np.median(absH_off), np.percentile(absH_off, 75),
            np.percentile(absH_off, 90), absH_off.max()))

    # Mann-Whitney: are |H(node)| < |H(off-node)| ?
    # log-transform for normality
    logH_at_nodes = np.log10(absH_nodes)
    logH_off = np.log10(absH_off)
    log('  median log|H| at nodes = %.3f' % np.median(logH_at_nodes))
    log('  median log|H| off-node = %.3f' % np.median(logH_off))
    out['d6_node_test'] = dict(
        n_nodes=int(mask.sum()),
        n_offnode=int(absH_off.size),
        absH_nodes_stats=dict(
            min=float(absH_nodes.min()), p10=float(np.percentile(absH_nodes, 10)),
            p25=float(np.percentile(absH_nodes, 25)),
            median=float(np.median(absH_nodes)),
            p75=float(np.percentile(absH_nodes, 75)),
            p90=float(np.percentile(absH_nodes, 90)),
            max=float(absH_nodes.max())),
        absH_offnode_stats=dict(
            min=float(absH_off.min()), p10=float(np.percentile(absH_off, 10)),
            p25=float(np.percentile(absH_off, 25)),
            median=float(np.median(absH_off)),
            p75=float(np.percentile(absH_off, 75)),
            p90=float(np.percentile(absH_off, 90)),
            max=float(absH_off.max())),
        median_logH_nodes=float(np.median(logH_at_nodes)),
        median_logH_offnode=float(np.median(logH_off)),
    )

    # ---- d=4 ----
    n_random = 200
    while True:
        k_re_int = rng.integers(-25 * 16, 25 * 16, n_random)
        k_im_int = rng.integers(-25 * 16, 25 * 16, n_random)
        nodes_re = k_re_int / 16.0
        nodes_im = k_im_int / 16.0
        u = (nodes_re + 1j * nodes_im).astype(np.complex128)
        Habs, valid, _ = H_at(u, 4, g4, h4, hz4['AU0_4'], fp4['delta4'])
        d_to_cat = np.array([nearest_cat_dist(r, i, cat4_re, cat4_im) for r, i in zip(nodes_re, nodes_im)])
        mask = valid & (d_to_cat > 0.05)
        if mask.sum() >= 20:
            break
    absH_nodes = Habs[mask]
    log('random k/16 nodes test d=4: n=%d far_from_cat=%d'
        % (n_random, mask.sum()))
    log('  |H| at k/16 nodes (no cat nearby) d=4: median=%.3e p10=%.3e p90=%.3e'
        % (np.median(absH_nodes), np.percentile(absH_nodes, 10),
            np.percentile(absH_nodes, 90)))

    # off-node
    n_rand2 = 500
    re_rand = rng.uniform(-20, 20, n_rand2)
    im_rand = rng.uniform(-20, 20, n_rand2)
    fr = 16 * re_rand - np.floor(16 * re_rand)
    fi = 16 * im_rand - np.floor(16 * im_rand)
    fr = np.minimum(fr, 1 - fr); fi = np.minimum(fi, 1 - fi)
    dist_to_node = np.sqrt((fr / 16) ** 2 + (fi / 16) ** 2)
    mask_off = dist_to_node > 0.01
    re_rand = re_rand[mask_off]; im_rand = im_rand[mask_off]
    u = (re_rand + 1j * im_rand).astype(np.complex128)
    Habs, valid, _ = H_at(u, 4, g4, h4, hz4['AU0_4'], fp4['delta4'])
    absH_off = Habs[valid]
    log('  |H| off-node d=4: median=%.3e p10=%.3e p90=%.3e'
        % (np.median(absH_off), np.percentile(absH_off, 10),
            np.percentile(absH_off, 90)))
    out['d4_node_test'] = dict(
        n_nodes=int(mask.sum()),
        n_offnode=int(absH_off.size),
        absH_nodes_stats=dict(median=float(np.median(absH_nodes)),
                                p10=float(np.percentile(absH_nodes, 10)),
                                p90=float(np.percentile(absH_nodes, 90))),
        absH_offnode_stats=dict(median=float(np.median(absH_off)),
                                 p10=float(np.percentile(absH_off, 10)),
                                 p90=float(np.percentile(absH_off, 90))),
    )

    with open(os.path.join(OUT_DIR, 'local_nodes.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('wrote local_nodes.json  elapsed=%.1fs' % (time.time() - T0))


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '2')
    os.environ.setdefault('MKL_NUM_THREADS', '2')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
    main()
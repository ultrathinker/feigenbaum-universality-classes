#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F8 figures: Rayleigh statistic plot + drift-by-m plot."""
import os
import json
import time
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT_IN = r'D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип'
OUT_DIR = r'.'
T0 = time.time()
LOG_PATH = os.path.join(OUT_DIR, 'progress.log')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    with open(LOG_PATH, 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def main():
    log('--- figs: load ---')
    R = json.load(open(os.path.join(OUT_DIR, 'rayleigh_results.json'), encoding='utf-8'))
    Z = json.load(open(os.path.join(OUT_DIR, 'zeros_vs_m.json'), encoding='utf-8'))

    # Figure 1: Rayleigh coarse + fine
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    ax = axes[0, 0]
    qc = np.array(R['coarse']['q'])
    ax.plot(qc, R['coarse']['R6_re'], label='d=6 Re', lw=1)
    ax.plot(qc, R['coarse']['R6_im'], label='d=6 Im', lw=1)
    ax.plot(qc, R['coarse']['R4_re'], label='d=4 Re', lw=1)
    ax.plot(qc, R['coarse']['R4_im'], label='d=4 Im', lw=1)
    ax.axvline(16.0, color='r', ls='--', lw=0.7, label='q=16')
    ax.set_xlabel('q')
    ax.set_ylabel('R(q)')
    ax.set_title('Rayleigh R(q) — coarse q in [0.5, 40] step 0.05')
    ax.set_xlim(0.5, 40)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[0, 1]
    qf = np.array(R['fine']['q'])
    ax.plot(qf, R['fine']['R6_re'], label='d=6 Re', lw=1)
    ax.plot(qf, R['fine']['R6_im'], label='d=6 Im', lw=1)
    ax.plot(qf, R['fine']['R4_re'], label='d=4 Re', lw=1)
    ax.plot(qf, R['fine']['R4_im'], label='d=4 Im', lw=1)
    ax.axvline(16.0, color='r', ls='--', lw=0.7, label='q=16')
    ax.set_xlabel('q')
    ax.set_ylabel('R(q)')
    ax.set_title('Fine q in [15.5, 16.5] step 0.001 — peak exactly at 16.000')
    ax.set_xlim(15.5, 16.5)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # Figure 2: drift-vs-m
    ax = axes[1, 0]
    for tr in Z['targets']:
        bms = tr['by_m']
        ms = [b['m'] for b in bms]
        ds = [b['dist_to_node'] for b in bms]
        ax.semilogy(ms, ds, '-o', label=f'target ({tr["target"][0]:.2f},{tr["target"][1]:.2f})')
    ax.set_xlabel('m')
    ax.set_ylabel('distance to nearest k/16 node')
    ax.set_title('Distance of refined zero from lattice node vs m (d=6)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which='both')

    # Figure 3: |H| at lattice nodes vs off-node
    L = json.load(open(os.path.join(OUT_DIR, 'local_nodes.json'), encoding='utf-8'))
    ax = axes[1, 1]
    nodes = L['d6_node_test']['absH_nodes_stats']
    off = L['d6_node_test']['absH_offnode_stats']
    labels = ['p10', 'p25', 'median', 'p75', 'p90']
    x = np.arange(len(labels))
    nodes_vals = [nodes['p10'], nodes['p25'], nodes['median'], nodes['p75'], nodes['p90']]
    off_vals = [off['p10'], off['p25'], off['median'], off['p75'], off['p90']]
    w = 0.35
    ax.bar(x - w/2, nodes_vals, w, label='at k/16 nodes (no cat zero nearby)')
    ax.bar(x + w/2, off_vals, w, label='off-node random points')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('|H|')
    ax.set_title('|H| distribution at k/16 nodes vs off-node (d=6, m=8)')
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis='y')

    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, 'fig_f8_summary.png')
    plt.savefig(out_png, dpi=120)
    log('wrote %s' % out_png)

    # Additional: Rayleigh fine (zoomed) standalone
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(qf, R['fine']['R6_re'], label='d=6 Re', lw=1.2)
    ax.plot(qf, R['fine']['R6_im'], label='d=6 Im', lw=1.2)
    ax.plot(qf, R['fine']['R4_re'], label='d=4 Re', lw=1.2)
    ax.plot(qf, R['fine']['R4_im'], label='d=4 Im', lw=1.2)
    ax.axvline(16.0, color='r', ls='--', lw=0.8, label='q=16.000')
    ax.set_xlabel('q')
    ax.set_ylabel('R(q) = |⟨exp(2πi q x)⟩|')
    ax.set_title('Rayleigh statistic — fine scan q∈[15.5,16.5] step 0.001')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, 'fig_f8_rayleigh_fine.png')
    plt.savefig(out_png, dpi=120)
    log('wrote %s' % out_png)


if __name__ == '__main__':
    main()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build patch figure: scatter of zeros found by scans A,B,C,D + lattice 1/16 +
S104 marker.  Two side-by-side panels: (1) raw scan-unique zeros, (2) only those
passing genuine-zero validation (drift 7->8->9)."""
import json
import os
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def fig_patch(rs):
    patch = rs['meta']['patch']
    re_lo, re_hi = patch['re_lo'], patch['re_hi']
    im_lo, im_hi = patch['im_lo'], patch['im_hi']

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5))

    # panel 1: all unique zeros, color by which scans found them
    ax = axes[0]
    ax.set_xlim(re_lo, re_hi)
    ax.set_ylim(im_lo, im_hi)
    ax.set_aspect('equal')
    ax.set_xlabel('Re u')
    ax.set_ylabel('Im u')
    ax.set_title('All unique zeros (clustered 1e-8)')

    # background: lattice 1/16 points in patch
    L = 16
    rr = []
    for k in range(int(re_lo * L), int(re_hi * L) + 1):
        if re_lo - 1e-9 <= k / L <= re_hi + 1e-9:
            rr.append(k / L)
    for k in range(int(im_lo * L), int(im_hi * L) + 1):
        pass  # we'll just plot on axes below
    # draw lattice 1/16 nodes (full AND half-integer intermediate)
    for Lq in (16, 8):
        for re in [k / Lq for k in range(int(re_lo * Lq), int(re_hi * Lq) + 1)
                    if re_lo - 1e-9 <= k / Lq <= re_hi + 1e-9]:
            ax.axvline(re, color='gray', linewidth=0.5, alpha=0.4, zorder=1)
        for im in [k / Lq for k in range(int(im_lo * Lq), int(im_hi * Lq) + 1)
                    if im_lo - 1e-9 <= k / Lq <= im_hi + 1e-9]:
            ax.axhline(im, color='gray', linewidth=0.5, alpha=0.4, zorder=1)

    # catalog zeros (purple stars)
    for cz in rs['catalog_patch']:
        ax.plot(cz['re'], cz['im'], '*', color='purple', markersize=14,
                markeredgecolor='black', zorder=5, label='catalog 3')

    # scan A: thick black circle if also in catalog, else thin
    for z in rs['scans']['A']['clustered']:
        ax.plot(z['re'], z['im'], 'o', color='black', markersize=9, zorder=4)
        ax.text(z['re'] + 0.05, z['im'] + 0.05, 'A', fontsize=8,
                color='black', zorder=4)
    for z in rs['scans']['B']['clustered']:
        ax.plot(z['re'], z['im'], 's', color='blue', markersize=8, zorder=4)
        ax.text(z['re'] + 0.05, z['im'] + 0.05, 'B', fontsize=8,
                color='blue', zorder=4)
    for z in rs['scans']['C']['clustered']:
        ax.plot(z['re'], z['im'], 'D', color='green', markersize=7, zorder=4)
    # D lots of points - small grey x
    dx = [z['re'] for z in rs['scans']['D']['clustered']]
    dy = [z['im'] for z in rs['scans']['D']['clustered']]
    ax.scatter(dx, dy, marker='x', color='red', s=15, alpha=0.55, zorder=3,
               label='D (random)')

    # S104 reference
    ax.plot(20.125013, 19.250015, 'P', color='orange', markersize=15,
            markeredgecolor='black', zorder=6, label='S104 (20.125, 19.250)')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.85)

    # panel 2: only validated genuine zeros (subset of D, etc.)
    ax = axes[1]
    ax.set_xlim(re_lo, re_hi)
    ax.set_ylim(im_lo, im_hi)
    ax.set_aspect('equal')
    ax.set_xlabel('Re u')
    ax.set_ylabel('Im u')
    ax.set_title('Genuine (drift 7->8->9 + domain OK)')
    for Lq in (16, 8):
        for re in [k / Lq for k in range(int(re_lo * Lq), int(re_hi * Lq) + 1)
                    if re_lo - 1e-9 <= k / Lq <= re_hi + 1e-9]:
            ax.axvline(re, color='gray', linewidth=0.5, alpha=0.4, zorder=1)
        for im in [k / Lq for k in range(int(im_lo * Lq), int(im_hi * Lq) + 1)
                    if im_lo - 1e-9 <= k / Lq <= im_hi + 1e-9]:
            ax.axhline(im, color='gray', linewidth=0.5, alpha=0.4, zorder=1)
    for cz in rs['catalog_patch']:
        ax.plot(cz['re'], cz['im'], '*', color='purple', markersize=14,
                markeredgecolor='black', zorder=5)

    # if validation done, show genuine/ambiguous/fake differently
    if 'validation' in rs:
        for group, color, marker in [('genuine', 'green', 'o'),
                                     ('ambiguous', 'orange', 's'),
                                     ('fake', 'red', 'x')]:
            xs = [g['re'] for g in rs['validation'][group]]
            ys = [g['im'] for g in rs['validation'][group]]
            ax.scatter(xs, ys, marker=marker, color=color, s=45, alpha=0.85,
                       zorder=4, label=f'{group} ({len(xs)})')
    else:
        for z in rs['scans']['A']['clustered']:
            ax.plot(z['re'], z['im'], 'o', color='black', markersize=9, zorder=4)
        for z in rs['scans']['B']['clustered']:
            ax.plot(z['re'], z['im'], 's', color='blue', markersize=8, zorder=4)

    ax.plot(20.125013, 19.250015, 'P', color='orange', markersize=15,
            markeredgecolor='black', zorder=6)
    ax.legend(loc='upper left', fontsize=8, framealpha=0.85)

    fig.suptitle('Rescan of patch [18,22]x[17,21] (lattice 1/16 lines grey): brief S')
    fig.tight_layout()
    out = os.path.join(HERE, 'patch_rescan.png')
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print('wrote', out)


if __name__ == '__main__':
    rs = json.load(open(os.path.join(HERE, 'rescan_results.json'), encoding='utf-8'))
    fig_patch(rs)

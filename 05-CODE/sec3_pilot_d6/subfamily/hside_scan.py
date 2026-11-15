#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hside_scan.py -- BRIEF-G, cheap H-side leg: full-resolution 2D grid of
|H_m(u)| (numpy complex128, no Newton) over the window u_cat +/- 1 local
spacing, step = spacing/64, for ALL 9 anomalous-subfamily zeros at levels
m=9 and m=10.

Purpose: test hypothesis (b) [search-grid artifact -- does a cleaner/deeper
H-zero sit near the catalogued position that the coarse 1/16 pilot scan
could have mislocated?] and give a first, cheap look at hypothesis (a)
[interference / near-double zero -- would show as a second comparably-deep
local minimum inside the window].  Single-threaded, fast (calibrated
~3-5s per zero per level), NOT backgrounded.
[oid:deg6G] NOT RIGOROUS -- float64, no interval arithmetic.
"""
import json
import os
import time

import numpy as np

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)


A_np = np.array([float(x) for x in C.A_MP])
Cc_np = np.array([float(x) for x in C.C_MP])
AU0f = float(C.AU0_6)
DELf = float(C.DELTA6)


def poly6(z, co):
    t = z * z
    t = t * t * t
    r = np.full_like(z, co[-1])
    for c in reversed(co[:-1]):
        r = r * t + c
    return r


def scan_grid(u0, spacing, ppa, m):
    n = ppa  # +/- 1 spacing at step spacing/ppa -> n steps each side
    offs = np.linspace(-1.0, 1.0, 2 * n + 1) * spacing
    RE, IM = np.meshgrid(offs, offs)
    grid = (u0.real + RE) + 1j * (u0.imag + IM)
    s = complex(AU0f / DELf ** m) * grid
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    x2m = None
    with np.errstate(over='ignore', invalid='ignore'):
        for k in range(n2):
            x = poly6(x, A_np) + s * poly6(x, Cc_np)
            if k + 1 == n1:
                x2m = x.copy()
        Habs = np.abs(x) / np.abs(x2m)
    valid = np.isfinite(Habs)
    return grid, Habs, valid, offs


def local_minima(Habs, valid, grid, offs):
    """Crude local-minimum finder on the grid (4-neighbour comparison)."""
    H = np.where(valid, Habs, np.inf)
    mins = []
    ny, nx = H.shape
    for iy in range(1, ny - 1):
        for ix in range(1, nx - 1):
            v = H[iy, ix]
            if not np.isfinite(v):
                continue
            nb = H[iy - 1:iy + 2, ix - 1:ix + 2]
            if v == nb.min():
                mins.append((float(v), complex(grid[iy, ix]), float(offs[ix]), float(offs[iy])))
    mins.sort(key=lambda r: r[0])
    return mins


targets = json.load(open(os.path.join(HERE, 'targets.json'), encoding='utf-8'))['targets']
PPA = 64  # grid points per spacing unit (step = spacing/64), window = +/-1 spacing

out = {}
log('BRIEF-G hside_scan: PPA=%d window=+/-1 spacing, levels m=9,10, targets=%d'
    % (PPA, len(targets)))
for i, (tid, t) in enumerate(targets.items()):
    u0 = complex(t['u_re'], t['u_im'])
    spacing = t['nearest_neighbour_spacing']
    entry = {'id': tid, 'u_cat': [u0.real, u0.imag], 'spacing': spacing, 'levels': {}}
    for m in (9, 10):
        t0 = time.time()
        grid, Habs, valid, offs = scan_grid(u0, spacing, PPA, m)
        mins = local_minima(Habs, valid, grid, offs)
        dt = time.time() - t0
        # value & offset of |H| exactly at the catalog point (nearest grid node)
        cidx = PPA
        Hcat = float(Habs[cidx, cidx])
        top = mins[:5]
        entry['levels'][str(m)] = {
            'time_sec': dt,
            'H_at_catalog_point': Hcat,
            'n_local_minima': len(mins),
            'top5_minima': [
                {'absH': v, 'u_re': u.real, 'u_im': u.imag,
                 'offset_re_in_spacing': dre / spacing, 'offset_im_in_spacing': dim / spacing,
                 'dist_from_catalog_in_spacing': abs(u - u0) / spacing}
                for (v, u, dre, dim) in top
            ],
        }
        log('%s m=%d: H@catalog=%.3e, deepest min |H|=%.3e at dist=%.4f*spacing, '
            '#minima=%d (%.2fs)'
            % (tid, m, Hcat, top[0][0] if top else float('nan'),
               (abs(top[0][1] - u0) / spacing) if top else float('nan'),
               len(mins), dt))
    out[tid] = entry
    with open(os.path.join(HERE, 'hside_scan.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('progress: %d/%d zeros done, wrote hside_scan.json' % (i + 1, len(targets)))

log('hside_scan DONE, total %.1fs' % (time.time() - T0))

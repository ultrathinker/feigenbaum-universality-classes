#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calibrate.py -- BRIEF-G step 2: time the cost primitives before committing to
a grid size.  Single-threaded (budget: <=2 CPU threads total across the whole
investigation).
[oid:deg6G] NOT RIGOROUS.
"""
import json
import os
import time

from mpmath import mp, mpf, mpc, fabs

mp.dps = 80

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
targets = json.load(open(os.path.join(HERE, 'targets.json'), encoding='utf-8'))['targets']

results = {}

# --- 1. raw orbit_x_dx cost (single call, no Newton) at n=1024/2048 --------
for n, tag in ((1024, 'n1024_m9'), (2048, 'n2048_m10')):
    a = C.A_F  # any finite a works for timing
    t0 = time.time()
    C.orbit_x_dx(a, n)
    dt = time.time() - t0
    results['orbit_x_dx_' + tag] = dt
    print('orbit_x_dx %s: %.4fs' % (tag, dt), flush=True)

# --- 2. best_du (full Newton, direct + cascade fallback) at m=9 and m=10 ---
sample_ids = ['S036', 'S104']  # 2 of the priority quartet
for m in (9, 10):
    for tid in sample_ids:
        t = targets[tid]
        u = mpc(mpf(t['u_re']), mpf(t['u_im']))
        t0 = time.time()
        d, a_ss, tag = C.best_du(u, m)
        dt = time.time() - t0
        key = '%s_m%d' % (tid, m)
        results['best_du_' + key] = {'sec': dt, 'd_u': d, 'tag': tag}
        print('best_du %s: %.3fs  d_u=%s  tag=%s' % (key, dt, d, tag), flush=True)

# --- 3. H-side orbit (cheap, mpmath, m=9) single evaluation ----------------
mp.dps = 50
for m in (9, 10):
    u = mpc(mpf(targets['S036']['u_re']), mpf(targets['S036']['u_im']))
    t0 = time.time()
    C.H_abs(u, m)
    dt = time.time() - t0
    results['H_abs_mpmath_m%d' % m] = dt
    print('H_abs mpmath(dps50) m=%d: %.4fs' % (m, dt), flush=True)
mp.dps = 80

# --- 4. numpy vectorized H_m grid scan cost (float64, cheap end) -----------
import numpy as np

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


def scan_grid_numpy(u0, spacing, halfwidth_in_spacing, ppa, m):
    """halfwidth_in_spacing * spacing on each side, ppa points per spacing
    unit (so step = spacing/ppa); returns |H_m| array."""
    n = int(round(halfwidth_in_spacing * ppa))
    offs = np.linspace(-halfwidth_in_spacing, halfwidth_in_spacing, 2 * n + 1) * spacing
    RE, IM = np.meshgrid(offs, offs)
    grid = (u0.real + RE) + 1j * (u0.imag + IM)
    s = complex(AU0f / DELf ** m) * grid
    n1, n2 = 2 ** m, 2 ** (m + 1)
    x = np.zeros(grid.shape, dtype=np.complex128)
    x2m = None
    for k in range(n2):
        x = poly6(x, A_np) + s * poly6(x, Cc_np)
        if k + 1 == n1:
            x2m = x.copy()
    with np.errstate(divide='ignore', invalid='ignore'):
        Habs = np.abs(x) / np.abs(x2m)
    return Habs, grid.shape


for m in (9, 10):
    t0 = time.time()
    Habs, shp = scan_grid_numpy(complex(3.62497, 16.12503), 2.885, 1.0, 64, m)
    dt = time.time() - t0
    results['numpy_grid_m%d' % m] = {'sec': dt, 'shape': list(shp), 'npts': int(shp[0] * shp[1])}
    print('numpy grid scan m=%d: shape=%s pts=%d time=%.3fs (%.2e s/pt)'
          % (m, shp, shp[0] * shp[1], dt, dt / (shp[0] * shp[1])), flush=True)

with open(os.path.join(HERE, 'calibration.json'), 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2, ensure_ascii=False, default=str)
print('wrote calibration.json')

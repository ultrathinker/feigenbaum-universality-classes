#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1b_common.py -- shared helpers for the R1B d=6 forensics package (analogue
of r1a_common.py for d=4).  Imports p2_H6_zeros.py functions AS-IS (H,
newton_refine, domain_check, scan_level, orbit_H_mpmath); copies the small,
already-reviewed statistics helpers (frac/near_stats/rayleigh) verbatim.

d=6-specific differences from d=4 handled here:
  - p2_H6_zeros.newton_refine's OWN default is rmax=70.0 (NOT 22.0 like d=4).
    r1_d6_rescan.py calls it with rmax=R_max+5=53.0 for STEP24 -- this is
    SMALLER than p2's own default (opposite-direction deviation from the d=4
    case, where r1_d4_rescan.py's rmax=25 was LARGER than p2's 22).
  - domain_ok criterion for d=6 (both in p2_H6_zeros.main AND r1_d6_rescan.py)
    is `maxt6 < TMAX6(=1.15) and den > 1e-12` -- NO truncation-error (merr)
    threshold, unlike d=4's `merr < 1e-11` extra condition.
  - mp_polish uses iters=16 (not 14) per p2_H6_zeros' own mp_polish_job calls.
  - catalog 'str' field format is 're_str' + '+' + 'im_str' + 'i' (NO spaces,
    NO parens) -- always exactly one literal '+' injected as separator
    regardless of im's own sign (so splitting on the FIRST '+' is safe and
    unambiguous for this dataset's magnitude range, no exponential notation).
"""
import json
import math
import os
import time

import numpy as np
from mpmath import mp, mpf, mpc, nstr, re as mpre, im as mpim
from scipy import stats as sstats

mp.dps = 50

import p2_H6_zeros as p2  # noqa: E402

CWD = os.path.dirname(os.path.abspath(__file__))
R_MAX = 48.0
STEP16 = 0.0625
TMAX6 = p2.TMAX6  # 1.15
ORIG_PATH = os.path.join(CWD, 'H6_zeros_ORIG.json')

P2_DEFAULT_RMAX = 70.0  # p2_H6_zeros.newton_refine's own default (main() calls
# it WITHOUT rmax).  r1_d6_rescan.py instead passes rmax=R_max+5=53.0 for
# STEP24 -- SMALLER than p2's own reach (opposite direction from the d=4 case).


def mk_logger(path):
    t0 = time.time()
    f = open(path, 'a', encoding='utf-8')

    def log(msg):
        line = "[%7.1fs] %s" % (time.time() - t0, msg)
        try:
            print(line, flush=True)
        except UnicodeEncodeError:
            print(line.encode('ascii', 'replace').decode('ascii'), flush=True)
        f.write(line + "\n")
        f.flush()
    return log


def load_au0(cache_path=os.path.join(CWD, 'H6_AU0_cache.json')):
    cache = json.load(open(cache_path, encoding='utf-8'))
    AU0 = mpf(cache['AU0_6'])
    Ls = mpf(cache['Ls'])
    return AU0, Ls


def mp_polish(u_start, m, AU0, iters=16, tol_step=mpf('1e-24')):
    u_mp = mpc(mpf(u_start.real), mpf(u_start.imag))
    for _ in range(iters):
        x2m, x2m1, dx2m, dx2m1 = p2.orbit_H_mpmath(u_mp, m, AU0, want_grad=True)
        H = x2m1 / x2m
        dH = (dx2m1 * x2m - x2m1 * dx2m) / x2m ** 2
        step = H / dH
        u_mp = u_mp - step
        if abs(step) < tol_step:
            break
    x2m, x2m1 = p2.orbit_H_mpmath(u_mp, m, AU0)
    return u_mp, abs(x2m1 / x2m)


def H8_mp_at(u_mp, AU0, m=8):
    x2m, x2m1 = p2.orbit_H_mpmath(u_mp, m, AU0)
    if abs(x2m) == 0:
        return None
    return abs(x2m1 / x2m)


def grid_orig(R=R_MAX, step=STEP16, shift_re=0.0, shift_im=0.0):
    xs = np.arange(-R + shift_re, R + step / 2, step)
    ys = np.arange(0.0 + shift_im, R + step / 2, step)
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    return UU + 1j * VV


def local_min_patch(u0, AU0, step, shift_re, shift_im, m=8, half=2):
    i0 = round((u0.real - shift_re) / step)
    j0 = round((u0.imag - shift_im) / step)
    idx = np.arange(-half, half + 1)
    xs = shift_re + (i0 + idx) * step
    ys = shift_im + (j0 + idx) * step
    UU, VV = np.meshgrid(xs, ys, indexing='ij')
    grid = UU + 1j * VV
    Habs, valid, maxt, errmax = p2.scan_level(m, AU0, grid)
    Hm = np.where(valid, Habs, np.inf)
    c = half
    core = Hm[c, c]
    neigh = [Hm[c - 1, c - 1], Hm[c - 1, c], Hm[c - 1, c + 1],
             Hm[c, c - 1], Hm[c, c + 1],
             Hm[c + 1, c - 1], Hm[c + 1, c], Hm[c + 1, c + 1]]
    is_local_min = bool(core <= min(neigh) and core < 0.2 and np.isfinite(core))
    return dict(node_re=float(grid[c, c].real), node_im=float(grid[c, c].imag),
                center_Habs=float(core) if np.isfinite(core) else None,
                center_valid=bool(valid[c, c]),
                is_local_min=is_local_min,
                neighbour_Habs=[float(x) if np.isfinite(x) else None for x in neigh])


def newton_from(u0, m, AU0, tol=1e-9, rmax=P2_DEFAULT_RMAX):
    try:
        u, res = p2.newton_refine(u0, m, AU0, tol=tol, rmax=rmax)
    except (OverflowError, ZeroDivisionError, ValueError):
        return None, float('inf')
    return u, res


def domain_ok(u, m, AU0):
    mt6, merr, den = p2.domain_check(u, m, AU0)
    return bool(mt6 < TMAX6 and den > 1e-12), dict(maxt6=mt6, maxerr=merr, den=den)


def full_candidate_pipeline(u0, AU0, m=8):
    out = {'u0': complex(u0)}
    u, res = newton_from(u0, m, AU0)
    out['newton_u'] = complex(u) if u is not None else None
    out['newton_absH'] = res
    if u is None or res > 1e-5:
        out['stage_failed'] = 'newton(|H|<=1e-5)'
        return out
    ok, dom = domain_ok(u, m, AU0)
    out['domain'] = dict(**dom, ok=ok)
    if not ok:
        out['stage_failed'] = 'domain_check'
        return out
    u_mp, r = mp_polish(u, m, AU0)
    out['polish_absH'] = float(r)
    out['polish_u'] = complex(mpre(u_mp), mpim(u_mp))
    out['polish_u_mp'] = nstr(u_mp, 25)
    if float(r) >= 1e-20:
        out['stage_failed'] = 'mp_polish(|H|<1e-20)'
        return out
    out['stage_failed'] = None
    out['accepted'] = True
    return out


def parse_u_mp(zrec):
    """Parse d=6 catalog 'str' field: 're_str+im_stri' (single literal '+'
    separator injected regardless of sign; split on the FIRST '+')."""
    s = zrec.get('str')
    if s:
        body = s.strip()
        assert body.endswith('i'), body
        body = body[:-1]
        if '+' not in body:
            raise ValueError('unrecognized d6 mpc string (no separator): %r' % s)
        re_s, im_s = body.split('+', 1)
        return mpc(mpf(re_s), mpf(im_s if im_s else '0'))
    return mpc(mpf(zrec['re']), mpf(zrec['im']))


# --------------------------------------------------------------------------- #
#  lattice statistics -- verbatim copy (frac/near_stats/rayleigh)
# --------------------------------------------------------------------------- #
def frac(x, q):
    f = (q * np.asarray(x, float)) % 1.0
    return np.minimum(f, 1.0 - f)


def near_stats(re, im, q=16, tol=1.0 / 16):
    re = np.asarray(re, float)
    im = np.asarray(im, float)
    n = len(re)
    if n == 0:
        return dict(n=0, re_near=None, im_near=None, joint=None, p_re=None, p_im=None, p_joint=None)
    fr, fi = frac(re, q), frac(im, q)
    nr, ni = fr < tol, fi < tol
    p0 = 2 * tol
    return dict(n=n, re_near=float(nr.mean()), im_near=float(ni.mean()),
                joint=float((nr & ni).mean()),
                p_re=float(sstats.binomtest(int(nr.sum()), n, p0, alternative='greater').pvalue),
                p_im=float(sstats.binomtest(int(ni.sum()), n, p0, alternative='greater').pvalue),
                p_joint=float(sstats.binomtest(int((nr & ni).sum()), n, p0 * p0, alternative='greater').pvalue))


def near_stats_shifted(re, im, q, shift_re, shift_im, tol=1.0 / 16):
    re = np.asarray(re, float) - shift_re
    im = np.asarray(im, float) - shift_im
    return near_stats(re, im, q=q, tol=tol)


def rayleigh(x, q):
    x = np.asarray(x, float)
    z = np.exp(2j * np.pi * q * x)
    R = abs(z.sum()) / len(x)
    Z = len(x) * R * R
    return float(R), float(math.exp(-Z))


def rayleigh_scan(re, im, q_range=range(4, 81)):
    out = []
    for q in q_range:
        Rr, pr = rayleigh(re, q)
        Ri, pi = rayleigh(im, q)
        out.append((q, Rr, pr, Ri, pi))
    return out


def rayleigh_peaks(scan, p_thresh=1e-4):
    return [(q, Rr, pr, Ri, pi) for (q, Rr, pr, Ri, pi) in scan if pr < p_thresh or pi < p_thresh]


def stats_block(re, im, q_list=(16, 24, 32, 48)):
    rec = {'n': len(re)}
    rec['near'] = {str(q): near_stats(re, im, q=q) for q in q_list}
    sc = rayleigh_scan(re, im)
    rec['rayleigh_peaks_p1e-4'] = [
        dict(q=q, R_re=Rr, p_re=pr, R_im=Ri, p_im=pi)
        for (q, Rr, pr, Ri, pi) in rayleigh_peaks(sc)]
    return rec


# --------------------------------------------------------------------------- #
#  catalog IO
# --------------------------------------------------------------------------- #
def load_catalog(path):
    j = json.load(open(path, encoding='utf-8'))
    zl = j['zeros_level8']
    zs = [(float(z['re']), float(z['im'])) for z in zl]
    return zs, zl, j


def match_indices(orig_zs, new_zs, tol=1e-6):
    orig_arr = np.array(orig_zs) if orig_zs else np.empty((0, 2))
    both_o, both_n, only_new = [], [], []
    for i, n in enumerate(new_zs):
        if len(orig_arr) == 0:
            only_new.append(i)
            continue
        d = np.abs(orig_arr - np.array(n)).max(axis=1)
        j = int(d.argmin())
        if d[j] < tol:
            if j not in both_o:
                both_o.append(j)
                both_n.append(i)
            else:
                only_new.append(i)
        else:
            only_new.append(i)
    only_orig = [j for j in range(len(orig_zs)) if j not in both_o]
    return both_o, both_n, only_new, only_orig


def to_json(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, mpf):
        return float(o)
    if isinstance(o, mpc):
        return {'re': float(mpre(o)), 'im': float(mpim(o))}
    if isinstance(o, complex):
        return {'re': o.real, 'im': o.imag}
    if isinstance(o, dict):
        return {k: to_json(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_json(v) for v in o]
    return o


def dump_json(obj, path):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(to_json(obj), f, indent=2, ensure_ascii=False)

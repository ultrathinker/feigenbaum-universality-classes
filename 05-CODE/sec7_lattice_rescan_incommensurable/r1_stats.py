#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1_stats.py -- lattice-rescan statistics for d=4 (62+ catalog) and d=6 (342+ catalog).

Inputs (all in this folder):
  H4_zeros.json          (orig catalog, 62)  -- NOT created here; we read the
                       canonical from quartic_pilot/H4_zeros.json
  H4_sanity.json         (rescan with orig 1/16 grid; control)
  H4_rescan_step24.json  (rescan with 1/24 + shift)
  H4_rescan_random.json  (rescan with random starts)
  H6_zeros_ORIG.json     (orig catalog, 342; copy of canonical)
  H6_rescan_random.json
  H6_rescan_step24.json

Outputs:
  r1_stats.json          (all numbers)
  pretty-printed table slice for the report

Metrics (frac<n> from reviewer_R5_lattice_stats.py):
  near(re, q, tol=1/16)  -> fraction of re with frac(q*re) < 1/16
  near stats for q in (16, 24, 32, 48, 64)
  Rayleigh test for q in [4..80]
  Catalog overlap (|deltau|<1e-6) -- 3 subsets: a/b/c

Late-9 check (d=6 only): list pre-defined late zeros and check if the new scan
recovers them within 1e-6.
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy import stats as sstats

# --------------------------------------------------------------------------- #
#  metrics (copy of reviewer_R5_lattice_stats.py -- NOT modified)
# --------------------------------------------------------------------------- #
def frac(x, q):
    f = (q * np.asarray(x, float)) % 1.0
    return np.minimum(f, 1.0 - f)


def near_stats(re, im, q=16, tol=1.0 / 16):
    fr, fi = frac(re, q), frac(im, q)
    nr, ni = fr < tol, fi < tol
    n = len(re)
    p0 = 2 * tol
    out = dict(n=n, re_near=float(nr.mean()), im_near=float(ni.mean()),
               joint=float((nr & ni).mean()),
               p_re=float(sstats.binomtest(int(nr.sum()), n, p0, alternative='greater').pvalue),
               p_im=float(sstats.binomtest(int(ni.sum()), n, p0, alternative='greater').pvalue),
               p_joint=float(sstats.binomtest(int((nr & ni).sum()), n, p0 * p0, alternative='greater').pvalue))
    return out


def rayleigh(x, q):
    x = np.asarray(x, float)
    z = np.exp(2j * np.pi * q * x)
    R = abs(z.sum()) / len(x)
    Z = len(x) * R * R
    return float(R), float(math.exp(-Z))


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #
def load_zeros(path):
    """Return list of (re, im) from a catalog JSON."""
    j = json.load(open(path, encoding='utf-8'))
    if 'zeros_level8' in j:
        zl = j['zeros_level8']
    elif 'zeros_limit_extrapolated' in j:
        zl = j['zeros_limit_extrapolated']
    else:
        raise ValueError('no zeros list in ' + path)
    return [(float(z['re']), float(z['im'])) for z in zl], j


def match_subsets(orig, new, tol=1e-6):
    """Return (both, only_new, only_orig) zero lists, |deltau|<tol."""
    orig_arr = np.array(orig)
    new_arr = np.array(new)
    if len(orig_arr) == 0 or len(new_arr) == 0:
        return [], list(new), list(orig)
    both_idx_o, both_idx_n = [], []
    only_new = []
    for i, n in enumerate(new):
        d = np.abs(orig_arr - n).max(axis=1)
        m = d.min()
        if m < tol:
            j = int(d.argmin())
            if j not in both_idx_o:
                both_idx_o.append(j)
                both_idx_n.append(i)
        else:
            only_new.append(n)
    both = [orig[j] for j in both_idx_o]
    only_orig = [o for i, o in enumerate(orig) if i not in both_idx_o]
    return both, only_new, only_orig


def rayleigh_scan(re, im, q_range=range(4, 81)):
    """Return list of (q, R_re, p_re, R_im, p_im) for each q in q_range."""
    out = []
    for q in q_range:
        Rr, pr = rayleigh(re, q)
        Ri, pi = rayleigh(im, q)
        out.append((q, Rr, pr, Ri, pi))
    return out


def rayleigh_peaks(scan, p_thresh=1e-4):
    """Return list of (q, R_re, p_re, R_im, p_im) where p<p_thresh."""
    return [(q, Rr, pr, Ri, pi) for (q, Rr, pr, Ri, pi) in scan if pr < p_thresh or pi < p_thresh]


# --------------------------------------------------------------------------- #
#  per-catalog stats
# --------------------------------------------------------------------------- #
def stats_for_label(name, re, im, q_list=(16, 24, 32, 48, 64)):
    """Return dict with near_stats for each q, and Rayleigh peaks."""
    rec = {'n': len(re)}
    rec['near'] = {str(q): near_stats(re, im, q=q) for q in q_list}
    rec['rayleigh_full'] = rayleigh_scan(re, im)
    rec['rayleigh_peaks_p1e-4'] = rayleigh_peaks(rec['rayleigh_full'])
    return rec


def meta_for(catalog_name, deg, j, n_list):
    out = {'catalog': catalog_name, 'deg': deg, 'source_n': len(j.get('zeros_level8', []))}
    return out


# --------------------------------------------------------------------------- #
#  main
# --------------------------------------------------------------------------- #
def main():
    cwd = os.path.dirname(os.path.abspath(__file__))
    out = {'date': time.strftime('%Y-%m-%d %H:%M:%S'),
           'meta': {'tolerance_matches': 1e-6, 'tolerance_dedup': 1e-8,
                    'seed': 20260816, 'q_list': [16, 24, 32, 48, 64]}}

    # ---- d=4 catalogs ----
    h4_orig, j_h4_orig = load_zeros('../sec3_pilot_d4/H4_zeros.json')
    cats_d4 = {
        'H4_orig62': h4_orig,
    }
    for nm in ('H4_sanity', 'H4_rescan_step24', 'H4_rescan_random'):
        path = os.path.join(cwd, nm + '.json')
        if os.path.exists(path):
            zs, _ = load_zeros(path)
            cats_d4[nm] = zs

    out['d4'] = {}
    for nm, zs in cats_d4.items():
        re = [z[0] for z in zs]
        im = [z[1] for z in zs]
        out['d4'][nm] = stats_for_label(nm, re, im)
        out['d4'][nm]['source_n'] = len(zs)

    # match ORIG vs each new catalog
    for nm in cats_d4:
        if nm == 'H4_orig62':
            continue
        new = cats_d4[nm]
        both, only_new, only_orig = match_subsets(h4_orig, new, tol=1e-6)
        out['d4'][nm + '_both'] = stats_for_label(nm + '_both', [z[0] for z in both], [z[1] for z in both])
        out['d4'][nm + '_only_new'] = stats_for_label(nm + '_only_new', [z[0] for z in only_new], [z[1] for z in only_new])
        out['d4'][nm + '_only_orig'] = stats_for_label(nm + '_only_orig', [z[0] for z in only_orig], [z[1] for z in only_orig])
        out['d4'][nm + '_overlap'] = {'n_both': len(both), 'n_only_new': len(only_new),
                                       'n_only_orig': len(only_orig),
                                       'frac_both': len(both) / len(h4_orig) if h4_orig else 0,
                                       'frac_only_new': len(only_new) / len(new) if new else 0}

    # ---- d=6 catalogs ----
    h6_orig, j_h6_orig = load_zeros(os.path.join(cwd, 'H6_zeros_ORIG.json'))
    cats_d6 = {'H6_orig342': h6_orig}
    for nm in ('H6_rescan_random', 'H6_rescan_step24'):
        path = os.path.join(cwd, nm + '.json')
        if os.path.exists(path):
            zs, _ = load_zeros(path)
            cats_d6[nm] = zs

    out['d6'] = {}
    for nm, zs in cats_d6.items():
        re = [z[0] for z in zs]
        im = [z[1] for z in zs]
        out['d6'][nm] = stats_for_label(nm, re, im)
        out['d6'][nm]['source_n'] = len(zs)

    for nm in cats_d6:
        if nm == 'H6_orig342':
            continue
        new = cats_d6[nm]
        both, only_new, only_orig = match_subsets(h6_orig, new, tol=1e-6)
        out['d6'][nm + '_both'] = stats_for_label(nm + '_both', [z[0] for z in both], [z[1] for z in both])
        out['d6'][nm + '_only_new'] = stats_for_label(nm + '_only_new', [z[0] for z in only_new], [z[1] for z in only_new])
        out['d6'][nm + '_only_orig'] = stats_for_label(nm + '_only_orig', [z[0] for z in only_orig], [z[1] for z in only_orig])
        out['d6'][nm + '_overlap'] = {'n_both': len(both), 'n_only_new': len(only_new),
                                       'n_only_orig': len(only_orig),
                                       'frac_both': len(both) / len(h6_orig) if h6_orig else 0,
                                       'frac_only_new': len(only_new) / len(new) if new else 0}

    # ---- late-9 check (d=6) ----
    # list pre-defined late zeros; they're stored in the canonical match_results.json
    try:
        late = json.load(open('../sec3_pilot_d6/match_results.json', encoding='utf-8'))
        late_ids = {'S031', 'S036', 'S070', 'S104', 'S187', 'S235', 'S279', 'S319', 'S339'}
        late9 = [r for r in late['catalog_results'] if r['id'] in late_ids]
        out['d6_late9'] = {'list': [{'id': r['id'], 're': r['u_cat_re'], 'im': r['u_cat_im']}
                                    for r in late9]}
        for nm in cats_d6:
            if nm == 'H6_orig342':
                continue
            new = cats_d6[nm]
            new_arr = np.array(new) if new else np.empty((0, 2))
            recovered = []
            for r in late9:
                if len(new_arr) == 0:
                    recovered.append({'id': r['id'], 'matched': False, 'min_d': None})
                    continue
                d = np.abs(new_arr - np.array([r['u_cat_re'], r['u_cat_im']])).max(axis=1)
                m = d.min()
                recovered.append({'id': r['id'], 'matched': bool(m < 1e-6), 'min_d': float(m)})
            out['d6_late9'][nm] = recovered
    except Exception as e:
        out['d6_late9'] = {'error': str(e)}

    # ---- write JSON ----
    def to_json(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, dict):
            return {k: to_json(v) for k, v in o.items()}
        if isinstance(o, list):
            return [to_json(v) for v in o]
        return o

    with open('r1_stats.json', 'w', encoding='utf-8') as f:
        json.dump(to_json(out), f, indent=2, ensure_ascii=False)
    print('wrote r1_stats.json')

    # ---- print summary table ----
    print_summary(out)


def print_summary(out):
    q_list = [16, 24, 32, 48, 64]
    for deg, sec in (('d4', out['d4']), ('d6', out['d6'])):
        print()
        print('=== %s lattice stats ===' % deg)
        hdr = '{:30s} {:>5s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s}'.format(
            'catalog', 'n', 're16', 'im16', 'j16', 're24', 'im24', 're32', 're48', 're64', 'p16')
        print(hdr)
        for nm, rec in sec.items():
            if 'near' not in rec:
                continue
            n = rec['n']
            def v(q, k):
                r = rec['near'][str(q)]
                return '%.1f%%' % (100 * r[k]) if k in r else '-'
            p16 = '-'
            if 'p_re' in rec['near']['16']:
                p16 = '%.1e' % rec['near']['16']['p_re']
            print('{:30s} {:5d} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s} {:>7s}'.format(
                nm, n,
                v(16, 're_near'), v(16, 'im_near'), v(16, 'joint'),
                v(24, 're_near'), v(24, 'im_near'),
                v(32, 're_near'), v(48, 're_near'), v(64, 're_near'),
                p16))
        print()
        print('--- Rayleigh peaks (p<1e-4) ---')
        for nm, rec in sec.items():
            if 'rayleigh_peaks_p1e-4' not in rec:
                continue
            pk = rec['rayleigh_peaks_p1e-4']
            q_peaks_re = [p[0] for p in pk if p[2] < 1e-4]
            q_peaks_im = [p[0] for p in pk if p[4] < 1e-4]
            print('  %30s: re_peaks=%s im_peaks=%s' % (nm, q_peaks_re[:10], q_peaks_im[:10]))


if __name__ == '__main__':
    main()

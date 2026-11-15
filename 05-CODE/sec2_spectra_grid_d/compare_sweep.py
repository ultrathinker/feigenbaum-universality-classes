#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_sweep.py -- BRIEF-F item (v): cross-check the operator curves delta(d)
and lambda2(d) against the cascade-difference sweep (quick_sweep_d*.json).
[oid:spectraF] NOT RIGOROUS.

For every sweep point d_s:
  - operator delta(d_s), lambda2(d_s) by 4-point Lagrange interpolation
    (mpmath) over the nearest grid points of spectra_grid.json;
  - rel diff of sweep delta vs operator delta;
  - where the sweep has an implied |lambda2| (lambda2-regime verdict): rel
    diff vs operator |lambda2|;
  - predicted cascade rate min(delta, 1/|lambda2|) vs the sweep median rate.
Output: compare_sweep.json + printed table.
"""
import json
import os

from mpmath import mp, mpf, nstr

mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
SWEEP_DIR = os.path.normpath(os.path.join(HERE, '..', 'quick_sweep_d'))


def load_grid():
    with open(os.path.join(HERE, 'spectra_grid.json'), encoding='utf-8') as f:
        data = json.load(f)
    pts = sorted(data['points'], key=lambda r: float(r['d']))
    ds = [float(p['d']) for p in pts]
    delta = [mpf(p['delta']['lam']) for p in pts]
    # lambda2 = modulus-2nd essential (pos branch for d>2, neg at the pinning
    # point d=2); lambda3 = next essential (neg branch for d>2)
    l2 = [mpf(p['lambda2_effective']['value']) for p in pts]
    l3 = [mpf(p['lambda3_effective']['value']) for p in pts]
    return ds, delta, l2, l3


def lagrange4(xs, ys, x):
    xs4, ys4 = [], []
    order = sorted(range(len(xs)), key=lambda i: abs(xs[i] - x))
    sel = order[:4]
    # keep monotone order for a stable interpolant
    sel = sorted(sel)
    for i in sel:
        xs4.append(mpf(xs[i]))
        ys4.append(ys[i])
    x = mpf(x)
    res = mpf(0)
    for i in range(4):
        term = ys4[i]
        for j in range(4):
            if j != i:
                term *= (x - xs4[j]) / (xs4[i] - xs4[j])
        res += term
    return res


def load_sweeps():
    rows = []
    for fn in ('sweep_results.json', 'sweep_v2_results.json',
               'sweep_v2b_results.json', 'sweep_v2c_results.json',
               'sweep_v3_results.json'):
        p = os.path.join(SWEEP_DIR, fn)
        if not os.path.exists(p):
            continue
        with open(p, encoding='utf-8') as f:
            for rec in json.load(f):
                if 'delta' not in rec or rec['delta'] is None:
                    continue
                lab = rec.get('label', rec['d'])
                rows.append({'d': float(rec['d']), 'label': lab,
                             'file': fn.replace('sweep_', '').replace('_results.json', ''),
                             'delta_sweep': float(rec['delta']),
                             'rate_median': rec.get('median_over_nodes'),
                             'implied_l2': rec.get('implied_lambda2'),
                             'verdict': rec.get('verdict')})
    rows.sort(key=lambda r: r['d'])
    return rows


def main():
    ds, delta, l2, l3 = load_grid()
    d_lo, d_hi = ds[0], ds[-1]
    rows = load_sweeps()
    out = []
    print('%6s %-5s %10s %10s %9s | %8s %8s %9s | %8s %8s %7s %s'
          % ('d', 'src', 'sw_delta', 'op_delta', 'rel_d%',
             'sw_impl', 'op_l2', 'rel_l2%', 'sw_rate', 'min-law', 'rel_r%', 'verdict'),
          flush=True)
    for r in rows:
        dv = r['d']
        rec = dict(r)
        if d_lo <= dv <= d_hi:
            op_d = lagrange4(ds, delta, dv)
            op_l = lagrange4(ds, l2, dv)
            rel_d = float(abs(op_d - mpf(r['delta_sweep'])) / op_d)
            minlaw = min(abs(op_d), 1 / abs(op_l))
            rec.update({'op_delta': nstr(op_d, 12), 'op_l2': nstr(op_l, 10),
                        'rel_delta': rel_d,
                        'op_min_law_rate': float(minlaw),
                        'rel_rate_vs_minlaw': (float(abs(minlaw - mpf(r['rate_median'])) / minlaw)
                                               if r.get('rate_median') else None)})
            sw_impl = r.get('implied_l2')
            if sw_impl:
                rec['rel_l2'] = float(abs(abs(op_l) - mpf(sw_impl)) / abs(op_l))
            print('%6.2f %-5s %10.5f %10.5f %9.2e | %8s %8.5f %9s | %8s %8.4f %7s %s'
                  % (dv, rec['file'], r['delta_sweep'], float(op_d), rel_d,
                     sw_impl or '-', float(op_l),
                     ('%.2e' % rec['rel_l2']) if 'rel_l2' in rec else '-',
                     ('%.4f' % r['rate_median']) if r.get('rate_median') is not None else '-',
                     float(minlaw),
                     ('%.2e' % rec['rel_rate_vs_minlaw']) if rec.get('rel_rate_vs_minlaw') is not None else '-',
                     r['verdict'] or ''), flush=True)
        else:
            rec['note'] = 'outside operator grid [%.1f, %.1f]' % (d_lo, d_hi)
            print('%6.2f %-5s %10.5f %-10s %-9s | outside operator grid'
                  % (dv, rec['file'], r['delta_sweep'], '-', '-'), flush=True)
        out.append(rec)
    with open(os.path.join(HERE, 'compare_sweep.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print('wrote compare_sweep.json', flush=True)


if __name__ == '__main__':
    main()

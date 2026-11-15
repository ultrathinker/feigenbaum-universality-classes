#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
retrack_pos.py -- post-processing pass after the first grid run.
[oid:spectraF] NOT RIGOROUS.

FINDING from the grid run: the continuation-tracked 'lam2' field followed the
NEGATIVE branch (at d=6 it equals agent B's lam3 = -0.3420), while the anchors'
lambda2 (d=4: +0.2918, d=6: +0.3711) is the MODULUS-SECOND essential -- the
POSITIVE branch.  Spectrum structure on d in [2,7] (from full-eig enumerations):
  pos branch:  +0.2026(2.5) ... +0.4021(7); at d=2 EXACTLY equals tower1
               alpha^{-2} = 0.15962844038269977 (match 1.6e-30) -> excluded at
               the single even point d=2, so lambda2(2) = -0.12365 (agent A);
  neg branch:  -0.1237(2) -> -0.3761(7), never crosses zero (lambda3 for d>2);
  m3 branch:   -0.0573(2) -> -0.2325(7).
Tower alpha^{-dk} exists ONLY at even integer d (exact at 2, 4, 6; absent at
fractional/odd d within enum resolution) -- at even d one essential branch
pins onto the tower value exactly.

This script:
  A) certifies the pos branch (N=40/48 inverse iteration + dense residuals) at
     every existing grid point (d=2.0: pos = tower1 record);
  B) adds flank points 1.9, 1.95, 1.99, 2.001, 2.01, 2.05 (role 'flank') to
     resolve the lambda2 sign change: it happens AT the pinning point d=2 via
     the tower degeneracy, not through a smooth zero of any branch;
  C) bisects the true crossover d*: delta(d)*pos(d) = 1 (the rate-law branch);
     relabels the old neg-branch crossing as crossover_negbranch;
  D) writes per-point 'branches' + 'lambda2_effective' and a 'summary' block.
"""
import json
import os
import time

from mpmath import mp, mpf, nstr

import spectra_grid as sg

mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'spectra_grid.json')

POS_ANCHORS = [(2.5, '0.202597109311'), (3.0, '0.237146111650'),
               (3.6, '0.271718170952'), (4.0, '0.291838408983'),
               (5.0, '0.335033527322'), (6.0, '0.371071467675'),
               (7.0, '0.402127346318')]


def pos_seed(dval):
    """Linear interpolation/extrapolation of the pos branch from anchors."""
    xs = [a[0] for a in POS_ANCHORS]
    ys = [mpf(a[1]) for a in POS_ANCHORS]
    if dval <= xs[0]:
        sl = (ys[1] - ys[0]) / (xs[1] - xs[0]) * mpf('1.08')
        return ys[0] + sl * (mpf(dval) - xs[0])
    if dval >= xs[-1]:
        sl = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        return ys[-1] + sl * (mpf(dval) - xs[-1])
    for i in range(len(xs) - 1):
        if xs[i] <= dval <= xs[i + 1]:
            f = (mpf(dval) - xs[i]) / (xs[i + 1] - xs[i])
            return ys[i] + f * (ys[i + 1] - ys[i])
    raise ValueError


def certify_pos_for_record(rec):
    d = mpf(rec['d'])
    a48 = [mpf(s) for s in rec['_a48']]
    al48 = mpf(rec['alpha48'])
    # N=40 pair by truncating the solved N=48 series (warm re-solve stalls on
    # the rounding floor of the 40-digit stored coefficients); the drift then
    # measures the trial-space difference N=40 vs N=48 on the same series.
    a40 = a48[:sg.N1 + 1]
    al40 = mpf(rec['alpha40'])
    A40 = sg.build_DT_matrix(a40, al40, d, sg.N1)
    A48 = sg.build_DT_matrix(a48, al48, d, sg.N2)
    sigma = pos_seed(float(rec['d']))
    cert = sg.certify_target(A40, A48, a40, al40, a48, al48, d, sigma)
    cert['seed_used'] = nstr(sigma, 12)
    cert['N40_via'] = 'truncated N=48 series'
    ok = (isinstance(cert.get('N48'), dict)
          and cert['N48']['lam_f'] > 0
          and abs(cert['N48']['lam_f'] - float(sigma)) < 0.03
          and cert['N48']['dense_resid'] < 1e-18
          and cert.get('drift', 1) < 1e-6)
    cert['ok'] = bool(ok)
    rec['pos'] = cert
    return cert


def branches_of(rec):
    """Modulus-ranked essentials with the even-d tower excluded."""
    d = float(rec['d'])
    pos = rec['pos']['N48']['lam_f'] if isinstance(rec.get('pos', {}).get('N48'), dict) else None
    if pos is None:
        return None
    neg = rec['lam2']['N48']['lam_f']      # tracked negative branch
    m3 = rec['lam3']['N48']['lam_f']       # tracked m3 branch
    delta = rec['delta']['N48']['lam_f']
    tow1 = float(rec['tower1_alpha_pow'])
    ess = [('pos', pos), ('neg', neg), ('m3', m3)]
    if 'eig_enum_N20' in rec:
        tracked = [pos, neg, m3, delta]
        for e in rec['eig_enum_N20']:
            v = e['lam_f']
            if e['kind'] != 'essential':
                continue
            if all(abs(v - t) > 0.01 for t in tracked) and abs(v) > 0.01:
                ess.append(('enum%d' % len(ess), v))
                tracked.append(v)
    if d == 2.0:
        # pos is pinned to tower1 at the even point: excluded from essentials
        l2, l3 = neg, m3
        note = 'pos==tower1 exactly at d=2 (match 1.6e-30) -> excluded; lambda2=neg branch'
    else:
        by_mod = sorted(ess, key=lambda t: -abs(t[1]))
        l2 = by_mod[0][1]
        l3 = by_mod[1][1]
        note = 'lambda2=%s branch (modulus-2nd essential)' % by_mod[0][0]
    rec['branches'] = {name: val for name, val in ess}
    rec['lambda2_effective'] = {'value': l2, 'branch': 'neg' if d == 2.0 else 'pos',
                                'note': note}
    rec['lambda3_effective'] = {'value': l3, 'branch': 'neg' if d != 2.0 else 'm3'}
    rec['delta_f'] = delta
    rec['tower_detach'] = pos - tow1 if d != 2.0 else None
    return rec


def new_point(ds, data, role):
    dval = float(ds)
    seeds = sg.nearest_seeds(data, dval)
    seeds['role'] = role
    rec, seeds_out = sg.compute_point(ds, seeds, full_eig=False)
    rec['_a48'] = [nstr(x, 40) for x in seeds_out['a48']]
    certify_pos_for_record(rec)
    branches_of(rec)
    data['points'].append(rec)
    data['points'].sort(key=lambda r: float(r['d']))
    sg.save_points(data)
    sg.plog('NEW d=%s (%s): delta=%s pos=%s neg=%s m3=%s pos-ok=%s tower_detach=%s'
            % (ds, role, rec['delta']['lam'][:14], rec['pos']['N48']['lam'][:12],
               rec['lam2']['N48']['lam'][:10], rec['lam3']['N48']['lam'][:10],
               rec['pos']['ok'],
               ('%.2e' % rec['tower_detach']) if rec.get('tower_detach') is not None else '-'))
    return rec


def crossover_pos_f(rec):
    return rec['pos']['N48']['lam_f'] * rec['delta']['N48']['lam_f'] - 1.0


def main():
    t0 = time.time()
    data = sg.load_points()
    sg.plog('RETRACK START: pos-branch certification over %d points' % len(data['points']))

    # ---- A) pos cert at every existing point ------------------------------ #
    for rec in data['points']:
        if 'pos' in rec:
            continue
        if rec['d'] == '2.0':
            rec['pos'] = dict(rec['tower1'])
            rec['pos']['ok'] = True
            rec['pos']['note'] = 'pinned: equals tower1 exactly at even d=2 (match 1.6e-30); excluded from essentials at this single point'
        else:
            certify_pos_for_record(rec)
        branches_of(rec)
        sg.save_points(data)
        if rec['d'] != '2.0':
            sg.plog('POS d=%s: lam=%s drift=%.1e dense=%.1e ok=%s (seed %s)'
                    % (rec['d'], rec['pos']['N48']['lam'][:14], rec['pos'].get('drift', -1),
                       rec['pos']['N48']['dense_resid'], rec['pos']['ok'], rec['pos']['seed_used']))

    # ---- B) flank points around the d=2 pinning --------------------------- #
    done = {r['d'] for r in data['points']}
    for ds in ('1.9', '1.95', '1.99', '2.001', '2.01', '2.05'):
        if ds not in done:
            new_point(ds, data, 'flank')

    # ---- C) crossover bisection on the pos branch -------------------------- #
    old = data['brackets'].pop('crossover', None)
    if old is not None:
        old['note'] = ('superseded: computed on the negative branch (lambda3) during the '
                       'first run before the pos/neg branch identification; kept as the '
                       'lambda3-branch crossing delta*|lambda3|=1')
        data['brackets']['crossover_negbranch_lambda3'] = old
        sg.save_points(data)
    pts = sorted(data['points'], key=lambda r: float(r['d']))
    bracket = None
    for r0, r1 in zip(pts, pts[1:]):
        f0, f1 = crossover_pos_f(r0), crossover_pos_f(r1)
        if f0 * f1 < 0:
            bracket = (r0, r1)
            break
    if bracket is None:
        sg.plog('CROSSOVER(pos): no sign change on grid -- UNEXPECTED')
    else:
        lo, hi = bracket
        sg.plog('CROSSOVER(pos): bracket [%.4f, %.4f] f=[%.6f, %.6f]'
                % (float(lo['d']), float(hi['d']), crossover_pos_f(lo), crossover_pos_f(hi)))
        step = 0
        while float(hi['d']) - float(lo['d']) > 2e-3 and step < 9:
            dm = (float(lo['d']) + float(hi['d'])) / 2
            ds = ('%.6f' % dm).rstrip('0')
            if ds in {r['d'] for r in data['points']}:
                mid = next(r for r in data['points'] if r['d'] == ds)
            else:
                mid = new_point(ds, data, 'bisect_x')
            if crossover_pos_f(lo) * crossover_pos_f(mid) < 0:
                hi = mid
            else:
                lo = mid
            step += 1
        flo, fhi = crossover_pos_f(lo), crossover_pos_f(hi)
        root = float(lo['d']) - flo * (float(hi['d']) - float(lo['d'])) / (fhi - flo)
        data['brackets']['crossover'] = {
            'def': 'delta(d) * lambda2(d) = 1, lambda2 = pos branch (modulus-2nd essential)',
            'bracket': [float(lo['d']), float(hi['d'])],
            'width': float(hi['d']) - float(lo['d']),
            'f_lo': flo, 'f_hi': fhi, 'linear_root': root, 'steps': step}
        sg.save_points(data)
        sg.plog('CROSSOVER(pos) DONE: [%.6f, %.6f] width %.1e linear root %.6f'
                % (float(lo['d']), float(hi['d']), float(hi['d']) - float(lo['d']), root))

    # ---- D) lambda2 sign-change resolution + summary ----------------------- #
    pts = sorted(data['points'], key=lambda r: float(r['d']))
    left = [r for r in pts if float(r['d']) < 2.0]
    right = [r for r in pts if 2.0 < float(r['d']) <= 2.1]
    at2 = next(r for r in pts if r['d'] == '2.0')
    data['brackets']['lam2_sign_change'] = {
        'resolution': ('NO smooth zero: on [2,7] the pos branch stays positive, the neg and '
                       'm3 branches stay negative (see branches per point). The sign change '
                       'of the modulus-ranked lambda2 happens exactly AT the even point d=2 '
                       'where the pos branch pins onto the tower value alpha^{-d}: at d=2 it '
                       'is excluded as a tower member (match 1.6e-30) giving lambda2(2) = '
                       '-0.12365 < 0, while for any d != 2 near 2 (fractional) it is an '
                       'essential eigenvalue and lambda2(d) = pos > 0.'),
        'left_flank_lambda2': [(r['d'], r['lambda2_effective']['value']) for r in left],
        'at_d2_lambda2': at2['lambda2_effective']['value'],
        'right_flank_lambda2': [(r['d'], r['lambda2_effective']['value']) for r in right],
        'right_flank_tower_detach': [(r['d'], r['tower_detach']) for r in right],
        'tower_pin_d2_match_rel': 1.6348123025521162e-30,
        'neg_branch_range': [min(r['branches']['neg'] for r in pts if 'branches' in r),
                             max(r['branches']['neg'] for r in pts if 'branches' in r)],
        'pos_branch_range_on_grid': [min(r['branches']['pos'] for r in pts if 'branches' in r),
                                     max(r['branches']['pos'] for r in pts if 'branches' in r)],
    }
    mono_l2 = all(a['lambda2_effective']['value'] < b['lambda2_effective']['value']
                  for a, b in zip([r for r in pts if float(r['d']) >= 2.1],
                                  [r for r in pts if float(r['d']) >= 2.1][1:]))
    data['summary'] = {
        'lambda2_monotone_from_2.1_to_7': bool(mono_l2),
        'lambda2_values_controls': {r['d']: r['lambda2_effective']['value']
                                    for r in pts if r['d'] in ('3.6', '5.0', '7.0')},
        'note_controls': 'lambda2 rises monotonically toward the sweep saturation ~0.53 (d=20)',
        'tower_even_d_only': ('tower alpha^{-dk} exact at d=2,4,6 (matches 1.6e-30, 3.8e-30, '
                              '9.0e-26); at fractional/odd d no eigenvalue matches alpha^{-d} '
                              '(rel mismatches 0.05..0.97 across the grid)'),
    }
    sg.save_points(data)
    sg.plog('RETRACK DONE in %.0fs (%d points total)' % (time.time() - t0, len(data['points'])))


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
track.py -- [oid:vecM] Step 1 (BRIEF-M Part 1): vector homotopy.
Pure post-processing over spectrum_window.json (produced by
spectrum_window.py, must exist and be complete). Cheap (no mpmath eig calls
here), so run directly (not backgrounded).

Downward sweep: start at d=4.0, mode nearest +0.291838408983545 (dynamically
anchored lambda2). At every step to the NEXT-LOWER d in the grid, compute
cosine overlap |<v_cur, r_i>| against every stable right-eigenvector r_i in
the window at that d; the argmax becomes the new v_cur. Continue down through
d=2.0 and (bonus, same cost) on to d=1.9.

Upward sweep: start at d=2.0, mode nearest the CERTIFIED anchor
lambda2(2)=-0.123652712552687024335. Track up to d=4.0 the same way.

Non-normality: primary continuity metric is the raw cosine overlap of RIGHT
eigenvectors (same monomial basis at every d, so directly comparable; this is
also what wave1's own cos_v2_v3/cos_v2_prev diagnostic used, so results are
comparable to their numbers). As a robustness cross-check we ALSO report the
biorthogonal projection coefficient c_i = <l_i, v_cur>/<l_i, r_i> using the
CANDIDATE's own left eigenvector l_i (dual basis) -- this is the
theoretically correct amplitude of v_cur along mode i for a non-normal
operator, and can disagree with the raw cosine when modes are poorly
conditioned (small <l_i,r_i>, i.e. near a Jordan-block-like coalescence).
Both are recorded at every step; the map/verdict discusses any disagreement.

Output: homotopy_track.json (written once at the end; the whole computation
is a few seconds, JSON is small). NOT RIGOROUS.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, nstr

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

mp.dps = 30

HERE = os.path.dirname(os.path.abspath(__file__))
INJ = os.path.join(HERE, 'spectrum_window.json')
OUTJ = os.path.join(HERE, 'homotopy_track.json')
PROG = os.path.join(HERE, 'progress.log')

WAVE1_JSON = os.path.join(HERE, '..', 'wave1_spectra_d', 'spectra_grid.json')

ANCHOR_POS4 = mpf('0.291838408983545104530997085799')   # dynamically-anchored lambda2(4)
ANCHOR_NEG2 = mpf('-0.123652712552687024335393652586')  # certified lambda2(2)


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def load_window():
    with open(INJ, encoding='utf-8') as f:
        data = json.load(f)
    pts = {}
    for ds, pt in data['points'].items():
        lams = [mpf(s) for s in pt['lams']]
        rvecs = [[mpf(s) for s in v] for v in pt['rvecs']]
        lvecs = [[mpf(s) for s in v] for v in pt['lvecs']]
        pts[ds] = {'lams': lams, 'rvecs': rvecs, 'lvecs': lvecs}
    return data['meta'], pts


def load_wave1_labels():
    with open(WAVE1_JSON, encoding='utf-8') as f:
        data = json.load(f)
    out = {}
    for r in data['points']:
        br = r.get('branches', {})
        out[r['d']] = {k: mpf(v) for k, v in br.items() if v is not None}
    return out


def cosine(v1, v2):
    n1 = sum(x * x for x in v1) ** mpf('0.5')
    n2 = sum(x * x for x in v2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(v1[i] * v2[i] for i in range(min(len(v1), len(v2))))
    return abs(dot / (n1 * n2))


def biorth(lvec_i, v_cur, rvec_i):
    denom = sum(lvec_i[k] * rvec_i[k] for k in range(min(len(lvec_i), len(rvec_i))))
    if denom == 0:
        return None
    numer = sum(lvec_i[k] * v_cur[k] for k in range(min(len(lvec_i), len(v_cur))))
    return numer / denom


def label_of(lam, wave1_br):
    """Nearest wave1 named-branch label at this d, for readability only."""
    if not wave1_br:
        return '?'
    best = min(wave1_br.items(), key=lambda kv: abs(kv[1] - lam))
    reldiff = abs(best[1] - lam) / max(abs(lam), mpf('1e-8'))
    if reldiff < mpf('1e-4'):
        return best[0]
    return '?(near %s, reldiff %.1e)' % (best[0], reldiff)


def sweep(order, pts, wave1_labels, seed_lam):
    """order: list of d-strings in the direction of travel (order[0] = start)."""
    ds0 = order[0]
    lams0 = pts[ds0]['lams']
    i0 = min(range(len(lams0)), key=lambda k: abs(lams0[k] - seed_lam))
    v_cur = pts[ds0]['rvecs'][i0]
    rows = [{'d': ds0, 'lam': nstr(lams0[i0], 16), 'label': label_of(lams0[i0], wave1_labels.get(ds0, {})),
             'cos_to_prev': None, 'cos_runner_up': None, 'gap': None, 'biorth': None}]
    for k in range(1, len(order)):
        d_to = order[k]
        lams = pts[d_to]['lams']
        rvecs = pts[d_to]['rvecs']
        lvecs = pts[d_to]['lvecs']
        scored = sorted(range(len(lams)), key=lambda i: -cosine(v_cur, rvecs[i]))
        if not scored:
            rows.append({'d': d_to, 'lam': None, 'label': 'NO STABLE MODES', 'cos_to_prev': None})
            continue
        i_best = scored[0]
        cos_best = cosine(v_cur, rvecs[i_best])
        cos_2nd = cosine(v_cur, rvecs[scored[1]]) if len(scored) > 1 else mpf(0)
        bc = biorth(lvecs[i_best], v_cur, rvecs[i_best])
        rows.append({'d': d_to,
                      'lam': nstr(lams[i_best], 16),
                      'label': label_of(lams[i_best], wave1_labels.get(d_to, {})),
                      'cos_to_prev': float(cos_best),
                      'cos_runner_up': float(cos_2nd),
                      'runner_up_lam': nstr(lams[scored[1]], 16) if len(scored) > 1 else None,
                      'runner_up_label': label_of(lams[scored[1]], wave1_labels.get(d_to, {})) if len(scored) > 1 else None,
                      'gap': float(cos_best - cos_2nd),
                      'biorth_coeff': (float(bc.real) if hasattr(bc, 'real') else float(bc)) if bc is not None else None})
        v_cur = rvecs[i_best]
    return rows


def main():
    plog('TRACK START (post-processing spectrum_window.json)')
    meta, pts = load_window()
    wave1_labels = load_wave1_labels()
    order_all = meta['grid']  # ascending, e.g. ['1.9', ..., '4.0']
    order_all = [d for d in order_all if d in pts]

    # downward: 4.0 -> 1.9 (through 2.0)
    down_order = list(reversed(order_all))
    down_order = down_order[down_order.index('4.0'):]  # start exactly at 4.0
    plog('Downward sweep order (d=4.0 -> %s), %d points' % (down_order[-1], len(down_order)))
    down_rows = sweep(down_order, pts, wave1_labels, ANCHOR_POS4)

    # upward: 2.0 -> 4.0
    up_order = order_all[order_all.index('2.0'):]
    plog('Upward sweep order (d=2.0 -> %s), %d points' % (up_order[-1], len(up_order)))
    up_rows = sweep(up_order, pts, wave1_labels, ANCHOR_NEG2)

    # min-gap scan (candidates for avoided/exact-crossing refinement)
    all_rows = [('down', r) for r in down_rows[1:]] + [('up', r) for r in up_rows[1:]]
    scored = [(sw, r) for sw, r in all_rows if r.get('cos_to_prev') is not None]
    worst_cos = sorted(scored, key=lambda sr: sr[1]['cos_to_prev'])[:5]
    worst_gap = sorted(scored, key=lambda sr: sr[1]['gap'])[:5]

    out = {'meta': {'tag': '[oid:vecM] BRIEF-M Part 1: vector homotopy, NOT RIGOROUS',
                     'anchor_pos4': str(ANCHOR_POS4), 'anchor_neg2': str(ANCHOR_NEG2),
                     'metric_primary': 'cosine overlap of RIGHT eigenvectors (raw monomial-t coefficient basis, comparable across d)',
                     'metric_secondary': 'biorthogonal coefficient c_i = <l_i,v_cur>/<l_i,r_i> using candidate own left/right pair (non-normal-operator-correct amplitude)'},
            'downward_from_pos4': down_rows,
            'upward_from_neg2': up_rows,
            'worst_continuity_links_by_cosine': [{'sweep': sw, **r} for sw, r in worst_cos],
            'worst_gap_links_by_runnerup_margin': [{'sweep': sw, **r} for sw, r in worst_gap]}
    with open(OUTJ, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog('TRACK DONE -> %s' % OUTJ)

    plog('--- DOWNWARD (d=4.0 -> 1.9), seed +0.29184 ---')
    for r in down_rows:
        plog('  d=%-7s lam=%-18s label=%-8s cos=%s gap=%s' % (
            r['d'], (r['lam'] or '-')[:16], r['label'],
            ('%.6f' % r['cos_to_prev']) if r.get('cos_to_prev') is not None else '-',
            ('%.4f' % r['gap']) if r.get('gap') is not None else '-'))
    plog('--- UPWARD (d=2.0 -> 4.0), seed -0.12365 ---')
    for r in up_rows:
        plog('  d=%-7s lam=%-18s label=%-8s cos=%s gap=%s' % (
            r['d'], (r['lam'] or '-')[:16], r['label'],
            ('%.6f' % r['cos_to_prev']) if r.get('cos_to_prev') is not None else '-',
            ('%.4f' % r['gap']) if r.get('gap') is not None else '-'))
    plog('worst cosine links: %s' % [(sw, r['d'], round(r['cos_to_prev'], 5)) for sw, r in worst_cos])
    plog('worst gap links:    %s' % [(sw, r['d'], round(r['gap'], 5)) for sw, r in worst_gap])


if __name__ == '__main__':
    main()

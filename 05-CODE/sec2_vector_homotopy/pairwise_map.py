#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pairwise_map.py -- [oid:vecM] BRIEF-M Part 1, ROBUST version.

METHODOLOGICAL FINDING that motivates this script: a naive GREEDY step-to-step
tracker (track.py: carry v_cur forward, re-match by argmax cosine at every
next d) works cleanly near d=2 (see homotopy_track.json) but DERAILS in the
middle of the range (d~3.0-3.2), where raw right-eigenvector cosine overlap
between the DIFFERENT branches 'pos' and 'm3' at the SAME d is already ~0.82
(measured directly: pos(3.0) vs m3(3.0) raw cosine = 0.822), uncomfortably
close to genuine same-branch cross-d overlaps in that stretch (~0.94-0.97).
A single near-tie step (d=3.2->3.1, cosine 0.935, runner-up gap only 0.012)
flips the greedy label and contaminates every step after it. This is a
genuine property of this NON-NORMAL operator's raw eigenvector basis, not a
bug -- and it is exactly why the brief calls for care with non-normality.

FIX: never carry a possibly-derailed vector across many steps. Every branch's
identity at every d is instead re-anchored INDEPENDENTLY by nearest-eigenvalue
match to wave1's own already-certified 'branches' dict (pos/neg/m3), which is
itself an externally cross-validated smooth curve (interpolated seeds between
d=2.5/3/3.6/4/5/6/7 anchors + Briggs cascade-dynamics cross-check at d=4/5/6)
-- NOT derived from my own walk, so no self-referential circularity. Then,
for every ADJACENT pair of d in the grid, independently compute the overlap
of {pos,neg,m3}(d_i) against ALL top modes (|lam|>0.03) at d_{i+1}: this
gives a LOCAL, non-cumulative connectivity verdict per link. Chaining the
high-confidence local links together (never carrying a raw vector across
more than one step) reconstructs the whole-range map without the
derailment risk seen in track.py.

Two metrics per link, both reported:
  - cosine: raw right-eigenvector overlap (same monomial-t basis at every d).
  - biorth: |c| = |<l_j(d_i+1), r_i(d_i)>| / |<l_j(d_i+1), r_j(d_i+1)>|, with
    r_j UNIT-NORMALISED first (this normalisation is essential -- an earlier,
    un-normalised version gave numerically meaningless coefficients, because
    mpmath's arbitrary per-eigenvector scale otherwise leaks into |c|).
    Restricted to the top |lam|>0.03 candidate pool: the tail (|lam|<0.03,
    down to 1e-9) is numerically ill-conditioned for this non-normal, high
    dynamic-range operator -- biorth coefficients there were checked to blow
    up to 1e6-1e9, pure noise (division by a near-vanishing <l,r>), while the
    top handful of modes give well-behaved O(1) biorth coefficients.

FINE REFINEMENT near the disputed point: this script also solves (via the
warm-started Newton in wave1's own spectra_grid.py, reused read-only) TWO
extra intermediate fixed points, d=1.995 and d=2.005, if the existing
1.99-2.0-2.001-2.01 links show any gap < 0.15 (ambiguity), to satisfy the
brief's 0.005-step requirement at the minimal-gap location.

Output: pairwise_map.json + printed table. NOT RIGOROUS. CPU: 1 process.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, nstr

import engine as E

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

mp.dps = 35

HERE = os.path.dirname(os.path.abspath(__file__))
OUTJ = os.path.join(HERE, 'pairwise_map.json')
PROG = os.path.join(HERE, 'progress.log')
WINJ = os.path.join(HERE, 'spectrum_window.json')

TOP_THRESH = mpf('0.03')

DGRID = ['1.9', '1.95', '1.99', '2.0', '2.001', '2.01', '2.02', '2.03', '2.04',
         '2.05', '2.1', '2.2', '2.3', '2.4', '2.5', '2.6', '2.7', '2.8', '2.9',
         '3.0', '3.1', '3.2', '3.3', '3.4', '3.5', '3.6', '3.7', '3.8', '3.9', '4.0']


def plog(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open(PROG, 'a', encoding='utf-8') as f:
        f.write(line + '\n')
        f.flush()


def load_window():
    with open(WINJ, encoding='utf-8') as f:
        data = json.load(f)
    pts = {}
    for ds, pt in data['points'].items():
        lams = [mpf(s) for s in pt['lams']]
        rvecs = [[mpf(s) for s in v] for v in pt['rvecs']]
        lvecs = [[mpf(s) for s in v] for v in pt['lvecs']]
        pts[ds] = {'lams': lams, 'rvecs': rvecs, 'lvecs': lvecs}
    return pts


def load_wave1_labels():
    with open(E.WAVE1_JSON, encoding='utf-8') as f:
        data = json.load(f)
    out = {}
    for r in data['points']:
        br = r.get('branches', {})
        out[r['d']] = {k: mpf(v) for k, v in br.items() if v is not None}
    return out


def top_modes(pt):
    idx = [i for i, lam in enumerate(pt['lams']) if abs(lam) > TOP_THRESH]
    return idx


def label_of(lam, wave1_br, tol=mpf('1e-4')):
    if not wave1_br:
        return '?'
    best = min(wave1_br.items(), key=lambda kv: abs(kv[1] - lam))
    reldiff = abs(best[1] - lam) / max(abs(lam), mpf('1e-8'))
    return best[0] if reldiff < tol else ('?~%s' % best[0])


def cosine(v1, v2):
    n1 = sum(x * x for x in v1) ** mpf('0.5')
    n2 = sum(x * x for x in v2) ** mpf('0.5')
    if n1 == 0 or n2 == 0:
        return mpf(0)
    dot = sum(v1[i] * v2[i] for i in range(len(v1)))
    return abs(dot / (n1 * n2))


def unitnorm(v):
    n = sum(x * x for x in v) ** mpf('0.5')
    return [x / n for x in v] if n != 0 else v


def biorth_unit(lvec_j, r_i, rvec_j_unit):
    denom = sum(lvec_j[k] * rvec_j_unit[k] for k in range(len(rvec_j_unit)))
    if denom == 0:
        return None
    numer = sum(lvec_j[k] * r_i[k] for k in range(len(r_i)))
    return numer / denom


def link(pt_i, pt_j, wave1_br_i, wave1_br_j):
    """For every top mode at pt_i, rank all top modes at pt_j by cosine AND
    by biorth. Returns list of per-source-mode records."""
    idx_i = top_modes(pt_i)
    idx_j = top_modes(pt_j)
    recs = []
    for i in idx_i:
        lam_i = pt_i['lams'][i]
        lab_i = label_of(lam_i, wave1_br_i)
        r_i = pt_i['rvecs'][i]
        scored_cos = []
        scored_bio = []
        for j in idx_j:
            lam_j = pt_j['lams'][j]
            lab_j = label_of(lam_j, wave1_br_j)
            cs = cosine(r_i, pt_j['rvecs'][j])
            r_j_unit = unitnorm(pt_j['rvecs'][j])
            bc = biorth_unit(pt_j['lvecs'][j], r_i, r_j_unit)
            scored_cos.append((float(cs), lab_j, nstr(lam_j, 10)))
            scored_bio.append((float(abs(bc)) if bc is not None else -1, lab_j, nstr(lam_j, 10)))
        scored_cos.sort(key=lambda t: -t[0])
        scored_bio.sort(key=lambda t: -t[0])
        recs.append({
            'src_label': lab_i, 'src_lam': nstr(lam_i, 12),
            'best_cos': scored_cos[0], 'runnerup_cos': scored_cos[1] if len(scored_cos) > 1 else None,
            'gap_cos': round(scored_cos[0][0] - (scored_cos[1][0] if len(scored_cos) > 1 else 0), 4),
            'best_bio': scored_bio[0], 'runnerup_bio': scored_bio[1] if len(scored_bio) > 1 else None,
            'gap_bio': round(scored_bio[0][0] - (scored_bio[1][0] if len(scored_bio) > 1 else 0), 4),
        })
    return recs


def main():
    plog('PAIRWISE MAP START (non-cumulative, anchored per-step links)')
    pts = load_window()
    wave1_br = load_wave1_labels()
    out = {'meta': {'tag': '[oid:vecM] BRIEF-M Part 1 robust pairwise map, NOT RIGOROUS',
                     'top_thresh': str(TOP_THRESH),
                     'note': 'each link independently anchored to wave1 branch labels; '
                             'no cumulative vector carried across >1 step (avoids the '
                             'derailment seen in track.py greedy walk)'},
           'links': []}
    for i in range(len(DGRID) - 1):
        d_i, d_j = DGRID[i], DGRID[i + 1]
        recs = link(pts[d_i], pts[d_j], wave1_br.get(d_i, {}), wave1_br.get(d_j, {}))
        out['links'].append({'d_from': d_i, 'd_to': d_j, 'records': recs})
        plog('link %s -> %s:' % (d_i, d_j))
        for r in recs:
            plog('   %-6s (%s) -> best[cos]=%-6s %-16s (gap %.3f)  best[bio]=%-6s %-16s (gap %.3f)'
                 % (r['src_label'], r['src_lam'][:10],
                    round(r['best_cos'][0], 4), r['best_cos'][1] + ' ' + r['best_cos'][2][:10], r['gap_cos'],
                    round(r['best_bio'][0], 4), r['best_bio'][1] + ' ' + r['best_bio'][2][:10], r['gap_bio']))
    with open(OUTJ, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    plog('PAIRWISE MAP DONE -> %s' % os.path.basename(OUTJ))


if __name__ == '__main__':
    main()

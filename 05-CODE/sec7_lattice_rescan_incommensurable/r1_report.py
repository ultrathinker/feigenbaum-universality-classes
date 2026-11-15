#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r1_report.py -- generate REPORT-R1-LATTICE-RESCAN.md from r1_stats.json.

Reads:
  r1_stats.json
  /tmp/.../H4_*.json, H6_*.json catalogs (already-loaded by r1_stats.py)
  match_results.json (for late-9)
Writes:
  REPORT-R1-LATTICE-RESCAN.md
"""
import json
import os
import sys
import time

import numpy as np


def fmt_pct(x):
    return '%.1f%%' % (100 * float(x)) if x is not None else 'n/a'


def fmt_p(p):
    if p is None:
        return 'n/a'
    if p < 1e-4:
        return '%.1e' % p
    return '%.3f' % p


def table_row(name, rec, q_subset=(16, 24, 32, 48, 64)):
    n = rec.get('n', 0)
    s = rec.get('near', {}).get('16', {})
    out = ['| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |'
           % (name, n,
              fmt_pct(s.get('re_near')),
              fmt_pct(s.get('im_near')),
              fmt_pct(s.get('joint')),
              fmt_pct(rec.get('near', {}).get('24', {}).get('re_near')),
              fmt_pct(rec.get('near', {}).get('24', {}).get('im_near')),
              fmt_pct(rec.get('near', {}).get('32', {}).get('re_near')),
              fmt_pct(rec.get('near', {}).get('48', {}).get('re_near')),
              fmt_pct(rec.get('near', {}).get('64', {}).get('re_near')),
              fmt_p(s.get('p_re')))]
    return out


def peaks_str(rec):
    if 'rayleigh_peaks_p1e-4' not in rec:
        return '-'
    pk = rec['rayleigh_peaks_p1e-4']
    re_peaks = [int(p[0]) for p in pk if p[2] < 1e-4]
    im_peaks = [int(p[0]) for p in pk if p[4] < 1e-4]
    return 're=%s, im=%s' % (re_peaks[:6], im_peaks[:6])


def main():
    stats = json.load(open('r1_stats.json', encoding='utf-8'))
    out = []
    out.append('# R1 — Решающий тест: «решётка 1/16» нулей H_d(u) — реальность или артефакт?\n')
    out.append('Дата: %s\n' % stats.get('date'))
    out.append('Допуск матчинга: |Δu|<1e-6; дедуп: |Δu|<1e-8; seed=20260816\n')
    out.append('\n')

    # ----- d=4 table -----
    out.append('## d=4: статистика по каталогам и подмножествам\n')
    out.append('| каталог | n | Re16 | Im16 | joint16 | Re24 | Im24 | Re32 | Re48 | Re64 | p_re16 |\n')
    out.append('|---|---|---|---|---|---|---|---|---|---|---|\n')
    for nm, rec in stats['d4'].items():
        if 'near' not in rec:
            continue
        out.append('| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n'
                   % (nm, rec.get('n', 0),
                      fmt_pct(rec['near']['16'].get('re_near')),
                      fmt_pct(rec['near']['16'].get('im_near')),
                      fmt_pct(rec['near']['16'].get('joint')),
                      fmt_pct(rec['near']['24'].get('re_near')),
                      fmt_pct(rec['near']['24'].get('im_near')),
                      fmt_pct(rec['near']['32'].get('re_near')),
                      fmt_pct(rec['near']['48'].get('re_near')),
                      fmt_pct(rec['near']['64'].get('re_near')),
                      fmt_p(rec['near']['16'].get('p_re'))))
    out.append('\nПики Рэлея (p<1e-4):\n')
    for nm, rec in stats['d4'].items():
        if 'rayleigh_peaks_p1e-4' in rec:
            out.append('* `%s`: %s\n' % (nm, peaks_str(rec)))
    out.append('\n')

    # ----- d=6 table -----
    out.append('## d=6: статистика по каталогам и подмножествам\n')
    out.append('| каталог | n | Re16 | Im16 | joint16 | Re24 | Im24 | Re32 | Re48 | Re64 | p_re16 |\n')
    out.append('|---|---|---|---|---|---|---|---|---|---|---|\n')
    for nm, rec in stats['d6'].items():
        if 'near' not in rec:
            continue
        out.append('| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n'
                   % (nm, rec.get('n', 0),
                      fmt_pct(rec['near']['16'].get('re_near')),
                      fmt_pct(rec['near']['16'].get('im_near')),
                      fmt_pct(rec['near']['16'].get('joint')),
                      fmt_pct(rec['near']['24'].get('re_near')),
                      fmt_pct(rec['near']['24'].get('im_near')),
                      fmt_pct(rec['near']['32'].get('re_near')),
                      fmt_pct(rec['near']['48'].get('re_near')),
                      fmt_pct(rec['near']['64'].get('re_near')),
                      fmt_p(rec['near']['16'].get('p_re'))))
    out.append('\nПики Рэлея (p<1e-4):\n')
    for nm, rec in stats['d6'].items():
        if 'rayleigh_peaks_p1e-4' in rec:
            out.append('* `%s`: %s\n' % (nm, peaks_str(rec)))
    out.append('\n')

    # ----- overlap table -----
    out.append('## Совпадения ORIG ↔ new (|Δu|<1e-6)\n')
    out.append('| deg | new | n_new | n_orig | n_both | n_only_new | n_only_orig | frac_both_of_orig | frac_only_new |\n')
    out.append('|---|---|---|---|---|---|---|---|---|\n')
    for nm, ov in stats['d4'].items():
        if '_overlap' not in nm:
            continue
        out.append('| d=4 | %s | %s | %s | %s | %s | %s | %s | %s |\n'
                   % (nm.replace('_overlap', ''),
                      stats['d4'].get(nm.replace('_overlap', ''), {}).get('n', '-'),
                      62, ov['n_both'], ov['n_only_new'], ov['n_only_orig'],
                      '%.1f%%' % (100 * ov['frac_both']),
                      '%.1f%%' % (100 * ov['frac_only_new'])))
    for nm, ov in stats['d6'].items():
        if '_overlap' not in nm:
            continue
        out.append('| d=6 | %s | %s | %s | %s | %s | %s | %s | %s |\n'
                   % (nm.replace('_overlap', ''),
                      stats['d6'].get(nm.replace('_overlap', ''), {}).get('n', '-'),
                      342, ov['n_both'], ov['n_only_new'], ov['n_only_orig'],
                      '%.1f%%' % (100 * ov['frac_both']),
                      '%.1f%%' % (100 * ov['frac_only_new'])))
    out.append('\n')

    # ----- late-9 -----
    if 'd6_late9' in stats and 'list' in stats['d6_late9']:
        out.append('## d=6: «поздние» 9 нулей — найдены ли новым сканом?\n')
        out.append('| id | u_cat (re, im) |')
        for k in stats['d6_late9'].keys():
            if k == 'list':
                continue
            out.append(' %s |' % k.replace('H6_', '').replace('_rescan_', '_'))
        out.append('\n')
        out.append('|---|---|')
        for k in stats['d6_late9'].keys():
            if k == 'list':
                continue
            out.append('---|')
        out.append('\n')
        ids = [r['id'] for r in stats['d6_late9']['list']]
        for i, lid in enumerate(ids):
            row = stats['d6_late9']['list'][i]
            out.append('| %s | %.6f, %.6f |' % (lid, row['re'], row['im']))
            for k in stats['d6_late9'].keys():
                if k == 'list':
                    continue
                r = stats['d6_late9'][k][i]
                if r['matched']:
                    out.append(' yes (d=%.1e) |' % r['min_d'])
                else:
                    out.append(' no (d=%s) |' % ('%.2e' % r['min_d'] if r['min_d'] is not None else '-'))
            out.append('\n')
        out.append('\n')

    # ----- timings -----
    if os.path.exists('timing_steps.txt'):
        out.append('## Время шагов\n')
        out.append('```\n')
        out.append(open('timing_steps.txt', encoding='utf-8').read())
        out.append('```\n')

    # ----- answer -----
    out.append('## Ответ на прямой вопрос (п.7 BRIEF)\n')
    out.append('Вписывается вручную после анализа таблиц выше.\n\n')
    out.append('### Вердикт (REAL / ARTIFACT / AMBIGUOUS)\n')
    out.append('*\n\n')
    out.append('### Объяснение механизма\n')
    out.append('*\n')

    out_text = ''.join(out)
    with open('REPORT-R1-LATTICE-RESCAN.md', 'w', encoding='utf-8') as f:
        f.write(out_text)
    print('wrote REPORT-R1-LATTICE-RESCAN.md (%d chars)' % len(out_text))


if __name__ == '__main__':
    main()

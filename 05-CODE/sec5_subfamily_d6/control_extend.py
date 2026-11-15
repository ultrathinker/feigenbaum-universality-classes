#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
control_extend.py -- BRIEF-G item 5 support: extend the 13 normal-control
zeros' d_u(m) table (match_results.json rate_study has m=8,9,10 only) with
m=7 and m=11, so the control fit uses the same level range as the anomalous
targets. Cheap: these are well-behaved catalog zeros (direct-seed Newton
converges almost always in the main pilot, 332/333 tier-1).
[oid:deg6G] NOT RIGOROUS.
"""
import json
import os
import time

from mpmath import mp, mpf, mpc

mp.dps = 80

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
T0 = time.time()


def log(msg):
    print('[%7.1fs] %s' % (time.time() - T0, msg), flush=True)


MR = json.load(open(os.path.join(BASE, 'match_results.json'), encoding='utf-8'))
rate_study = MR['rate_study']

out_path = os.path.join(HERE, 'control_extra.json')
out = json.load(open(out_path, encoding='utf-8')) if os.path.exists(out_path) else {}

for r in rate_study:
    tid = r['id']
    if tid in out and '7' in out[tid] and '11' in out[tid]:
        continue
    u = complex(0, 0)
    # recover u from id via H6_zeros catalog (rate_study doesn't store re/im directly here)
    HZ = json.load(open(os.path.join(BASE, 'H6_zeros.json'), encoding='utf-8'))
    idx = int(tid[1:]) - 1
    z = HZ['zeros_limit_extrapolated'][idx]
    u_mp = mpc(mpf(z['re']), mpf(z['im']))
    row = out.get(tid, {})
    for m in (7, 11):
        if str(m) in row:
            continue
        d, a_ss, tag = C.best_du(u_mp, m)
        row[str(m)] = d
        log('%s m=%d d_u=%s tag=%s' % (tid, m, d, tag))
    out[tid] = row
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)

log('DONE control_extend')

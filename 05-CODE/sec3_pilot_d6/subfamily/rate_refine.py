#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rate_refine.py -- BRIEF-G item 5 support: fill the rate-curve gaps left by
best_du's cheap 2-seed search (mside_scan.py) using the FULL tiered seed
search (common.full_search_du, ported from p4_match.py: direct + 3x3 a-grid
+ radial fan + cascade).  Only 9 zeros x 5 levels (m=7..11) -- bounded cost.
Run AFTER the two mside_scan.py background processes finish (keeps total
CPU usage at <=2 concurrent processes per budget: this one runs alone).
[oid:deg6G] NOT RIGOROUS.
"""
import json
import os
import time

from mpmath import mp, mpf, mpc

mp.dps = 80

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
plog = open(os.path.join(HERE, 'progress_raterefine.log'), 'a', encoding='utf-8')


def log(msg):
    line = '[%7.1fs] %s' % (time.time() - T0, msg)
    print(line, flush=True)
    plog.write(line + '\n')
    plog.flush()
    os.fsync(plog.fileno())


targets = json.load(open(os.path.join(HERE, 'targets.json'), encoding='utf-8'))['targets']
ALL_IDS = list(targets.keys())

out_path = os.path.join(HERE, 'rate_full.json')
out = {}
if os.path.exists(out_path):
    out = json.load(open(out_path, encoding='utf-8'))

log('START full-tier rate refinement, %d zeros x levels 7..11' % len(ALL_IDS))
for tid in ALL_IDS:
    if tid in out and all(str(m) in out[tid]['d_u'] for m in (7, 8, 9, 10, 11)):
        log('%s already complete, skip' % tid)
        continue
    t = targets[tid]
    u = mpc(mpf(t['u_re']), mpf(t['u_im']))
    row = out.get(tid, {'id': tid, 'd_u': {}, 'tag': {}, 'sec': {}})
    for m in (7, 8, 9, 10, 11):
        if str(m) in row['d_u']:
            continue
        t0 = time.time()
        d, a_ss, tag = C.full_search_du(u, m)
        dt = time.time() - t0
        row['d_u'][str(m)] = d
        row['tag'][str(m)] = tag
        row['sec'][str(m)] = dt
        log('%-6s m=%2d d_u=%s tag=%s (%.2fs)'
            % (tid, m, ('%.4e' % d) if d is not None else 'None', tag, dt))
    out[tid] = row
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log('%s done, wrote rate_full.json' % tid)

log('DONE rate_refine, total %.1fs' % (time.time() - T0))

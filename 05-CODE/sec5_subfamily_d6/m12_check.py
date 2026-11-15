#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""m12_check.py -- BRIEF-G supplementary: push the priority quartet's rate
curve to m=12 to see whether S036's m=10/m=11 plateau (d_u ~0.022-0.023,
essentially flat) continues, and whether S104/S187/S339's clean m=11
recovery (d_u~1e-4, direct tier) holds up one level further.
[oid:deg6G] NOT RIGOROUS."""
import json, os, time
from mpmath import mp, mpf, mpc
mp.dps = 80
import common as C
HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
targets = json.load(open(os.path.join(HERE,'targets.json'), encoding='utf-8'))['targets']
out_path = os.path.join(HERE,'rate_full.json')
out = json.load(open(out_path, encoding='utf-8'))
for tid in ['S036','S104','S187','S339']:
    t = targets[tid]
    u = mpc(mpf(t['u_re']), mpf(t['u_im']))
    t0=time.time()
    d, a_ss, tag = C.full_search_du(u, 12)
    dt = time.time()-t0
    out[tid]['d_u']['12'] = d
    out[tid]['tag']['12'] = tag
    out[tid]['sec']['12'] = dt
    print('%-6s m=12 d_u=%s tag=%s (%.2fs)' % (tid, ('%.4e'%d) if d is not None else 'None', tag, dt), flush=True)
    json.dump(out, open(out_path,'w',encoding='utf-8'), indent=2, ensure_ascii=False)
print('DONE m12_check, total %.1fs' % (time.time()-T0))

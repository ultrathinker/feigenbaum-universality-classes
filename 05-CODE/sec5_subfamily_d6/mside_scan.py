#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mside_scan.py -- BRIEF-G, expensive M-side leg (single-threaded background
process; run twice -- once per group -- so total CPU usage stays at 2
processes/threads, per budget).

For each target zero:
  (A) rate curve: d_u(m) at m = 7,8,9,10,11 via mpmath Newton (best_du,
      direct + cascade-from-m-1 seeds, exact-period-checked residual-based
      convergence -- brief item 4, never trust step size).
  (B) cross-hair scan: 1D Newton scan along Re(u) and Im(u) through the
      catalog point, window +/-1 local inter-zero spacing, step =
      spacing/STEP_DIV, at levels m=9 and m=10.  Tests hypothesis (c):
      is there ANOTHER nearby u where the M-side correspondence survives
      with comparable or smaller d_u (a rival/second branch), or is the
      catalog point the unique local optimum?

Usage: python mside_scan.py <group_name> <STEP_DIV> <id1,id2,...>
  group_name = 'priority' | 'secondary' (used for filenames only)
  STEP_DIV   = grid points per spacing unit (64 = full density, 32 = half)

Output (incremental, flushed after every point):
  progress_<group>.log   -- human progress with block/point/time/ETA
  scan_<group>.jsonl      -- cross-hair scan, one JSON record per point
  rate_<group>.json       -- rate-curve table (per zero, all m at once)
[oid:deg6G] NOT RIGOROUS.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf, mpc, fabs, nstr

mp.dps = 80

import common as C

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    group = sys.argv[1]
    step_div = int(sys.argv[2])
    ids = sys.argv[3].split(',')

    T0 = time.time()
    plog_path = os.path.join(HERE, 'progress_%s.log' % group)
    plog = open(plog_path, 'a', encoding='utf-8')

    def log(msg):
        line = '[%7.1fs] [%s] %s' % (time.time() - T0, group, msg)
        print(line, flush=True)
        plog.write(line + '\n')
        plog.flush()
        os.fsync(plog.fileno())

    targets = json.load(open(os.path.join(HERE, 'targets.json'), encoding='utf-8'))['targets']

    log('START group=%s step_div=%d ids=%s' % (group, step_div, ids))

    # ---------------- (A) rate curve m=7..11 -------------------------------
    rate_out_path = os.path.join(HERE, 'rate_%s.json' % group)
    rate_out = {}
    if os.path.exists(rate_out_path):
        rate_out = json.load(open(rate_out_path, encoding='utf-8'))
    for tid in ids:
        if tid in rate_out:
            continue
        t = targets[tid]
        u = mpc(mpf(t['u_re']), mpf(t['u_im']))
        row = {'id': tid, 'd_u': {}, 'sec': {}}
        for m in (7, 8, 9, 10, 11):
            t0 = time.time()
            d, a_ss, tag = C.best_du(u, m)
            dt = time.time() - t0
            row['d_u'][str(m)] = d
            row['sec'][str(m)] = dt
            log('rate %-6s m=%2d d_u=%s tag=%s (%.2fs)'
                % (tid, m, ('%.4e' % d) if d is not None else 'None', tag, dt))
        rate_out[tid] = row
        with open(rate_out_path, 'w', encoding='utf-8') as f:
            json.dump(rate_out, f, indent=2, ensure_ascii=False)
        log('rate curve done for %s, wrote %s' % (tid, os.path.basename(rate_out_path)))
    log('=== rate curves complete (%d zeros) ===' % len(ids))

    # ---------------- (B) cross-hair scan -----------------------------------
    jsonl_path = os.path.join(HERE, 'scan_%s.jsonl' % group)
    done_keys = set()
    if os.path.exists(jsonl_path):
        with open(jsonl_path, encoding='utf-8') as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    done_keys.add((rec['id'], rec['m'], rec['axis'], rec['k']))
                except Exception:
                    pass
        log('resume: %d points already done in %s' % (len(done_keys), jsonl_path))

    jf = open(jsonl_path, 'a', encoding='utf-8')

    # build full point list first (for ETA)
    plan = []
    for tid in ids:
        t = targets[tid]
        u0 = mpc(mpf(t['u_re']), mpf(t['u_im']))
        spacing = mpf(t['nearest_neighbour_spacing'])
        step = spacing / step_div
        for m in (9, 10):
            for axis in ('re', 'im'):
                for k in range(-step_div, step_div + 1):
                    if axis == 'im' and k == 0:
                        continue  # k=0 shared with 're' axis (the catalog point itself)
                    plan.append((tid, m, axis, k, u0, step))
    total = len(plan)
    log('plan: %d cross-hair points total (step_div=%d)' % (total, step_div))

    n_done = 0
    n_conv = 0
    t_start = time.time()
    for (tid, m, axis, k, u0, step) in plan:
        key = (tid, m, axis, k)
        if key in done_keys:
            n_done += 1
            continue
        if axis == 're':
            u = u0 + step * k
        else:
            u = u0 + mpc(0, 1) * step * k
        t0 = time.time()
        d, a_ss, tag = C.best_du(u, m)
        dt = time.time() - t0
        rec = {'id': tid, 'm': m, 'axis': axis, 'k': k,
               'u_re': float(u.real), 'u_im': float(u.imag),
               'd_u': d, 'tier': tag,
               'a_ss_re': (nstr(a_ss.real, 25) if a_ss is not None else None),
               'a_ss_im': (float(a_ss.imag) if a_ss is not None else None)}
        jf.write(json.dumps(rec) + '\n')
        jf.flush()
        os.fsync(jf.fileno())
        n_done += 1
        if d is not None:
            n_conv += 1
        if n_done % 25 == 0 or n_done == total:
            el = time.time() - t_start
            rate = el / max(1, (n_done - len(done_keys)))
            remaining = total - n_done
            eta = rate * remaining
            log('block %d/%d done (conv=%d) elapsed=%.0fs ETA=%.0fs (%.1fmin)'
                % (n_done, total, n_conv, el, eta, eta / 60.0))

    log('=== cross-hair scan complete: %d/%d points, %d converged+period-ok ==='
        % (n_done, total, n_conv))
    log('DONE group=%s total_elapsed=%.1fs' % (group, time.time() - T0))
    plog.close()
    jf.close()


if __name__ == '__main__':
    main()

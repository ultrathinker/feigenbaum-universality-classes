# -*- coding: utf-8 -*-
# [oid:kraw6L] extract_targets_late.py -- extract the 9 "late" subfamily
# superstable targets (seeds + levels m) from the agentG subfamily artifacts.
#
# Source of truth for WHICH m each target needs: REPORT_SUBFAMILY.md (the
# "first clean tier-1 hit" table) cross-checked against rate_full.json d_u(m).
# Per REPORT_SUBFAMILY.md:
#   S031 S070 S235 S279 S319 -> m = 10   (period 2^11 = 2048)
#   S104 S187 S339           -> m = 11   (period 2^12 = 4096)
#   S036                     -> m = 12   (period 2^13 = 8192)
# The seed is the CATALOG prediction a(u_cat, m) = A_F + u_cat/(A6*delta6^m)
# with the blind-run constants a_F/A6/delta6 from targets.json (strings, kept
# at full precision here).  The refined seeds (Newton at the true period, high
# precision) are produced separately by refine_late.py.
#
# This script only reads (read-only rule) and writes targets_late.json.
import json
import os
import sys
import time
from pathlib import Path

from mpmath import mp, mpc, mpf, re as mpre, im as mpim, nstr

HERE = Path(__file__).resolve().parent
SUB = Path(r"../sec5_subfamily_d6")
PROGRESS = HERE / "progress.log"

# frozen assignment (from REPORT_SUBFAMILY.md first-clean-tier-1-hit table):
#   label -> (m, period, priority/edge tags from targets.json)
TARGETS_M = {
    "S031": 10, "S070": 10, "S235": 10, "S279": 10, "S319": 10,
    "S104": 11, "S187": 11, "S339": 11,
    "S036": 12,
}

MP_DPS = 100


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    mp.dps = MP_DPS
    tj = json.loads((SUB / "targets.json").read_text(encoding="utf-8"))
    rf = json.loads((SUB / "rate_full.json").read_text(encoding="utf-8"))

    a_F = mpf(tj["a_F"])
    A6 = mpf(tj["A6"])
    d6 = mpf(tj["delta6"])

    # sanity: the REPORT table must match rate_full.json d_u(m): at its
    # assigned m each target has its FIRST d_u < 0.01 (the blind-run accept
    # threshold) at-or-after m=9, and is "direct" tier-1.
    targets = {}
    for lbl, m in sorted(TARGETS_M.items()):
        t = tj["targets"][lbl]
        u = mpc(mpf(t["u_re"]), mpf(t["u_im"]))
        a_seed = a_F + u / (A6 * d6 ** m)
        d_u = rf[lbl]["d_u"].get(str(m))
        tag = rf[lbl]["tag"].get(str(m))
        if d_u is None or d_u >= 0.01:
            print("WARN %s: d_u(m=%d)=%s not < 0.01" % (lbl, m, d_u))
        if d_u is not None and d_u >= 0.01:
            raise RuntimeError("%s: assigned m=%d has d_u=%s >= 0.01"
                               % (lbl, m, d_u))
        period = 1 << (m + 1)
        targets[lbl] = {
            "id": lbl,
            "m": m,
            "period": period,
            "n": period,
            "priority": t["priority"],
            "edge": t["edge"],
            "u_cat_re": str(t["u_re"]),
            "u_cat_im": str(t["u_im"]),
            # seed = catalog prediction at the TRUE level m (NOT m=9): the
            # direct (tier-1) starting point that the subfamily report's
            # Newton used.
            "seed_re": nstr(mpre(a_seed), 50),
            "seed_im": nstr(mpim(a_seed), 50),
            "seed_policy": "predicted_catalog_at_m",
            "rate_d_u": d_u,
            "rate_tag": tag,
        }
        print("%s: m=%d P=%d d_u=%.3e tag=%s seed=(%s, %s)"
              % (lbl, m, period, d_u, tag,
                 targets[lbl]["seed_re"][:28], targets[lbl]["seed_im"][:28]))

    out = {
        "run": "extract_targets_late",
        "oid": "kraw6L",
        "note": "9 late subfamily targets; m assigned from "
                "REPORT_SUBFAMILY.md first-clean-tier-1-hit table, "
                "cross-checked against rate_full.json d_u(m)<0.01; seed = "
                "catalog prediction a(u_cat,m)=A_F+u/(A6*delta6^m)",
        "a_F": str(a_F), "A6": str(A6), "delta6": str(d6),
        "targets": targets,
    }
    out_path = HERE / "targets_late.json"
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    log_progress("STEP extract_targets_late: %d targets written to %s "
                 "(m=[10,10,10,10,10,11,11,11,12])"
                 % (len(targets), out_path.name))
    print("wrote %s (%d targets)" % (out_path, len(targets)))
    return 0


if __name__ == "__main__":
    sys.exit(main())

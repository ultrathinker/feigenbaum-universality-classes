# -*- coding: utf-8 -*-
# [oid:kraw6H] refine6.py -- deterministic Newton refinement of the two
# blind-run seeds whose float64-im storage cannot be certified directly
# (S078, S090: |F'(a*)| ~ 8e14 / 1.1e18, the two largest in the catalog; the
# F'(X)-interval blows up across 0 for r >= 1e-25 / 1e-29, while the stored
# seed sits 1.7e-24 / 1.9e-24 from the root -- the float64 imaginary part
# carries only ~7 significant decimals at im ~ 3e-8).
#
# NOT the rigorous layer: rigor comes from the Krawczyk certificate at the
# refined center (krawczyk6_certificate.py).  What this script fixes is a
# DETERMINISTIC center: exactly REFINE_ITERS damped-free Newton steps on
# x_1024(a) = 0 at mp.dps = REFINE_DPS, starting from the stored seed strings,
# then truncation to OUT_DIGITS significant decimals.  The certificate pins
# sha256(refined_seeds6.json) and proves in Arb that the refined center lies
# within SEED_DIST_LIMIT of the stored seed (identity of roots is then forced
# by pairwise isolation: nearest foreign box is 2.6e-10 away).
import json
import sys
import time
from pathlib import Path

from mpmath import mp, mpf, mpc, fabs, nstr, re as mpre, im as mpim

HERE = Path(__file__).resolve().parent
PILOT = Path(r"../sec3_pilot_d6/match_results.json")
OUT = HERE / "refined_seeds6.json"
PROGRESS = HERE / "progress.log"

REFINE_LABELS = ("S078", "S090")
REFINE_DPS = 120
REFINE_ITERS = 8          # fixed: deterministic output
OUT_DIGITS = 60           # significant decimals kept for the box center
SEED_DIST_LIMIT = "5e-24" # Arb pin: |refined - stored| must stay below this


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def orbit(a, n):
    x, dx = mpc(0), mpc(0)
    for _ in range(n):
        dx = -(x ** 6) - 6 * a * x ** 5 * dx
        x = 1 - a * x ** 6
    return x, dx


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    mp.dps = REFINE_DPS
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    rows = {r["id"]: r for r in d["catalog_results"]}
    out = {
        "run": "refine6", "oid": "kraw6H",
        "note": "deterministic Newton refinement (REFINE_ITERS=%d plain steps, "
                "dps=%d, centers truncated to %d significant decimals) of the "
                "float64-im-limited seeds; rigor lives in the Krawczyk "
                "certificate, not here" % (REFINE_ITERS, REFINE_DPS, OUT_DIGITS),
        "seed_dist_limit": SEED_DIST_LIMIT,
        "seeds": {},
    }
    for lbl in REFINE_LABELS:
        rec = rows[lbl]
        if not (rec.get("accepted") and not rec.get("edge")):
            raise RuntimeError("%s is not an interior accepted row" % lbl)
        a = mpc(mpf(rec["a_ss_re"]), mpf(rec["a_ss_im"]))
        a_stored = a
        trace = []
        for it in range(REFINE_ITERS):
            x, dx = orbit(a, 1024)
            step = x / dx if dx != 0 else mpc(mp.inf, mp.inf)
            a = a - step
            trace.append({"iter": it + 1,
                          "residual_abs_x1024": float(fabs(x)),
                          "abs_step": float(fabs(step))})
        x, dx = orbit(a, 1024)
        re_s = nstr(mpre(a), OUT_DIGITS)
        im_s = nstr(mpim(a), OUT_DIGITS)
        dist = fabs(a - a_stored)
        out["seeds"][lbl] = {
            "re": re_s, "im": im_s,
            "dist_to_stored_seed": nstr(dist, 30),
            "final_residual_abs_x1024": nstr(fabs(x), 10),
            "trace": trace,
        }
        print("%s: refined, |a_ref - a_stored| = %s (limit %s), final "
              "|x_1024| = %s" % (lbl, nstr(dist, 8), SEED_DIST_LIMIT,
                                 nstr(fabs(x), 5)))
        for t in trace[:3] + trace[-1:]:
            print("   it%d |step|=%.3e" % (t["iter"], t["abs_step"]))
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    log_progress("STEP refine6: deterministic refinement of %s -> %s "
                 "(iters=%d, dps=%d, %d digits, dist limit %s)"
                 % (", ".join(REFINE_LABELS), OUT.name, REFINE_ITERS,
                    REFINE_DPS, OUT_DIGITS, SEED_DIST_LIMIT))
    print("wrote %s" % OUT)
    return 0


if __name__ == "__main__":
    sys.exit(main())

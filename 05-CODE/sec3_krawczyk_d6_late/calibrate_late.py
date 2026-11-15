# -*- coding: utf-8 -*-
# [oid:kraw6L] calibrate_late.py -- calibrate the Krawczyk ladder for the 9
# late boxes BEFORE freezing (brief pt. 3 of BRIEF-H, same discipline as wave-2
# calibrate6.py).  Tries a grid of radii x precisions on ALL 9 targets (one per
# m level is enough, but the set is small so we do all 9) and records, for each
# config, pass / df_lb / residual / kappa.  Output: calibration_late.json.
#
# The FROZEN selection rule (used by krawczyk6_late_certificate.py) is:
#   smallest ladder radius certifying with kappa_ub < 0.5; fallback = smallest
#   passing; else first config (audit).
# This script only measures; the certificate generator applies the frozen rule.
import json
import sys
import time
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
PROGRESS = HERE / "progress.log"

REFINE = HERE / "refined_seeds_late.json"

# candidate ladder (radii, tightest-first), same shape as wave-2.  The late
# roots sit at |F'(a*)| ~ 1e13..1e16 (bigger |F'| than most of wave-2's 287),
# and the refined centers are ~1e-15..1e-13 away from the true root -- so
# 1e-24 is likely already too tight for the biggest-|F'| ones and the real
# working band is 1e-23..1e-19.
RADII_LADDER = ("1e-26", "1e-25", "1e-24", "1e-23", "1e-22", "1e-21", "1e-20")
PREC_LADDER = (256, 512, 1024)


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def attempt(spec, r_str, prec):
    """Krawczyk attempt (identical operator assembly to the generator)."""
    flint.ctx.prec = prec
    A = flint.arb
    C = flint.acb
    c0 = C(A(spec["re"]), A(spec["im"]))
    W = C(A("[+/- %s]" % r_str), A("[+/- %s]" % r_str))
    r_re = W.real.rad().upper()
    r_im = W.imag.rad().upper()
    r = r_re if bool(r_re > r_im) else r_im
    X = c0 + W
    for comp in (X.real, X.imag):
        if not bool(comp.rad() >= r):
            raise RuntimeError("degenerate box: component radius < r")
    n = spec["period"]
    m = n.bit_length() - 2
    half = 1 << m

    z = C(0)
    dz = C(0)
    for _ in range(n):
        dz = -(z ** 6) - 6 * c0 * z ** 5 * dz
        z = 1 - c0 * z ** 6
    F0, Df0 = z, dz

    zb = C(0)
    dzb = C(0)
    lowers = []
    j = 0
    for k in range(1, n + 1):
        dzb = -(zb ** 6) - 6 * X * zb ** 5 * dzb
        zb = 1 - X * zb ** 6
        if k <= half and k == (1 << j):
            lowers.append(zb.abs_lower())
            j += 1
    if len(lowers) != m + 1:
        raise RuntimeError("period marks miscount")
    D = dzb

    df_lb = D.abs_lower()
    df_nonzero = bool(df_lb > 0)
    N = c0 - F0 / D
    n_dist = (N - c0).abs_upper()
    newton_inside = bool(n_dist < r)
    Y = 1 / Df0
    E = 1 - Y * D
    K = c0 - Y * F0 + E * W
    k_dist = (K - c0).abs_upper()
    k_inside = bool(k_dist < r)
    kappa = k_dist / r
    period_ok = [bool(v > 0) for v in lowers]
    return {
        "box_radius": str(r.upper()),
        "prec_bits": prec,
        "f0_abs_ub": str(F0.abs_upper().upper()),
        "df_abs_lb": str(df_lb.lower()),
        "newton_dist_ub": str(n_dist.upper()),
        "krawczyk_dist_ub": str(k_dist.upper()),
        "krawczyk_kappa_ub": str(kappa.upper()),
        "newton_inside": newton_inside,
        "krawczyk_inside": k_inside,
        "df_nonzero": df_nonzero,
        "period_all_ok": all(period_ok),
        "period_min_lb": str(min(lowers).lower()),
        "_pass": bool(df_nonzero and newton_inside and k_inside and all(period_ok)),
    }


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    t0 = time.time()
    ref = json.loads(REFINE.read_text(encoding="utf-8"))
    seeds = ref["seeds"]
    print("=== calibrate_late.py [oid:kraw6L]: ladder calibration for 9 late "
          "boxes ===")
    out = {"run": "calibrate_late", "oid": "kraw6L",
           "radii_ladder": list(RADII_LADDER), "prec_ladder": list(PREC_LADDER),
           "per_target": {}}
    for lbl in sorted(seeds):
        s = seeds[lbl]
        spec = {"re": s["re"], "im": s["im"], "period": s["period"]}
        per = {"m": s["m"], "period": s["period"], "configs": {}}
        for prec in PREC_LADDER:
            for r_str in RADII_LADDER:
                rec = attempt(spec, r_str, prec)
                key = "r=%s prec=%d" % (r_str, prec)
                per["configs"][key] = rec
                print("%s %-18s pass=%s kappa=%s df_lb=%s" % (
                    lbl, key, rec["_pass"],
                    rec["krawczyk_kappa_ub"][:14], rec["df_abs_lb"][:14]))
        out["per_target"][lbl] = per
        out_path = HERE / "calibration_late.json"
        out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                            encoding="utf-8")
    out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    log_progress("STEP calibrate_late: calibration grid %d radii x %d precs x "
                 "%d targets measured -> calibration_late.json (%.1fs)"
                 % (len(RADII_LADDER), len(PREC_LADDER), len(seeds),
                    time.time() - t0))
    print("wrote %s (%.1fs)" % (out_path, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

# -*- coding: utf-8 -*-
# [oid:kraw6L] krawczyk6_late_certificate.py -- strict Krawczyk boxes for the
# 9 "late" superstable parameters of the subfamily f_a(x)=1-a*x^6
# (BRIEF-P-LATE-BOXES).
#
# The 9 targets were REJECTED by the blind run at the frozen level m=9; the
# subfamily investigation (REPORT_SUBFAMILY.md) showed they are genuine
# superstable parameters that first appear cleanly at m=10 (S031/S070/S235/
# S279/S319), m=11 (S104/S187/S339) or m=12 (S036).  Here we certify each at
# its TRUE level: period P = 2^(m+1) in {2048, 4096, 8192}, i.e. we prove
# existence + uniqueness of a zero of F(a) := x_P(a) in a box AND that its
# exact period is P (all lower marks 2^j, j<m+1, stay away from 0).
#
# Adaptation of the wave-2 generator (krawczyk6_certificate.py, oid kraw6H):
#   - seeds from refined_seeds_late.json (refine_late.py; model lesson: the
#     catalog-predicted a(u_cat, m) must be damped-Newton-refined at high
#     precision BEFORE Krawczyk -- plain fixed-step Newton overshoots for two
#     of the nine because |F'(a*)| ~ 1e13..1e16).
#   - per-target m and period (NOT a single global m).
#   - ladders FROZEN from calibrate_late.py -> calibration_late.json.
#   - workers = 3 (BRIEF: <=3 cores).
#   - period criterion: STRICT POSITIVITY of inf|x_{2^j}| (no fixed numeric
#     cutoffs; the cascade here compresses the marks only mildly, they are
#     still ~1e-2 .. 5e-2 at radius 1e-26).
#   - isolation: pairwise among the 9 new boxes AND against all 287 boxes of
#     the wave-2 certificate (read-only, exact decimal arithmetic).
#
# METHOD (per seed, all inclusions via Arb ball arithmetic -- python-flint):
#   x_0 = 0,  x_{k+1} = 1 - a*x_k^6;
#   F(a) := x_P(a),  F'(a) := dx_P(a),  dx_{k+1} = -x_k^6 - 6*a*x_k^5*dx_k.
#   X := c0 + W, W = [+/- r]^2 (square superset of the disk B(0,r)).
#   Interval Newton N(X) := c0 - F0/D;  Krawczyk K(X) := c0 - Y*F0 +
#       (1 - Y*D)*(X - c0), Y = 1/Df0.
#   PASS (strict ball comparisons):
#     (i) 0 not in D; (ii) |N-c0|_ub < r; (iii) |K-c0|_ub < r;
#     (iv) inf_{a in X} |x_{2^j}(a)| > 0 for j = 0..m  (x_1==1 excludes p=1).
#   THEOREM per box: existence + uniqueness of the zero of F in B(c0,r), and
#   exact period P (return times p*Z, p|P; (iv) kills all p=2^j).
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
REFINE = HERE / "refined_seeds_late.json"
TJL = HERE / "targets_late.json"
CALIB = HERE / "calibration_late.json"
W2CERT = Path(r"../sec3_krawczyk_d6/"
              r"krawczyk6_certificates.json")
OUT_JSON = HERE / "krawczyk6_late_certificates.json"
PROGRESS = HERE / "progress.log"

MAX_WORKERS = 3               # hard cap (BRIEF: <= 3 cores)
# FROZEN from calibrate_late.py (all 9 pass at the tightest radius 1e-26 with
# kappa < 0.5; the frame covers the refined centers, whose residuals are
# ~1e-64..1e-100 while |F'(X)| ~ 1e12..1e16, with huge headroom).
RADII_LADDER = ("1e-26", "1e-25", "1e-24", "1e-23", "1e-22")
PREC_LADDER = (256, 512, 1024)
KAPPA_MARGIN = "0.5"

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic for f_a(x)=1-a*x^6; "
    "9 late subfamily targets each at its TRUE level m (10/11/12), "
    "F(a)=x_P(a), P=2^(m+1) in {2048,4096,8192}, x_0=0; X=c0+[+/-r]^2 "
    "(square superset of B(c0,r)); N(X)=c0-F(c0)/F'(X); "
    "K(X)=c0-Y*F(c0)+(1-Y*F'(X))*(X-c0), Y=1/F'(c0); pass iff 0 notin F'(X), "
    "N(X) and K(X) strictly inside X, and inf_{a in X}|x_{2^j}(a)|>0 for "
    "j=0..m; implies existence and uniqueness of a superstable parameter of "
    "exact period P in B(c0,r)"
)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_specs():
    ref = json.loads(REFINE.read_text(encoding="utf-8"))
    tj = json.loads(TJL.read_text(encoding="utf-8"))
    specs = []
    for lbl in sorted(ref["seeds"]):
        s = ref["seeds"][lbl]
        t = tj["targets"][lbl]
        specs.append({
            "label": lbl,
            "source": "refined_seeds_late.json (subfamily targets.json, "
                      "catalog prediction refined by damped Newton)",
            "n": s["period"],
            "m": s["m"],
            "u_cat_re": t["u_cat_re"],
            "u_cat_im": t["u_cat_im"],
            "re": s["re"],
            "im": s["im"],
            "seed_policy": "refined",
            "rate_d_u": t["rate_d_u"],
            "rate_tag": t["rate_tag"],
        })
    return specs


def attempt(spec, r_str, prec, n_override=None):
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
    n = n_override if n_override is not None else spec["n"]
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
        "box_rad_re": str(r_re.upper()),
        "box_rad_im": str(r_im.upper()),
        "prec_bits": prec,
        "f0_abs_ub": str(F0.abs_upper().upper()),
        "df_abs_lb": str(df_lb.lower()),
        "newton_dist_ub": str(n_dist.upper()),
        "newton_inside": newton_inside,
        "krawczyk_dist_ub": str(k_dist.upper()),
        "krawczyk_kappa_ub": str(kappa.upper()),
        "krawczyk_inside": k_inside,
        "df_nonzero": df_nonzero,
        "period_checks": [
            {"j": jj, "d": 1 << jj, "abs_lb": str(lowers[jj].lower())}
            for jj in range(m + 1)
        ],
        "period_all_ok": all(period_ok),
        "period_min_lb": str(min(lowers).lower()),
        "_pass": bool(df_nonzero and newton_inside and k_inside and all(period_ok)),
        "_n_used": n, "_m_used": m,
    }


def certify_root(spec):
    chosen = None
    fallback = None
    for prec in PREC_LADDER:
        for r_str in RADII_LADDER:
            rec = attempt(spec, r_str, prec)
            if rec["_pass"]:
                if fallback is None:
                    fallback = rec
                if bool(flint.arb(rec["krawczyk_kappa_ub"]) <
                        flint.arb(KAPPA_MARGIN)):
                    chosen = rec
                    break
        if chosen is not None:
            break
    if chosen is None:
        chosen = fallback if fallback is not None else attempt(
            spec, RADII_LADDER[0], PREC_LADDER[0])
    out = {
        "label": spec["label"],
        "source": spec["source"],
        "n": spec["n"], "m": spec["m"],
        "u_cat_re": spec["u_cat_re"], "u_cat_im": spec["u_cat_im"],
        "box_center_re": spec["re"], "box_center_im": spec["im"],
        "seed_policy": spec["seed_policy"],
    }
    out.update(chosen)
    out["pass"] = out.pop("_pass")
    return out


def pairwise_isolation(records):
    """Rigorous disjointness of the certified SQUARES within the set, from
    exact centers + reported radii (component separation > r_i + r_j)."""
    flint.ctx.prec = 512
    A = flint.arb
    worst = None
    worst_pair = None
    all_disjoint = True
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            a, b = records[i], records[j]
            dre = A(a["box_center_re"]) - A(b["box_center_re"])
            dim = A(a["box_center_im"]) - A(b["box_center_im"])
            rsum = A(a["box_radius"]) + A(b["box_radius"])
            sep_re = abs(dre).lower()
            sep_im = abs(dim).lower()
            sep = sep_re if bool(sep_re > sep_im) else sep_im
            margin = sep / rsum.upper()
            disjoint = bool(sep > rsum.upper())
            all_disjoint = all_disjoint and disjoint
            if worst is None or bool(margin < worst):
                worst, worst_pair = margin, (a["label"], b["label"])
    return {
        "n_boxes": len(records),
        "n_pairs": len(records) * (len(records) - 1) // 2,
        "all_disjoint": bool(all_disjoint),
        "worst_separation_margin_lb": str(worst.lower()),
        "worst_pair": list(worst_pair),
    }


def isolation_vs_287(new_records):
    """Disjointness of each new box from each of the 287 wave-2 boxes
    (exact decimal arithmetic; the wave-2 file is read-only)."""
    flint.ctx.prec = 512
    A = flint.arb
    w2 = json.loads(W2CERT.read_text(encoding="utf-8"))
    old = w2["roots"]
    worst = None
    worst_pair = None
    all_disjoint = True
    n_pairs = 0
    for a in new_records:
        for b in old:
            n_pairs += 1
            dre = A(a["box_center_re"]) - A(b["box_center_re"])
            dim = A(a["box_center_im"]) - A(b["box_center_im"])
            rsum = A(a["box_radius"]) + A(b["box_radius"])
            sep_re = abs(dre).lower()
            sep_im = abs(dim).lower()
            sep = sep_re if bool(sep_re > sep_im) else sep_im
            margin = sep / rsum.upper()
            disjoint = bool(sep > rsum.upper())
            all_disjoint = all_disjoint and disjoint
            if worst is None or bool(margin < worst):
                worst, worst_pair = margin, (a["label"], b["label"])
    return {
        "n_new": len(new_records),
        "n_old": len(old),
        "n_pairs": n_pairs,
        "all_disjoint": bool(all_disjoint),
        "worst_separation_margin_lb": str(worst.lower()),
        "worst_pair": list(worst_pair),
    }


def correspondence_block(records):
    """NOT RIGOROUS (blind-run constants uncertified): scaled u-position of
    each certified center vs the subfamily catalog u, using u = A6*(a-a_F)
    *delta6^m at the box's OWN level m."""
    tj = json.loads(TJL.read_text(encoding="utf-8"))
    a_F = float(tj["a_F"])
    A6 = float(tj["A6"])
    d6 = float(tj["delta6"])
    d_us = []
    box_u_sizes = []
    for rec in records:
        scale = A6 * d6 ** rec["m"]
        u_re = (float(flint.arb(rec["box_center_re"]).mid()) - a_F) * scale
        u_im = float(flint.arb(rec["box_center_im"]).mid()) * scale
        du = ((u_re - float(rec["u_cat_re"])) ** 2
              + (u_im - float(rec["u_cat_im"])) ** 2) ** 0.5
        d_us.append(du)
        box_u_sizes.append(float(flint.arb(rec["box_radius"]).mid()) * scale)
    return {
        "_not_rigorous": True,
        "note": "u-positions use uncertified float constants a_F/A6/delta6 of "
                "the blind run; the rigorous statement lives in a-space",
        "d_u_min": min(d_us), "d_u_max": max(d_us),
        "box_u_size_max": max(box_u_sizes),
    }


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== krawczyk6_late_certificate.py [oid:kraw6L]: strict boxes for "
          "the 9 late subfamily parameters (f_a(x)=1-a*x^6, P=2048/4096/8192) "
          "===")
    specs = load_specs()
    print("loaded %d refined seeds (m=%s, P=%s)" % (
        len(specs),
        sorted(set(s["m"] for s in specs)),
        sorted(set(s["n"] for s in specs))))
    log_progress("STEP generator start: %d seeds, radii=%s prec=%s pool=%d"
                 % (len(specs), RADII_LADDER, PREC_LADDER, MAX_WORKERS))

    records = []

    def flush_cert():
        payload = {
            "run": "krawczyk6_late_certificate", "oid": "kraw6L",
            "method": METHOD_STRING,
            "sha256_pins": {"refined_seeds_late.json": sha256(REFINE),
                            "targets_late.json": sha256(TJL),
                            "calibration_late.json": sha256(CALIB),
                            "krawczyk6_certificates.json": sha256(W2CERT)},
            "family": "1 - a*x^6",
            "frozen": {"radii_ladder": list(RADII_LADDER),
                       "prec_ladder": list(PREC_LADDER),
                       "max_workers": MAX_WORKERS,
                       "selection_rule": "smallest ladder radius passing "
                                         "with kappa_ub < 0.5; fallback: "
                                         "smallest passing; else first config "
                                         "(audit)",
                       "kappa_margin": KAPPA_MARGIN,
                       "period_criterion": "STRICT POSITIVITY of "
                                           "inf|x_{2^j}|, j=0..m (no fixed "
                                           "numeric thresholds; resolved via "
                                           "Arb ball arithmetic at same prec)"},
            "n_roots": len(records),
            "roots": records,
            "is_complete": False,
        }
        OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                            encoding="utf-8")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i, rec in enumerate(ex.map(certify_root, specs), start=1):
            records.append(rec)
            flush_cert()
            dt = time.time() - t0
            per = dt / i
            eta = per * (len(specs) - i)
            line = ("PARAM %3d/%d id=%s P=%d r<=%s kappa<=%s pass=%s "
                    "t=%.1fs ETA=%.0fs"
                    % (i, len(specs), rec["label"], rec["n"],
                       rec["box_radius"][:14], rec["krawczyk_kappa_ub"][:12],
                       rec["pass"], dt, eta))
            print("  " + line, flush=True)
            log_progress(line)

    log_progress("STEP pairwise isolation among %d new boxes"
                 % len(records))
    iso = pairwise_isolation(records)
    log_progress("STEP isolation vs 287 wave-2 boxes")
    iso287 = isolation_vs_287(records)
    n_pass = sum(1 for r in records if r["pass"])
    all_pass = bool(n_pass == len(records) and iso["all_disjoint"]
                    and iso287["all_disjoint"])

    for r in records:
        r.pop("_n_used", None)
        r.pop("_m_used", None)

    payload = {
        "run": "krawczyk6_late_certificate", "oid": "kraw6L",
        "method": METHOD_STRING,
        "sha256_pins": {"refined_seeds_late.json": sha256(REFINE),
                        "targets_late.json": sha256(TJL),
                        "calibration_late.json": sha256(CALIB),
                        "krawczyk6_certificates.json": sha256(W2CERT)},
        "family": "1 - a*x^6",
        "frozen": {"radii_ladder": list(RADII_LADDER),
                   "prec_ladder": list(PREC_LADDER),
                   "max_workers": MAX_WORKERS,
                   "selection_rule": "smallest ladder radius passing "
                                     "with kappa_ub < 0.5; fallback: "
                                     "smallest passing; else first config "
                                     "(audit)",
                   "kappa_margin": KAPPA_MARGIN,
                   "period_criterion": "STRICT POSITIVITY of "
                                       "inf|x_{2^j}|, j=0..m (no fixed "
                                       "numeric thresholds; resolved via "
                                       "Arb ball arithmetic at same prec)"},
        "n_roots": len(records),
        "roots": records,
        "pairwise_isolation": iso,
        "isolation_vs_wave2": iso287,
        "degenerate_boxes": [],
        "correspondence": correspondence_block(records),
        "all_pass": all_pass,
        "is_complete": True,
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    log_progress("STEP generator done: %d/%d pass, iso@9 all_disjoint=%s, "
                 "iso@287 all_disjoint=%s (worst margin %s, pair %s), "
                 "all_pass=%s, total %.1fs"
                 % (n_pass, len(records), iso["all_disjoint"],
                    iso287["all_disjoint"],
                    iso287["worst_separation_margin_lb"][:10],
                    iso287["worst_pair"], all_pass, time.time() - t0))
    print("=== SUMMARY: %d/%d pass; isolation@9=%s isolation@287=%s; "
          "all_pass=%s (%.1fs) -> %s" % (
              n_pass, len(records), iso["all_disjoint"],
              iso287["all_disjoint"], all_pass, time.time() - t0, OUT_JSON))
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
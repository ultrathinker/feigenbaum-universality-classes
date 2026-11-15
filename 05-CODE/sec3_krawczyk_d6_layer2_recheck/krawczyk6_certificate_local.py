# -*- coding: utf-8 -*-
# [oid:kraw6H] krawczyk6_certificate.py -- strict Krawczyk boxes for the 287
# interior accepted superstable parameters of the sextic family
# f_a(x) = 1 - a*x^6 (BRIEF-H-KRAWCZYK6; seeds: agentB_deg6/match_results.json
# blind run, verdict set |u| <= EDGE_U = 44, m = 9, period P = 2^10 = 1024,
# F(a) := x_P(a), orbit of the critical point x = 0).
#
# Direct adaptation of the gate-2-audited wave-1 generator
# (wave1_krawczyk4/krawczyk4_certificate.py, oid kraw4E); changes per BRIEF-H:
# degree 6 (x^6; derivative -x^6 - 6a*x^5*dx), deg6 seed format (a_ss_re is an
# exact decimal string, a_ss_im is a float64 -> exact round-trip repr string),
# 287 seeds, ladders CALIBRATED on the first 10 boxes (calibrate6.py,
# calibration6.json) before freezing -- brief pt. 3.
#
# METHOD (per seed, all inclusions via Arb ball arithmetic -- python-flint,
# no float comparisons anywhere in the decision path):
#   x_0 = 0,  x_{k+1} = 1 - a*x_k^6;
#   F(a)  := x_n(a)  with n = 2^(m+1) = 1024,
#   F'(a) := dx_n(a), dx_{k+1} = -x_k^6 - 6*a*x_k^5*dx_k  (updated BEFORE x).
#   X := c0 + W, W = [+/- r]^2 (SQUARE, superset of the disk B(0,r); AUDIT_FIX
#   from the mside package AUDIT_FINDING.md: never a zero-width component);
#   r := outward upper bound of the MAXIMUM radius of BOTH components of W,
#   with a regression guard (both component radii of X >= r) BEFORE any orbit.
#   Point orbit: F0 = F(c0), Df0 = F'(c0).  Box orbit over X: D = F'(X) and
#   lower bounds of |x_{2^j}(X)|, j = 0..m (j=0 is x_1 == 1 identically).
#   Interval Newton (Moore):   N(X) := c0 - F0/D.
#   Krawczyk (Y := 1/Df0):     K(X) := c0 - Y*F0 + (1 - Y*D)*(X - c0).
#   PASS rule (strict ball comparisons):
#     (i)   0 not in D                (D.abs_lower() > 0);
#     (ii)  N(X) strictly inside X    ((N-c0).abs_upper() < r);
#     (iii) K(X) strictly inside X    ((K-c0).abs_upper() < r);
#     (iv)  inf_{a in X} |x_{2^j}(a)| > 0 for j = 0..m (marks at 2^j < n).
#   THEOREM: existence + uniqueness of a zero of F in the DISK B(c0, r) (see
#   REPORT-KRAWCZYK6.md for the statement; uniqueness holds on the square X),
#   of EXACT period n: return times of the critical orbit form p*Z, p | n;
#   (iv) excludes p in {2^j : j <= m} (x_1 == 1 excludes p = 1), hence p = n.
#
# PAIRWISE ISOLATION (after all boxes): squares X_i, X_j disjoint iff a
# component separation exceeds r_i + r_j (exact decimal arithmetic).
#
# PARALLELISM: ProcessPoolExecutor, max_workers=4 EXACTLY (BRIEF: CPU budget
# <= 4 processes/threads).  Incremental output: krawczyk6_certificates.json
# rewritten after EVERY root (is_complete=false until the end); progress.log
# line after every root, flushed.
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
PILOT = HERE / "match_results_dummy.json"      # not used by attempt()/certify_root()
CALIB = HERE / "calibration_dummy.json"
OUT_JSON = HERE / "krawczyk6_certificates_local.json"
PROGRESS = HERE / "progress_local.log"

MAX_WORKERS = 4               # hard cap (BRIEF: pool = 4)
# FROZEN after calibration on the first 10 boxes (calibrate6.py; brief pt. 3):
# tight-first; float64-im seeds have residual <= half-ulp(6.3e-8) = 6.6e-24,
# covered by 1e-23/1e-22; 1e-21 is the catch-all fallback, 1e-25 serves
# super-precise (real) seeds that pass even tighter than 1e-24.
RADII_LADDER = ("1e-24", "1e-23", "1e-22", "1e-21", "1e-25")
PREC_LADDER = (256, 512, 1024)
# Layer 2 (refined seeds, S078/S090 only): |F'(a*)| ~ 8e14/1.1e18 makes the
# F'(X)-interval straddle 0 for r >= 1e-25/1e-29 while the float64-im stored
# seed sits 1.7e-24/1.6e-24 from the root -- no source-string box can pass.
# Centers = deterministic Newton refinements (refine6.py -> refined_seeds6.json,
# sha-pinned), certified by the SAME pipeline; identity with the blind-run root
# pinned in Arb: |c_ref - stored_seed|_ub < 5e-24 (nearest foreign box: 2.6e-10).
LAYER2_RADII_LADDER = ("1e-30", "1e-29", "1e-28", "1e-27", "1e-31")
SEED_DIST_LIMIT = "5e-24"
REFINE = HERE / "refined_seeds6.json" if (HERE / "refined_seeds6.json").exists() else HERE / "refined_seeds_dummy.json"

# frozen inventory of the blind run (asserted in load_seeds):
EXPECTED = {"catalog_total": 342, "interior_total": 294, "interior_accepted": 287,
            "edge_total": 48, "edge_accepted": 46,
            "rejected_interior": ["S031", "S036", "S070", "S104", "S187",
                                  "S235", "S279"],
            "rejected_edge": ["S319", "S339"]}

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic for f_a(x)=1-a*x^6: "
    "X=c0+[+/-r]^2 (square superset of B(c0,r)); F(a)=x_n(a), n=2^(m+1)=1024, "
    "x_0=0; N(X)=c0-F(c0)/F'(X); K(X)=c0-Y*F(c0)+(1-Y*F'(X))*(X-c0), "
    "Y=1/F'(c0); pass iff 0 notin F'(X), N(X) and K(X) strictly inside X, "
    "and inf_{a in X}|x_{2^j}(a)|>0 for j=0..m; implies existence and "
    "uniqueness of a superstable parameter of exact period n in B(c0,r)"
)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def ub_str(x) -> str:
    return str(flint.arb(x).upper())


def lb_str(x) -> str:
    return str(flint.arb(x).lower())


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_seeds():
    """Seeds = the 287 interior accepted rows (verdict set |u| <= 44) of the
    deg6 blind run.  Seed = exact decimal a_ss_re string + exact round-trip
    repr of the a_ss_im float64 (the blind run stores the imaginary part as a
    float).  The blind run's row inventory is asserted against EXPECTED so the
    catalog cannot drift silently."""
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    fr = d["meta"]["frozen"]
    m_lvl = int(fr["M_MATCH"])
    n = 1 << (m_lvl + 1)
    if int(fr["PERIOD"]) != n:
        raise RuntimeError("pilot frozen PERIOD != 2^(M_MATCH+1)")
    if fr["family"] != "1 - a*x^6":
        raise RuntimeError("unexpected family %r" % fr["family"])
    seeds = []
    rej_int, rej_edge = [], []
    n_edge_acc = 0
    refined = json.loads(REFINE.read_text(encoding="utf-8"))["seeds"]
    for rec in d["catalog_results"]:
        if rec.get("accepted"):
            if rec.get("edge"):
                n_edge_acc += 1
            else:
                spec = {
                    "label": rec["id"],
                    "source": "match_results:catalog_results[accepted,interior]",
                    "n": n, "m": m_lvl,
                    "re": str(rec["a_ss_re"]), "im": repr(float(rec["a_ss_im"])),
                    "u_cat_re": rec["u_cat_re"], "u_cat_im": rec["u_cat_im"],
                    "seed_policy": "source_string", "ladder": RADII_LADDER}
                if rec["id"] in refined:
                    r6 = refined[rec["id"]]
                    spec.update({
                        "re": r6["re"], "im": r6["im"],
                        "seed_policy": "refined", "ladder": LAYER2_RADII_LADDER,
                        "stored_seed_re": str(rec["a_ss_re"]),
                        "stored_seed_im": repr(float(rec["a_ss_im"]))})
                seeds.append(spec)
        else:
            (rej_edge if rec.get("edge") else rej_int).append(rec["id"])
    got = {"catalog_total": len(d["catalog_results"]),
           "interior_accepted": len(seeds), "edge_accepted": n_edge_acc,
           "rejected_interior": rej_int, "rejected_edge": rej_edge}
    if (got["catalog_total"] != EXPECTED["catalog_total"]
            or got["interior_accepted"] != EXPECTED["interior_accepted"]
            or got["edge_accepted"] != EXPECTED["edge_accepted"]
            or got["rejected_interior"] != EXPECTED["rejected_interior"]
            or got["rejected_edge"] != EXPECTED["rejected_edge"]):
        raise RuntimeError("blind-run inventory drift: %r" % got)
    return seeds


def attempt(spec, r_str, prec, family_eps=0, n_override=None):
    """One certification attempt at (radius, precision). family_eps != 0 or
    n_override are ONLY used by the tamper battery (Engine A corruptions);
    the clean run always uses family_eps=0 and n from the spec."""
    flint.ctx.prec = prec
    A = flint.arb
    C = flint.acb

    c0 = C(A(spec["re"]), A(spec["im"]))
    W = C(A("[+/- %s]" % r_str), A("[+/- %s]" % r_str))   # AUDIT_FIX: square
    r_re = W.real.rad().upper()
    r_im = W.imag.rad().upper()
    r = r_re if bool(r_re > r_im) else r_im
    X = c0 + W
    for comp in (X.real, X.imag):                          # regression guard
        if not bool(comp.rad() >= r):
            raise RuntimeError("degenerate box: component radius < r "
                               "(segment regression, AUDIT_FINDING)")
    if family_eps:
        eps = A(1) + A(family_eps)
        c0 = C(c0.real * eps, c0.imag * eps)
        X = C(X.real * eps, X.imag * eps)
    n = n_override if n_override is not None else spec["n"]
    m = n.bit_length() - 2          # n = 2^(m+1)  ->  m
    half = 1 << m

    # ---- point orbit: F0 = F(c0), Df0 = F'(c0)
    z = C(0)
    dz = C(0)
    for _ in range(n):
        dz = -(z ** 6) - 6 * c0 * z ** 5 * dz
        z = 1 - c0 * z ** 6
    F0, Df0 = z, dz

    # ---- box orbit: D = F'(X); marks |x_{2^j}(X)| at j = 0..m
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
    D = dzb
    if len(lowers) != m + 1:
        raise RuntimeError("period marks miscount")

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
        "box_radius": ub_str(r),
        "box_rad_re": ub_str(r_re),
        "box_rad_im": ub_str(r_im),
        "prec_bits": prec,
        "f0_abs_ub": ub_str(F0.abs_upper()),
        "df_abs_lb": lb_str(df_lb),
        "newton_dist_ub": ub_str(n_dist),
        "newton_inside": newton_inside,
        "krawczyk_dist_ub": ub_str(k_dist),
        "krawczyk_kappa_ub": ub_str(kappa),
        "krawczyk_inside": k_inside,
        "df_nonzero": df_nonzero,
        "period_checks": [
            {"j": jj, "d": 1 << jj, "abs_lb": lb_str(lowers[jj])}
            for jj in range(m + 1)
        ],
        "period_all_ok": all(period_ok),
        "_pass": bool(df_nonzero and newton_inside and k_inside and all(period_ok)),
        "_n_used": n, "_m_used": m,
    }


KAPPA_MARGIN = "0.5"        # selection rule: smallest radius certifying with
                            # kappa_ub < 0.5 (2x contraction headroom; the
                            # float64-im seeds have offsets up to ~6.6e-24, so
                            # tight-first-without-margin lands some roots at
                            # kappa ~ 0.9998 -- valid but fragile to reproduce)


def certify_root(spec):
    """Ladder search with a contraction-margin rule (FROZEN):
    (1) the SMALLEST ladder radius whose attempt passes with kappa_ub < 0.5;
    (2) fallback: the smallest radius that passes at all;
    (3) on total failure freeze the first configuration so the failure is
    fully reported and auditable.
    Layer-2 (refined) seeds use their own frozen ladder and additionally get a
    rigorous Arb pin |c_ref - stored_seed|_ub < SEED_DIST_LIMIT."""
    ladder = spec.get("ladder", RADII_LADDER)
    A = flint.arb
    chosen = None
    fallback = None
    for prec in PREC_LADDER:
        for r_str in ladder:
            rec = attempt(spec, r_str, prec)
            if rec["_pass"]:
                if fallback is None:
                    fallback = rec
                if bool(flint.arb(rec["krawczyk_kappa_ub"]) < A(KAPPA_MARGIN)):
                    chosen = rec
                    break
        if chosen is not None:
            break
    if chosen is None:
        chosen = fallback if fallback is not None else attempt(
            spec, ladder[0], PREC_LADDER[0])
    out = {
        "label": spec["label"], "source": spec["source"],
        "n": spec["n"], "m": spec["m"],
        "u_cat_re": spec["u_cat_re"], "u_cat_im": spec["u_cat_im"],
        "box_center_re": spec["re"], "box_center_im": spec["im"],
        "seed_policy": spec["seed_policy"],
    }
    if spec["seed_policy"] == "refined":
        flint.ctx.prec = 512
        A, C = flint.arb, flint.acb
        c_ref = C(A(spec["re"]), A(spec["im"]))
        c_sto = C(A(spec["stored_seed_re"]), A(spec["stored_seed_im"]))
        dist = (c_ref - c_sto).abs_upper()
        if not bool(dist < A(SEED_DIST_LIMIT)):
            raise RuntimeError("%s: refined center drifted from the stored "
                               "seed: %s >= %s" % (spec["label"],
                                                   ub_str(dist), SEED_DIST_LIMIT))
        out.update({
            "stored_seed_re": spec["stored_seed_re"],
            "stored_seed_im": spec["stored_seed_im"],
            "seed_dist_ub": ub_str(dist),
            "refine_iters": 8,
        })
    out.update(chosen)
    out["pass"] = out.pop("_pass")
    return out


def pairwise_isolation(records):
    """Rigorous disjointness of the certified SQUARES from exact decimal
    centers + reported radii (component separation > r_i + r_j)."""
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
            sep_re = abs(dre).lower()          # exact: decimals, no rounding
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
        "worst_separation_margin_lb": lb_str(worst),
        "worst_pair": list(worst_pair),
    }


def correspondence_block(records):
    """NOT RIGOROUS (float constants a_F, A6, delta6 from the blind run are
    uncertified): scaled u-positions of the certified boxes vs catalog zeros.
    u(a) = A6*(a - a_F)*delta6^m with m = 9; box u-radius = r*A6*delta6^9."""
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    a_F = float(d["a_F"])
    A6 = float(d["A6"])
    d6 = float(d["delta6"])
    scale = A6 * d6 ** 9
    d_us, box_u_sizes = [], []
    for rec in records:
        u_re = (float(flint.arb(rec["box_center_re"]).mid()) - a_F) * scale
        u_im = float(flint.arb(rec["box_center_im"]).mid()) * scale
        du = ((u_re - rec["u_cat_re"]) ** 2 + (u_im - rec["u_cat_im"]) ** 2) ** 0.5
        d_us.append(du)
        box_u_sizes.append(float(flint.arb(rec["box_radius"]).mid()) * scale)
    return {
        "_not_rigorous": True,
        "note": "u-positions use uncertified float constants a_F/A6/delta6 "
                "from the blind run; the rigorous statement lives in a-space",
        "d_u_min": min(d_us), "d_u_max": max(d_us),
        "box_u_size_max": max(box_u_sizes),
    }


def main():
    t0 = time.time()
    try:                                   # Cyrillic paths on cp1251/cp1252 consoles
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== krawczyk6_certificate.py [oid:kraw6H]: Krawczyk boxes for 287 "
          "sextic superstable parameters (f_a(x)=1-a*x^6, n=1024) ===")
    seeds = load_seeds()
    print("loaded %d interior accepted seeds from match_results.json "
          "(rejected interior: %s; edge rows excluded from the verdict set: "
          "46 accepted + %s rejected)" % (len(seeds),
                                          EXPECTED["rejected_interior"],
                                          EXPECTED["rejected_edge"]))
    log_progress("STEP generator start: %d seeds, ladders radii=%s prec=%s, "
                 "pool=%d" % (len(seeds), RADII_LADDER, PREC_LADDER, MAX_WORKERS))

    records = []

    def flush_cert():
        payload = {
            "run": "krawczyk6_certificate", "oid": "kraw6H",
            "method": METHOD_STRING,
            "sha256_pins": {"match_results.json": sha256(PILOT),
                            "calibration6.json": sha256(CALIB),
                            "refined_seeds6.json": sha256(REFINE)},
            "family": "1 - a*x^6",
            "frozen": {"radii_ladder": list(RADII_LADDER),
                       "layer2_radii_ladder": list(LAYER2_RADII_LADDER),
                       "prec_ladder": list(PREC_LADDER),
                       "max_workers": MAX_WORKERS,
                       "verdict_set": "interior |u|<=EDGE_U=44 (287 rows); "
                                      "46 accepted edge rows are NOT certified",
                       "ladder_origin": "calibrated on the first 10 boxes "
                                        "(calibrate6.py -> calibration6.json, "
                                        "pinned), then frozen",
                       "selection_rule": "smallest ladder radius passing "
                                        "with kappa_ub < 0.5; fallback: "
                                        "smallest passing; else first config "
                                        "(audit)",
                       "layer2_policy": "S078/S090: |F'(a*)|~8e14/1.1e18 "
                                        "blows the F'(X)-interval across 0 "
                                        "before the float64-im seed offset "
                                        "(1.7e-24/1.6e-24) is inside the box; "
                                        "centers = deterministic refinements "
                                        "(refine6.py), pinned to the stored "
                                        "seeds in Arb: seed_dist_ub < "
                                        + SEED_DIST_LIMIT,
                       "seed_dist_limit": SEED_DIST_LIMIT},
            "n_roots": len(records),
            "roots": records,
            "is_complete": False,
        }
        OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                            encoding="utf-8")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i, rec in enumerate(ex.map(certify_root, seeds), start=1):
            records.append(rec)
            flush_cert()                       # incremental, after EVERY box
            dt = time.time() - t0
            per = dt / i
            eta = per * (len(seeds) - i)
            line = ("PARAM %3d/%d id=%s P=%d r<=%s kappa<=%s pass=%s "
                    "t=%.1fs ETA=%.0fs"
                    % (i, len(seeds), rec["label"], rec["_n_used"],
                       rec["box_radius"][:12], rec["krawczyk_kappa_ub"][:10],
                       rec["pass"], dt, eta))
            print("  " + line, flush=True)
            log_progress(line)

    # ---- pairwise isolation + degenerate-box audit + final flush
    log_progress("STEP pairwise isolation check (%d pairs)"
                 % (len(records) * (len(records) - 1) // 2))
    iso = pairwise_isolation(records)
    # degenerate-box audit (redundant with the in-attempt regression guard):
    # all three stored radius strings come from the same parse, so they must
    # be identical, and neither component radius may be CERTAINLY below r
    # (flint >= on balls means "certainly", so identical balls compare False).
    flint.ctx.prec = 512
    degen = []
    for r in records:
        same = (r["box_rad_re"] == r["box_radius"] == r["box_rad_im"])
        re_less = bool(flint.arb(r["box_rad_re"]) < flint.arb(r["box_radius"]))
        im_less = bool(flint.arb(r["box_rad_im"]) < flint.arb(r["box_radius"]))
        if not same or re_less or im_less:
            degen.append(r["label"])
    if degen:
        raise RuntimeError("degenerate boxes survived: %s" % degen)
    n_pass = sum(1 for r in records if r["pass"])
    all_pass = bool(n_pass == len(records) and iso["all_disjoint"])

    for r in records:
        r.pop("_n_used", None)
        r.pop("_m_used", None)

    payload = {
        "run": "krawczyk6_certificate", "oid": "kraw6H",
        "method": METHOD_STRING,
        "sha256_pins": {"match_results.json": sha256(PILOT),
                        "calibration6.json": sha256(CALIB),
                        "refined_seeds6.json": sha256(REFINE)},
        "family": "1 - a*x^6",
        "frozen": {"radii_ladder": list(RADII_LADDER),
                   "layer2_radii_ladder": list(LAYER2_RADII_LADDER),
                   "prec_ladder": list(PREC_LADDER),
                   "max_workers": MAX_WORKERS,
                   "verdict_set": "interior |u|<=EDGE_U=44 (287 rows); "
                                  "46 accepted edge rows are NOT certified",
                   "ladder_origin": "calibrated on the first 10 boxes "
                                    "(calibrate6.py -> calibration6.json, "
                                    "pinned), then frozen",
                   "selection_rule": "smallest ladder radius passing "
                                    "with kappa_ub < 0.5; fallback: "
                                    "smallest passing; else first config "
                                    "(audit)",
                   "layer2_policy": "S078/S090: |F'(a*)|~8e14/1.1e18 "
                                    "blows the F'(X)-interval across 0 "
                                    "before the float64-im seed offset "
                                    "(1.7e-24/1.6e-24) is inside the box; "
                                    "centers = deterministic refinements "
                                    "(refine6.py), pinned to the stored "
                                    "seeds in Arb: seed_dist_ub < "
                                    + SEED_DIST_LIMIT,
                   "seed_dist_limit": SEED_DIST_LIMIT},
        "n_roots": len(records),
        "roots": records,
        "pairwise_isolation": iso,
        "degenerate_boxes": [],
        "correspondence": correspondence_block(records),
        "all_pass": all_pass,
        "is_complete": True,
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    log_progress("STEP generator done: %d/%d pass, isolation all_disjoint=%s "
                 "(worst margin %s, pair %s), all_pass=%s, total %.1fs"
                 % (n_pass, len(records), iso["all_disjoint"],
                    iso["worst_separation_margin_lb"][:10], iso["worst_pair"],
                    all_pass, time.time() - t0))
    print("=== SUMMARY: %d/%d pass; isolation %s; all_pass=%s (%.1fs) -> %s"
          % (n_pass, len(records), iso["all_disjoint"], all_pass,
             time.time() - t0, OUT_JSON))
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())

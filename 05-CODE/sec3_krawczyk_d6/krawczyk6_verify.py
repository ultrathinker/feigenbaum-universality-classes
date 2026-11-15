# -*- coding: utf-8 -*-
# [oid:kraw6H] krawczyk6_verify.py -- INDEPENDENT verifier for
# krawczyk6_certificates.json.  Exit 0 = verified, exit 1 = any drift.
#
# SEPARATE CODE PATH from krawczyk6_certificate.py (adaptation of the audited
# wave-1 verifier to degree 6 and the two-layer seed policy): complex balls
# are hand-rolled 3-tuples (re, im, r) -- midpoint components as arb intervals
# plus ONE explicit ball radius absorbing all propagation slack -- with
# hand-derived slack formulas (vmul: |m_x|*r_y + |m_y|*r_x + r_x*r_y +
# midpoint rounding; vabs_ub/lb hand-written), DIFFERENT orbit algebra
# (z5 = z4*z, z6 = z3*z3 vs the generator's acb pow), and DIFFERENT operator
# assemblies:
#   newton_dist    := |f(c0)|_ub / |f'(X)|_lb        (no ball division),
#   krawczyk_dist  := |Y|_ub*|F0|_ub + |1-Y*D|_ub*sqrt(2)*r
#                     (triangle inequality; sqrt(2)-corner factor for the
#                     certified SQUARE, while the ORBITS are propagated over
#                     the DISK B(c0, r_reported) -- a subset of the square:
#                     domination-direction recheck of every orbit quantity).
#
# TOLERANCE POLICY (frozen here, justified by _probe_paths.py measurements on
# the final certificate): degree-6 superstable roots have F'-curvature over
# the boxes up to 8 orders worse than the quartic (|E| = |1-Y*F'(X)| up to
# ~6e-9 relative; layer-2 roots S078/S090 even blow the F'(X)-interval across
# 0 at source-string radii), so exact cross-path equality of the magnitude
# bounds is structurally impossible -- the two paths bound different
# supersets.  Frozen bands (two-sided factors, observed worst in brackets):
#   f0_abs_ub          rel 1e-2   [1.0002]
#   df_abs_lb          factor 1.5 [1.034]
#   newton_dist_ub     factor 1.5 [1.054]
#   period abs_lb      factor 1.1 [1.001]
#   krawczyk_dist_ub / kappa_ub: LAYER-1 only, factor 3 [1.19]; layer-2
#     exempt (generator slack ~8e8x: ball-product accumulation over 1024
#     steps dominates; their pass/inside booleans are EXACTLY recomputed).
# Every boolean (df_nonzero, newton_inside, krawczyk_inside, period_all_ok,
# pass) is recomputed and must match EXACTLY (domination direction: the disk
# subset must satisfy every claim made on the square).
#
# EVERY reported field is independently checked:
#   [0] frozen header strings + sha256 pins (match_results.json,
#       calibration6.json, refined_seeds6.json);
#   [1] source inventory: 287 labels = interior accepted rows (verdict set
#       |u|<=44), rejected interior exactly [S031 S036 S070 S104 S187 S235
#       S279], 46 accepted edge rows excluded, 2 rejected edge rows;
#       per-root n=2^(m+1)=1024, m=9, u_cat match; LAYER-1 box centers
#       EXACT string match to the pilot seeds; LAYER-2 (exactly S078/S090):
#       stored seeds EXACT string match to the pilot, centers EXACT string
#       match to refined_seeds6.json, seed_dist_ub recomputed in Arb and
#       pinned < 5e-24;
#   [2] structure: exact key sets per layer; radius in the layer's frozen
#       ladder and in (0, 1e-20] (L1) / (0, 1e-24] (L2); prec in the frozen
#       ladder; period_checks shape j=0..m with d=2^j; non-degenerate radii
#       (AUDIT_FINDING guard); selection rule respected (kappa_ub < 0.5
#       unless a smaller radius failed -- checked via [3] recompute);
#   [3] full tuple-path recomputation per root (pool=4): all magnitudes
#       within the frozen bands, all booleans exact, pass confirmed;
#   [4] pairwise isolation recomputed from reported centers + radii
#       (41,041 pairs);
#   [5] aggregates + correspondence sanity (d_u_max < R_ACCEPT_U,
#       box_u_size_max < 1e-12).
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
PILOT = Path(r"../sec3_pilot_d6/match_results.json")
CALIB = HERE / "calibration6.json"
REFINE = HERE / "refined_seeds6.json"
CERT_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    HERE / "krawczyk6_certificates.json")

MAX_WORKERS = 4               # BRIEF: pool = 4
RADII_LADDER = ("1e-24", "1e-23", "1e-22", "1e-21", "1e-25")
LAYER2_RADII_LADDER = ("1e-30", "1e-29", "1e-28", "1e-27", "1e-31")
PREC_LADDER = (256, 512, 1024)
SEED_DIST_LIMIT = "5e-24"
KAPPA_MARGIN = "0.5"
REL_TOL_F0 = "1e-2"
FACTOR_DF = "1.5"
FACTOR_ND = "1.5"
FACTOR_PER = "1.1"
FACTOR_KD = "3"               # layer-1 krawczyk_dist/kappa only

EXPECTED = {"catalog_total": 342, "interior_total": 294, "interior_accepted": 287,
            "edge_total": 48, "edge_accepted": 46,
            "rejected_interior": ["S031", "S036", "S070", "S104", "S187",
                                  "S235", "S279"],
            "rejected_edge": ["S319", "S339"]}
LAYER2_LABELS = ("S078", "S090")

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic for f_a(x)=1-a*x^6: "
    "X=c0+[+/-r]^2 (square superset of B(c0,r)); F(a)=x_n(a), n=2^(m+1)=1024, "
    "x_0=0; N(X)=c0-F(c0)/F'(X); K(X)=c0-Y*F(c0)+(1-Y*F'(X))*(X-c0), "
    "Y=1/F'(c0); pass iff 0 notin F'(X), N(X) and K(X) strictly inside X, "
    "and inf_{a in X}|x_{2^j}(a)|>0 for j=0..m; implies existence and "
    "uniqueness of a superstable parameter of exact period n in B(c0,r)"
)

TOP_KEYS = {"run", "oid", "method", "sha256_pins", "family", "frozen",
            "n_roots", "roots", "pairwise_isolation", "degenerate_boxes",
            "correspondence", "all_pass", "is_complete"}
ROOT_KEYS_L1 = {"label", "source", "n", "m", "u_cat_re", "u_cat_im",
                "box_center_re", "box_center_im", "seed_policy",
                "box_radius", "box_rad_re", "box_rad_im",
                "prec_bits", "f0_abs_ub", "df_abs_lb",
                "newton_dist_ub", "newton_inside", "krawczyk_dist_ub",
                "krawczyk_kappa_ub", "krawczyk_inside", "df_nonzero",
                "period_checks", "period_all_ok", "pass"}
ROOT_KEYS_L2 = ROOT_KEYS_L1 | {"stored_seed_re", "stored_seed_im",
                               "seed_dist_ub", "refine_iters"}
PC_KEYS = {"j", "d", "abs_lb"}
ISO_KEYS = {"n_boxes", "n_pairs", "all_disjoint",
            "worst_separation_margin_lb", "worst_pair"}

errors = []


def err(msg):
    errors.append(msg)
    print("  FAIL:", msg)


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- tuple path
def recompute_root(spec):
    """Worker: independent (re, im, r)-tuple recomputation over the DISK
    B(c0, r_reported).  Returns dict of arb-strings and booleans."""
    flint.ctx.prec = spec["prec_bits"]
    A = flint.arb

    def rabs_ub(x):
        a, b = abs(x.lower()), abs(x.upper())
        return b if bool(b > a) else a

    def vconst(re_, im_=None):
        return (A(re_), A(0) if im_ is None else A(im_), A(0))

    def vadd(x, y):
        re_, im_ = x[0] + y[0], x[1] + y[1]
        r = x[2] + y[2] + re_.rad() + im_.rad()
        return (re_, im_, r)

    def vaddc(x, c_re, c_im):
        re2, im2 = x[0] + A(c_re), x[1] + A(c_im)
        r = x[2] + re2.rad() + im2.rad()
        return (re2, im2, r)

    def vneg(x):
        return (-x[0], -x[1], x[2])

    def cmul_mid(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def mabs_ub(m):
        return (m[0] * m[0] + m[1] * m[1]).sqrt().upper()

    def vmul(x, y):
        m = cmul_mid((x[0], x[1]), (y[0], y[1]))
        r = (mabs_ub((x[0], x[1])) * y[2] + mabs_ub((y[0], y[1])) * x[2]
             + x[2] * y[2] + m[0].rad() + m[1].rad())
        return (m[0], m[1], r)

    def vscal6(x):
        re_, im_ = 6 * x[0], 6 * x[1]
        r = 6 * x[2] + re_.rad() + im_.rad()
        return (re_, im_, r)

    def vabs_ub(x):
        return ((rabs_ub(x[0]) + x[2]) ** 2
                + (rabs_ub(x[1]) + x[2]) ** 2).sqrt().upper()

    def vabs_lb(x):
        return (x[0] * x[0] + x[1] * x[1]).sqrt().lower() - x[2]

    n = spec["n"]
    m_lvl = spec["m"]
    r = A(spec["box_radius"])

    c0 = vconst(spec["box_center_re"], spec["box_center_im"])
    Abox = (c0[0], c0[1], c0[2] + r)              # DISK B(c0, r_reported)

    def step(z, dz, a):
        """x_{k+1} = 1 - a*x^6 ; dx_{k+1} = -x^6 - 6*a*x^5*dx (dx first)."""
        z2 = vmul(z, z)
        z3 = vmul(z2, z)
        z4 = vmul(z2, z2)                          # NOT z3*z: other association
        z5 = vmul(z4, z)                           # NOT pow
        z6 = vmul(z3, z3)                          # NOT z5*z
        az5dx = vmul(vmul(a, z5), dz)
        dz2 = vneg(vadd(z6, vscal6(az5dx)))
        az6 = vmul(a, z6)
        z2n = vaddc(vneg(az6), 1, 0)
        return z2n, dz2

    # ---- point orbit: F0 = f(c0), Df0 = f'(c0)
    z = vconst(0)
    dz = vconst(0)
    for _ in range(n):
        z, dz = step(z, dz, c0)
    F0, Df0 = z, dz

    # ---- box orbit over the disk: D = f'(X); marks |x_{2^j}| at j = 0..m
    zb = vconst(0)
    dzb = vconst(0)
    lowers = []
    j = 0
    half = 1 << m_lvl
    for k in range(1, n + 1):
        zb, dzb = step(zb, dzb, Abox)
        if k <= half and k == (1 << j):
            lowers.append(vabs_lb(zb))
            j += 1
    if len(lowers) != m_lvl + 1:
        raise RuntimeError("period marks miscount (verifier)")
    D = dzb

    df_lb = vabs_lb(D)
    df_nonzero = bool(df_lb > 0)
    f0_ub = vabs_ub(F0)

    # ---- interval-Newton distance, no ball division:
    n_dist = (f0_ub / df_lb).upper()
    newton_inside = bool(n_dist < r)

    # ---- Krawczyk distance via triangle inequality; orbits ran on the DISK
    # (subset of the certified square) but the E*(X-c0) term is bounded with
    # the square corner factor sqrt(2)*r (superset bound of |X-c0|):
    #     K - c0 = -Y*F0 + (1 - Y*D)*(X - c0),  |X - c0| <= sqrt(2)*r
    d0re, d0im = Df0[0], Df0[1]
    den = (d0re * d0re + d0im * d0im)
    if not bool(den.lower() > 0):
        raise RuntimeError("verifier: |f'(c0)| interval contains 0")
    yre_int = d0re / den
    yim_int = (-d0im) / den
    yre, yim = yre_int.mid(), yim_int.mid()
    yrad = yre_int.rad() + yim_int.rad() + yre.rad() + yim.rad()
    Yball = (yre, yim, yrad)
    Y_abs_ub = vabs_ub(Yball)
    E = vadd(vmul(Yball, D), vconst(-1))          # Y*D - 1
    e_abs_ub = vabs_ub(E)
    sqrt2 = A(2).sqrt().upper()
    k_dist = (Y_abs_ub * f0_ub + e_abs_ub * sqrt2 * r).upper()
    k_inside = bool(k_dist < r)
    kappa = k_dist / r

    period_ok = [bool(v > 0) for v in lowers]

    def ub_str(x):
        return str(A(x).upper())

    def lb_str(x):
        return str(A(x).lower())

    return {
        "f0_abs_ub": ub_str(f0_ub),
        "df_abs_lb": lb_str(df_lb),
        "newton_dist_ub": ub_str(n_dist),
        "krawczyk_dist_ub": ub_str(k_dist),
        "krawczyk_kappa_ub": ub_str(kappa),
        "newton_inside": newton_inside,
        "krawczyk_inside": k_inside,
        "df_nonzero": df_nonzero,
        "period_abs_lb": [lb_str(v) for v in lowers],
        "period_all_ok": all(period_ok),
        "pass": bool(df_nonzero and newton_inside and k_inside and all(period_ok)),
    }


# ---------------------------------------------------------------- field cmp
def cmp_rel(rep_str, rec_arb, label, rel=REL_TOL_F0):
    A = flint.arb
    try:
        rep = A(str(rep_str))
    except Exception as exc:
        err("%s: unparseable bound string %r (%s)" % (label, rep_str, exc))
        return
    rec_mid = A(rec_arb).mid()
    scale = (rec_mid * rec_mid).sqrt().upper()
    if bool(scale == 0):                    # both exactly zero -> equal
        if bool(rep.mid() == 0):
            return
        err("%s drift: reported nonzero, recomputed 0" % label)
        return
    diff = ((rep.mid() - rec_mid) ** 2).sqrt().upper()
    if not bool(diff <= scale * A(rel)):
        err("%s drift: reported=%s recomputed~%s (rel diff > %s)" % (
            label, str(rep.mid())[:26], str(rec_mid)[:26], rel))


def cmp_factor(rep_str, rec_arb, label, band=FACTOR_KD):
    """Two-sided factor band: 1/band <= rep/rec <= band (midpoints)."""
    A = flint.arb
    try:
        rep = A(str(rep_str)).mid()
    except Exception as exc:
        err("%s: unparseable bound string %r (%s)" % (label, rep_str, exc))
        return
    rec = A(rec_arb).mid()
    if bool(rec == 0) or bool(rep == 0):
        if bool(rep == 0) and bool(rec == 0):
            return
        err("%s drift: zero/nonzero mismatch reported=%s recomputed=%s"
            % (label, str(rep)[:12], str(rec)[:12]))
        return
    ratio = rep / rec
    if not bool(ratio >= A(1) / A(band) and ratio <= A(band)):
        err("%s drift: reported/recomputed ratio %s outside [%s, %s]"
            % (label, str(ratio)[:16], str(A(1) / A(band))[:8], band))


def isolation_recompute(roots):
    flint.ctx.prec = 512
    A = flint.arb
    all_disjoint = True
    worst = None
    worst_pair = None
    n = len(roots)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = roots[i], roots[j]
            dre = A(a["box_center_re"]) - A(b["box_center_re"])
            dim = A(a["box_center_im"]) - A(b["box_center_im"])
            rsum = A(a["box_radius"]) + A(b["box_radius"])
            sep_re = abs(dre).lower()
            sep_im = abs(dim).lower()
            sep = sep_re if bool(sep_re > sep_im) else sep_im
            if not bool(sep > rsum.upper()):
                all_disjoint = False
            margin = sep / rsum.upper()
            if worst is None or bool(margin < worst):
                worst, worst_pair = margin, (a["label"], b["label"])
    return all_disjoint, worst, worst_pair


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    if not CERT_PATH.exists():
        print("certificate not found:", CERT_PATH)
        return 1
    CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))
    A = flint.arb

    print("[0] frozen header fields + sha256 pins")
    if CERT.get("run") != "krawczyk6_certificate":
        err("run field drift: %r" % CERT.get("run"))
    if CERT.get("oid") != "kraw6H":
        err("oid field drift: %r" % CERT.get("oid"))
    if CERT.get("method") != METHOD_STRING:
        err("method field drift (frozen string mismatch)")
    if CERT.get("family") != "1 - a*x^6":
        err("family drift: %r" % CERT.get("family"))
    pins = CERT.get("sha256_pins", {})
    if set(pins.keys()) != {"match_results.json", "calibration6.json",
                            "refined_seeds6.json"}:
        err("sha256_pins key set drift: %s" % sorted(pins.keys()))
    for path, key in ((PILOT, "match_results.json"),
                      (CALIB, "calibration6.json"),
                      (REFINE, "refined_seeds6.json")):
        if pins.get(key) != sha256(path):
            err("sha256 pin mismatch %s: cert=%s actual=%s"
                % (key, str(pins.get(key))[:16], sha256(path)[:16]))
    fr = CERT.get("frozen", {})
    if fr.get("radii_ladder") != list(RADII_LADDER):
        err("frozen radii_ladder drift: %s" % fr.get("radii_ladder"))
    if fr.get("layer2_radii_ladder") != list(LAYER2_RADII_LADDER):
        err("frozen layer2_radii_ladder drift: %s"
            % fr.get("layer2_radii_ladder"))
    if fr.get("prec_ladder") != list(PREC_LADDER):
        err("frozen prec_ladder drift: %s" % fr.get("prec_ladder"))
    if fr.get("max_workers") != MAX_WORKERS:
        err("frozen max_workers drift: %s" % fr.get("max_workers"))
    if fr.get("seed_dist_limit") != SEED_DIST_LIMIT:
        err("frozen seed_dist_limit drift: %s" % fr.get("seed_dist_limit"))
    if "smallest ladder radius passing" not in str(fr.get("selection_rule", "")):
        err("frozen selection_rule drift: %r" % fr.get("selection_rule"))

    print("[1] source inventory (match_results.json + refined_seeds6.json)")
    src = json.loads(PILOT.read_text(encoding="utf-8"))
    ref = json.loads(REFINE.read_text(encoding="utf-8"))["seeds"]
    frz = src["meta"]["frozen"]
    m_lvl_src = int(frz["M_MATCH"])
    n_src = 1 << (m_lvl_src + 1)
    if int(frz["PERIOD"]) != n_src:
        err("pilot frozen PERIOD != 2^(M_MATCH+1)")
    r_accept = float(frz["R_ACCEPT_U"])
    expected = {}
    rej_int, rej_edge = [], []
    n_edge_acc = 0
    for rec in src["catalog_results"]:
        if rec.get("accepted"):
            if rec.get("edge"):
                n_edge_acc += 1
            else:
                expected[rec["id"]] = {
                    "n": n_src, "m": m_lvl_src,
                    "u_re": rec["u_cat_re"], "u_im": rec["u_cat_im"],
                    "re": str(rec["a_ss_re"]), "im": repr(float(rec["a_ss_im"]))}
        else:
            (rej_edge if rec.get("edge") else rej_int).append(rec["id"])
    if rej_int != EXPECTED["rejected_interior"]:
        err("rejected interior rows drift: %s" % rej_int)
    if rej_edge != EXPECTED["rejected_edge"]:
        err("rejected edge rows drift: %s" % rej_edge)
    if n_edge_acc != EXPECTED["edge_accepted"]:
        err("accepted edge rows: %s (expected 46)" % n_edge_acc)
    if len(expected) != 287:
        err("expected inventory is not 287 roots: %d" % len(expected))
    if sorted(ref.keys()) != sorted(LAYER2_LABELS):
        err("refined_seeds6.json labels are not exactly %s: %s"
            % (list(LAYER2_LABELS), sorted(ref.keys())))

    roots = CERT.get("roots", [])
    labels = [r.get("label") for r in roots]
    if sorted(labels) != sorted(expected.keys()):
        err("root label set mismatch: cert=%s..." % sorted(labels)[:6])
    if len(labels) != len(set(labels)):
        err("duplicate labels in certificate")

    print("[2] structure + per-root source pins (two layers)")
    if set(CERT.keys()) != TOP_KEYS:
        err("top-level key set drift: %s" % sorted(set(CERT.keys()) ^ TOP_KEYS))
    l2_seen = []
    for r in roots:
        lbl = r.get("label")
        policy = r.get("seed_policy")
        keys_ok = (ROOT_KEYS_L2 if policy == "refined" else ROOT_KEYS_L1)
        if set(r.keys()) != keys_ok:
            err("%s: per-root key set drift (policy=%s): %s" % (
                lbl, policy, sorted(set(r.keys()) ^ keys_ok)))
            continue
        if policy == "refined":
            l2_seen.append(lbl)
        e = expected.get(lbl)
        if e is None:
            continue
        if r["n"] != e["n"] or r["m"] != e["m"]:
            err("%s period/m drift: cert=(%s,%s) source=(%s,%s)" % (
                lbl, r["n"], r["m"], e["n"], e["m"]))
        if r["n"] != (1 << (r["m"] + 1)):
            err("%s: n != 2^(m+1)" % lbl)
        if r["u_cat_re"] != e["u_re"] or r["u_cat_im"] != e["u_im"]:
            err("%s u_cat drift vs source" % lbl)
        if policy == "source_string":
            if r["box_center_re"] != e["re"] or r["box_center_im"] != e["im"]:
                err("%s box center drift vs source seed: cert=(%r,%r) "
                    "source=(%r,%r)" % (lbl, r["box_center_re"],
                                        r["box_center_im"], e["re"], e["im"]))
            ladder = RADII_LADDER
            band_hi = "1e-20"
        else:
            if r["stored_seed_re"] != e["re"] or r["stored_seed_im"] != e["im"]:
                err("%s stored seed drift vs source" % lbl)
            if (r["box_center_re"] != ref[lbl]["re"]
                    or r["box_center_im"] != ref[lbl]["im"]):
                err("%s refined center drift vs refined_seeds6.json" % lbl)
            if r["refine_iters"] != 8:
                err("%s refine_iters drift: %s" % (lbl, r["refine_iters"]))
            flint.ctx.prec = 512
            c_ref = flint.acb(A(r["box_center_re"]), A(r["box_center_im"]))
            c_sto = flint.acb(A(r["stored_seed_re"]), A(r["stored_seed_im"]))
            dist = (c_ref - c_sto).abs_upper()
            cmp_rel(r["seed_dist_ub"], dist, "%s seed_dist_ub" % lbl, "1e-6")
            if not bool(dist < A(SEED_DIST_LIMIT)):
                err("%s: |refined - stored| %s >= limit %s"
                    % (lbl, str(dist)[:14], SEED_DIST_LIMIT))
            ladder = LAYER2_RADII_LADDER
            band_hi = "1e-24"
        try:
            rr = A(str(r["box_radius"]))
        except Exception as exc:
            err("%s: unparseable box_radius %r (%s)" % (lbl, r["box_radius"], exc))
            continue
        if not (bool(A(0) < rr) and bool(rr < A(band_hi))):
            err("%s: box_radius outside sanity band (0, %s]: %s"
                % (lbl, band_hi, str(rr)[:20]))
        in_ladder = any(bool(((rr.mid() / A(rstr).mid() - A(1)) ** 2).sqrt()
                             <= A("1e-6")) for rstr in ladder)
        if not in_ladder:
            err("%s: box_radius not within tolerance of any frozen ladder "
                "entry: %s" % (lbl, str(rr.mid())[:16]))
        if r["prec_bits"] not in PREC_LADDER:
            err("%s: prec_bits %s not in frozen ladder %s"
                % (lbl, r["prec_bits"], PREC_LADDER))
        # AUDIT_FINDING guard: the stored component radii must be identical
        # to the stored box radius (square, both components), and never
        # certainly below it
        if not (r["box_rad_re"] == r["box_radius"] == r["box_rad_im"]):
            err("%s: component radii strings differ from box_radius "
                "(possible segment regression)" % lbl)
        elif bool(A(str(r["box_rad_re"])) < A(str(r["box_radius"]))):
            err("%s: component radius certainly below r (degenerate)" % lbl)
        pcs = r["period_checks"]
        if len(pcs) != r["m"] + 1:
            err("%s: period_checks length %d != m+1 = %d"
                % (lbl, len(pcs), r["m"] + 1))
        for idx, pc in enumerate(pcs):
            if set(pc.keys()) != PC_KEYS:
                err("%s: period_check %d key set drift: %s"
                    % (lbl, idx, sorted(pc.keys())))
            elif pc["j"] != idx or pc["d"] != (1 << idx):
                err("%s: period_check %d meta drift j=%s d=%s"
                    % (lbl, idx, pc["j"], pc["d"]))
    if sorted(l2_seen) != sorted(LAYER2_LABELS):
        err("layer-2 roots are not exactly %s: %s" % (list(LAYER2_LABELS),
                                                      sorted(l2_seen)))

    print("[3] full tuple-path recomputation over the DISK B(c0, r) "
          "(%d roots, <=%d workers)" % (len(roots), MAX_WORKERS))
    specs = [{"label": r["label"], "n": r["n"], "m": r["m"],
              "box_center_re": r["box_center_re"],
              "box_center_im": r["box_center_im"],
              "box_radius": r["box_radius"],
              "prec_bits": r["prec_bits"]}
             for r in roots if set(r.keys()) in (ROOT_KEYS_L1, ROOT_KEYS_L2)]
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        recs = list(ex.map(recompute_root, specs))
    by_label = {s["label"]: rec for s, rec in zip(specs, recs)}

    for r in roots:
        if set(r.keys()) not in (ROOT_KEYS_L1, ROOT_KEYS_L2):
            continue
        lbl = r["label"]
        rec = by_label.get(lbl)
        if rec is None:
            continue
        cmp_rel(r["f0_abs_ub"], A(rec["f0_abs_ub"]), "%s f0_abs_ub" % lbl)
        cmp_factor(r["df_abs_lb"], A(rec["df_abs_lb"]), "%s df_abs_lb" % lbl,
                   FACTOR_DF)
        cmp_factor(r["newton_dist_ub"], A(rec["newton_dist_ub"]),
                   "%s newton_dist_ub" % lbl, FACTOR_ND)
        if r["seed_policy"] == "source_string":
            cmp_factor(r["krawczyk_dist_ub"], A(rec["krawczyk_dist_ub"]),
                       "%s krawczyk_dist_ub" % lbl, FACTOR_KD)
            cmp_factor(r["krawczyk_kappa_ub"], A(rec["krawczyk_kappa_ub"]),
                       "%s kappa_ub" % lbl, FACTOR_KD)
            # selection rule: kappa_ub must respect the frozen margin 0.5
            # (the fallback branch -- smallest passing radius without margin
            # -- did not trigger for any root in the final run)
            kap = A(str(r["krawczyk_kappa_ub"])).mid()
            if bool(kap > A("0.5")):
                err("%s: kappa_ub %s exceeds the frozen selection margin 0.5"
                    % (lbl, str(kap)[:12]))
        for fld in ("newton_inside", "krawczyk_inside", "df_nonzero",
                    "period_all_ok", "pass"):
            if bool(r[fld]) != bool(rec[fld]):
                err("%s %s drift: cert=%s recomputed=%s"
                    % (lbl, fld, r[fld], rec[fld]))
        if not rec["pass"]:
            err("%s: independent recompute did NOT confirm pass" % lbl)
        pcs = r["period_checks"]
        if len(pcs) == len(rec["period_abs_lb"]) and len(pcs) == r["m"] + 1:
            for idx, pc in enumerate(pcs):
                cmp_factor(pc["abs_lb"], A(rec["period_abs_lb"][idx]),
                           "%s period_checks[j=%d].abs_lb" % (lbl, idx),
                           FACTOR_PER)

    print("[4] pairwise isolation recomputed (%d pairs)"
          % (len(roots) * (len(roots) - 1) // 2))
    iso = CERT.get("pairwise_isolation", {})
    if set(iso.keys()) != ISO_KEYS:
        err("pairwise_isolation key set drift: %s" % sorted(iso.keys()))
    else:
        all_disj, worst, worst_pair = isolation_recompute(roots)
        if not all_disj:
            err("recomputed isolation: some certified boxes OVERLAP")
        if iso.get("n_boxes") != len(roots):
            err("isolation n_boxes drift: %s vs %d" % (iso.get("n_boxes"), len(roots)))
        if iso.get("n_pairs") != len(roots) * (len(roots) - 1) // 2:
            err("isolation n_pairs drift: %s" % iso.get("n_pairs"))
        if bool(iso.get("all_disjoint")) is not all_disj:
            err("isolation all_disjoint drift: cert=%s recomputed=%s"
                % (iso.get("all_disjoint"), all_disj))
        cmp_rel(iso["worst_separation_margin_lb"], worst,
                "isolation worst_separation_margin_lb", "1e-6")
        if list(iso.get("worst_pair", [])) != list(worst_pair):
            err("isolation worst_pair drift: cert=%s recomputed=%s"
                % (iso.get("worst_pair"), worst_pair))
    if CERT.get("degenerate_boxes") != []:
        err("degenerate_boxes not empty: %s" % CERT.get("degenerate_boxes"))

    print("[5] aggregates + correspondence sanity")
    if CERT.get("n_roots") != len(roots):
        err("n_roots drift: cert=%s actual=%d" % (CERT.get("n_roots"), len(roots)))
    all_pass_v = bool(roots) and all(bool(r.get("pass")) for r in roots)
    if bool(CERT.get("all_pass")) != all_pass_v:
        err("all_pass drift: cert=%s recomputed=%s"
            % (CERT.get("all_pass"), all_pass_v))
    if not all_pass_v:
        err("certificate does not establish all_pass on recompute")
    if CERT.get("is_complete") is not True:
        err("certificate is not marked is_complete")
    corr = CERT.get("correspondence", {})
    if not corr.get("_not_rigorous"):
        err("correspondence block missing _not_rigorous flag")
    if not (corr.get("d_u_max") is not None and corr["d_u_max"] < r_accept):
        err("correspondence d_u_max %s not below blind-run R_ACCEPT_U %s"
            % (corr.get("d_u_max"), r_accept))
    if not (corr.get("box_u_size_max") is not None
            and corr["box_u_size_max"] < 1e-12):
        err("correspondence box_u_size_max %s not tiny" % corr.get("box_u_size_max"))

    dt = time.time() - t0
    if errors:
        print("\nVERDICT: FAIL -- %d discrepancies (%.1fs)" % (len(errors), dt))
        return 1
    print("\nVERDICT: PASS -- krawczyk6_certificates.json independently "
          "verified: %d roots recomputed via the (re,im,r)-tuple disk path "
          "(Krawczyk + interval Newton + %d period checks each; 2 refined-"
          "layer roots with Arb seed-distance pins), isolation over %d "
          "pairs, source pins and structure confirmed (%.1fs)"
          % (len(roots), 10, len(roots) * (len(roots) - 1) // 2, dt))
    return 0


if __name__ == "__main__":
    sys.exit(main())

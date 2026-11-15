# -*- coding: utf-8 -*-
# [oid:kraw4E] krawczyk4_verify.py -- INDEPENDENT verifier for
# krawczyk4_certificates.json.  Exit 0 = verified, exit 1 = any drift.
#
# SEPARATE CODE PATH from krawczyk4_certificate.py (same independence pair
# as the audited mside package): complex balls are hand-rolled 3-tuples
# (re, im, r) -- midpoint components as arb intervals plus ONE explicit ball
# radius absorbing all propagation slack -- with hand-derived slack formulas
# (vmul: |m_x|*r_y + |m_y|*r_x + r_x*r_y + midpoint rounding; vabs_ub/lb
# hand-written), DIFFERENT orbit algebra (x^3 = x^2*x, x^4 = x^2*x^2 vs the
# generator's acb pow), and DIFFERENT operator assemblies:
#   newton_dist    := |f(c0)|_ub / |f'(X)|_lb        (no ball division),
#   krawczyk_dist  := |Y|_ub*|F0|_ub + |1-Y*D|_ub*r  (triangle inequality),
# and the box is propagated as the DISK B(c0, r_reported) (subset of the
# generator's certified square): a domination-direction recheck.
#
# EVERY reported field is independently checked:
#   [0] frozen header strings + sha256 pin of match_results.json;
#   [1] source inventory: 61 labels = accepted pilot rows, source rejected
#       exactly [Q57]; per-root n=2^(m+1)=1024=pilot frozen PERIOD, m=9,
#       box centers EXACT string match to the pilot a_ss seeds, u_cat match;
#   [2] structure: exact key sets; radius in the frozen ladder and in
#       (0, 1e-10]; prec in the frozen ladder; period_checks shape
#       j=0..m with d=2^j; non-degenerate radii (AUDIT_FINDING guard);
#   [3] full tuple-path recomputation per root (pool=4): f0_abs_ub,
#       df_abs_lb, newton_dist_ub, krawczyk_dist_ub, krawczyk_kappa_ub,
#       every period abs_lb (all within 1e-6 relative), all booleans, pass;
#   [4] pairwise isolation recomputed from reported centers + radii;
#   [5] aggregates + correspondence sanity (d_u_max < pilot R_ACCEPT_U).
import hashlib
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

NSTR_RE = re.compile(r'^\(([-+0-9.eE]+)\s*([+-])\s*([-+0-9.eE]+)j\)$')


def parse_nstr(s):
    m = NSTR_RE.match(s.strip())
    if m is None:
        raise ValueError("unparseable nstr: %r" % s)
    return m.group(1), ("-" if m.group(2) == "-" else "") + m.group(3)

HERE = Path(__file__).resolve().parent
PILOT = Path(r"../sec3_pilot_d4/match_results.json")
CERT_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else (
    HERE / "krawczyk4_certificates.json")

MAX_WORKERS = 4               # BRIEF: pool = 4
RADII_LADDER = ("1e-26", "1e-27", "1e-25", "1e-24", "1e-28")
PREC_LADDER = (256, 512, 1024)
REL_TOL = "1e-6"
# krawczyk_dist_ub / kappa_ub: the generator propagates E*W over the SQUARE
# (corner factor up to sqrt(2) on the E*W term) while the verifier bounds the
# DISK with r (and a slightly tighter D); both are valid upper bounds of the
# same true |K-c0|, so their midpoints may differ by a few 1e-6 relative on
# roots where the E*W term is visible.  1e-4 leaves 100x margin to that and
# still catches any plausible forgery (e.g. kappa scaled by 0.01).
REL_TOL_KDIST = "1e-4"

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic for f_a(x)=1-a*x^4: "
    "X=c0+[+/-r]^2 (square superset of B(c0,r)); F(a)=x_n(a), n=2^(m+1)=1024, "
    "x_0=0; N(X)=c0-F(c0)/F'(X); K(X)=c0-Y*F(c0)+(1-Y*F'(X))*(X-c0), "
    "Y=1/F'(c0); pass iff 0 notin F'(X), N(X) and K(X) strictly inside X, "
    "and inf_{a in X}|x_{2^j}(a)|>0 for j=0..m; implies existence and "
    "uniqueness of a superstable parameter of exact period n in B(c0,r)"
)

TOP_KEYS = {"run", "oid", "method", "sha256_pins", "family", "frozen",
            "n_roots", "roots", "pairwise_isolation", "degenerate_boxes",
            "correspondence", "all_pass", "is_complete"}
ROOT_KEYS = {"label", "source", "n", "m", "u_cat_re", "u_cat_im",
             "box_center_re", "box_center_im", "box_radius", "box_rad_re",
             "box_rad_im", "prec_bits", "f0_abs_ub", "df_abs_lb",
             "newton_dist_ub", "newton_inside", "krawczyk_dist_ub",
             "krawczyk_kappa_ub", "krawczyk_inside", "df_nonzero",
             "period_checks", "period_all_ok", "pass"}
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

    def vscal4(x):
        re_, im_ = 4 * x[0], 4 * x[1]
        r = 4 * x[2] + re_.rad() + im_.rad()
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
        """x_{k+1} = 1 - a*x^4 ; dx_{k+1} = -x^4 - 4*a*x^3*dx (dx first)."""
        z2 = vmul(z, z)
        z3 = vmul(z2, z)
        z4 = vmul(z2, z2)                          # NOT z3*z: other association
        az3dx = vmul(vmul(a, z3), dz)
        dz2 = vneg(vadd(z4, vscal4(az3dx)))
        az4 = vmul(a, z4)
        z2n = vaddc(vneg(az4), 1, 0)
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

    # ---- Krawczyk distance via triangle inequality over the disk:
    #     K - c0 = -Y*F0 + (1 - Y*D)*(X - c0),  |X - c0| <= r
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
    k_dist = (Y_abs_ub * f0_ub + e_abs_ub * r).upper()
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
def cmp_bound(rep_str, rec_arb, label, rel=REL_TOL):
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

    print("[0] frozen header fields + sha256 pin")
    if CERT.get("run") != "krawczyk4_certificate":
        err("run field drift: %r" % CERT.get("run"))
    if CERT.get("oid") != "kraw4E":
        err("oid field drift: %r" % CERT.get("oid"))
    if CERT.get("method") != METHOD_STRING:
        err("method field drift (frozen string mismatch)")
    if CERT.get("family") != "1 - a*x^4":
        err("family drift: %r" % CERT.get("family"))
    pins = CERT.get("sha256_pins", {})
    if set(pins.keys()) != {"match_results.json"}:
        err("sha256_pins key set drift: %s" % sorted(pins.keys()))
    if pins.get("match_results.json") != sha256(PILOT):
        err("sha256 pin mismatch match_results.json: cert=%s actual=%s" % (
            pins.get("match_results.json"), sha256(PILOT)))
    fr = CERT.get("frozen", {})
    if fr.get("radii_ladder") != list(RADII_LADDER):
        err("frozen radii_ladder drift: %s" % fr.get("radii_ladder"))
    if fr.get("prec_ladder") != list(PREC_LADDER):
        err("frozen prec_ladder drift: %s" % fr.get("prec_ladder"))
    if fr.get("max_workers") != MAX_WORKERS:
        err("frozen max_workers drift: %s" % fr.get("max_workers"))

    print("[1] source inventory (match_results.json)")
    src = json.loads(PILOT.read_text(encoding="utf-8"))
    m_lvl_src = int(src["meta"]["frozen"]["M_MATCH"])
    n_src = 1 << (m_lvl_src + 1)
    if int(src["meta"]["frozen"]["PERIOD"]) != n_src:
        err("pilot frozen PERIOD != 2^(M_MATCH+1)")
    r_accept = float(src["meta"]["frozen"]["R_ACCEPT_U"])
    expected = {}
    rejected = []
    for rec in src["catalog_results"]:
        if rec.get("accepted"):
            expected[rec["id"]] = {"n": n_src, "m": m_lvl_src,
                                   "u_re": rec["u_cat_re"], "u_im": rec["u_cat_im"]}
        else:
            rejected.append(rec["id"])
    if rejected != ["Q57"]:
        err("rejected pilot rows are not exactly [Q57]: %s" % rejected)
    if len(expected) != 61:
        err("expected inventory is not 61 roots: %d" % len(expected))

    roots = CERT.get("roots", [])
    labels = [r.get("label") for r in roots]
    if sorted(labels) != sorted(expected.keys()):
        err("root label set mismatch: cert=%s..." % sorted(labels)[:6])
    if len(labels) != len(set(labels)):
        err("duplicate labels in certificate")

    print("[2] structure + per-root source pins")
    if set(CERT.keys()) != TOP_KEYS:
        err("top-level key set drift: %s" % sorted(set(CERT.keys()) ^ TOP_KEYS))
    for r in roots:
        lbl = r.get("label")
        if set(r.keys()) != ROOT_KEYS:
            err("%s: per-root key set drift: %s" % (
                lbl, sorted(set(r.keys()) ^ ROOT_KEYS)))
            continue
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
        # box centers must EXACTLY equal the pilot seed strings (via the
        # generator's nstr parser): rebuild source strings and compare
        src_rec = next(s for s in src["catalog_results"] if s["id"] == lbl)
        s_re, s_im = parse_nstr(src_rec["a_ss"])
        if r["box_center_re"] != s_re or r["box_center_im"] != s_im:
            err("%s box center drift vs source seed: cert=(%r,%r) "
                "source=(%r,%r)" % (lbl, r["box_center_re"],
                                    r["box_center_im"], s_re, s_im))
        try:
            rr = A(str(r["box_radius"]))
        except Exception as exc:
            err("%s: unparseable box_radius %r (%s)" % (lbl, r["box_radius"], exc))
            continue
        if not (bool(A(0) < rr) and bool(rr < A("1e-10"))):
            err("%s: box_radius outside sanity band (0, 1e-10]: %s"
                % (lbl, str(rr)[:20]))
        in_ladder = any(bool(((rr.mid() / A(rstr).mid() - A(1)) ** 2).sqrt()
                             <= A(REL_TOL)) for rstr in RADII_LADDER)
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

    print("[3] full tuple-path recomputation over the DISK B(c0, r) "
          "(%d roots, <=%d workers)" % (len(roots), MAX_WORKERS))
    specs = [{"label": r["label"], "n": r["n"], "m": r["m"],
              "box_center_re": r["box_center_re"],
              "box_center_im": r["box_center_im"],
              "box_radius": r["box_radius"],
              "prec_bits": r["prec_bits"]}
             for r in roots if set(r.keys()) == ROOT_KEYS]
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        recs = list(ex.map(recompute_root, specs))
    by_label = {s["label"]: rec for s, rec in zip(specs, recs)}

    for r in roots:
        if set(r.keys()) != ROOT_KEYS:
            continue
        lbl = r["label"]
        rec = by_label.get(lbl)
        if rec is None:
            continue
        cmp_bound(r["f0_abs_ub"], A(rec["f0_abs_ub"]), "%s f0_abs_ub" % lbl)
        cmp_bound(r["df_abs_lb"], A(rec["df_abs_lb"]), "%s df_abs_lb" % lbl)
        cmp_bound(r["newton_dist_ub"], A(rec["newton_dist_ub"]),
                  "%s newton_dist_ub" % lbl)
        cmp_bound(r["krawczyk_dist_ub"], A(rec["krawczyk_dist_ub"]),
                  "%s krawczyk_dist_ub" % lbl, rel=REL_TOL_KDIST)
        cmp_bound(r["krawczyk_kappa_ub"], A(rec["krawczyk_kappa_ub"]),
                  "%s kappa_ub" % lbl, rel=REL_TOL_KDIST)
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
                cmp_bound(pc["abs_lb"], A(rec["period_abs_lb"][idx]),
                          "%s period_checks[j=%d].abs_lb" % (lbl, idx))

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
        cmp_bound(iso["worst_separation_margin_lb"], worst,
                  "isolation worst_separation_margin_lb")
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
        err("correspondence d_u_max %s not below pilot R_ACCEPT_U %s"
            % (corr.get("d_u_max"), r_accept))
    if not (corr.get("box_u_size_max") is not None
            and corr["box_u_size_max"] < 1e-15):
        err("correspondence box_u_size_max %s not tiny" % corr.get("box_u_size_max"))

    dt = time.time() - t0
    if errors:
        print("\nVERDICT: FAIL -- %d discrepancies (%.1fs)" % (len(errors), dt))
        return 1
    print("\nVERDICT: PASS -- krawczyk4_certificates.json independently "
          "verified: %d roots recomputed via the (re,im,r)-tuple disk path "
          "(Krawczyk + interval Newton + %d period checks each), isolation "
          "over %d pairs, source pins and structure confirmed (%.1fs)"
          % (len(roots), 10, len(roots) * (len(roots) - 1) // 2, dt))
    return 0


if __name__ == "__main__":
    sys.exit(main())

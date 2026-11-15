# -*- coding: utf-8 -*-
# [oid:msideK1] mside_verify.py -- INDEPENDENT verifier for
# mside_certificate.json. Exit 0 = certificate verified, exit 1 = any drift.
#
# SEPARATE CODE PATH from mside_certificate.py (same independence pair as the
# audited winding_bridge package): complex balls are hand-rolled 3-tuples
# (re, im, r) -- midpoint components as arb intervals plus ONE explicit ball
# radius that absorbs all propagation slack -- with hand-derived propagation
# formulas (vsqu/vmul slack terms 2|M|r + r^2 etc.), hand-written modulus
# bounds (|z|_ub = sqrt(|re|_u^2+|im|_u^2)+r, |z|_lb via component distance
# to the origin), and DIFFERENT operator assemblies:
#   newton_dist := |f(c0)|_ub / |f'(X)|_lb          (no ball division),
#   krawczyk_dist := |Y|_ub*|F0|_ub + |1-Y*D|_ub*r  (triangle inequality,
#                                                    Y = point 1/f'(c0)),
# whereas the generator uses native acb ball division and ball products.
# The verifier re-derives the box radius from the REPORTED value (a superset
# of the certified disk, hence a domination-direction recheck).
#
# EVERY reported field is independently checked (BRIEF acceptance rule,
# lesson K1-K4 of the winding_bridge revision):
#   [0] run/oid/method exact frozen-string match; sha256 pin of results.json;
#   [1] root inventory: 23 labels; Z72 must come from stage3b_z72_corrected;
#       the 22 others from accepted stage3_catalog rows (and results.json's
#       rejected catalog row must be exactly Z72); per-root n, m and the box
#       center strings must EXACTLY equal the source values;
#   [2] structure: exact key sets (top level, per root, per period check);
#       radius in the frozen ladder (up to outward-rounding slack), radius in
#       (0, 1e-10]; prec_bits in the frozen precision ladder;
#       period_checks has exactly m+1 entries j=0..m with d=2^j;
#   [3] full recomputation per root (tuple path, parallel, max_workers=10):
#       f0_abs_ub, df_abs_lb, newton_dist_ub, krawczyk_dist_ub,
#       krawczyk_kappa_ub (all within 1e-6 relative of the recomputation),
#       newton_inside, krawczyk_inside, df_nonzero, period_all_ok, every
#       period_checks abs_lb, pass;
#   [4] aggregates: n_roots, all_pass.
#
# PARALLELISM: ProcessPoolExecutor, max_workers=10 EXACTLY (owner directive).
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
RESULTS_JSON = Path(r"D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/"
                    r"article_v1/v2_work/unified_mandelbrot/results.json")
CERT_PATH = Path(sys.argv[1]) if len(sys.argv) > 1 else (HERE / "mside_certificate.json")

MAX_WORKERS = 10          # hard cap per owner directive (NOT more than 10)
RADII_LADDER = ("1e-27", "1e-28", "1e-26", "1e-25", "1e-29")
PREC_LADDER = (512, 1024, 2048)
REL_TOL = "1e-6"

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic: X=B(c0,r); "
    "N(X)=c0-f(c0)/f'(X); K(X)=c0-Y*f(c0)+(1-Y*f'(X))*(X-c0), Y=1/f'(c0); "
    "pass iff 0 notin f'(X), N(X) and K(X) strictly inside X, and "
    "inf_{c in X}|P_c^{2^j}(0)|>0 for all j=0..m; implies existence and "
    "uniqueness of a zero of P_c^n(0) in X of exact period n"
)

TOP_KEYS = {"run", "oid", "method", "sha256_pins", "n_roots", "roots", "all_pass"}
ROOT_KEYS = {"label", "source", "n", "m", "box_center_re", "box_center_im",
             "box_radius", "prec_bits", "f0_abs_ub", "df_abs_lb",
             "newton_dist_ub", "newton_inside", "krawczyk_dist_ub",
             "krawczyk_kappa_ub", "krawczyk_inside", "df_nonzero",
             "period_checks", "period_all_ok", "pass"}
PC_KEYS = {"j", "d", "abs_lb"}

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
    """Worker: independent (re, im, r)-tuple recomputation of one root's
    certificate fields. Returns dict with arb-strings and booleans."""
    flint.ctx.prec = spec["prec_bits"]
    A = flint.arb

    def rabs_ub(x):
        """|x| upper for a real arb interval: max(|lower|, |upper|)."""
        a, b = abs(x.lower()), abs(x.upper())
        return b if bool(b > a) else a

    def vconst(re_, im_=None):
        return (A(re_), A(0) if im_ is None else A(im_), A(0))

    def vadd(x, y):
        re_, im_ = x[0] + y[0], x[1] + y[1]
        r = x[2] + y[2] + re_.rad() + im_.rad()
        return (re_, im_, r)

    def vaddc(x, re_, im_):
        re2, im2 = x[0] + A(re_), x[1] + A(im_)
        r = x[2] + re2.rad() + im2.rad()
        return (re2, im2, r)

    def cmul_mid(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def mabs_ub(m):
        return (m[0] * m[0] + m[1] * m[1]).sqrt().upper()

    def mabs_lb(m):
        return (m[0] * m[0] + m[1] * m[1]).sqrt().lower()

    def vmul(x, y):
        m = cmul_mid((x[0], x[1]), (y[0], y[1]))
        r = (mabs_ub((x[0], x[1])) * y[2] + mabs_ub((y[0], y[1])) * x[2]
             + x[2] * y[2] + m[0].rad() + m[1].rad())
        return (m[0], m[1], r)

    def vscal2(x):
        re_, im_ = 2 * x[0], 2 * x[1]
        r = 2 * x[2] + re_.rad() + im_.rad()
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
    X = (c0[0], c0[1], c0[2] + r)          # disk B(c0, r_reported)

    # ---- point orbit: F0 = f(c0), Df0 = f'(c0)  (r-parts stay ~0)
    z = vconst(0)
    dz = vconst(0)
    for _ in range(n):
        dz = vaddc(vscal2(vmul(z, dz)), 1, 0)          # 2*z*dz + 1
        z = vaddc(vmul(z, z), spec["box_center_re"], spec["box_center_im"])
    F0 = z
    Df0 = dz

    # ---- box orbit: D = f'(X); marks |P^d(0)| at d = 2^j, j = 0..m
    zb = vconst(0)
    dzb = vconst(0)
    lowers = []
    j = 0
    half = 1 << m_lvl
    for k in range(1, n + 1):
        dzb = vaddc(vscal2(vmul(zb, dzb)), 1, 0)       # 2*zb*dzb + 1
        # zb <- zb^2 + X   (square then add the disk)
        sq = vmul(zb, zb)
        zb = (sq[0] + X[0], sq[1] + X[1], sq[2] + X[2]
              + (sq[0] + X[0]).rad() + (sq[1] + X[1]).rad())
        if k <= half and k == (1 << j):
            lowers.append(vabs_lb(zb))
            j += 1
    if len(lowers) != m_lvl + 1:
        raise RuntimeError("period marks miscount (verifier)")

    D = dzb
    df_lb = vabs_lb(D)                       # |f'(X)| lower
    df_nonzero = bool(df_lb > 0)

    f0_ub = vabs_ub(F0)

    # ---- interval-Newton distance, no ball division:
    #     |N(X) - c0| = |F0/D| <= |F0|_ub / |D|_lb
    n_dist = (f0_ub / df_lb).upper()
    newton_inside = bool(n_dist < r)

    # ---- Krawczyk distance via triangle inequality:
    #     K - c0 = -Y*F0 + (1 - Y*D)*W,  |W| <= r
    Y_abs_ub = (A(1) / vabs_lb(Df0)).upper()          # |1/f'(c0)|
    # Y ball = conj(Df0)/|Df0|^2 on midpoint intervals; slack from the
    # spread of the interval division itself
    d0re, d0im = Df0[0], Df0[1]
    den = (d0re * d0re + d0im * d0im)                  # interval |Df0|^2
    if not bool(den.lower() > 0):
        raise RuntimeError("verifier: |f'(c0)| interval contains 0")
    yre_int = d0re / den
    yim_int = (-d0im) / den
    yre, yim = yre_int.mid(), yim_int.mid()
    yrad = (yre_int.rad() + yim_int.rad() + yre.rad() + yim.rad())
    Yball = (yre, yim, yrad)
    E = vaddc(vmul(Yball, D), -1, 0)                   # Y*D - 1 (|E| symmetric)
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
    """Reported bound string vs recomputed arb ball: midpoints must agree to
    `rel` relative tolerance."""
    A = flint.arb
    try:
        rep = A(str(rep_str))
    except Exception as exc:
        err("%s: unparseable bound string %r (%s)" % (label, rep_str, exc))
        return
    rec_mid = A(rec_arb).mid()
    scale = (rec_mid * rec_mid).sqrt().upper()
    diff = ((rep.mid() - rec_mid) ** 2).sqrt().upper()
    if not bool(diff <= scale * A(rel)):
        err("%s drift: reported=%s recomputed~%s (rel diff > %s)" % (
            label, str(rep.mid())[:26], str(rec_mid)[:26], rel))


def main():
    t0 = time.time()
    if not CERT_PATH.exists():
        print("certificate not found:", CERT_PATH)
        return 1
    CERT = json.loads(CERT_PATH.read_text(encoding="utf-8"))
    A = flint.arb

    print("[0] frozen header fields + sha256 pin")
    if CERT.get("run") != "mside_certificate":
        err("run field drift: %r" % CERT.get("run"))
    if CERT.get("oid") != "msideK1":
        err("oid field drift: %r" % CERT.get("oid"))
    if CERT.get("method") != METHOD_STRING:
        err("method field drift (frozen string mismatch)")
    pins = CERT.get("sha256_pins", {})
    if set(pins.keys()) != {"results.json"}:
        err("sha256_pins key set drift: %s" % sorted(pins.keys()))
    actual_pin = sha256(RESULTS_JSON)
    if pins.get("results.json") != actual_pin:
        err("sha256 pin mismatch results.json: cert=%s actual=%s" % (
            pins.get("results.json"), actual_pin))

    print("[1] source inventory (results.json)")
    src = json.loads(RESULTS_JSON.read_text(encoding="utf-8"))
    expected = {}
    rejected = []
    for rec in src["stage3_catalog"]:
        if rec.get("accepted") and rec.get("c_ss") is not None:
            expected[rec["id"]] = {
                "source": "stage3_catalog", "n": int(rec["period"]),
                "m": int(rec["m"]), "re": str(rec["c_ss"]["re_str"]),
                "im": str(rec["c_ss"]["im_str"])}
        else:
            rejected.append(rec["id"])
    if rejected != ["Z72"]:
        err("rejected catalog rows in source are not exactly [Z72]: %s" % rejected)
    b = src["stage3b_z72_corrected"]
    if not b.get("accepted"):
        err("stage3b_z72_corrected no longer accepted in source")
    expected["Z72"] = {
        "source": "stage3b_z72_corrected", "n": int(b["period"]), "m": int(b["m"]),
        "re": str(b["c_ss"]["re_str"]), "im": str(b["c_ss"]["im_str"])}
    if len(expected) != 23:
        err("expected inventory is not 23 roots: %d" % len(expected))

    roots = CERT.get("roots", [])
    labels = [r.get("label") for r in roots]
    if sorted(labels) != sorted(expected.keys()):
        err("root label set mismatch: cert=%s expected=%s" % (
            sorted(labels), sorted(expected.keys())))
    if len(labels) != len(set(labels)):
        err("duplicate labels in certificate")

    print("[2] structure + per-root source pins")
    if set(CERT.keys()) != TOP_KEYS:
        err("top-level key set drift: %s" % sorted(CERT.keys()))
    for r in roots:
        lbl = r.get("label")
        if set(r.keys()) != ROOT_KEYS:
            err("%s: per-root key set drift: %s" % (lbl, sorted(set(r.keys()) ^ ROOT_KEYS)))
            continue
        e = expected.get(lbl)
        if e is None:
            continue
        if r["source"] != e["source"]:
            err("%s source drift: cert=%s expected=%s" % (lbl, r["source"], e["source"]))
        if r["n"] != e["n"]:
            err("%s period n drift: cert=%s source=%s" % (lbl, r["n"], e["n"]))
        if r["m"] != e["m"]:
            err("%s m drift: cert=%s source=%s" % (lbl, r["m"], e["m"]))
        if r["box_center_re"] != e["re"]:
            err("%s box_center_re drift vs source: cert=%r source=%r" % (
                lbl, r["box_center_re"], e["re"]))
        if r["box_center_im"] != e["im"]:
            err("%s box_center_im drift vs source: cert=%r source=%r" % (
                lbl, r["box_center_im"], e["im"]))
        if r["n"] != (1 << (r["m"] + 1)):
            err("%s: n != 2^(m+1)" % lbl)
        # radius: in ladder, in (0, 1e-10]
        try:
            rr = A(str(r["box_radius"]))
        except Exception as exc:
            err("%s: unparseable box_radius %r (%s)" % (lbl, r["box_radius"], exc))
            continue
        if not (bool(A(0) < rr) and bool(rr < A("1e-10"))):
            err("%s: box_radius outside sanity band (0, 1e-10]: %s" % (lbl, str(rr)[:20]))
        in_ladder = any(bool(((rr.mid() / A(rstr).mid() - A(1)) ** 2).sqrt()
                             <= A(REL_TOL)) for rstr in RADII_LADDER)
        if not in_ladder:
            err("%s: box_radius not within tolerance of any frozen ladder entry" % lbl)
        if r["prec_bits"] not in PREC_LADDER:
            err("%s: prec_bits %s not in frozen ladder %s" % (lbl, r["prec_bits"], PREC_LADDER))
        pcs = r["period_checks"]
        if len(pcs) != r["m"] + 1:
            err("%s: period_checks length %d != m+1 = %d" % (lbl, len(pcs), r["m"] + 1))
        for idx, pc in enumerate(pcs):
            if set(pc.keys()) != PC_KEYS:
                err("%s: period_check %d key set drift: %s" % (lbl, idx, sorted(pc.keys())))

    print("[3] full tuple-path recomputation (%d roots, <=%d workers)" % (
        len(roots), MAX_WORKERS))
    specs = [{"label": r["label"], "n": r["n"], "m": r["m"],
              "box_center_re": r["box_center_re"], "box_center_im": r["box_center_im"],
              "box_radius": r["box_radius"], "prec_bits": r["prec_bits"]}
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
        cmp_bound(r["newton_dist_ub"], A(rec["newton_dist_ub"]), "%s newton_dist_ub" % lbl)
        cmp_bound(r["krawczyk_dist_ub"], A(rec["krawczyk_dist_ub"]), "%s krawczyk_dist_ub" % lbl)
        cmp_bound(r["krawczyk_kappa_ub"], A(rec["krawczyk_kappa_ub"]), "%s kappa_ub" % lbl)
        for fld in ("newton_inside", "krawczyk_inside", "df_nonzero",
                    "period_all_ok", "pass"):
            if bool(r[fld]) != bool(rec[fld]):
                err("%s %s drift: cert=%s recomputed=%s" % (lbl, fld, r[fld], rec[fld]))
        if not rec["pass"]:
            err("%s: independent recompute did NOT confirm pass" % lbl)
        pcs = r["period_checks"]
        if len(pcs) == len(rec["period_abs_lb"]) and len(pcs) == r["m"] + 1:
            for idx, pc in enumerate(pcs):
                if pc["j"] != idx:
                    err("%s period_checks[%d].j=%s != %d" % (lbl, idx, pc["j"], idx))
                if pc["d"] != (1 << idx):
                    err("%s period_checks[%d].d=%s != %d" % (lbl, idx, pc["d"], 1 << idx))
                cmp_bound(pc["abs_lb"], A(rec["period_abs_lb"][idx]),
                          "%s period_checks[j=%d].abs_lb" % (lbl, idx))

    print("[4] aggregates")
    if CERT.get("n_roots") != len(roots):
        err("n_roots drift: cert=%s actual=%d" % (CERT.get("n_roots"), len(roots)))
    all_pass_v = bool(roots) and all(bool(r.get("pass")) for r in roots)
    if bool(CERT.get("all_pass")) != all_pass_v:
        err("all_pass drift: cert=%s recomputed=%s" % (CERT.get("all_pass"), all_pass_v))
    if not all_pass_v:
        err("certificate does not establish all_pass on recompute")

    dt = time.time() - t0
    if errors:
        print("\nVERDICT: FAIL -- %d discrepancies (%.1fs)" % (len(errors), dt))
        return 1
    print("\nVERDICT: PASS -- mside_certificate.json independently verified: "
          "%d roots recomputed via the (re,im,r)-tuple path (Krawczyk, interval "
          "Newton, all %d period checks each), source pins and structure "
          "confirmed (%.1fs)" % (len(roots), 11, dt))
    return 0


if __name__ == "__main__":
    sys.exit(main())

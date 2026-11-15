# -*- coding: utf-8 -*-
# [oid:msideK1] mside_certificate.py -- strict certification of the 23
# M-side superstable parameters (Mandelbrot family P_c(z)=z^2+c).
#
# TASK (BRIEF.md, mside_superstable): turn each of the 23 numerical Newton
# roots in v2_work/unified_mandelbrot/results.json (mpmath, 150 dps) into a
# THEOREM: "a superstable parameter of exact period n=2^(m+1) exists and is
# unique in the complex ball X_k", via interval Newton / Krawczyk in ball
# arithmetic (python-flint / Arb). Z72 uses the corrected root from
# stage3b_z72_corrected, NOT the legacy catalog row.
#
# METHOD (per root, all inclusions via Arb ball operations, no float
# comparisons anywhere in the decision path):
#   f(c)   := P_c^n(0),  z_0=0,  z_{k+1}=z_k^2+c;
#   f'(c)  := dz_n(c),   dz_0=0, dz_{k+1}=2*z_k*dz_k+1.
#   X := c0 + W, where W = [-r_str, r_str] x [-r_str, r_str] is the SQUARE,
#   a superset of the closed disk B(0, r_str); c0 = the numerical root
#   (exact decimal string from results.json), r_str from a frozen
#   deterministic ladder. AUDIT_FIX (AUDIT_FINDING.md): the previous version
#   constructed W with A(0) as the imaginary part, which degenerated the
#   enclosure to the horizontal SEGMENT [Re c0 +/- r] x {Im c0}, so f'(X)
#   and the period checks covered a segment, not the claimed disk.
#   r := outward upper bound of the MAXIMUM of the radii of BOTH components
#   of W (acb .rad() is not trusted for a two-component ball), and a
#   regression guard before any orbit asserts rad(Re X) >= r and
#   rad(Im X) >= r, which also yields B(c0, r) subset X (disk in square).
#   Point evaluation: F0 = f(c0), Df0 = f'(c0) (balls at c0).
#   Box inclusion:    D  = f'(X)   (ball propagation of (z,dz) over X).
#   Interval Newton operator (Moore):
#       N(X) := c0 - F0/D.
#   Krawczyk operator (Krawczyk 1969; Y := 1/Df0):
#       K(X) := c0 - Y*F0 + (1 - Y*D)*(X - c0).
#   PASS rule (all must hold, each a strict Arb ball comparison):
#     (i)   0 not in D            (D.abs_lower() > 0);
#     (ii)  N(X) strictly inside X ((N-c0).abs_upper() < r);
#     (iii) K(X) strictly inside X ((K-c0).abs_upper() < r);
#     (iv)  for every j=0..m: inf_{c in X} |P_c^{2^j}(0)| > 0 (ball
#           propagation of the orbit over X, lower bound of the modulus).
#   THEOREM (see REPORT_MSIDE.md for the statement and proof sketch).
#   All enclosures (D, N, K, period checks) are propagated over the
#   superset square X with B(c0, r) subset X (regression guard). (i)+(iii):
#   |K-c0|_ub < r puts the ball K inside the open disk B(c0, r) subset X;
#   k(x) = x - Y*f(x) maps every x in the convex square X into K (complex
#   mean-value form over the segment [c0, x] subset X, whose f'-image lies
#   in D), so Brouwer's fixed point theorem gives a zero of f in
#   K subset B(c0, r); (i) plus the segment-integral mean-value argument
#   gives uniqueness in all of X. (ii) is the Moore interval-Newton
#   containment: every zero of f in X lies in N(X) with |N-c0|_ub < r, i.e.
#   inside the disk B(c0, r). With (iv): for every c in X, the critical
#   orbit does not return to 0 before iteration n, hence the unique zero
#   c* in B(c0, r) is a superstable parameter of EXACT period n. The
#   theorem statement for the DISK B(c0, r) is preserved verbatim.
#
# PARALLELISM: ProcessPoolExecutor, max_workers=10 EXACTLY (hard owner
# directive: the machine is shared; never raise this).
#
# INPUTS (read-only, sha256-pinned in the output):
#   article_v1/v2_work/unified_mandelbrot/results.json
#
# OUTPUT: mside_certificate.json (this folder only).
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
OUT_JSON = HERE / "mside_certificate.json"

MAX_WORKERS = 10          # hard cap per owner directive (NOT more than 10)
RADII_LADDER = ("1e-27", "1e-28", "1e-26", "1e-25", "1e-29")
PREC_LADDER = (512, 1024, 2048)

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic: X=B(c0,r); "
    "N(X)=c0-f(c0)/f'(X); K(X)=c0-Y*f(c0)+(1-Y*f'(X))*(X-c0), Y=1/f'(c0); "
    "pass iff 0 notin f'(X), N(X) and K(X) strictly inside X, and "
    "inf_{c in X}|P_c^{2^j}(0)|>0 for all j=0..m; implies existence and "
    "uniqueness of a zero of P_c^n(0) in X of exact period n"
)

REJECTED_IDS_EXPECTED = ["Z72"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def ub_str(x) -> str:
    """Outer (outward-rounded-up) decimal ball string for an upper bound."""
    return str(flint.arb(x).upper())


def lb_str(x) -> str:
    """Outer (outward-rounded-down) decimal ball string for a lower bound."""
    return str(flint.arb(x).lower())


def load_roots():
    d = json.loads(RESULTS_JSON.read_text(encoding="utf-8"))
    roots = []
    rejected = []
    for rec in d["stage3_catalog"]:
        if rec.get("accepted") and rec.get("c_ss") is not None:
            roots.append({
                "label": rec["id"], "source": "stage3_catalog",
                "n": int(rec["period"]), "m": int(rec["m"]),
                "re": str(rec["c_ss"]["re_str"]), "im": str(rec["c_ss"]["im_str"]),
            })
        else:
            rejected.append(rec["id"])
    if rejected != REJECTED_IDS_EXPECTED:
        raise RuntimeError("unexpected rejected catalog rows: %s" % rejected)
    b = d["stage3b_z72_corrected"]
    if not b.get("accepted"):
        raise RuntimeError("stage3b_z72_corrected is not accepted")
    roots.append({
        "label": "Z72", "source": "stage3b_z72_corrected",
        "n": int(b["period"]), "m": int(b["m"]),
        "re": str(b["c_ss"]["re_str"]), "im": str(b["c_ss"]["im_str"]),
    })
    roots.sort(key=lambda r: r["label"])
    if len(roots) != 23:
        raise RuntimeError("expected 23 roots, got %d" % len(roots))
    for r in roots:
        if r["n"] != (1 << (r["m"] + 1)):
            raise RuntimeError("n != 2^(m+1) for %s" % r["label"])
    return roots


def attempt(spec, r_str, prec):
    """One certification attempt at (radius, precision). Returns record or
    None-flags structure; all decision comparisons are strict ball tests."""
    flint.ctx.prec = prec
    A = flint.arb
    C = flint.acb

    c0 = C(A(spec["re"]), A(spec["im"]))
    # AUDIT_FIX: enclosure W is the SQUARE [+/- r_str] x [+/- r_str],
    # a superset of the closed disk B(0, r_str) (the old code had A(0) as
    # the imaginary part -> degenerate horizontal segment, see
    # AUDIT_FINDING.md).
    W = C(A("[+/- %s]" % r_str), A("[+/- %s]" % r_str))
    # r: outward upper bound of the MAXIMUM of the radii of BOTH components
    # of W (do not rely on acb .rad() for a two-component ball).
    r_re = W.real.rad().upper()
    r_im = W.imag.rad().upper()
    r = r_re if bool(r_re > r_im) else r_im
    X = c0 + W                                # X = square c0 + [+/- r]^2 >= B(c0, r)
    # AUDIT_FIX regression guard, BEFORE any orbit propagation: BOTH
    # components of X must carry radius >= r. This catches a recurrence of
    # the segment bug AND guarantees B(c0, r) subset X, which is what makes
    # the theorem valid for the DISK B(c0, r) while every enclosure (D, N,
    # K, period checks) is propagated over the superset square X.
    for comp in (X.real, X.imag):
        if not bool(comp.rad() >= r):
            raise RuntimeError("degenerate box: a component of X has radius "
                               "< r (segment regression, AUDIT_FINDING.md)")
    n = spec["n"]
    m = spec["m"]

    # ---- point orbit: F0 = f(c0), Df0 = f'(c0)
    z = C(0)
    dz = C(0)
    for _ in range(n):
        dz = 2 * z * dz + 1
        z = z * z + c0
    F0 = z
    Df0 = dz

    # ---- box orbit: D = f'(X), and Z_d = P_c^d(0) over X for d=2^j, j<=m
    zb = C(0)
    dzb = C(0)
    half = 1 << m
    lowers = []
    j = 0
    for k in range(1, n + 1):
        dzb = 2 * zb * dzb + 1
        zb = zb * zb + X
        if k <= half and k == (1 << j):
            lowers.append(zb.abs_lower())
            j += 1
    D = dzb
    if len(lowers) != m + 1:
        raise RuntimeError("period marks miscount")

    df_lb = D.abs_lower()
    df_nonzero = bool(df_lb > 0)

    # ---- interval Newton operator
    N = c0 - F0 / D
    n_dist = (N - c0).abs_upper()
    newton_inside = bool(n_dist < r)

    # ---- Krawczyk operator
    Y = 1 / Df0
    E = 1 - Y * D
    K = c0 - Y * F0 + E * W
    k_dist = (K - c0).abs_upper()
    k_inside = bool(k_dist < r)
    kappa = k_dist / r

    period_ok = [bool(v > 0) for v in lowers]
    period_all_ok = all(period_ok)

    return {
        "box_radius": ub_str(r),
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
        "period_all_ok": period_all_ok,
        "_pass": bool(df_nonzero and newton_inside and k_inside and period_all_ok),
    }


def certify_root(spec):
    t0 = time.time()
    chosen = None
    for prec in PREC_LADDER:
        for r_str in RADII_LADDER:
            rec = attempt(spec, r_str, prec)
            if rec["_pass"]:
                chosen = rec
                break
        if chosen is not None:
            break
    if chosen is None:
        # freeze at the first (most conservative standard) configuration so
        # the failure is fully reported and auditable
        chosen = attempt(spec, RADII_LADDER[0], PREC_LADDER[0])
    out = {
        "label": spec["label"],
        "source": spec["source"],
        "n": spec["n"],
        "m": spec["m"],
        "box_center_re": spec["re"],
        "box_center_im": spec["im"],
    }
    out.update(chosen)
    out["pass"] = out.pop("_pass")
    # NOTE: no timestamps/elapsed fields -- every reported field must be
    # independently re-checkable by mside_verify.py (BRIEF acceptance rule).
    return out


def main():
    t0 = time.time()
    print("=== mside_certificate.py [oid:msideK1]: Krawczyk/Newton certification "
          "of 23 superstable parameters ===")
    roots = load_roots()
    print("loaded %d roots from results.json (Z72 <- stage3b_z72_corrected)" % len(roots))

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        records = list(ex.map(certify_root, roots))

    n_pass = sum(1 for r in records if r["pass"])
    for r in records:
        print("  %-4s n=%d r<=%s kappa<=%s N-dist<=%s min|P^d(0)|>=%s pass=%s" % (
            r["label"], r["n"], r["box_radius"][:12], r["krawczyk_kappa_ub"][:12],
            r["newton_dist_ub"][:12],
            min(pc["abs_lb"] for pc in r["period_checks"])[:12], r["pass"]))

    payload = {
        "run": "mside_certificate",
        "oid": "msideK1",
        "method": METHOD_STRING,
        "sha256_pins": {"results.json": sha256(RESULTS_JSON)},
        "n_roots": len(records),
        "roots": records,
        "all_pass": bool(n_pass == len(records)),
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    print("=== SUMMARY: %d/%d pass, all_pass=%s (%.1fs) -> %s" % (
        n_pass, len(records), payload["all_pass"], time.time() - t0, OUT_JSON))
    return 0 if payload["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())

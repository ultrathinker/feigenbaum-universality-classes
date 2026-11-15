# -*- coding: utf-8 -*-
# [oid:kraw4E] krawczyk4_certificate.py -- strict Krawczyk boxes for the 61
# accepted superstable parameters of the quartic family f_a(x) = 1 - a*x^4
# (BRIEF-E-KRAWCZYK4; seeds: quartic_pilot/match_results.json, m=9, period
# P = 2^10 = 1024, F(a) := x_P(a), orbit of the critical point x=0).
#
# METHOD (per seed, all inclusions via Arb ball arithmetic -- python-flint,
# no float comparisons anywhere in the decision path):
#   x_0 = 0,  x_{k+1} = 1 - a*x_k^4;
#   F(a)  := x_n(a)  with n = 2^(m+1) = 1024,
#   F'(a) := dx_n(a), dx_{k+1} = -x_k^4 - 4*a*x_k^3*dx_k  (updated BEFORE x).
#   X := c0 + W, W = [+/- r]^2 (SQUARE, superset of the disk B(0,r); AUDIT_FIX
#   from D:\Islanded-1 ... mside AUDIT_FINDING.md: never a zero-width component);
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
#   REPORT-KRAWCZYK4.md for the statement; uniqueness holds on the square X),
#   of EXACT period n: return times of the critical orbit form p*Z, p | n;
#   (iv) excludes p in {2^j : j <= m} (x_1 == 1 excludes p = 1), hence p = n.
#
# PAIRWISE ISOLATION (after all boxes): squares X_i, X_j disjoint iff a
# component separation exceeds r_i + r_j (exact decimal arithmetic).
#
# PARALLELISM: ProcessPoolExecutor, max_workers=4 EXACTLY (BRIEF: CPU budget
# <= 4 processes/threads).  Incremental output: krawczyk4_certificates.json
# rewritten after EVERY root (is_complete=false until the end); progress.log
# line after every root, flushed.
import hashlib
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
PILOT = Path(r"../sec3_pilot_d4/match_results.json")
OUT_JSON = HERE / "krawczyk4_certificates.json"
PROGRESS = HERE / "progress.log"

MAX_WORKERS = 4               # hard cap (BRIEF: pool = 4)
RADII_LADDER = ("1e-26", "1e-27", "1e-25", "1e-24", "1e-28")
PREC_LADDER = (256, 512, 1024)

METHOD_STRING = (
    "Krawczyk + interval Newton in Arb ball arithmetic for f_a(x)=1-a*x^4: "
    "X=c0+[+/-r]^2 (square superset of B(c0,r)); F(a)=x_n(a), n=2^(m+1)=1024, "
    "x_0=0; N(X)=c0-F(c0)/F'(X); K(X)=c0-Y*F(c0)+(1-Y*F'(X))*(X-c0), "
    "Y=1/F'(c0); pass iff 0 notin F'(X), N(X) and K(X) strictly inside X, "
    "and inf_{a in X}|x_{2^j}(a)|>0 for j=0..m; implies existence and "
    "uniqueness of a superstable parameter of exact period n in B(c0,r)"
)

NSTR_RE = re.compile(r'^\(([-+0-9.eE]+)\s*([+-])\s*([-+0-9.eE]+)j\)$')


def parse_nstr(s):
    """mpmath nstr '(re + imj)' -> (re_str, im_str) exact decimal strings."""
    m = NSTR_RE.match(s.strip())
    if m is None:
        raise ValueError("unparseable nstr: %r" % s)
    return m.group(1), ("-" if m.group(2) == "-" else "") + m.group(3)


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
    """Seeds = accepted a_ss strings from match_results.json (frozen pilot).
    Period link taken from the pilot's own frozen fields, not from memory."""
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    fr = d["meta"]["frozen"]
    m_lvl = int(fr["M_MATCH"])
    n = 1 << (m_lvl + 1)
    if int(fr["PERIOD"]) != n:
        raise RuntimeError("pilot frozen PERIOD != 2^(M_MATCH+1)")
    if fr["family"] != "1 - a*x^4":
        raise RuntimeError("unexpected family %r" % fr["family"])
    seeds, rejected = [], []
    for rec in d["catalog_results"]:
        if rec.get("accepted"):
            re_s, im_s = parse_nstr(rec["a_ss"])
            seeds.append({"label": rec["id"], "source": "match_results:catalog_results[accepted]",
                          "n": n, "m": m_lvl, "re": re_s, "im": im_s,
                          "u_cat_re": rec["u_cat_re"], "u_cat_im": rec["u_cat_im"]})
        else:
            rejected.append(rec["id"])
    if len(seeds) != 61 or rejected != ["Q57"]:
        raise RuntimeError("expected 61 accepted seeds and rejected==[Q57], got %d/%r"
                           % (len(seeds), rejected))
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
        dz = -(z ** 4) - 4 * c0 * z ** 3 * dz
        z = 1 - c0 * z ** 4
    F0, Df0 = z, dz

    # ---- box orbit: D = F'(X); marks |x_{2^j}(X)| at j = 0..m
    zb = C(0)
    dzb = C(0)
    lowers = []
    j = 0
    for k in range(1, n + 1):
        dzb = -(zb ** 4) - 4 * X * zb ** 3 * dzb
        zb = 1 - X * zb ** 4
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


def certify_root(spec):
    """Ladder search: first (prec, radius) that passes; on total failure freeze
    the first configuration so the failure is fully reported and auditable."""
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
        chosen = attempt(spec, RADII_LADDER[0], PREC_LADDER[0])
    out = {
        "label": spec["label"], "source": spec["source"],
        "n": spec["n"], "m": spec["m"],
        "u_cat_re": spec["u_cat_re"], "u_cat_im": spec["u_cat_im"],
        "box_center_re": spec["re"], "box_center_im": spec["im"],
    }
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
    """NOT RIGOROUS (float constants a_F, A4, delta4 from the pilot are
    uncertified): scaled u-positions of the certified boxes vs catalog zeros.
    u(a) = A4*(a - a_F)*delta4^m with m = 9; box u-radius = r*A4*delta4^9."""
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    a_F = float(d["a_F"])
    A4 = float(d["A4"])
    d4 = float(d["delta4"])
    scale = A4 * d4 ** 9
    d_us, box_u_sizes = [], []
    for rec in records:
        u_re = (float(flint.arb(rec["box_center_re"]).mid()) - a_F) * scale
        u_im = float(flint.arb(rec["box_center_im"]).mid()) * scale
        du = ((u_re - rec["u_cat_re"]) ** 2 + (u_im - rec["u_cat_im"]) ** 2) ** 0.5
        d_us.append(du)
        box_u_sizes.append(float(flint.arb(rec["box_radius"]).mid()) * scale)
    return {
        "_not_rigorous": True,
        "note": "u-positions use uncertified float constants a_F/A4/delta4 "
                "from the pilot; the rigorous statement lives in a-space",
        "d_u_min": min(d_us), "d_u_max": max(d_us),
        "box_u_size_max": max(box_u_sizes),
    }


def main():
    t0 = time.time()
    try:                                   # Cyrillic paths on cp1251/cp1252 consoles
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== krawczyk4_certificate.py [oid:kraw4E]: Krawczyk boxes for 61 "
          "quartic superstable parameters (f_a(x)=1-a*x^4, n=1024) ===")
    seeds = load_seeds()
    print("loaded %d seeds from match_results.json (rejected: Q57)" % len(seeds))
    log_progress("STEP generator start: %d seeds, ladders radii=%s prec=%s, pool=%d"
                 % (len(seeds), RADII_LADDER, PREC_LADDER, MAX_WORKERS))

    records = []

    def flush_cert():
        payload = {
            "run": "krawczyk4_certificate", "oid": "kraw4E",
            "method": METHOD_STRING,
            "sha256_pins": {"match_results.json": sha256(PILOT)},
            "family": "1 - a*x^4",
            "frozen": {"radii_ladder": list(RADII_LADDER),
                       "prec_ladder": list(PREC_LADDER),
                       "max_workers": MAX_WORKERS},
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
            line = ("PARAM %2d/%d id=%s P=%d r<=%s kappa<=%s pass=%s "
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
        "run": "krawczyk4_certificate", "oid": "kraw4E",
        "method": METHOD_STRING,
        "sha256_pins": {"match_results.json": sha256(PILOT)},
        "family": "1 - a*x^4",
        "frozen": {"radii_ladder": list(RADII_LADDER),
                   "prec_ladder": list(PREC_LADDER),
                   "max_workers": MAX_WORKERS},
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

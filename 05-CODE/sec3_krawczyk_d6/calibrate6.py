# -*- coding: utf-8 -*-
# [oid:kraw6H] calibrate6.py -- PROTOCOL STEP (brief H, pt. 3): calibrate the
# radius/precision ladders on the FIRST 10 boxes BEFORE freezing them for the
# full 287-box run.  NOT RIGOROUS (protocol exploration): the choice of ladder
# is data-driven; every certified box afterwards is rigorous regardless of how
# the ladder was chosen.
#
# Same arithmetic core as the wave-1 audited generator (krawczyk4_certificate.py),
# adapted to the sextic family f_a(x) = 1 - a*x^6:
#   x_0 = 0,  x_{k+1} = 1 - a*x_k^6;
#   F(a)  := x_n(a),  n = 2^(m+1) = 1024,
#   F'(a) := dx_n(a), dx_{k+1} = -x_k^6 - 6*a*x_k^5*dx_k  (updated BEFORE x).
# Box X = c0 + [+/-r]^2 (square, superset of the disk B(c0,r)); pass rule:
# 0 notin F'(X), N(X) and K(X) strictly inside X, inf_X |x_{2^j}| > 0 (j=0..m).
#
# Grid: radii x prec below; seeds = first 10 of the 287 interior accepted rows
# of agentB_deg6/match_results.json (catalog order S001..).  Output:
# calibration6.json + progress.log lines.  Pool = 4 (CPU budget).
import json
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
PILOT = Path(r"../sec3_pilot_d6/match_results.json")
PROGRESS = HERE / "progress.log"
OUT = HERE / "calibration6.json"

MAX_WORKERS = 4
CAL_RADII = ("1e-21", "1e-22", "1e-23", "1e-24", "1e-25", "1e-26")
CAL_PREC = (256, 512)
N_CAL = 10

NSTR_RE = re.compile(r'^\(([-+0-9.eE]+)\s*([+-])\s*([-0-9.eE]+)j\)$')


def log_progress(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_seeds():
    """287 interior accepted rows (verdict set |u|<=44) of the deg6 blind run.
    Seed = exact decimal a_ss_re string + exact round-trip repr of the a_ss_im
    float (the blind run stores the imaginary part as float64 only)."""
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    fr = d["meta"]["frozen"]
    m_lvl = int(fr["M_MATCH"])
    n = 1 << (m_lvl + 1)
    if int(fr["PERIOD"]) != n:
        raise RuntimeError("pilot frozen PERIOD != 2^(M_MATCH+1)")
    if fr["family"] != "1 - a*x^6":
        raise RuntimeError("unexpected family %r" % fr["family"])
    seeds = []
    for rec in d["catalog_results"]:
        if rec.get("accepted") and not rec.get("edge"):
            seeds.append({"label": rec["id"],
                          "source": "match_results:catalog_results[accepted,interior]",
                          "n": n, "m": m_lvl,
                          "re": str(rec["a_ss_re"]), "im": repr(float(rec["a_ss_im"])),
                          "u_cat_re": rec["u_cat_re"], "u_cat_im": rec["u_cat_im"]})
    if len(seeds) != 287:
        raise RuntimeError("expected 287 interior accepted seeds, got %d" % len(seeds))
    return seeds


def attempt(spec, r_str, prec):
    """One certification attempt at (radius, precision). Core identical to the
    wave-1 generator (adapted to x^6); returns pass/fail diagnostics."""
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
    n = spec["n"]
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
    ok = bool(df_nonzero and newton_inside and k_inside and all(period_ok))

    def ub(x):
        return str(flint.arb(x).upper())

    def lb(x):
        return str(flint.arb(x).lower())

    return {
        "r": r_str, "prec": prec,
        "pass": ok,
        "f0_abs_ub": ub(F0.abs_upper()),
        "df_abs_lb": lb(df_lb),
        "newton_dist_ub": ub(n_dist),
        "kappa_ub": ub(kappa),
        "krawczyk_inside": k_inside, "newton_inside": newton_inside,
        "df_nonzero": df_nonzero, "period_all_ok": all(period_ok),
        "min_period_lb": lb(min(lowers[1:], key=lambda v: v.mid())),
    }


def job(args):
    spec, r_str, prec = args
    t0 = time.time()
    try:
        rec = attempt(spec, r_str, prec)
        rec["err"] = None
    except Exception as exc:              # e.g. degenerate box guard
        rec = {"r": r_str, "prec": prec, "pass": False, "err": "%s: %s"
               % (type(exc).__name__, exc)}
    rec["dt"] = time.time() - t0
    rec["label"] = spec["label"]
    return rec


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== calibrate6.py [oid:kraw6H]: ladder calibration on first %d "
          "boxes (NOT RIGOROUS protocol step) ===" % N_CAL)
    seeds = load_seeds()[:N_CAL]
    jobs = [(s, r, p) for s in seeds for p in CAL_PREC for r in CAL_RADII]
    print("grid: %d seeds x %d prec x %d radii = %d attempts, pool=%d"
          % (len(seeds), len(CAL_PREC), len(CAL_RADII), len(jobs), MAX_WORKERS))
    log_progress("STEP calibration start: %d attempts (first %d seeds, radii=%s "
                 "prec=%s), pool=%d"
                 % (len(jobs), len(seeds), CAL_RADII, CAL_PREC, MAX_WORKERS))
    recs = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as ex:
        for i, rec in enumerate(ex.map(job, jobs), start=1):
            recs.append(rec)
            if i % 12 == 0 or i == len(jobs):
                print("  %3d/%d (%.0fs)" % (i, len(jobs), time.time() - t0),
                      flush=True)

    # pivot: per (label, prec) the set of passing radii
    summary = {}
    for rec in recs:
        s = summary.setdefault(rec["label"], {})
        key = "prec%d" % rec["prec"]
        s.setdefault(key, []).append(rec)
    print("\n%-6s | %s" % ("label", "per-prec passing radii (kappa_ub)"))
    for lbl, s in summary.items():
        parts = []
        for key in sorted(s):
            okr = [(r["r"], r["kappa_ub"][:6]) for r in s[key] if r["pass"]]
            parts.append("%s: %s" % (key, okr if okr else "NONE"))
        print("%-6s | %s" % (lbl, " ; ".join(parts)))

    payload = {
        "run": "calibrate6", "oid": "kraw6H",
        "note": "NOT RIGOROUS protocol step (brief H pt.3): ladder calibration "
                "on the first 10 boxes; ladders frozen afterwards in "
                "krawczyk6_certificate.py",
        "grid": {"radii": list(CAL_RADII), "prec": list(CAL_PREC), "n_seeds": len(seeds)},
        "records": recs,
        "elapsed_sec": time.time() - t0,
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    log_progress("STEP calibration done: results in calibration6.json (%.1fs)"
                 % (time.time() - t0))
    print("wrote %s (%.1fs)" % (OUT, time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())

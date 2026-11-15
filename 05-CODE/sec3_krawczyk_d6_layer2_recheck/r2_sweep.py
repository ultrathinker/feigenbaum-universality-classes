# -*- coding: utf-8 -*-
# R2 sweep: replicate the Krawczyk certificate from
# krawczyk6_certificate_local.py (an exact copy of wave2_krawczyk6/
# krawczyk6_certificate.py with paths re-pointed at this folder so the
# import works offline) for the two layer-2 seeds S078/S090 and three
# layer-1 controls S001/S172/S180. The brief requires a sweep over
# prec in {256, 512, 1024} x radius in {1e-30, 1e-31, 1e-32, 1e-33, 1e-34}
# x pow_order in {pow, mul}. F, the orbit, and the PASS criteria are NOT
# modified; we only add a parallel `attempt_mul` that builds x^6 as x^5*x
# (the mathematically equivalent but differently rounded path used by the
# reviewer in R1_krawczyk_check.py).  All outputs go to r2_table.json.
import json, time, sys, os, math
import flint

# Import the original attempt() verbatim from the local copy
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from krawczyk6_certificate_local import attempt as attempt_pow   # x^6 = x**6 (author path)

PROGRESS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "progress.log")
TABLE_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "r2_table.json")


def log(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    with open(PROGRESS, "a", encoding="utf-8") as f:
        f.write(line + "\n")
    print(line, flush=True)


def fprime_ub_pow(spec, r_str, prec):
    """Recompute F'(X) over the box for the pow variant and return the
    upper bound on |F'(X)|. The original attempt() only stores the lower
    bound; for the relative-width of the magnitude ball we need both.
    Algorithmically identical to the box-orbit loop in attempt() (pow),
    but we throw away everything but D.abs_upper()."""
    flint.ctx.prec = prec
    A = flint.arb; C = flint.acb
    c0 = C(A(spec["re"]), A(spec["im"]))
    W = C(A("[+/- %s]" % r_str), A("[+/- %s]" % r_str))
    X = c0 + W
    n = spec["n"]
    zb = C(0); dzb = C(0)
    for _ in range(n):
        dzb = -(zb ** 6) - 6 * X * zb ** 5 * dzb
        zb = 1 - X * zb ** 6
    return float(dzb.abs_upper())


def attempt_pow_with_ub(spec, r_str, prec):
    """Wrapper around the author's attempt() that also computes df_abs_ub
    in the SAME variant for a consistent relative-width."""
    rec = attempt_pow(spec, r_str, prec)
    rec["df_abs_ub"] = fprime_ub_pow(spec, r_str, prec)
    return rec


def attempt_mul(spec, r_str, prec):
    """Mirror of attempt() but with x^6 = x^5*x instead of x**6.  Same
    definitions of F, F', N, K, PASS.  Kept tiny to make the structural
    difference vs the author's path obvious."""
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
            raise RuntimeError("degenerate box")
    n = spec["n"]; m = n.bit_length() - 2; half = 1 << m

    z = C(0); dz = C(0)
    for _ in range(n):
        x5 = z ** 5
        x6 = x5 * z
        dz = -x6 - 6 * c0 * x5 * dz
        z = 1 - c0 * x6
    F0, Df0 = z, dz

    zb = C(0); dzb = C(0); lowers = []; j = 0
    for k in range(1, n + 1):
        x5 = zb ** 5
        x6 = x5 * zb
        dzb = -x6 - 6 * X * x5 * dzb
        zb = 1 - X * x6
        if k <= half and k == (1 << j):
            lowers.append(zb.abs_lower()); j += 1
    D = dzb
    df_lb = D.abs_lower(); df_ub = D.abs_upper(); df_nonzero = bool(df_lb > 0)
    N = c0 - F0 / D
    n_dist = (N - c0).abs_upper(); newton_inside = bool(n_dist < r)
    Y = 1 / Df0; E = 1 - Y * D
    K = c0 - Y * F0 + E * W
    k_dist = (K - c0).abs_upper(); k_inside = bool(k_dist < r); kappa = k_dist / r
    period_ok = [bool(v > 0) for v in lowers]
    return {
        "box_radius": str(r),
        "prec_bits": prec,
        "df_abs_lb": float(df_lb),
        "df_abs_ub": float(df_ub),
        "df_nonzero": df_nonzero,
        "newton_inside": newton_inside,
        "krawczyk_inside": k_inside,
        "krawczyk_dist_ub": float(k_dist),
        "krawczyk_kappa_ub": float(kappa),
        "newton_dist_ub": float(n_dist),
        "period_all_ok": all(period_ok),
        "min_abs_lb_x2j": min(float(v) for v in lowers) if lowers else None,
        "_pass": bool(df_nonzero and newton_inside and k_inside and all(period_ok)),
        "_D": D,
    }


def ub(x):
    return float(flint.arb(x).upper())


def mid(x):
    return float(flint.arb(x).mid())


def sweep(spec, radii, precs, controls_name):
    rows = []
    n_attempts = 0
    t0 = time.time()
    for prec in precs:
        for r_str in radii:
            for variant, fn in (("pow", attempt_pow_with_ub), ("mul", attempt_mul)):
                n_attempts += 1
                t1 = time.time()
                rec = fn(spec, r_str, prec)
                dt = time.time() - t1
                # Normalize: pow path returns stringy flint-str fields;
                # mul path returns plain floats. Cast to arb so both go
                # through the same flt-path.
                def _arb(v):
                    if isinstance(v, str):
                        return flint.arb(v)
                    if isinstance(v, flint.arb):
                        return v
                    return flint.arb(float(v))
                kappa_arb = _arb(rec["krawczyk_kappa_ub"])
                k_dist_arb = _arb(rec["krawczyk_dist_ub"])
                n_dist_arb = _arb(rec["newton_dist_ub"])
                # df_abs_lb: lower bound of |F'(X)| over the box.
                df_abs_lb_arb = _arb(rec["df_abs_lb"])
                # Relative width of the |F'(X)| magnitude ball:
                # need both lb and ub. The pow path stores only lb in
                # `df_abs_lb`; for the relative-width column we recover
                # an upper bound from a fresh computation over X in the
                # SAME variant if available, else default to the lb
                # (relwidth=0 then means the upper-bound computation was
                # not done in the same call; report it as 0 and flag).
                if "df_abs_ub" in rec:
                    df_abs_ub = float(rec["df_abs_ub"])
                else:
                    # pow path: re-evaluate F'(X) with the same variant
                    # to get the magnitude ub cheaply. Use the stored
                    # complex ball D if the call site exposed it; else
                    # approximate: relwidth = df_abs_lb_rad/df_abs_lb_mid
                    # where df_abs_lb_rad is the half-width of the lb
                    # string (effectively zero in the pow path -- the
                    # ball on the magnitude is dominated by its lb).
                    df_abs_ub = float(df_abs_lb_arb)  # conservative:
                                                       # relwidth reports
                                                       # a lower bound on
                                                       # width, not an UB
                df_abs_lb_mid = float(df_abs_lb_arb)
                relwidth = (df_abs_ub - df_abs_lb_mid) / df_abs_lb_mid if df_abs_lb_mid > 0 else float("inf")
                r_arb = flint.arb("[+/- %s]" % r_str)
                r_mid = mid(r_arb); r_rad = ub(r_arb) - r_mid
                # min |x_{2^j}|
                if "min_abs_lb_x2j" in rec:
                    min_lb = rec["min_abs_lb_x2j"]
                else:
                    lbs = [pc["abs_lb"] for pc in rec["period_checks"]]
                    min_lb = min(mid(flint.arb(s)) for s in lbs[1:])
                rows.append({
                    "label": spec["label"], "variant": variant, "prec_bits": prec,
                    "r_str": r_str, "r_mid": r_mid, "r_rad": r_rad,
                    "krawczyk_inside": rec["krawczyk_inside"],
                    "newton_inside": rec["newton_inside"],
                    "df_nonzero": rec["df_nonzero"], "period_all_ok": rec["period_all_ok"],
                    "kappa_ub": ub(kappa_arb),
                    "krawczyk_dist_ub": ub(k_dist_arb),
                    "newton_dist_ub": ub(n_dist_arb),
                    "df_abs_lb_mid": df_abs_lb_mid,
                    "df_abs_ub": df_abs_ub,
                    "relwidth_Fprime": relwidth,
                    "min_abs_lb_x2j": min_lb,
                    "pass": rec["_pass"],
                    "secs": round(dt, 3),
                })
                log("ATTEMPT %s variant=%s prec=%d r=%s -> pass=%s kappa=%.3e "
                    "relw=%.3e dt=%.2fs"
                    % (spec["label"], variant, prec, r_str, rec["_pass"],
                       ub(kappa_arb), relwidth, dt))
    log("SUBJECT %s done: %d attempts in %.1fs"
        % (controls_name, n_attempts, time.time() - t0))
    return rows


# Box centers taken FROM the wave-2 certificates (verbatim), so the
# sweep uses the EXACT same starting center as the author's pipeline:
CONTROLS = [
    # Layer-2 (refined centers from refined_seeds6.json via original certs)
    {"label": "S078", "re": "1.68326021471585695676789162994044437928369871845204035174003",
     "im": "0.0000000304688386034566582645639286555284827209706220879207270316018",
     "layer": 2, "n": 1024},
    {"label": "S090", "re": "1.68326021407733028941891203851869493069161775556088452095435",
     "im": "0.0000000333875456590396721462828842474765862804275526757952217564935",
     "layer": 2, "n": 1024},
    # Layer-1 controls (centers = exact source_string from certificates)
    {"label": "S001", "re": "1.68326019674506167709076643255",
     "im": "0.0", "layer": 1, "n": 1024},
    {"label": "S172", "re": "1.68326019238794326990144482838",
     "im": "5.054698896122019e-08", "layer": 1, "n": 1024},
    {"label": "S180", "re": "1.68326019225830984909817035003",
     "im": "5.091402932792374e-08", "layer": 1, "n": 1024},
]

R_LAYER2 = ("1e-30", "1e-31", "1e-32", "1e-33", "1e-34")
R_LAYER1 = ("1e-22", "1e-23", "1e-24", "1e-25", "1e-30", "1e-31")
PRECS = (256, 512, 1024)


def main():
    t0 = time.time()
    # Sanity time on one cheap box (S001, r=1e-24, prec=256, mul): only
    # to gauge wall-clock before fanning out the full grid. Brief: "prec
    # 1024 может быть медленным — измерить на одном боксе, потом решать".
    log("CALIBRATION timing one S001 r=1e-24 prec=1024 mul attempt")
    t1 = time.time()
    attempt_mul(CONTROLS[2], "1e-24", 1024)
    log("CALIBRATION S001 r=1e-24 prec=1024 mul: %.2fs" % (time.time() - t1))
    log("CALIBRATION timing one S090 r=1e-32 prec=512 mul attempt")
    t1 = time.time()
    attempt_mul(CONTROLS[1], "1e-32", 512)
    log("CALIBRATION S090 r=1e-32 prec=512 mul: %.2fs" % (time.time() - t1))

    table = {"meta": {"family": "1 - a*x^6", "n": 1024, "m": 9,
                      "precs": list(PRECS),
                      "r_layer2": list(R_LAYER2),
                      "r_layer1": list(R_LAYER1),
                      "methods": {"pow": "author: x**6, x**5",
                                  "mul": "x^5 then x^6 = x^5*x (equivalent path)"}},
             "rows": []}
    for ctrl in CONTROLS:
        radii = R_LAYER2 if ctrl["layer"] == 2 else R_LAYER1
        rs = sweep(ctrl, radii, PRECS, ctrl["label"])
        table["rows"].extend(rs)
    table["meta"]["total_secs"] = round(time.time() - t0, 2)
    with open(TABLE_OUT, "w", encoding="utf-8") as f:
        json.dump(table, f, ensure_ascii=False, indent=2)
    log("WROTE %s rows=%d total=%.1fs" % (TABLE_OUT, len(table["rows"]),
                                            table["meta"]["total_secs"]))


if __name__ == "__main__":
    main()
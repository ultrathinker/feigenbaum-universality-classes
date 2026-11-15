# -*- coding: utf-8 -*-
# [oid:kraw4E] run_tampers4.py -- deliberate-corruption battery for
# krawczyk4_certificates.json -> TAMPER_LOG.md.
#
# Principle "lie plausibly": every tamper keeps the artifact INTERNALLY
# consistent (correct formats, plausible magnitudes, self-consistent flags)
# so that only RECOMPUTATION or source-pin reconciliation can catch it.
#
# Two engines:
#   A (rebuild)  -- krawczyk4_certificate.attempt() with a corrupted input
#                   (shifted seed, tampered family, wrong period, fake
#                   parameter, inflated radius); caught = NO ladder
#                   configuration certifies (pass-any is False).
#   B (verify)   -- mutated certificate JSON; caught = krawczyk4_verify.py
#                   exits 1.
#
# Tampers (the six BRIEF-mandated ones first):
#   T1  seed_shift          A: seed a_ss Re + 1e-20 (root leaves the box ->
#                           Krawczyk fails at every radius)
#                           B: Q05 box_center_re last decimal digit +1
#                           (shift ~1e-30, box still contains the true root;
#                           only the exact source pin catches it)
#   T2  family_coeff        A: orbit of g_a(x)=1-(1+1e-9)*a*x^4 on Q01,Q28,
#                           Q39,Q62 (seed is not a root of the tampered
#                           family)
#   T3  wrong_period        A: n=512 (down: x_512(c0) ~ 1e-2, Newton/Krawczyk
#                           distance ~ 1e-7 >> r) and n=2048 (up: the box
#                           contains the period-1024 root, so the exact-
#                           period check at d=1024 gives lower bound 0)
#                           B: Q19 cert row downgraded n->512, m->8, period
#                           checks trimmed (self-consistent)
#   T4  inflated_radius     A: Q51/Q52 forced r=1e-8 (nearest pair: boxes
#                           overlap -> pairwise isolation violated)
#                           B: all radii -> "1e-8" (off-ladder AND overlapping)
#   T5  fake_parameter      A: rebuild at fake seeds: Q57 predicted position
#                           (the pilot's rejected zero) and CTRL1/CTRL2
#                           control points (u -> a via the pilot constants);
#                           none is a period-1024 superstable root
#   T6  box_swap            B: Q02 <-> Q03 centers swapped
#   T7  pin_zeroed          B: sha256 pin -> 64 zeros
#   T8  kappa_forged        B: Q06 krawczyk_kappa_ub * 0.01 (forged tighter)
#   T9  period_check_dropped B: Q61 j=9 (d=512) check removed
#   T10 pass_flag_flip      B: Q30 pass -> false while all_pass stays true
import copy
import json
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
CERT = HERE / "krawczyk4_certificates.json"
VERIFIER = HERE / "krawczyk4_verify.py"
TMP = HERE / "_tampered_krawczyk4.json"
LOG = HERE / "TAMPER_LOG.md"
PILOT = Path(r"../sec3_pilot_d4/match_results.json")

sys.path.insert(0, str(HERE))
import krawczyk4_certificate as GEN  # Engine A: reuse attempt() only

mp.mp.dps = 50


# ---------------------------------------------------------------- engine A
def ladder_pass_any(spec, family_eps=0, n_override=None):
    """Walk the FULL frozen ladder; return (pass_any, first_fail_reason)."""
    reason = None
    for prec in GEN.PREC_LADDER:
        for r_str in GEN.RADII_LADDER:
            try:
                rec = GEN.attempt(spec, r_str, prec,
                                  family_eps=family_eps, n_override=n_override)
            except Exception as exc:
                return False, "exception: %s" % exc
            if rec["_pass"]:
                return True, "certified at r=%s prec=%d (!)" % (r_str, prec)
            fails = []
            if not rec["df_nonzero"]:
                fails.append("0 in F'(X)")
            if not rec["newton_inside"]:
                fails.append("N(X) not inside X")
            if not rec["krawczyk_inside"]:
                fails.append("K(X) not inside X")
            if not rec["period_all_ok"]:
                fails.append("exact-period lower bound = 0")
            if reason is None:
                reason = "; ".join(fails) + " (r=%s)" % r_str
    return False, reason


def shifted_spec(label, shift="1e-20"):
    seeds = {s["label"]: s for s in GEN.load_seeds()}
    s = seeds[label]
    re_shifted = mp.nstr(mp.mpf(s["re"]) + mp.mpf(shift), 30)
    spec = dict(s)
    spec["re"] = re_shifted
    return spec


def fake_spec_from_u(label, u_re, u_im):
    d = json.loads(PILOT.read_text(encoding="utf-8"))
    a = (mp.mpf(d["a_F"])
         + mp.mpc(mp.mpf(u_re), mp.mpf(u_im))
         / (mp.mpf(d["A4"]) * mp.mpf(d["delta4"]) ** 9))
    a_s = mp.nstr(a, 35)
    s_re, s_im = GEN.parse_nstr(a_s)
    return {"label": label, "source": "fake", "n": 1024, "m": 9,
            "re": s_re, "im": s_im, "u_cat_re": u_re, "u_cat_im": u_im}


def isolation_ok_two(re1, im1, r1, re2, im2, r2):
    import flint
    flint.ctx.prec = 512
    A = flint.arb
    dre = abs(A(re1) - A(re2)).lower()
    dim = abs(A(im1) - A(im2)).lower()
    sep = dre if bool(dre > dim) else dim
    rsum = (A(r1) + A(r2)).upper()
    return bool(sep > rsum)


# ---------------------------------------------------------------- engine B
def get_root(cert, label):
    return next(r for r in cert["roots"] if r["label"] == label)


def scale_bound_str(s, factor):
    import flint
    val = float(flint.arb(str(s)).mid())
    return "%.30e" % (val * factor)


def run_verifier(path):
    p = subprocess.run([sys.executable, "-X", "utf8", str(VERIFIER), str(path)],
                       capture_output=True, text=True, timeout=1800)
    fails = [ln.strip() for ln in p.stdout.splitlines() if "FAIL:" in ln]
    return p.returncode, fails


def make_tampered_B():
    base = json.loads(CERT.read_text(encoding="utf-8"))
    src = json.loads(PILOT.read_text(encoding="utf-8"))
    out = []

    c = copy.deepcopy(base)
    r = get_root(c, "Q05")
    r["box_center_re"] = r["box_center_re"][:-1] + str(
        (int(r["box_center_re"][-1]) + 1) % 10)
    out.append(("T1b seed_shift(pin)", c))

    c = copy.deepcopy(base)
    r = get_root(c, "Q19")
    r["n"] = 512
    r["m"] = 8
    r["period_checks"] = [pc for pc in r["period_checks"] if pc["j"] <= 8]
    out.append(("T3c period_downgrade", c))

    c = copy.deepcopy(base)
    for r in c["roots"]:
        r["box_radius"] = "1e-8"
        r["box_rad_re"] = "1e-8"
        r["box_rad_im"] = "1e-8"
    out.append(("T4b radius_inflate", c))

    c = copy.deepcopy(base)
    a, b = get_root(c, "Q02"), get_root(c, "Q03")
    a["box_center_re"], b["box_center_re"] = b["box_center_re"], a["box_center_re"]
    a["box_center_im"], b["box_center_im"] = b["box_center_im"], a["box_center_im"]
    out.append(("T6 box_swap", c))

    c = copy.deepcopy(base)
    c["sha256_pins"]["match_results.json"] = "0" * 64
    out.append(("T7 pin_zeroed", c))

    c = copy.deepcopy(base)
    r = get_root(c, "Q06")
    r["krawczyk_kappa_ub"] = scale_bound_str(r["krawczyk_kappa_ub"], 0.01)
    out.append(("T8 kappa_forged", c))

    c = copy.deepcopy(base)
    r = get_root(c, "Q61")
    r["period_checks"].pop()
    out.append(("T9 period_check_dropped", c))

    c = copy.deepcopy(base)
    get_root(c, "Q30")["pass"] = False
    out.append(("T10 pass_flag_flip", c))

    return out


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== run_tampers4.py [oid:kraw4E]: corruption battery ===")

    seeds = {s["label"]: s for s in GEN.load_seeds()}
    src = json.loads(PILOT.read_text(encoding="utf-8"))
    rows = []          # (tag, engine, desc, caught, detail)

    # ---- positive controls (must NOT fail) ------------------------------
    print("[positive controls]")
    ok, why = ladder_pass_any(seeds["Q01"])
    rows.append(("C1 clean_rebuild", "A",
                 "uncorrupted Q01 rebuilds and certifies on the ladder",
                 ok, why if ok else "EXPECTED cert, got none"))
    print("  C1 clean rebuild pass_any=%s" % ok)
    rc, _ = run_verifier(CERT)
    print("  C2 clean verifier exit=%d" % rc)
    rows.append(("C2 clean_verify", "B", "verifier exit 0 on the intact "
                 "certificate", rc == 0, "exit=%d" % rc))

    # ---- T1a seed shift (engine A) ---------------------------------------
    print("[T1a] seed shift +1e-20 (rebuild)")
    ok, why = ladder_pass_any(shifted_spec("Q01"))
    rows.append(("T1a seed_shift(rebuild)", "A",
                 "seed a_ss Re + 1e-20: true root leaves every box, "
                 "K(X) not inside X at all radii", not ok,
                 "certified(!)" if ok else why))
    print("  pass_any=%s | %s" % (ok, why))

    # ---- T2 family coefficient (engine A) ---------------------------------
    print("[T2] family coefficient x(1+1e-9) (rebuild)")
    allok = True
    detail = []
    for lbl in ("Q01", "Q28", "Q39", "Q62"):
        ok, why = ladder_pass_any(seeds[lbl], family_eps=1e-9)
        allok = allok and (not ok)
        detail.append("%s:%s" % (lbl, "fail" if not ok else "CERT(!)"))
    rows.append(("T2 family_coeff(rebuild)", "A",
                 "orbit of g_a(x)=1-(1+1e-9)a*x^4 at 4 seeds: no ladder "
                 "configuration certifies", allok, "; ".join(detail)))
    print("  %s" % "; ".join(detail))

    # ---- T3 wrong period (engine A) ---------------------------------------
    print("[T3a/T3b] wrong period 512 / 2048 (rebuild)")
    ok512, why512 = ladder_pass_any(seeds["Q01"], n_override=512)
    ok2048, why2048 = ladder_pass_any(seeds["Q01"], n_override=2048)
    rows.append(("T3a period_512(rebuild)", "A",
                 "claim period 512 for a period-1024 root: x_512(c0)~1e-2, "
                 "Newton/Krawczyk distance >> r", not ok512,
                 "certified(!)" if ok512 else why512))
    rows.append(("T3b period_2048(rebuild)", "A",
                 "claim period 2048: the box contains the period-1024 root, "
                 "exact-period lower bound at d=1024 collapses to 0",
                 not ok2048, "certified(!)" if ok2048 else why2048))
    print("  n=512 pass_any=%s | %s" % (ok512, why512))
    print("  n=2048 pass_any=%s | %s" % (ok2048, why2048))

    # ---- T4 inflated radius (engine A: nearest pair isolation) ------------
    print("[T4a] inflated radius 1e-8 on the nearest pair Q51/Q52 (rebuild)")
    import flint
    flint.ctx.prec = 512
    sep_isolated = isolation_ok_two(seeds["Q51"]["re"], seeds["Q51"]["im"],
                                    "1e-8", seeds["Q52"]["re"],
                                    seeds["Q52"]["im"], "1e-8")
    ok51, why51 = ladder_pass_any(seeds["Q51"], n_override=None)
    # force r=1e-8 attempts explicitly (off-ladder value):
    rec51 = GEN.attempt(seeds["Q51"], "1e-8", 512)
    rec52 = GEN.attempt(seeds["Q52"], "1e-8", 512)
    caught = (not sep_isolated) or (not rec51["_pass"]) or (not rec52["_pass"])
    detail = ("isolation(Q51,Q52)=%s, K-inside=%s/%s"
              % (sep_isolated, rec51["krawczyk_inside"], rec52["krawczyk_inside"]))
    rows.append(("T4a radius_1e-8(rebuild)", "A",
                 "nearest pair Q51/Q52 inflated to r=1e-8: boxes overlap "
                 "(pairwise isolation violated) and/or K(X) leaves X",
                 caught, detail))
    print("  %s" % detail)

    # ---- T5 fake parameters (engine A) ------------------------------------
    print("[T5] fake parameters: Q57 predicted, CTRL1, CTRL2 (rebuild)")
    q57 = next(r for r in src["catalog_results"] if r["id"] == "Q57")
    fakes = [("Q57pred", q57["u_cat_re"], q57["u_cat_im"])]
    for cr in src["controls"][:2]:
        fakes.append((cr["id"], cr["u_re"], cr["u_im"]))
    allok = True
    detail = []
    for lbl, ur, ui in fakes:
        spec = fake_spec_from_u(lbl, ur, ui)
        ok, why = ladder_pass_any(spec)
        allok = allok and (not ok)
        detail.append("%s:%s" % (lbl, "fail" if not ok else "CERT(!)"))
    rows.append(("T5 fake_parameters(rebuild)", "A",
                 "rebuild at the rejected Q57 predicted position and at "
                 "CTRL1/CTRL2 control points: none certifies", allok,
                 "; ".join(detail)))
    print("  %s" % "; ".join(detail))

    # ---- engine B ----------------------------------------------------------
    print("[engine B] mutated certificates -> independent verifier")
    for name, cert in make_tampered_B():
        TMP.write_text(json.dumps(cert, ensure_ascii=False, indent=2),
                       encoding="utf-8")
        rc, fails = run_verifier(TMP)
        first = fails[0] if fails else "(no FAIL line)"
        rows.append((name, "B", "mutated certificate JSON", rc == 1, first[:110]))
        print("  %-24s exit=%d caught=%s | %s" % (name, rc, rc == 1, first[:80]))
    if TMP.exists():
        TMP.unlink()

    # ---- log ----------------------------------------------------------------
    n_caught = sum(1 for r in rows if r[3] and not r[0].startswith("C"))
    n_tamp = sum(1 for r in rows if not r[0].startswith("C"))
    ctrl_ok = all(r[3] for r in rows if r[0].startswith("C"))

    lines = []
    lines.append("# TAMPER_LOG -- corruption battery for "
                 "krawczyk4_certificates.json [oid:kraw4E]")
    lines.append("")
    lines.append("Date: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("Subject: 61 Krawczyk certificates (quartic f_a(x)=1-a*x^4, "
                 "period 1024) produced by `krawczyk4_certificate.py` (Arb "
                 "acb ball path) and verified by the independent "
                 "`krawczyk4_verify.py` (hand-rolled (re,im,r) tuple path). "
                 "Principle: every tamper lies plausibly -- internally "
                 "consistent formats and magnitudes -- so only recomputation "
                 "or exact source-pin reconciliation can catch it.")
    lines.append("")
    lines.append("Engines: **A** = rebuild via the generator with corrupted "
                 "input (caught = no ladder configuration certifies); "
                 "**B** = mutated certificate JSON (caught = verifier exit 1).")
    lines.append("")
    lines.append("Controls (must pass): %s" % (
        "OK" if ctrl_ok else "**FAILED**"))
    for r in rows:
        if r[0].startswith("C"):
            lines.append("- %s: %s -- %s" % (r[0], "OK" if r[3] else "FAIL", r[4]))
    lines.append("")
    lines.append("| # | tamper | engine | what was corrupted | caught | detail / first verifier complaint |")
    lines.append("|---|--------|--------|--------------------|--------|----------------------------------|")
    for tag, engine, desc, caught, detail in rows:
        if tag.startswith("C"):
            continue
        esc = str(detail).replace("|", "/").strip()
        lines.append("| %s | %s | %s | %s | %s | %s |" % (
            tag.split()[0], tag.split(None, 1)[1], engine, desc,
            "YES" if caught else "**NO**", esc))
    lines.append("")
    lines.append("Summary: **%d/%d tampers caught**; controls %s. %s" % (
        n_caught, n_tamp, "OK" if ctrl_ok else "FAILED",
        "Battery PASSED." if (n_caught == n_tamp and ctrl_ok)
        else "BATTERY FAILED -- a hole remains."))
    lines.append("")
    LOG.write_text("\n".join(lines), encoding="utf-8")
    print("=== SUMMARY: %d/%d caught, controls %s -> %s (%.1fs)"
          % (n_caught, n_tamp, "OK" if ctrl_ok else "FAILED", LOG,
             time.time() - t0))
    return 0 if (n_caught == n_tamp and ctrl_ok) else 1


if __name__ == "__main__":
    sys.exit(main())

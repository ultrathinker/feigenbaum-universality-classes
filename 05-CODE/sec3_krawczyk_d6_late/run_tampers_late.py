# -*- coding: utf-8 -*-
# [oid:kraw6L] run_tampers_late.py -- deliberate-corruption battery for
# krawczyk6_late_certificates.json -> TAMPER_LATE.md.
#
# Principle "lie plausibly" (same as wave-2): every tamper keeps the artifact
# INTERNALLY consistent so only RECOMPUTATION or source-pin reconciliation can
# catch it.  Two engines:
#   A (rebuild) -- krawczyk6_late_certificate.attempt() with a corrupted input
#                  (shifted center, swapped-in foreign seed, wrong period,
#                  inflated radius, flipped mantissa digit); caught = NO
#                  ladder configuration certifies.
#   B (verify)  -- mutated certificate JSON; caught = krawczyk6_late_verify.py
#                  exits 1.
#
# BRIEF-mandated tampers (>= 8):
#   T1  center_shift        A: refined center Re + 1e-30 (root leaves the box)
#                           B: S036 box_center_re last decimal digit +1
#   T2  foreign_seed        A: swap S031's center into S104 (the S104 box then
#                           does not contain its root)
#   T3  wrong_period        A: S031 at P=1024 (2^10, m=9) -- the true root of
#                           x_2048 is not a root of x_1024; S036 at P=4096
#                           (2^12, m=11) -- the exact-period check at d=4096
#                           collapses to 0
#                           B: S187 period_checks trimmed (self-consistent)
#   T4  radius_inflate      A: S104 radius 1e-8 (overlaps S187 + K leaves X)
#                           B: all radii -> 1e-8 (off-ladder AND overlapping)
#   T5  mantissa_digit_flip B: S235 box_center_im last decimal digit +1
#                           (~1e-70 shift, below box radius -- the exact
#                           string pin on refined_seeds_late.json catches it)
#   T6  pin_zeroed          B: sha256 pin -> 64 zeros
#   T7  period_mark_forged  B: S031 period_checks[j=3].abs_lb * 0.001 (forged
#                           tighter; factor-1.1 period band catches it)
#   T8  pass_flag_flip      B: S339 pass -> false while all_pass stays true
#   T9  period_check_dropped B: S279 j=10 (d=1024) check removed
#   T10 center_swap_new     B: S104 <-> S187 centers swapped (nearest pair)
#   T11 u_cat_drift         B: S031 u_cat_re drifted (source pin reconciliation)
#   T12 wrong_ladder        B: frozen radii_ladder edited (header drift)
import copy
import json
import subprocess
import sys
import time
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
CERT = HERE / "krawczyk6_late_certificates.json"
VERIFIER = HERE / "krawczyk6_late_verify.py"
TMP = HERE / "_tampered_late.json"
LOG = HERE / "TAMPER_LATE.md"

sys.path.insert(0, str(HERE))
import krawczyk6_late_certificate as GEN  # Engine A: reuse attempt() only


# ---------------------------------------------------------------- engine A
def ladder_pass_any(spec, n_override=None):
    reason = None
    for prec in GEN.PREC_LADDER:
        for r_str in GEN.RADII_LADDER:
            try:
                rec = GEN.attempt(spec, r_str, prec, n_override=n_override)
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


def load_specs_by_label():
    return {s["label"]: s for s in GEN.load_specs()}


def shifted_spec(label, shift="1e-30", shift_im=False):
    s = dict(load_specs_by_label()[label])
    import mpmath as mpm
    mpm.mp.dps = 80
    if shift_im:
        s["im"] = mpm.nstr(mpm.mpf(s["im"]) + mpm.mpf(shift), 40)
    else:
        s["re"] = mpm.nstr(mpm.mpf(s["re"]) + mpm.mpf(shift), 40)
    return s


def foreign_spec(label, foreign_from):
    s = dict(load_specs_by_label()[label])
    f = load_specs_by_label()[foreign_from]
    s["re"], s["im"] = f["re"], f["im"]
    s["note"] = "center of %s injected into %s" % (foreign_from, label)
    return s


def flip_mantissa_spec(label, last_digit):
    s = dict(load_specs_by_label()[label])
    s["re"] = s["re"][:-1] + str((int(s["re"][-1]) + last_digit) % 10)
    return s


# ---------------------------------------------------------------- engine B
def get_root(cert, label):
    return next(r for r in cert["roots"] if r["label"] == label)


def scale_bound_str(s, factor):
    val = float(flint.arb(str(s)).mid())
    return "%.30e" % (val * factor)


def run_verifier(path):
    p = subprocess.run([sys.executable, "-X", "utf8", str(VERIFIER), str(path)],
                       capture_output=True, text=True, timeout=1800)
    fails = [ln.strip() for ln in p.stdout.splitlines() if "FAIL:" in ln]
    return p.returncode, fails


def make_tampered_B():
    base = json.loads(CERT.read_text(encoding="utf-8"))
    out = []

    c = copy.deepcopy(base)
    r = get_root(c, "S036")
    r["box_center_re"] = r["box_center_re"][:-1] + str(
        (int(r["box_center_re"][-1]) + 1) % 10)
    out.append(("T1b center_shift(pin)", c))

    c = copy.deepcopy(base)
    r = get_root(c, "S187")
    r["period_checks"] = [pc for pc in r["period_checks"] if pc["j"] <= 10]
    out.append(("T3c period_checks_trimmed", c))

    c = copy.deepcopy(base)
    for r in c["roots"]:
        r["box_radius"] = "1e-8"
        r["box_rad_re"] = "1e-8"
        r["box_rad_im"] = "1e-8"
    out.append(("T4b radius_inflate", c))

    c = copy.deepcopy(base)
    c["sha256_pins"]["refined_seeds_late.json"] = "0" * 64
    out.append(("T6 pin_zeroed", c))

    c = copy.deepcopy(base)
    r = get_root(c, "S031")
    r["period_checks"][3]["abs_lb"] = scale_bound_str(
        r["period_checks"][3]["abs_lb"], 0.001)
    out.append(("T7 period_mark_forged", c))

    c = copy.deepcopy(base)
    get_root(c, "S339")["pass"] = False
    out.append(("T8 pass_flag_flip", c))

    c = copy.deepcopy(base)
    get_root(c, "S279")["period_checks"].pop()
    out.append(("T9 period_check_dropped", c))

    c = copy.deepcopy(base)
    a, b = get_root(c, "S104"), get_root(c, "S187")
    a["box_center_re"], b["box_center_re"] = b["box_center_re"], a["box_center_re"]
    a["box_center_im"], b["box_center_im"] = b["box_center_im"], a["box_center_im"]
    out.append(("T10 center_swap_S104_S187", c))

    c = copy.deepcopy(base)
    r = get_root(c, "S031")
    r["u_cat_re"] = str(float(r["u_cat_re"]) + 1.0)
    out.append(("T11 u_cat_drift", c))

    c = copy.deepcopy(base)
    c["frozen"]["radii_ladder"] = ["1e-30", "1e-29", "1e-28"]
    out.append(("T12 ladder_forged", c))

    return out


def main():
    t0 = time.time()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    print("=== run_tampers_late.py [oid:kraw6L]: corruption battery ===")

    rows = []

    # ---- positive controls (must NOT fail) ------------------------------
    print("[positive controls]")
    ok, why = ladder_pass_any(load_specs_by_label()["S036"])
    rows.append(("C1 clean_rebuild", "A",
                 "uncorrupted S036 (P=8192, m=12) rebuilds and certifies on "
                 "the ladder", ok, why if ok else "EXPECTED cert, got none"))
    print("  C1 clean rebuild pass_any=%s" % ok)
    rc, _ = run_verifier(CERT)
    print("  C2 clean verifier exit=%d" % rc)
    rows.append(("C2 clean_verify", "B", "verifier exit 0 on the intact "
                 "certificate", rc == 0, "exit=%d" % rc))

    # ---- T1a center shift (engine A) --------------------------------------
    print("[T1a] center shift +1e-20 (rebuild; beyond the 1e-26 radius)")
    ok, why = ladder_pass_any(shifted_spec("S031", shift="1e-20"))
    rows.append(("T1a center_shift(rebuild)", "A",
                 "refined center Re + 1e-20 (>> box radius 1e-26): true root "
                 "leaves the box, K(X) not inside X at all radii", not ok,
                 "certified(!)" if ok else why))
    print("  pass_any=%s | %s" % (ok, why))

    # ---- T2 foreign seed (engine A) ----------------------------------------
    print("[T2] foreign seed swap S031 -> S104 (rebuild)")
    ok, why = ladder_pass_any(foreign_spec("S104", "S031"))
    rows.append(("T2 foreign_seed(rebuild)", "A",
                 "S104 certified with S031's center (distance ~2e-9 vs "
                 "radius 1e-26): root outside every box", not ok,
                 "certified(!)" if ok else why))
    print("  pass_any=%s | %s" % (ok, why))

    # ---- T3 wrong period (engine A) -----------------------------------------
    print("[T3a/T3b] wrong period 1024 (S031) / 4096 (S036) (rebuild)")
    ok1024, why1024 = ladder_pass_any(load_specs_by_label()["S031"],
                                      n_override=1024)
    ok4096, why4096 = ladder_pass_any(load_specs_by_label()["S036"],
                                      n_override=4096)
    rows.append(("T3a period_1024(rebuild)", "A",
                 "claim period 1024 for the S031 period-2048 root: "
                 "Newton/Krawczyk distance >> r", not ok1024,
                 "certified(!)" if ok1024 else why1024))
    rows.append(("T3b period_4096(rebuild)", "A",
                 "claim period 4096 for the S036 period-8192 root: the box "
                 "contains the period-8192 root, exact-period lower bound at "
                 "d=4096 collapses to 0", not ok4096,
                 "certified(!)" if ok4096 else why4096))
    print("  P=1024 pass_any=%s | %s" % (ok1024, why1024))
    print("  P=4096 pass_any=%s | %s" % (ok4096, why4096))

    # ---- T4a inflated radius (engine A: nearest new pair isolation) ---------
    print("[T4a] inflated radius 1e-8 on the nearest new pair S104/S187 "
          "(rebuild)")
    A = flint.arb
    flint.ctx.prec = 512
    s104 = load_specs_by_label()["S104"]
    s187 = load_specs_by_label()["S187"]
    dre = abs(A(s104["re"]) - A(s187["re"])).lower()
    dim = abs(A(s104["im"]) - A(s187["im"])).lower()
    sep = dre if bool(dre > dim) else dim
    sep_isolated = bool(sep > A("1e-8") + A("1e-8"))
    rec104 = GEN.attempt(s104, "1e-8", 512)
    rec187 = GEN.attempt(s187, "1e-8", 512)
    caught = (not sep_isolated) or (not rec104["_pass"]) or (not rec187["_pass"])
    detail = ("isolation(S104,S187)@1e-8=%s, K-inside=%s/%s"
              % (sep_isolated, rec104["krawczyk_inside"],
                 rec187["krawczyk_inside"]))
    rows.append(("T4a radius_1e-8(rebuild)", "A",
                 "nearest new pair S104/S187 (sep ~2e-10) inflated to r=1e-8: "
                 "boxes overlap and/or K(X) leaves X", caught, detail))
    print("  %s" % detail)

    # ---- T5 mantissa digit flip (engine B: exact string pin) ----------------
    print("[T5] mantissa digit flip on stored center (pin reconciliation)")
    c = copy.deepcopy(json.loads(CERT.read_text(encoding="utf-8")))
    r = get_root(c, "S235")
    r["box_center_im"] = r["box_center_im"][:-1] + str(
        (int(r["box_center_im"][-1]) + 1) % 10)
    TMP.write_text(json.dumps(c, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    rc, fails = run_verifier(TMP)
    first = fails[0] if fails else "(no FAIL line)"
    rows.append(("T5 mantissa_flip(pin)", "B",
                 "S235 box_center_im last decimal digit +1 (~1e-70 shift, "
                 "below box radius -- the EXACT string pin on refined_seeds_"
                 "late.json catches it)", rc == 1, first[:110]))
    print("  exit=%d caught=%s | %s" % (rc, rc == 1, first[:80]))

    # ---- engine B ----------------------------------------------------------
    print("[engine B] mutated certificates -> independent verifier")
    for name, cert in make_tampered_B():
        TMP.write_text(json.dumps(cert, ensure_ascii=False, indent=2),
                       encoding="utf-8")
        rc, fails = run_verifier(TMP)
        first = fails[0] if fails else "(no FAIL line)"
        rows.append((name, "B", "mutated certificate JSON", rc == 1,
                     first[:110]))
        print("  %-30s exit=%d caught=%s | %s" % (name, rc, rc == 1,
                                                  first[:80]))
    if TMP.exists():
        TMP.unlink()

    # ---- log ---------------------------------------------------------------
    n_caught = sum(1 for r in rows if r[3] and not r[0].startswith("C"))
    n_tamp = sum(1 for r in rows if not r[0].startswith("C"))
    ctrl_ok = all(r[3] for r in rows if r[0].startswith("C"))

    lines = []
    lines.append("# TAMPER_LATE -- corruption battery for "
                 "krawczyk6_late_certificates.json [oid:kraw6L]")
    lines.append("")
    lines.append("Date: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("Subject: 9 late-subfamily Krawczyk certificates (sextic "
                 "f_a(x)=1-a*x^6, periods 2048/4096/8192; refined seeds from "
                 "refine_late.py) produced by `krawczyk6_late_certificate.py` "
                 "(Arb acb ball path) and verified by the independent "
                 "`krawczyk6_late_verify.py` (hand-rolled (re,im,r) tuple "
                 "disk path). Principle: every tamper lies plausibly -- "
                 "internally consistent formats and magnitudes -- so only "
                 "recomputation or exact source-pin reconciliation can catch "
                 "it.")
    lines.append("")
    lines.append("Engines: **A** = rebuild via the generator with corrupted "
                 "input (caught = no ladder configuration certifies); "
                 "**B** = mutated certificate JSON (caught = verifier exit 1).")
    lines.append("")
    lines.append("Controls (must pass): %s" % (
        "OK" if ctrl_ok else "**FAILED**"))
    for r in rows:
        if r[0].startswith("C"):
            lines.append("- %s: %s -- %s" % (r[0], "OK" if r[3] else "FAIL",
                                             r[4]))
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

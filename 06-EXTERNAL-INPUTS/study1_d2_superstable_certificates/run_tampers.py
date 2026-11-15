# -*- coding: utf-8 -*-
# [oid:msideK1] run_tampers.py -- deliberate-corruption battery for
# mside_certificate.json -> TAMPER_LOG.md.
#
# Project principle "lie plausibly": every tamper keeps the certificate
# INTERNALLY consistent (correct formats, plausible magnitudes, self-consistent
# pass flags) so that only RECOMPUTATION from primary sources or source-pin
# reconciliation can catch it. Each tamper must make mside_verify.py exit 1;
# the clean control must exit 0.
#
# Tampers:
#  T1 center_last_digit    Z05 box_center_re, last decimal digit bumped
#                          (shift ~1e-30; box would still contain the true
#                          root -- only the exact source pin catches it)
#  T2 z72_legacy_center    Z72 center replaced by the LEGACY value (c_pred at
#                          the certified disk centre, the known rejected
#                          provenance) -- internally consistent otherwise
#  T3 newton_dist_scaled   Z17 newton_dist_ub * 0.9 (forged tighter)
#  T4 kappa_scaled         Z40 krawczyk_kappa_ub * 0.01 (forged tighter)
#  T5 period_check_dropped Z63: the j=10 (d=1024) period check removed
#  T6 period_downgrade     Z19: n 2048->1024, m->9, period checks trimmed to
#                          j<=8 (fully self-consistent downgrade)
#  T7 pin_zeroed           sha256 pin of results.json -> 64 zeros
#  T8 pass_flag_flip       Z57 pass -> false while all_pass stays true
#  T9 radius_tamper        Z66 box_radius -> "1e-29" (an in-ladder value, so
#                          structure checks pass; only recomputation of the
#                          radius-dependent kappa/contraction catches it)
#  T10 label_center_swap   Z68 <-> Z69 centers swapped
import copy
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE / "mside_certificate.json"
VERIFIER = HERE / "mside_verify.py"
TMP = HERE / "_tampered_mside_certificate.json"
RESULTS_JSON = Path(r"D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/"
                    r"article_v1/v2_work/unified_mandelbrot/results.json")
LOG = HERE / "TAMPER_LOG.md"


def get_root(cert, label):
    return next(r for r in cert["roots"] if r["label"] == label)


def scale_bound_str(s, factor):
    """Parse an Arb bound string, scale its midpoint by factor, emit a
    plausible plain-decimal bound string of the same kind."""
    import flint
    val = float(flint.arb(str(s)).mid())
    return "%.30e" % (val * factor)


def make_tampered():
    """Returns list of (name, description, mutated-cert-dict)."""
    base = json.loads(CERT.read_text(encoding="utf-8"))
    src = json.loads(RESULTS_JSON.read_text(encoding="utf-8"))
    z72_legacy = next(r for r in src["stage3_catalog"] if r["id"] == "Z72")["c_pred"]
    out = []

    # T1
    c = copy.deepcopy(base)
    r = get_root(c, "Z05")
    r["box_center_re"] = r["box_center_re"][:-1] + str(
        (int(r["box_center_re"][-1]) + 1) % 10)
    out.append(("T1 center_last_digit",
                "Z05 box_center_re last digit +1 (shift ~1e-30, still inside "
                "the certified box; caught only by the exact source pin)", c))

    # T2
    c = copy.deepcopy(base)
    r = get_root(c, "Z72")
    r["box_center_re"] = str(z72_legacy["re_str"])
    r["box_center_im"] = str(z72_legacy["im_str"])
    out.append(("T2 z72_legacy_center",
                "Z72 center replaced by the legacy certified-centre prediction "
                "c_pred (rejected provenance), other fields untouched", c))

    # T3
    c = copy.deepcopy(base)
    r = get_root(c, "Z17")
    r["newton_dist_ub"] = scale_bound_str(r["newton_dist_ub"], 0.9)
    out.append(("T3 newton_dist_scaled",
                "Z17 newton_dist_ub forged 0.9x (plausibly tighter)", c))

    # T4
    c = copy.deepcopy(base)
    r = get_root(c, "Z40")
    r["krawczyk_kappa_ub"] = scale_bound_str(r["krawczyk_kappa_ub"], 0.01)
    out.append(("T4 kappa_scaled",
                "Z40 krawczyk_kappa_ub forged 0.01x (plausibly tighter)", c))

    # T5
    c = copy.deepcopy(base)
    r = get_root(c, "Z63")
    r["period_checks"].pop()
    out.append(("T5 period_check_dropped",
                "Z63: j=10 (d=1024) period check removed", c))

    # T6
    c = copy.deepcopy(base)
    r = get_root(c, "Z19")
    r["n"] = 1024
    r["m"] = 9
    r["period_checks"] = [pc for pc in r["period_checks"] if pc["j"] <= 8]
    out.append(("T6 period_downgrade",
                "Z19: n 2048->1024, m 10->9, period checks trimmed to j<=8 "
                "(self-consistent downgrade)", c))

    # T7
    c = copy.deepcopy(base)
    c["sha256_pins"]["results.json"] = "0" * 64
    out.append(("T7 pin_zeroed",
                "sha256 pin of results.json zeroed", c))

    # T8
    c = copy.deepcopy(base)
    get_root(c, "Z57")["pass"] = False
    out.append(("T8 pass_flag_flip",
                "Z57 pass -> false while all_pass stays true", c))

    # T9
    c = copy.deepcopy(base)
    get_root(c, "Z66")["box_radius"] = "1e-29"
    out.append(("T9 radius_tamper",
                'Z66 box_radius -> "1e-29" (in-ladder value; only the '
                "recomputed radius-dependent fields can catch it)", c))

    # T10
    c = copy.deepcopy(base)
    a, b = get_root(c, "Z68"), get_root(c, "Z69")
    a["box_center_re"], b["box_center_re"] = b["box_center_re"], a["box_center_re"]
    a["box_center_im"], b["box_center_im"] = b["box_center_im"], a["box_center_im"]
    out.append(("T10 label_center_swap",
                "Z68 <-> Z69 box centers swapped", c))

    return out


def run_verifier(path):
    p = subprocess.run([sys.executable, str(VERIFIER), str(path)],
                       capture_output=True, text=True, timeout=1200)
    fails = [ln.strip() for ln in p.stdout.splitlines() if "FAIL:" in ln]
    return p.returncode, fails


def main():
    t0 = time.time()
    print("=== run_tampers.py [oid:msideK1]: corruption battery ===")

    print("[0] clean control")
    rc_clean, fails_clean = run_verifier(CERT)
    print("    exit=%d %s" % (rc_clean, "OK" if rc_clean == 0 else "UNEXPECTED"))

    rows = []
    n_caught = 0
    for name, desc, cert in make_tampered():
        TMP.write_text(json.dumps(cert, ensure_ascii=False, indent=2),
                       encoding="utf-8")
        rc, fails = run_verifier(TMP)
        caught = (rc == 1)
        n_caught += caught
        first = fails[0] if fails else "(no FAIL line)"
        first2 = fails[1] if len(fails) > 1 else ""
        print("  %-24s exit=%d caught=%s | %s" % (name, rc, caught, first[:90]))
        rows.append((name, desc, rc, caught, first, first2))
    if TMP.exists():
        TMP.unlink()

    lines = []
    lines.append("# TAMPER_LOG -- corruption battery for mside_certificate.json [oid:msideK1]")
    lines.append("")
    lines.append("Date: %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    lines.append("")
    lines.append("Subject: `mside_certificate.json` (23 Krawczyk/interval-Newton")
    lines.append("certificates of superstable parameters) verified by the independent")
    lines.append("tuple-path verifier `mside_verify.py`. Principle: every tamper is")
    lines.append("internally consistent (\"lie plausibly\") and must be caught by")
    lines.append("recomputation from primary sources or by exact source-pin reconciliation.")
    lines.append("")
    lines.append("Clean control: exit code **%d** (expected 0) %s" % (
        rc_clean, "-- OK" if rc_clean == 0 else "-- UNEXPECTED"))
    lines.append("")
    lines.append("| # | tamper | what was corrupted | exit | caught | first verifier complaint |")
    lines.append("|---|--------|--------------------|------|--------|--------------------------|")
    for name, desc, rc, caught, first, first2 in rows:
        esc = first.replace("|", "/").strip()
        lines.append("| %s | %s | %s | %d | %s | %s |" % (
            name.split()[0], name.split(None, 1)[1], desc, rc,
            "YES" if caught else "**NO**", esc))
    lines.append("")
    lines.append("Follow-up complaints (2nd FAIL line), for the record:")
    for name, desc, rc, caught, first, first2 in rows:
        if first2:
            lines.append("- %s: %s" % (name.split()[0], first2.replace("|", "/").strip()))
    lines.append("")
    lines.append("Summary: **%d/%d tampers caught**; clean control exit %d. %s" % (
        n_caught, len(rows), rc_clean,
        "Battery PASSED." if (n_caught == len(rows) and rc_clean == 0)
        else "BATTERY FAILED -- verifier has a hole."))
    lines.append("")
    LOG.write_text("\n".join(lines), encoding="utf-8")
    print("=== SUMMARY: %d/%d caught, clean exit %d -> %s (%.1fs)" % (
        n_caught, len(rows), rc_clean, LOG, time.time() - t0))
    return 0 if (n_caught == len(rows) and rc_clean == 0) else 1


if __name__ == "__main__":
    sys.exit(main())

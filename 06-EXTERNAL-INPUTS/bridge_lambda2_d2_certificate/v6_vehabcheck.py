# -*- coding: utf-8 -*-
"""
v6_vehabcheck.py -- independent re-verification of the Z72 milestone-B
conditional theorem (REPORT_VEHA_B.md, parallel agent).

Parts:
  A  own float Newton cascade m=8..13 for Z72 (my seed chain, my u_inf),
     own d_u(8), d_u(9); m=10 circularity probe vs archived box.
  B  own flint/Arb Krawczyk certification at m=11 (n=4096) and m=12
     (n=8192): own seed, own precision/radius ladders, own code; center
     and radius compared with sonnet_veha_b step6c/step6e.
  C  error-budget arithmetic recomputed independently at dps=200 (Delta11
     as an interval from MY certified box), two C-calibrations (theirs
     from agent-A d_u, mine own), sensitivities, eps table, break-even.
  D  deliberate corruption: shift of the m=11 anchor center by s*e^{i.theta},
     budget(s) at eps=1e-10 and eps=0, bisection threshold s*.
  E  two-layer (delta + lambda2) tail model from my own u(11..13) fit;
     effect on the final budget.

All external files are read-only inputs. Constants verbatim from
Investigation-1 certificates/calibrations.
"""
import json
import time
import sys
import mpmath as mp
import flint
from flint import arb, acb, ctx

sys.stdout.reconfigure(encoding="utf-8")

OUT = r"."
RES1 = r"../study1_d2_frozen_experiment/results.json"
ZCERT = r"../study1_d2_z72_localization/z72_zero_certificate.json"
MCERT = r"../study1_d2_superstable_certificates/mside_certificate.json"
RATE = r"../../05-CODE/sec5_lambda2_deg2/rate_deg2.json"
S6C = r"D:\Исследования\И-4 — Мост (реализация)\Исследование-4-Реализация\sonnet_veha_b\step6c_m11_krawczyk_result.json"
S6E = r"D:\Исследования\И-4 — Мост (реализация)\Исследование-4-Реализация\sonnet_veha_b\step6e_m12_result.json"

LAMBDA2_STR = "-0.1236527125526870243353936525855166103901"  # my v1/v3 certified midpoint (kilo_verify)

T0 = time.time()


def log(msg=""):
    print(f"[{time.time()-T0:7.1f}s] {msg}")
    sys.stdout.flush()


# ---------------------------------------------------------------- constants
meta = json.load(open(RES1, encoding="utf-8"))["meta"]
mp.mp.dps = 200
c_F = mp.mpmathify(meta["c_F_str"].replace("(", "").replace(")", "").replace("+ 0.0j", ""))
A = mp.mpmathify(meta["A_str"].replace("(", "").replace(")", "").replace("+ 0.0j", ""))
delta = mp.mpf(meta["delta_str"])
DELTA_LOWER_LIT = ("4.6692016091029906718532038204662016172581855774757686327456513"
                   "4300413355586854842639071238126635124200000000000000000000000000000"
                   "00000000000000000000000000000000")
delta_lower = mp.mpf(DELTA_LOWER_LIT) - mp.mpf("8.52e-154")  # theta_certificate point5_recompute

zc = json.load(open(ZCERT, encoding="utf-8"))["small_disk"]
u0 = mp.mpc(mp.mpf(zc["centre_re"]), mp.mpf(zc["centre_im"]))
r0 = mp.mpf(zc["r0"])

mc = json.load(open(MCERT, encoding="utf-8"))
z72_m10 = [rr for rr in mc["roots"] if rr["label"] == "Z72"][0]
c10_arch = mp.mpc(mp.mpf(z72_m10["box_center_re"]), mp.mpf(z72_m10["box_center_im"]))

log("inputs loaded")
log(f"  |c_F|={mp.nstr(abs(c_F),20)}  |A|={mp.nstr(abs(A),20)}")
log(f"  delta={mp.nstr(delta,30)}")
log(f"  delta_lower={mp.nstr(delta_lower,30)}  (delta-delta_lower={mp.nstr(delta-delta_lower,5)})")
log(f"  u0={mp.nstr(u0,20)}  r0={r0}")

# ================================================================ PART A
# my own float cascade m=8..13, seed chain: m=8 from c_pred(u0), then
# renormalization extrapolation upward (no u0 input beyond the rough m=8 seed)
mp.mp.dps = 150


def F_and_DF(c, n):
    z = mp.mpc(0)
    dz = mp.mpc(0)
    for _ in range(n):
        dz = 2 * z * dz + 1
        z = z * z + c
    return z, dz


def newton_super(c_init, n, tol=mp.mpf("1e-135"), itmax=60):
    c = c_init
    for it in range(itmax):
        z, dz = F_and_DF(c, n)
        step = z / dz
        # damping for the first steps if wild
        k = 1
        while abs(step) > mp.mpf("0.5") and k < 40:
            step = step / 2
            k *= 2
        c = c - step
        if abs(step) < tol:
            return c, it + 1, True
    return c, itmax, False


def exact_period_min(c, m):
    """min over j=0..m of |z_{2^j}(c)| (float)."""
    z = mp.mpc(0)
    j = 0
    vals = []
    for k in range(1, (1 << m) + 1):
        z = z * z + c
        if k == (1 << j):
            vals.append(abs(z))
            j += 1
    return min(vals), vals


cascade = {}
c_prev = None
c_seed8 = c_F + u0 / (A * delta ** 8)
for m in range(7, 14):
    n = 1 << (m + 1)
    if m == 7:
        seed = c_F + u0 / (A * delta ** 7)
    elif m == 8:
        seed = c_seed8
    else:
        seed = cascade[m - 1]["c"] + (cascade[m - 1]["c"] - cascade[m - 2]["c"]) / delta
    cstar, iters, conv = newton_super(seed, n)
    pmin, _ = exact_period_min(cstar, m)
    u_m = A * (cstar - c_F) * delta ** m
    cascade[m] = {"c": cstar, "seed_err": abs(cstar - seed), "iters": iters,
                  "conv": conv, "period_min": pmin, "u": u_m}
    log(f"A: m={m} n={n}: conv={conv} it={iters} seed_err={mp.nstr(abs(cstar-seed),4)} "
        f"period_min={mp.nstr(pmin,4)} |u-u0|={mp.nstr(abs(u_m-u0),12)}")

# seed-independence: perturb the m=7 c_pred seed, redo the WHOLE cascade
# (with full history, so every extrapolation uses converged points), compare
mp.mp.dps = 150


def run_chain(seed7_pert=None, m_hi=11):
    ch = {}
    for m in range(6, m_hi + 1):
        n = 1 << (m + 1)
        if m <= 7:
            s = c_F + u0 / (A * delta ** m)
            if m == 7 and seed7_pert is not None:
                s = s + seed7_pert
        else:
            s = ch[m - 1] + (ch[m - 1] - ch[m - 2]) / delta
        cstar, _, conv = newton_super(s, n)
        assert conv, f"chain failed at m={m}"
        ch[m] = cstar
    return ch


chain_ref = run_chain()
pert_small = run_chain(mp.mpc(mp.mpf("1.9e-8"), mp.mpf("-1.1e-8")))
pert_wide = run_chain(mp.mpc(mp.mpf("7e-5"), mp.mpf("4e-5")))
seed_indep_small = {m: abs(pert_small[m] - chain_ref[m]) for m in range(6, 12)}
seed_indep_wide = {m: abs(pert_wide[m] - chain_ref[m]) for m in range(6, 12)}
log("A: seed-independence, pert 2.2e-8 @ m=7 seed: "
    + ", ".join(f"m{m}:{mp.nstr(seed_indep_small[m], 3)}" for m in range(7, 12)))
log("A: seed-independence, pert 8.1e-5 @ m=7 seed (neighboring branch): "
    + ", ".join(f"m{m}:{mp.nstr(seed_indep_wide[m], 3)}" for m in range(7, 12)))

# my u_inf (Richardson from pair 12,13) and my d_u(8), d_u(9)
u_inf_mine = (delta * cascade[13]["u"] - cascade[12]["u"]) / (delta - 1)
d_u = {m: abs(cascade[m]["u"] - u_inf_mine) for m in (8, 9, 10, 11, 12, 13)}
log(f"A: u_inf_mine = {mp.nstr(u_inf_mine, 30)}")
for m in (8, 9, 11, 12):
    log(f"A: my d_u({m}) = {mp.nstr(d_u[m], 14)}")

# agent A quoted d_u (external, NOT RIGOROUS) for comparison
rate = json.load(open(RATE, encoding="utf-8"))["rows"]
row72 = [r for r in rate if r["id"] == "Z72_true"][0]
d_uA = {int(k): mp.mpf(v) for k, v in row72["d_u"].items()}
log(f"A: agent-A d_u(8)={mp.nstr(d_uA[8],10)} d_u(9)={mp.nstr(d_uA[9],10)} d_u(11)={mp.nstr(d_uA[11],10)}")
log(f"A: my-vs-agentA d_u8 rel diff = {mp.nstr(abs(d_u[8]-d_uA[8])/d_uA[8],6)}, "
    f"d_u9 rel = {mp.nstr(abs(d_u[9]-d_uA[9])/d_uA[9],6)}, d_u11 rel = {mp.nstr(abs(d_u[11]-d_uA[11])/d_uA[11],6)}")

# m=10 circularity probe (independent confirmation of their section 0)
mp.mp.dps = 100
c_pred10 = c_F + u0 / (A * delta ** 10)
circ_diff = abs(c_pred10 - c10_arch)
log(f"A: circularity probe |c_pred(10) - archived Z72 box center| = {mp.nstr(circ_diff, 6)}")

# ================================================================ PART B
# my own Krawczyk in flint; own code, own ladders


def kraw_attempt(re_str, im_str, r_str, prec, n, m):
    ctx.prec = prec
    c0 = acb(arb(re_str), arb(im_str))
    W = acb(arb("[+/- %s]" % r_str), arb("[+/- %s]" % r_str))
    r = W.real.rad().upper()
    if not bool(W.imag.rad().upper() <= r):
        r = W.imag.rad().upper()
    X = c0 + W
    # point orbit at c0
    z = acb(0)
    dz = acb(0)
    for _ in range(n):
        dz = 2 * z * dz + 1
        z = z * z + c0
    F0, Df0 = z, dz
    # box orbit with exact-period marks at k = 2^j, j=0..m
    zb = acb(0)
    dzb = acb(0)
    marks = []
    j = 0
    for k in range(1, n + 1):
        dzb = 2 * zb * dzb + 1
        zb = zb * zb + X
        if k <= (n >> 1) and k == (1 << j):
            marks.append(zb.abs_lower())
            j += 1
    D = dzb
    if len(marks) != m + 1:
        raise RuntimeError("mark miscount")
    df_lb = D.abs_lower()
    df_nonzero = bool(df_lb > 0)
    # interval Newton: N = c0 - F0/D(X)
    Nop = c0 - F0 / D
    n_dist = (Nop - c0).abs_upper()
    newton_inside = bool(n_dist < r)
    # Krawczyk: K = c0 - Y F0 + (1 - Y D) W,  Y = 1/Df0
    Y = 1 / Df0
    E = 1 - Y * D
    K = c0 - Y * F0 + E * W
    k_dist = (K - c0).abs_upper()
    k_inside = bool(k_dist < r)
    period_all_ok = all(bool(v > 0) for v in marks)
    passed = bool(df_nonzero and newton_inside and k_inside and period_all_ok)
    return {
        "prec_bits": prec, "r_str": r_str, "r_ub": str(r),
        "df_abs_lb": str(df_lb), "df_nonzero": df_nonzero,
        "newton_dist_ub": str(n_dist), "newton_inside": newton_inside,
        "krawczyk_dist_ub": str(k_dist), "krawczyk_inside": k_inside,
        "kappa": str(k_dist / r) if bool(r > 0) else "n/a",
        "period_lowers": [str(v) for v in marks],
        "period_all_ok": period_all_ok, "pass": passed,
    }


def certify(re_str, im_str, radii, precs, n, m, tag):
    attempts = []
    chosen = None
    for prec in precs:
        for rs in radii:
            rec = kraw_attempt(re_str, im_str, rs, prec, n, m)
            attempts.append({"prec": prec, "r_str": rs, "pass": rec["pass"],
                             "newton_inside": rec["newton_inside"],
                             "krawczyk_inside": rec["krawczyk_inside"],
                             "period_all_ok": rec["period_all_ok"],
                             "df_nonzero": rec["df_nonzero"]})
            log(f"B[{tag}]: prec={prec} r={rs}: pass={rec['pass']} "
                f"(N={rec['newton_inside']}, K={rec['krawczyk_inside']}, "
                f"P={rec['period_all_ok']}, df!=0={rec['df_nonzero']})")
            if rec["pass"] and chosen is None:
                chosen = rec
        if chosen is not None:
            break
    return chosen, attempts


mp.mp.dps = 160
R11 = mp.nstr(cascade[11]["c"].real, 130)
I11 = mp.nstr(cascade[11]["c"].imag, 130)
R12 = mp.nstr(cascade[12]["c"].real, 130)
I12 = mp.nstr(cascade[12]["c"].imag, 130)
log("B: my m=11 center strings (130 digits) prepared from my own cascade")

RADII_11 = tuple(f"1e-{p}" for p in range(33, 19, -1))  # tight -> loose
PRECS = (1024, 1536, 2048)
chosen11, attempts11 = certify(R11, I11, RADII_11, PRECS, 4096, 11, "m11")

RADII_12 = tuple(f"1e-{p}" for p in range(30, 15, -1))
chosen12, attempts12 = certify(R12, I12, RADII_12, PRECS, 8192, 12, "m12")

# center agreement with sonnet's certified boxes
s6c = json.load(open(S6C, encoding="utf-8"))
s6e = json.load(open(S6E, encoding="utf-8"))
c11_sonnet = mp.mpc(mp.mpf(s6c["box_center_re"]), mp.mpf(s6c["box_center_im"]))
c12_sonnet = mp.mpc(mp.mpf(s6e["c_final_re"]), mp.mpf(s6e["c_final_im"]))
mp.mp.dps = 120
d11_centers = abs(cascade[11]["c"] - c11_sonnet)
d12_centers = abs(cascade[12]["c"] - c12_sonnet)
log(f"B: |my c*(11) - sonnet c*(11)| = {mp.nstr(d11_centers, 8)}")
log(f"B: |my c*(12) - sonnet c*(12)| = {mp.nstr(d12_centers, 8)}")

# ================================================================ PART C
# budget arithmetic, dps=200, my order of operations
mp.mp.dps = 200
c11_pt = mp.mpc(mp.mpf(R11), mp.mpf(I11))
c12_pt = mp.mpc(mp.mpf(R12), mp.mpf(I12))

# interval u(11), u(12) from MY certified boxes (A, c_F, delta as points)
ctx.prec = 1024


def u_ball_from_box(re_str, im_str, r_str, m):
    c = acb(arb(re_str), arb(im_str)) + acb(arb("[+/- %s]" % r_str), arb("[+/- %s]" % r_str))
    Ab = acb(arb(mp.nstr(A, 150)), arb(0))
    cFb = acb(arb(mp.nstr(c_F, 150)), arb(0))
    db = acb(arb(mp.nstr(delta, 150)), arb(0))
    u = Ab * (c - cFb) * db ** m
    return u


r11_str = chosen11["r_str"] if chosen11 else "1e-27"
r12_str = chosen12["r_str"] if chosen12 else "1e-24"
u11_ball = u_ball_from_box(R11, I11, r11_str, 11)
u12_ball = u_ball_from_box(R12, I12, r12_str, 12)
u0_ball = acb(arb(zc["centre_re"]), arb(zc["centre_im"]))
D11_ball = (u11_ball - u0_ball).abs_upper()
D12_ball = (u12_ball - u0_ball).abs_upper()
log(f"C: Delta11 ball = {str(D11_ball)[:120]}")
log(f"C: Delta12 ball = {str(D12_ball)[:120]}")

mp.mp.dps = 200
Delta11_ub = mp.mpf(str(D11_ball).split(" +/-")[0].strip("[]"))
Delta12_ub = mp.mpf(str(D12_ball).split(" +/-")[0].strip("[]"))
Delta11_rad = mp.mpf(str(D11_ball).split("+/- ")[1].split("]")[0])

# C calibrations
C8A = d_uA[8] * delta ** 8
C9A = d_uA[9] * delta ** 9
C_geoA = mp.sqrt(C8A * C9A)
C8m = d_u[8] * delta ** 8
C9m = d_u[9] * delta ** 9
C_geom = mp.sqrt(C8m * C9m)
log(f"C: C8(agentA)={mp.nstr(C8A,12)} C9(agentA)={mp.nstr(C9A,12)} geomean={mp.nstr(C_geoA,12)}")
log(f"C: C8(mine)  ={mp.nstr(C8m,12)} C9(mine)  ={mp.nstr(C9m,12)} geomean={mp.nstr(C_geom,12)}")


def tail_bound(C, m0, dl=delta_lower):
    return C * dl ** (-m0) / (1 - 1 / dl)


dudA_11 = abs((c11_pt - c_F) * delta ** 11)
dudcF_11 = abs(A) * delta ** 11
dudA_12 = abs((c12_pt - c_F) * delta ** 12)
dudcF_12 = abs(A) * delta ** 12
log(f"C: m=11 |du/dA|={mp.nstr(dudA_11,12)} |du/dcF|={mp.nstr(dudcF_11,12)}")
log(f"C: m=12 |du/dA|={mp.nstr(dudA_12,12)} |du/dcF|={mp.nstr(dudcF_12,12)}")

eps_list = [mp.mpf("1e-20"), mp.mpf("5.238268502643473e-19"), mp.mpf("1e-10"),
            mp.mpf("1e-9"), mp.mpf("1e-6")]


def budget(Delta_m, tail_m, dA, dcF, eps):
    return Delta_m + tail_m + (dA + dcF) * eps


tableA, tableM = {}, {}
for eps in eps_list:
    bA = budget(Delta11_ub, tail_bound(C_geoA, 11), dudA_11, dudcF_11, eps)
    bM = budget(Delta11_ub, tail_bound(C_geom, 11), dudA_11, dudcF_11, eps)
    tableA[mp.nstr(eps, 4)] = bA
    tableM[mp.nstr(eps, 4)] = bM
    log(f"C: eps={mp.nstr(eps,4)}: budget(C_theirs)={mp.nstr(bA,10)} margin={mp.nstr(r0/bA,5)}x | "
        f"budget(C_mine)={mp.nstr(bM,10)} margin={mp.nstr(r0/bM,5)}x")

eps_breakA = (r0 - Delta11_ub - tail_bound(C_geoA, 11)) / (dudA_11 + dudcF_11)
eps_breakM = (r0 - Delta11_ub - tail_bound(C_geom, 11)) / (dudA_11 + dudcF_11)
log(f"C: eps break-even: C_theirs {mp.nstr(eps_breakA,8)}  C_mine {mp.nstr(eps_breakM,8)}")

# m=12 anchor cross-check + m=10 (circular) anchor for reference
b12_e10 = budget(Delta12_ub, tail_bound(C_geoA, 12), dudA_12, dudcF_12, mp.mpf("1e-10"))
dudA_10 = abs((mp.mpc(mp.mpf(z72_m10["box_center_re"]), mp.mpf(z72_m10["box_center_im"])) - c_F) * delta ** 10)
dudcF_10 = abs(A) * delta ** 10
c10_pt_k = mp.mpc(mp.mpf(z72_m10["box_center_re"]), mp.mpf(z72_m10["box_center_im"]))
Delta10_pt = abs(A * (c10_pt_k - c_F) * delta ** 10 - u0)
b10_e10 = Delta10_pt + tail_bound(C_geoA, 10) + (dudA_10 + dudcF_10) * mp.mpf("1e-10")
log(f"C: m=12 anchor budget @1e-10 = {mp.nstr(b12_e10, 10)} (margin {mp.nstr(r0/b12_e10,4)}x) "
    f"[their claim: 2.504e-2, 0.40x]")
log(f"C: m=10 (circular) anchor budget @1e-10 = {mp.nstr(b10_e10, 10)} (margin {mp.nstr(r0/b10_e10,4)}x) "
    f"[their claim: 1.151e-3, 8.69x]; Delta10(point, archived box)={mp.nstr(Delta10_pt,8)}")

# ================================================================ PART D
# corruption: shift the m=11 anchor center
theta = mp.mpf("2.348672")
e_shift = mp.e ** (mp.mpc(0, 1) * theta)
amp = abs(A) * delta ** 11  # |du/ds| bound


def budget_shift(s, eps, C=C_geoA):
    c_shift = c11_pt + s * e_shift
    Dshift = abs(A * (c_shift - c_F) * delta ** 11 - u0)
    return Dshift + tail_bound(C, 11) + (dudA_11 + dudcF_11) * eps


corruption = {"theta": str(theta), "amp_absA_delta11": mp.nstr(amp, 12), "rows": []}
for s in [mp.mpf("1e-20"), mp.mpf("1e-18"), mp.mpf("1e-16"), mp.mpf("1e-14"),
          mp.mpf("1e-12"), mp.mpf("1e-11"), mp.mpf("3.16e-11"), mp.mpf("1e-10"),
          mp.mpf("3.16e-10"), mp.mpf("1e-9"), mp.mpf("1e-8")]:
    b0 = budget_shift(s, mp.mpf(0))
    b10 = budget_shift(s, mp.mpf("1e-10"))
    row = {"s": mp.nstr(s, 4), "budget_eps0": mp.nstr(b0, 8), "closed_eps0": bool(b0 < r0),
           "budget_eps1e_10": mp.nstr(b10, 8), "closed_eps1e_10": bool(b10 < r0)}
    corruption["rows"].append(row)
    log(f"D: s={row['s']}: eps0 -> {row['budget_eps0']} ({'OK' if row['closed_eps0'] else 'FAIL'}); "
        f"eps=1e-10 -> {row['budget_eps1e_10']} ({'OK' if row['closed_eps1e_10'] else 'FAIL'})")


def bisect_threshold(eps):
    lo, hi = mp.mpf("1e-14"), mp.mpf("1e-7")
    f_lo = budget_shift(lo, eps) - r0
    f_hi = budget_shift(hi, eps) - r0
    if not (f_lo < 0 < f_hi):
        return None
    for _ in range(80):
        mid = (lo + hi) / 2
        if budget_shift(mid, eps) < r0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


s_star0 = bisect_threshold(mp.mpf(0))
s_star10 = bisect_threshold(mp.mpf("1e-10"))
log(f"D: closure-loss threshold s*: eps=0 -> {mp.nstr(s_star0,8)}, eps=1e-10 -> {mp.nstr(s_star10,8)}")
log(f"D: analytic guide (0.01-base)/(|A|delta^11): "
    f"{mp.nstr((r0-Delta11_ub-tail_bound(C_geoA,11))/amp, 8)} (eps=0), "
    f"{mp.nstr((r0-Delta11_ub-tail_bound(C_geoA,11)-(dudA_11+dudcF_11)*mp.mpf('1e-10'))/amp, 8)} (eps=1e-10)")

# ================================================================ PART E
# two-layer model u(m) = u_inf + C*delta^-m + D*lambda2^m, fit on 11,12,13
mp.mp.dps = 100
lam2 = mp.mpf(LAMBDA2_STR)
U = {m: cascade[m]["u"] for m in (11, 12, 13)}
# linear system for (u_inf, C, D)
Mrows = []
for m in (11, 12, 13):
    Mrows.append([mp.mpf(1), delta ** (-m), lam2 ** m])
rhs = [U[11], U[12], U[13]]
Min = mp.matrix(Mrows)
rhs_v = mp.matrix(rhs)
sol = mp.lu_solve(Min, rhs_v)
u_inf2, C2, D2 = sol[0, 0], sol[1, 0], sol[2, 0]
log(f"E: two-layer fit: u_inf={mp.nstr(u_inf2, 24)}")
log(f"E:   |C|={mp.nstr(abs(C2),10)}  |D|={mp.nstr(abs(D2),10)}")
# residuals on independent levels
for m in (8, 9, 10):
    pred = u_inf2 + C2 * delta ** (-m) + D2 * lam2 ** m
    log(f"E:   resid m={m}: {mp.nstr(abs(cascade[m]['u']-pred), 8)}")

layer1_11 = abs(C2) * delta_lower ** (-11) / (1 - 1 / delta_lower)
layer2_11 = abs(D2) * abs(lam2) ** 11 / (1 - abs(lam2))
tail_2layer = layer1_11 + layer2_11
tail_1layer = tail_bound(C_geoA, 11)
log(f"E: tail(11) one-layer  = {mp.nstr(tail_1layer, 10)}")
log(f"E: tail(11) two-layer  = {mp.nstr(tail_2layer, 10)} "
    f"(L1={mp.nstr(layer1_11,8)}, L2={mp.nstr(layer2_11,8)})")
for eps in (mp.mpf("1e-20"), mp.mpf("1e-10")):
    b2 = Delta11_ub + tail_2layer + (dudA_11 + dudcF_11) * eps
    log(f"E: budget two-layer eps={mp.nstr(eps,3)}: {mp.nstr(b2, 10)} margin={mp.nstr(r0/b2,5)}x")

# ================================================================ dump
out = {
    "partA": {
        "cascade": {str(m): {"c_re": mp.nstr(cascade[m]["c"].real, 60),
                             "c_im": mp.nstr(cascade[m]["c"].imag, 60),
                             "conv": cascade[m]["conv"],
                             "iters": cascade[m]["iters"],
                             "period_min": mp.nstr(cascade[m]["period_min"], 8),
                             "u": mp.nstr(cascade[m]["u"], 40),
                             "dist_to_u0": mp.nstr(abs(cascade[m]["u"] - u0), 20)}
                    for m in cascade},
        "seed_indep_small_pert22e_8": {str(m): mp.nstr(seed_indep_small[m], 6) for m in range(6, 12)},
        "seed_indep_wide_pert81e_5": {str(m): mp.nstr(seed_indep_wide[m], 6) for m in range(6, 12)},
        "u_inf_mine": mp.nstr(u_inf_mine, 40),
        "d_u_mine": {str(m): mp.nstr(d_u[m], 16) for m in d_u},
        "d_u_agentA": {str(m): mp.nstr(d_uA[m], 16) for m in d_uA},
        "circularity_probe_c10": mp.nstr(circ_diff, 8),
    },
    "partB": {
        "m11": {"chosen": chosen11, "attempts": attempts11,
                "center_diff_vs_sonnet": mp.nstr(d11_centers, 10)},
        "m12": {"chosen": chosen12, "attempts": attempts12,
                "center_diff_vs_sonnet": mp.nstr(d12_centers, 10)},
    },
    "partC": {
        "Delta11_ub_ball": str(D11_ball), "Delta12_ub_ball": str(D12_ball),
        "Delta11_rad": mp.nstr(Delta11_rad, 6),
        "C_theirs": {"C8": mp.nstr(C8A, 14), "C9": mp.nstr(C9A, 14), "geomean": mp.nstr(C_geoA, 14)},
        "C_mine": {"C8": mp.nstr(C8m, 14), "C9": mp.nstr(C9m, 14), "geomean": mp.nstr(C_geom, 14)},
        "tail11_C_theirs": mp.nstr(tail_bound(C_geoA, 11), 14),
        "tail11_C_mine": mp.nstr(tail_bound(C_geom, 11), 14),
        "dudA_11": mp.nstr(dudA_11, 14), "dudcF_11": mp.nstr(dudcF_11, 14),
        "dudA_12": mp.nstr(dudA_12, 14), "dudcF_12": mp.nstr(dudcF_12, 14),
        "budget_table_C_theirs": {k: mp.nstr(v, 12) for k, v in tableA.items()},
        "budget_table_C_mine": {k: mp.nstr(v, 12) for k, v in tableM.items()},
        "eps_break_even_C_theirs": mp.nstr(eps_breakA, 10),
        "eps_break_even_C_mine": mp.nstr(eps_breakM, 10),
        "budget_m12_eps1e_10": mp.nstr(b12_e10, 12),
        "budget_m10_circular_eps1e_10": mp.nstr(b10_e10, 12),
        "Delta10_point_from_archived_box": mp.nstr(Delta10_pt, 12),
    },
    "partD": {
        **corruption,
        "s_star_eps0": mp.nstr(s_star0, 8) if s_star0 else None,
        "s_star_eps1e_10": mp.nstr(s_star10, 8) if s_star10 else None,
    },
    "partE": {
        "u_inf_2layer": mp.nstr(u_inf2, 30), "C": mp.nstr(C2, 14), "D": mp.nstr(D2, 14),
        "resid_m8": mp.nstr(abs(cascade[8]["u"] - (u_inf2 + C2 * delta ** (-8) + D2 * lam2 ** 8)), 8),
        "resid_m9": mp.nstr(abs(cascade[9]["u"] - (u_inf2 + C2 * delta ** (-9) + D2 * lam2 ** 9)), 8),
        "resid_m10": mp.nstr(abs(cascade[10]["u"] - (u_inf2 + C2 * delta ** (-10) + D2 * lam2 ** 10)), 8),
        "tail11_1layer": mp.nstr(tail_1layer, 12),
        "tail11_2layer": mp.nstr(tail_2layer, 12),
        "layer2_share": mp.nstr(layer2_11 / tail_1layer, 8),
    },
}
with open(OUT + r"\v6_vehabcheck.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
log("wrote v6_vehabcheck.json")
print("DONE")

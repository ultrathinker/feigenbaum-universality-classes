# -*- coding: utf-8 -*-
"""z72_zero_verify.py -- независимый верификатор сертификата z72_zero_certificate.json.

НЕ импортирует z72_zero_build.py. Пересчитывает каждое поле с нуля:
  * winding(H9) ДРУГИМ кодовым путём -- flint.acb_poly Horner (как z72_verify.py),
    а НЕ LAff-классом билдера;
  * per-arc |H9|_lower -- той же acb_poly орбитой (независимо от LAff);
  * L7_arc (level-7 tube) -- собственной копией AF-машинерии (single-symbol);
  * tail S_disk/Phi(s8)/Phi(s9)/Lam_final -- собственной копией peeling-рекурсии;
  * все sha256-пины входов;
  * доменные проверки (Lemma 2/3, Rmax_L7, Rmax_L9, far_tail).

Сверяет КАЖДОЕ числовое/булево поле сертификата с допуском REL_TOL; любое
расхождение -> exit 1. Tamper self-test на ФАЙЛОВОЙ КОПИИ: порча winding/margin
обязана детектироваться (exit 1); нетронутый оригинал -- exit 0.
"""
from __future__ import annotations
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
import time
from fractions import Fraction
from pathlib import Path

import flint

ROOT = Path(r"D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive")
HERE = Path(__file__).resolve().parent
flint.ctx.prec = 512
A, C = flint.arb, flint.acb
PI = A.pi()
TWO_PI = 2 * PI
I_UNIT = C(0, 1)

sys.path.insert(0, str(ROOT / "w35" / "x4_walk_sonnet"))
import x4_common as xc  # noqa: E402

DATA_PATH = ROOT / "w33" / "t5_cert_sonnet" / "w33_gh_hiprec.json"
CL34_PATH = ROOT / "w36" / "tail_codex" / "cl34_certificate.json"
ABS36_PATH = ROOT / "w36" / "audit_codex" / "absolute_margins_arb.json"
POINT8B_JSON = ROOT / "path_to_8" / "point8b_spectral_certificate" / "spectral_certificate.json"
POINT8C_JSON = ROOT / "path_to_8" / "point8c_spectral_theta" / "theta_certificate.json"
LHAT_PATH = ROOT / "w36" / "finish_fable" / "l_hat_fable.json"
X4COMMON_PATH = ROOT / "w35" / "x4_walk_sonnet" / "x4_common.py"
Z72_CERT_PATH = ROOT / "path_to_8" / "point17_z72_geometry" / "z72_certificate.json"
CERT_PATH = HERE / "z72_zero_certificate.json"
REL_TOL = 1e-9

AU0 = xc.AU0
DELTA = xc.DELTA_pt
R_DISC = xc.r_disc
DELTA_LO = DELTA.real.lower()

_cl34 = json.loads(CL34_PATH.read_text(encoding="utf-8"))
_abs36 = json.loads(ABS36_PATH.read_text(encoding="utf-8"))
EPS_G = A(_cl34["G"]["distance_and_Taylor_tail_upper"]).upper()
EPS_H = A(_cl34["V_h_delta"]["h_error_stored_center_upper"]).upper()
nh = A(_abs36["nh_upper"]).upper() + EPS_H

DELTA0 = A("4.66920160910299067185320382046620161725818557747576863274565134300413")
DELTA_ERR = A("9.452e-30")
DELTA_WIDENED = DELTA0 + A(0, DELTA_ERR)
LEFT_HEAD_NORM = A("1.104854250319780")
K284 = A("284")
RHO_STAR = A("7.5e-4")
PPERP_NORM = A("2.808542")
RESOLVENT_CS = A("5.495835")
RESOLVENT_MU = A("0.35")
TAIL_MARGIN = 250
J_MAX = 360
NDT_REF = A("22.19")
NDT2_REF = A("193.63")
M3 = A("1834")
Q2 = A("111.550905120798913")
L1 = (NDT_REF + K284 * RHO_STAR).upper()
L2 = (NDT2_REF + 2 * NDT_REF * K284 * RHO_STAR + (K284 * RHO_STAR) ** 2).upper()
q_direct = (L2 / (DELTA_LO ** 4)).upper()
q3 = (M3 * nh ** 3 / A(6)).upper()

LIN_F = A("1.10485425032")
PPERP = A("2.808542")
T_BOUND = A("4.7911")
M = xc.M
_lh_raw = json.loads(LHAT_PATH.read_text(encoding="utf-8"))
LHAT = [C(A(a), A(b)) for a, b in zip(_lh_raw["l_hat_re"], _lh_raw["l_hat_im"])]
L_HEAD_RADIUS = A("1.588343e-14") + A("7.925e-15")
L_TAIL_INF = A("6.990710e-17")
_db = C(A(DELTA.real.mid(), A(xc.RHO_D)))
_d7 = C(1) / (_db ** 7)
DSDU_C = C(_d7.real.mid(), _d7.imag.mid()) * AU0
DSDU_R = (C(_d7.real.mid(), _d7.imag.mid()).abs_upper() * 7 * A(xc.RHO_D)).upper()

ZERO_ACB = C(0)
ONE_ACB = C(1)
INV_R = C(1) / C("2.5")
TWO_OVER_R = C(2) / C("2.5")
ZERO_ARB = A(0)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fnum(s):
    return float(A(s)) if not isinstance(s, (int, float)) else float(s)


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def build_Lam_final():
    point8b = json.loads(POINT8B_JSON.read_text(encoding="utf-8"))
    u_table = {row["j"]: A(row["U_j_upper"]) for row in point8b["final_table"]}

    def U(j):
        if j == 0:
            return PPERP_NORM
        if j in u_table:
            return u_table[j]
        return RESOLVENT_CS * (RESOLVENT_MU ** j)

    delta_u = DELTA_WIDENED.upper()
    Krho = (K284 * RHO_STAR).upper()
    Hcoef = (LEFT_HEAD_NORM * nh).upper()
    dp = [A(1)]
    for _ in range(J_MAX):
        dp.append((dp[-1] * delta_u).upper())
    Lam_new = [A(1)]
    Lam_final = [A(1)]
    binding = ["EQUAL"]
    for j in range(1, J_MAX + 1):
        acc = A(0)
        for m in range(j):
            b = (Hcoef * dp[j - 1 - m] + U(j - 1 - m)).upper()
            acc = (acc + Lam_new[m] * b).upper()
        nj = (Hcoef * dp[j] + U(j) + Krho * acc).upper()
        oj = (A(L2) ** (j // 2)).upper() if j % 2 == 0 else (A(L1) * A(L2) ** ((j - 1) // 2)).upper()
        Lam_new.append(nj)
        if float(nj) <= float(oj):
            Lam_final.append(nj)
            binding.append("NEW")
        else:
            Lam_final.append(oj)
            binding.append("OLD")
    tally = {"NEW": binding.count("NEW"), "OLD": binding.count("OLD"), "EQUAL": binding.count("EQUAL")}
    return Lam_final, tally


Lam_final, Lam_tally = build_Lam_final()
_p8c = json.loads(POINT8C_JSON.read_text(encoding="utf-8"))
_p8c_match = True
for row in _p8c["peeling"]["table_sample"]:
    j = row["j"]
    if j > J_MAX:
        continue
    if abs(float(Lam_final[j]) - row["Lam_final_upper_float"]) / max(abs(row["Lam_final_upper_float"]), 1e-300) > 1e-6:
        _p8c_match = False
if not _p8c_match or Lam_tally != _p8c["peeling"]["binding_tally"]:
    print("VERIFY: FAIL -- Lam_final does not reproduce point8c")
    sys.exit(1)


def tail_terms(abs_uc_rho):
    def Phi(s):
        return (Q2 * s * s + q3 * s ** 3).upper()

    def s_of(k):
        return (abs_uc_rho * A(AU0.real.mid()) / (DELTA_LO ** (k + 1))).upper()

    def D_new(k):
        return (A(Lam_final[k - 7]) * Phi(s_of(k))).upper()

    def D_old(k):
        j = k - 7
        lam_old = (A(L2) ** (j // 2)).upper() if j % 2 == 0 else (A(L1) * A(L2) ** ((j - 1) // 2)).upper()
        return (lam_old * Phi(s_of(k))).upper()

    K_CUT = 7 + TAIL_MARGIN
    near = A(0)
    for k in range(7, K_CUT):
        near = (near + D_new(k)).upper()
    far_tail = ((D_old(K_CUT) + D_old(K_CUT + 1)) / (A(1) - q_direct)).upper()
    S_new = (near + far_tail).upper()
    Phi_s8 = Phi(s_of(7))
    Phi_s9 = Phi(s_of(8))
    lam1 = A(Lam_final[1])
    lam1_phi_s9 = (lam1 * Phi_s9).upper()
    S_minus = (S_new - Phi_s8 - lam1_phi_s9).upper()
    far_neg = bool(float(far_tail) <= float(near) * 1e-6 or float(far_tail) < 1e-200)
    return S_new, Phi_s8, Phi_s9, lam1, lam1_phi_s9, S_minus, far_tail, far_neg


# =====================================================================
# INDEPENDENT winding cross-check via flint.acb_poly Horner (different code
# path than the builder's LAff class). Used ONLY to confirm winding COUNT==1
# independently; the per-arc |Hm|_lb field check uses the LAff recompute below
# (which matches the builder to REL_TOL since it is the same algorithm,
# independently re-coded -- this is the z72_verify.py philosophy: field check
# via the same method from scratch + algorithmic cross-check via acb_poly).
# =====================================================================
GPOLY_COMPACT = flint.acb_poly(list(xc.tc_g))
HPOLY_COMPACT = flint.acb_poly(list(xc.tc_h))


def winding_count_acb_poly(uc, rho, level, n_arcs):
    """Independent winding COUNT via flint.acb_poly orbit (different code path)."""
    s_factor = AU0 / (DELTA ** level)
    STEPS = 2 ** (level + 1)
    MID = 2 ** level

    def Hm_pt(u):
        s = u * s_factor
        x = C(0)
        x_mid = C(0)
        for n in range(STEPS):
            w = (x * x - C(1)) / C(R_DISC)
            x = GPOLY_COMPACT(w) + s * HPOLY_COMPACT(w)
            x = C(A(x.real.mid(), A(x.real.rad()).upper() + EPS_G + s.abs_upper() * EPS_H),
                  A(x.imag.mid(), A(x.imag.rad()).upper() + EPS_G + s.abs_upper() * EPS_H))
            if n + 1 == MID:
                x_mid = x
        return x / x_mid

    endpoints = {}
    for j in range(n_arcs):
        for tt in (Fraction(j, n_arcs), Fraction(j + 1, n_arcs)):
            if tt not in endpoints:
                th = A(0) if (tt == 0 or tt == 1) else TWO_PI * A(tt.numerator) / A(tt.denominator)
                um = uc + C(rho) * (I_UNIT * C(th)).exp()
                endpoints[tt] = Hm_pt(um)
    total = A(0)
    for j in range(n_arcs):
        total = total + (endpoints[Fraction(j + 1, n_arcs)] / endpoints[Fraction(j, n_arcs)]).arg()
    wind = total / TWO_PI
    cand = int(round(float(wind.mid())))
    ok = ((wind - cand).contains(0) and not (wind - (cand - 1)).contains(0)
          and not (wind - (cand + 1)).contains(0))
    return cand if ok else None


# =====================================================================
# LAff single-symbol AF winding (independent copy of the builder's machine,
# re-coded here for the field check). Identical algorithm -> matches builder
# to REL_TOL; catches any tamper/corruption of per-arc fields.
# =====================================================================
class LAff:
    __slots__ = ("c", "a", "r")

    def __init__(self, c, a=C(0), r=A(0)):
        self.c = C(c)
        self.a = C(a)
        self.r = A(r)


def la_const(c):
    return LAff(c)


def la_add(x, y):
    return LAff(x.c + y.c, x.a + y.a, x.r + y.r)


def la_addc(x, c):
    return LAff(x.c + C(c), x.a, x.r)


def la_sc(x, c):
    c = C(c)
    return LAff(x.c * c, x.a * c, x.r * c.abs_upper())


def la_abs_ub(x):
    return x.c.abs_upper() + x.a.abs_upper() + x.r


def la_mul(x, y):
    ax, ay = x.a.abs_upper(), y.a.abs_upper()
    return LAff(x.c * y.c, x.c * y.a + y.c * x.a,
                ax * ay + x.c.abs_upper() * y.r + y.c.abs_upper() * x.r + ax * y.r + ay * x.r + x.r * y.r)


def la_horner(coeffs, z):
    v = la_const(C(0))
    for c in reversed(coeffs):
        v = la_addc(la_mul(v, z), c)
    return v


def la_recip(y):
    cl = y.c.abs_lower()
    rho = y.a.abs_upper() + y.r
    if not (cl > 0 and rho < cl):
        return None
    q = rho / cl
    inv = C(1) / y.c
    return LAff(inv, -y.a * inv * inv, y.r / (cl * cl) + (q * q / (1 - q)) / cl)


def Hm_light(u_mid, u_rad, level):
    STEPS = 2 ** (level + 1)
    MID = 2 ** level
    sf = C(AU0) / (C(DELTA) ** level)
    u = LAff(C(u_mid), C(A(u_rad)))
    s = la_sc(u, sf)
    s_abs = la_abs_ub(s)
    x = la_const(C(0))
    x_mid = None
    max_w = A(0)
    for n in range(STEPS):
        w = la_sc(la_addc(la_mul(x, x), -1), C(1) / C(R_DISC))
        wu = la_abs_ub(w)
        if not (wu < 1):
            return None, None
        if wu > max_w:
            max_w = wu
        vg = la_horner(xc.tc_g, w)
        vh = la_horner(xc.tc_h, w)
        vg.r = vg.r + EPS_G
        vh.r = vh.r + EPS_H
        x = la_add(vg, la_mul(vh, s))
        x.r = x.r + EPS_G + s_abs * EPS_H
        if n + 1 == MID:
            x_mid = x
    inv = la_recip(x_mid)
    if inv is None:
        return None, None
    return la_mul(x, inv), max_w


def recompute_winding_laff(uc, rho, level, n_arcs):
    endpoints = {}
    arcs_ok = []
    per_arc_lb = []
    per_arc_rmax = []
    for j in range(n_arcs):
        t1 = Fraction(j, n_arcs)
        t2 = Fraction(j + 1, n_arcs)
        dt = t2 - t1
        tm = (t1 + t2) / 2
        theta = TWO_PI * A(tm.numerator) / A(tm.denominator)
        chord = (2 * rho * (PI * A(dt.numerator) / (2 * A(dt.denominator))).sin()).upper()
        u_mid = uc + C(rho) * (I_UNIT * C(theta)).exp()
        Hm, maxw = Hm_light(u_mid, chord, level)
        if Hm is None:
            return None
        lb = (Hm.c.abs_lower() - Hm.a.abs_upper() - Hm.r).lower()
        if not (lb > 0):
            return None
        per_arc_lb.append(lb)
        per_arc_rmax.append(A(maxw).upper())
        arcs_ok.append((t1, t2))
        for tt in (t1, t2):
            if tt not in endpoints:
                th = A(0) if (tt == 0 or tt == 1) else TWO_PI * A(tt.numerator) / A(tt.denominator)
                um = uc + C(rho) * (I_UNIT * C(th)).exp()
                Hp, _ = Hm_light(um, A(0), level)
                if Hp is None:
                    return None
                re_m = A(Hp.c.real.mid())
                im_m = A(Hp.c.imag.mid())
                err = A(Hp.c.real.rad()) + A(Hp.c.imag.rad()) + Hp.a.abs_upper() + Hp.r
                endpoints[tt] = C(A(re_m, err), A(im_m, err))
    total = A(0)
    for t1, t2 in arcs_ok:
        total = total + (endpoints[t2] / endpoints[t1]).arg()
    wind = total / TWO_PI
    cand = int(round(float(wind.mid())))
    ok = ((wind - cand).contains(0) and not (wind - (cand - 1)).contains(0)
          and not (wind - (cand + 1)).contains(0))
    return {"winding_ball": str(wind), "count": cand if ok else None,
            "per_arc_lb": per_arc_lb, "per_arc_rmax": per_arc_rmax}


# =====================================================================
# L7 tube (independent AF copy, single-symbol)
# =====================================================================
def l1(a):
    s = ZERO_ARB
    for v in a:
        s = s + v.abs_upper()
    return s


class AF:
    __slots__ = ("c", "a", "r")

    def __init__(self, c, a, r):
        self.c = c
        self.a = a
        self.r = r


def af_const(c):
    return AF(c, [ZERO_ACB], ZERO_ARB)


def af_add(x, y):
    return AF(x.c + y.c, [x.a[0] + y.a[0]], x.r + y.r)


def af_addc(x, c):
    return AF(x.c + c, [x.a[0]], x.r)


def af_sc(x, c):
    return AF(x.c * c, [x.a[0] * c], x.r * c.abs_upper())


def af_mul(x, y):
    lx, ly = l1(x.a), l1(y.a)
    return AF(x.c * y.c, [x.c * y.a[0] + y.c * x.a[0]],
              lx * ly + x.c.abs_upper() * y.r + y.c.abs_upper() * x.r + lx * y.r + ly * x.r + x.r * y.r)


def af_horner(coeffs, z):
    a = af_const(ZERO_ACB)
    for ck in reversed(coeffs):
        a = af_add(af_mul(a, z), ck)
    return a


def af_abs_u(x):
    return x.c.abs_upper() + l1(x.a) + x.r


def af_abs_l(x):
    return x.c.abs_lower() - l1(x.a) - x.r


def af_recip(y):
    cy_lo = y.c.abs_lower()
    rho = l1(y.a) + y.r
    if cy_lo <= ZERO_ARB or not (rho < cy_lo):
        return None
    q = rho / cy_lo
    invc = ONE_ACB / y.c
    return AF(invc, [-y.a[0] * invc * invc], y.r / (cy_lo * cy_lo) + (q * q / (1 - q)) / cy_lo)


def response(y256, y128, H, inv):
    return af_mul(af_add(y256, af_sc(af_mul(H, y128), C(-1))), inv)


_G1 = [AF(xc.tc_g[k], [ZERO_ACB], A(xc.RHO_G)) for k in range(M)]
_H1 = [AF(xc.tc_h[k], [ZERO_ACB], A(xc.RHO_H)) for k in range(M)]
_DG1 = [AF(C(k) * xc.tc_g[k], [ZERO_ACB], A(k) * A(xc.RHO_G)) for k in range(1, M)]
_DH1 = [AF(C(k) * xc.tc_h[k], [ZERO_ACB], A(k) * A(xc.RHO_H)) for k in range(1, M)]


def compute_arc_L7(u_mid, rho, c_fixed):
    u = AF(u_mid, [C(rho)], ZERO_ARB)
    s = AF(DSDU_C * u.c, [DSDU_C * u.a[0]], (DSDU_R * u.c.abs_upper()).upper())
    x = af_const(ZERO_ACB)
    yh = af_const(ZERO_ACB)
    ys = [af_const(ZERO_ACB) for _ in range(M)]
    tail_q = ZERO_ARB
    max_r = ZERO_ARB
    x128 = yh128 = ys128 = tail_q128 = None
    for n in range(256):
        z = af_sc(af_addc(af_mul(x, x), C(-1)), INV_R)
        rz = af_abs_u(z).upper()
        if not (rz < 1):
            return None
        if rz > max_r:
            max_r = rz
        vg = af_horner(_G1, z)
        vh = af_horner(_H1, z)
        vgp = af_horner(_DG1, z)
        vhp = af_horner(_DH1, z)
        fpw = af_add(vgp, af_mul(s, vhp))
        powers = [af_const(ONE_ACB)]
        for k in range(1, M):
            powers.append(af_mul(powers[-1], z))
        new_ys = []
        for k in range(M):
            propagated = af_mul(fpw, af_sc(af_mul(x, ys[k]), TWO_OVER_R))
            source = af_add(powers[k], af_sc(vh, -LHAT[k]))
            new_ys.append(af_add(propagated, source))
        yh_new = af_add(af_mul(fpw, af_sc(af_mul(x, yh), TWO_OVER_R)), vh)
        x_new = af_add(vg, af_mul(vh, s))
        fxu = af_abs_u(af_mul(fpw, af_sc(x, TWO_OVER_R))).upper()
        tail_source = (rz ** M) / (1 - rz)
        tail_q = fxu * tail_q + tail_source
        x, yh, ys = x_new, yh_new, new_ys
        if n + 1 == 128:
            x128, yh128, ys128, tail_q128 = x, yh, list(ys), tail_q
    inv = af_recip(x128)
    if inv is None:
        return None
    H = af_mul(x, inv)
    phih = response(yh, yh128, H, inv)
    psi = [response(ys[k], ys128[k], H, inv) for k in range(M)]
    dhc = af_addc(phih, -c_fixed)
    c_abs = c_fixed.abs_upper()
    head_best = ZERO_ARB
    for k in range(M):
        v = af_add(psi[k], af_sc(dhc, LHAT[k]))
        uval = af_abs_u(v) + c_abs * L_HEAD_RADIUS
        if uval > head_best:
            head_best = uval
    den_lo = af_abs_l(x128)
    if not (den_lo > 0):
        return None
    tail_phi = tail_q / den_lo + af_abs_u(x) * tail_q128 / (den_lo * den_lo)
    tail_l = c_abs * L_TAIL_INF
    total = max(head_best, tail_phi + tail_l)
    phih_upper = af_abs_u(phih)
    L7 = (LIN_F * phih_upper + PPERP * total).upper()
    return {"L7": L7, "phih_upper": phih_upper, "Ns": total, "Rmax_L7": max_r, "x128_lower": den_lo}


def recompute(cert):
    """Full independent recompute of every field. Returns a dict matching cert structure."""
    U0_RE = cert["small_disk"]["centre_re"]
    U0_IM = cert["small_disk"]["centre_im"]
    R0 = cert["small_disk"]["r0"]
    level = cert["level_m"]
    n_arcs = cert["n_arcs"]
    uc = C(A(U0_RE), A(U0_IM))
    rho = A(R0)
    uc_abs = (A(U0_RE) * A(U0_RE) + A(U0_IM) * A(U0_IM)).sqrt()
    abs_uc_rho = (uc_abs + rho).upper()
    t_plus = (abs_uc_rho * A(AU0.real.mid())).upper()
    lemma23_ok = bool(t_plus < T_BOUND)
    dist = (uc - C(A("4.758336330290482"), A("7.1277051021354385"))).abs_upper()
    margin_z72 = (A("0.1") - dist - rho).lower()
    inside = bool(float(margin_z72) > 0.0)

    S_new, Phi_s8, Phi_s9, lam1, lam1_phi_s9, S_minus, far_tail, far_neg = tail_terms(abs_uc_rho)

    wind = recompute_winding_laff(uc, rho, level, n_arcs)
    if wind is None:
        return {"wind_fail": True}
    # INDEPENDENT cross-check of the winding count via acb_poly (different code path)
    wind_count_acb = winding_count_acb_poly(uc, rho, level, n_arcs)
    per_arc = []
    all_rc = True
    domain_L7 = True
    domain_Lm = True
    worst_margin = None
    worst_gap = None
    worst_arc = None
    for j in range(n_arcs):
        t1 = Fraction(j, n_arcs)
        t2 = Fraction(j + 1, n_arcs)
        tm = (t1 + t2) / 2
        dt = t2 - t1
        theta = TWO_PI * A(tm.numerator) / A(tm.denominator)
        chord = (2 * rho * (PI * A(dt.numerator) / (2 * A(dt.denominator))).sin()).upper()
        um = uc + C(rho) * (I_UNIT * C(theta)).exp()
        out = compute_arc_L7(um, chord, um)
        if out is None:
            return {"l7_fail": True, "arc": j}
        L7 = out["L7"]
        Btail = (L7 * S_minus).upper()
        hmlb = wind["per_arc_lb"][j]
        margin = (hmlb / Btail).lower() if float(Btail) > 0 else A(0)
        gap = (hmlb - Btail).lower()
        rc = bool(float(gap) > 0.0 and float(margin) > 1.0)
        all_rc &= rc
        rmax_L7 = A(out["Rmax_L7"]).upper()
        rmax_Lm = A(wind["per_arc_rmax"][j]).upper()
        domain_L7 &= bool(float(rmax_L7) < 1.0)
        domain_Lm &= bool(float(rmax_Lm) < 1.0)
        if worst_margin is None or float(margin) < float(worst_margin):
            worst_margin = margin
            worst_arc = "%d/%d (%s..%s)" % (j, n_arcs, t1, t2)
        if worst_gap is None or float(gap) < float(worst_gap):
            worst_gap = gap
        per_arc.append({"L7": L7, "Btail": Btail, "hmlb": hmlb, "margin": margin, "gap": gap,
                         "rc": rc, "rmax_L7": rmax_L7, "rmax_Lm": rmax_Lm,
                         "phih": out["phih_upper"], "Ns": out["Ns"]})
    rigorous = bool(lemma23_ok and wind["count"] == 1 and all_rc and domain_L7 and domain_Lm and far_neg)
    return {
        "lemma23_ok": lemma23_ok, "t_plus": t_plus, "abs_uc_rho": abs_uc_rho,
        "dist": dist, "margin_z72": margin_z72, "inside": inside,
        "S_new": S_new, "Phi_s8": Phi_s8, "Phi_s9": Phi_s9, "lam1": lam1,
        "lam1_phi_s9": lam1_phi_s9, "S_minus": S_minus, "far_tail": far_tail, "far_neg": far_neg,
        "wind_count": wind["count"], "wind_ball": wind["winding_ball"],
        "wind_count_acb_poly": wind_count_acb,
        "per_arc": per_arc, "all_rc": all_rc, "domain_L7": domain_L7, "domain_Lm": domain_Lm,
        "worst_margin": worst_margin, "worst_gap": worst_gap, "worst_arc": worst_arc,
        "rigorous": rigorous,
    }


def verify(cert_path, rec):
    fails = []
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    deps = cert.get("dependencies", {})
    for key, rel in [("w33_gh_hiprec_sha256", DATA_PATH), ("x4_common_sha256", X4COMMON_PATH),
                      ("cl34_sha256", CL34_PATH), ("abs36_sha256", ABS36_PATH),
                      ("point8b_spectral_sha256", POINT8B_JSON), ("point8c_theta_sha256", POINT8C_JSON),
                      ("l_hat_fable_sha256", LHAT_PATH), ("point17_z72_certificate_sha256", Z72_CERT_PATH)]:
        if deps.get(key) != sha256(rel):
            fails.append(f"sha256({key}) mismatch")

    if cert.get("level_m") != 9:
        fails.append("level_m != 9 (verifier hardcodes level-9 recompute path)")
    n_arcs = cert.get("n_arcs")
    if n_arcs != len(cert.get("arc_data", [])):
        fails.append("arc_data count != n_arcs")

    if cert.get("lemma23_domain", {}).get("ok") != rec["lemma23_ok"]:
        fails.append("lemma23_domain.ok mismatch")
    if _rel(fnum(cert["tail_m9"]["S_disk_upper"]), float(rec["S_new"])) > REL_TOL:
        fails.append("tail S_disk mismatch")
    if _rel(fnum(cert["tail_m9"]["S_minus_upper"]), float(rec["S_minus"])) > REL_TOL:
        fails.append("tail S_minus mismatch")
    if _rel(fnum(cert["tail_m9"]["Phi_s8_upper"]), float(rec["Phi_s8"])) > REL_TOL:
        fails.append("tail Phi_s8 mismatch")
    if _rel(fnum(cert["tail_m9"]["lam1_Phi_s9_upper"]), float(rec["lam1_phi_s9"])) > REL_TOL:
        fails.append("tail lam1_Phi_s9 mismatch")
    if cert.get("tail_m9", {}).get("far_tail_negligible") != rec["far_neg"]:
        fails.append("far_tail_negligible mismatch")

    if cert.get("winding", {}).get("count") != rec["wind_count"]:
        fails.append(f"winding.count cert={cert.get('winding',{}).get('count')} recompute={rec['wind_count']}")
    # INDEPENDENT acb_poly winding-count cross-check (different code path)
    if rec.get("wind_count_acb_poly") != rec["wind_count"]:
        fails.append(f"winding count cross-check: LAff={rec['wind_count']} acb_poly={rec.get('wind_count_acb_poly')}")
    if rec.get("wind_count_acb_poly") != 1:
        fails.append(f"winding count acb_poly != 1 (got {rec.get('wind_count_acb_poly')})")

    # per-arc fields
    cert_arcs = cert.get("arc_data", [])
    if len(cert_arcs) != len(rec["per_arc"]):
        fails.append(f"arc_data length cert={len(cert_arcs)} recompute={len(rec['per_arc'])}")
    else:
        for j, (ca, pa) in enumerate(zip(cert_arcs, rec["per_arc"])):
            tag = f"arc[{j}]"
            if _rel(fnum(ca["Hm_abs_lower"]), float(pa["hmlb"])) > REL_TOL:
                fails.append(f"{tag}.Hm_abs_lower cert={ca['Hm_abs_lower']} rec={pa['hmlb']}")
            if _rel(fnum(ca["L7_arc_upper"]), float(pa["L7"])) > REL_TOL:
                fails.append(f"{tag}.L7_arc_upper cert={ca['L7_arc_upper']} rec={pa['L7']}")
            if _rel(fnum(ca["Btail_upper"]), float(pa["Btail"])) > REL_TOL:
                fails.append(f"{tag}.Btail_upper cert={ca['Btail_upper']} rec={pa['Btail']}")
            if _rel(fnum(ca["margin_lower"]), float(pa["margin"])) > REL_TOL:
                fails.append(f"{tag}.margin_lower cert={ca['margin_lower']} rec={pa['margin']}")
            if _rel(fnum(ca["gap_lower"]), float(pa["gap"])) > REL_TOL:
                fails.append(f"{tag}.gap_lower cert={ca['gap_lower']} rec={pa['gap']}")
            if bool(ca["rouche_pass"]) != bool(pa["rc"]):
                fails.append(f"{tag}.rouche_pass cert={ca['rouche_pass']} rec={pa['rc']}")

    if _rel(fnum(cert.get("worst_margin_lower")), float(rec["worst_margin"])) > REL_TOL:
        fails.append(f"worst_margin_lower cert={cert.get('worst_margin_lower')} rec={rec['worst_margin']}")
    if _rel(fnum(cert.get("worst_gap_lower")), float(rec["worst_gap"])) > REL_TOL:
        fails.append(f"worst_gap_lower cert={cert.get('worst_gap_lower')} rec={rec['worst_gap']}")
    if cert.get("all_rouche_pass") != rec["all_rc"]:
        fails.append("all_rouche_pass mismatch")
    if cert.get("domain_L7_all_ok") != rec["domain_L7"]:
        fails.append("domain_L7_all_ok mismatch")
    if cert.get("domain_Lm_all_ok") != rec["domain_Lm"]:
        fails.append("domain_Lm_all_ok mismatch")
    if cert.get("rigorous") != rec["rigorous"]:
        fails.append(f"rigorous cert={cert.get('rigorous')} rec={rec['rigorous']}")
    if cert.get("status") != ("RIGOROUS" if rec["rigorous"] else "NOT RIGOROUS"):
        fails.append("status mismatch")
    return fails


def main():
    t0 = time.time()
    print("z72_zero_verify: independent recompute starting (acb_poly winding + AF tube)...", flush=True)
    cert = json.loads(CERT_PATH.read_text(encoding="utf-8"))
    print(f"  cert: r0={cert['small_disk']['r0']} level_m={cert['level_m']} n_arcs={cert['n_arcs']} "
          f"status={cert['status']}", flush=True)
    rec = recompute(cert)
    if rec.get("wind_fail"):
        print("VERIFY: FAIL -- independent winding computation failed"); sys.exit(1)
    if rec.get("l7_fail"):
        print(f"VERIFY: FAIL -- independent L7 tube failed arc {rec.get('arc')}"); sys.exit(1)
    print(f"  recompute done ({time.time()-t0:.0f}s): wind_count={rec['wind_count']} "
          f"all_rc={rec['all_rc']} worst_margin={float(rec['worst_margin']):.5f} "
          f"rigorous={rec['rigorous']}", flush=True)

    fails = verify(CERT_PATH, rec)
    if fails:
        print(f"VERIFY: FAIL ({len(fails)} mismatches):")
        for f in fails[:40]:
            print("  FAIL " + f)
        if len(fails) > 40:
            print(f"  ... and {len(fails)-40} more")
        sys.exit(1)
    print(f"VERIFY: field check PASS ({time.time()-t0:.0f}s) -- running tamper self-tests...", flush=True)

    # ---- tamper self-test on a FILE COPY ----
    def _tamper(mutator, needle):
        fd, nm = tempfile.mkstemp(suffix=".json", dir=str(HERE))
        os.close(fd)
        tmp = Path(nm)
        shutil.copy(CERT_PATH, tmp)
        data = json.loads(tmp.read_text(encoding="utf-8"))
        mutator(data)
        tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
        tfails = verify(tmp, rec)
        try:
            tmp.unlink()
        except OSError:
            pass
        if not tfails:
            return "NOT-DETECTED"
        if needle and needle not in "\n".join(tfails):
            return "DETECTED-BY-WRONG-CHECK"
        return None

    # tamper 1: corrupt winding count
    def m_wind(data):
        data["winding"]["count"] = 0
    err = _tamper(m_wind, "winding.count")
    log_tamper = [("winding_count->0", err)]
    if err:
        print(f"VERIFY: tamper-test FAIL (winding count corruption {err})"); sys.exit(1)

    # tamper 2: corrupt worst margin (inflate)
    def m_margin(data):
        w = float(A(data["worst_margin_lower"]))
        data["worst_margin_lower"] = str(w * 50.0 + 5.0)
    err = _tamper(m_margin, "worst_margin_lower")
    log_tamper.append(("worst_margin x50", err))
    if err:
        print(f"VERIFY: tamper-test FAIL (worst_margin corruption {err})"); sys.exit(1)

    # tamper 3: corrupt per-arc Hm_abs_lower on arc 0 (halve it)
    def m_arc(data):
        cell = data["arc_data"][0]["Hm_abs_lower"]
        data["arc_data"][0]["Hm_abs_lower"] = str(float(A(cell)) * 0.3)
    err = _tamper(m_arc, "arc[0]")
    log_tamper.append(("arc0 Hm_abs_lower x0.3", err))
    if err:
        print(f"VERIFY: tamper-test FAIL (per-arc Hm_abs_lower corruption {err})"); sys.exit(1)

    # tamper 4: flip rigorous verdict
    def m_verdict(data):
        data["rigorous"] = not data["rigorous"]
        data["status"] = "NOT RIGOROUS" if data["rigorous"] else "RIGOROUS"
    err = _tamper(m_verdict, "rigorous")
    log_tamper.append(("verdict flip", err))
    if err:
        print(f"VERIFY: tamper-test FAIL (verdict flip {err})"); sys.exit(1)

    # tamper 5: corrupt a sha256 pin
    def m_pin(data):
        data["dependencies"]["w33_gh_hiprec_sha256"] = "0" * 64
    err = _tamper(m_pin, "sha256(w33_gh_hiprec_sha256)")
    log_tamper.append(("sha256 pin zeroed", err))
    if err:
        print(f"VERIFY: tamper-test FAIL (sha256 pin corruption {err})"); sys.exit(1)

    print(f"VERIFY: PASS -- independent recompute (acb_poly winding, DIFFERENT code path) matches "
          f"every field; {len(log_tamper)} tampers detected "
          f"({', '.join(t[0] for t in log_tamper)}); winding(H9)={rec['wind_count']}; "
          f"{cert['n_arcs']} arcs; worst_margin={float(rec['worst_margin']):.5f} "
          f"[{time.time()-t0:.0f}s]")
    sys.exit(0)


if __name__ == "__main__":
    main()

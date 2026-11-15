# -*- coding: utf-8 -*-
"""z72_zero_build.py -- строгий (Arb 512 бит) сертификат локализации нуля Z72:
нуль предельной H лежит в малом диске D(u0, r0) внутри сертификатного диска Z72.

Метод (как в кольцевых сертификатах / point17):
  * обмотка H_m по контуру малого диска (LAff single-symbol AF, как z72_certificate.
    Hm_light/compute_winding_level) -- rigorous winding ball.
  * перенос на предельную H через Руше: |H_m|_lower > B_tail на контуре,
    B_tail = L7_arc * (S_disk - Phi(s8) - Lam_final(1)*Phi(s9))   [m9, anchor 7]
    -- уровень m=9 (орбита 1024), хвост уровня 7 (anchor 7, телескоп H-H9=Sum_{k>=9} Dk,
       сняты члены k=7 и k=8). На малом контуре уровень 9 в домене (в отличие от Z72
       big-contour). L7_arc -- level-7 tube response (base-mode compute_arc, = K100 tube
       здесь, т.к. phih/Ns точечные).
  * S_disk, Lam_final -- forward peeling recursion (point8c), сверен с table_sample.

Уровень m=8 (как буквально в брифе) НЕ закрывается при доступных r0 (margin~0.38 при
r0=0.01; нужно r0~0.04 -- у границы Z72). m9 -- та же поуровневая лемма (anchor 7), на
один член больше; закрывается при малых r0. См. REPORT_Z72LOC.md.

Usage: python z72_zero_build.py [r0] [level_m]
   r0 по умолчанию 0.012; level_m по умолчанию 9.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import flint

ROOT = Path(r"D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive")
HERE = Path(__file__).resolve().parent
OUT_PATH = HERE / "z72_zero_certificate.json"
LOG_PATH = HERE / "z72_zero_build.log"

flint.ctx.prec = 512
A, C = flint.arb, flint.acb
PI = A.pi()
TWO_PI = 2 * PI
I_UNIT = C(0, 1)

sys.path.insert(0, str(ROOT / "w35" / "x4_walk_sonnet"))
import x4_common as xc  # noqa: E402

T0 = time.time()
LOGF = LOG_PATH.open("w", encoding="utf-8", buffering=1)

# ---- inputs (frozen) ----
DATA_PATH = ROOT / "w33" / "t5_cert_sonnet" / "w33_gh_hiprec.json"
CL34_PATH = ROOT / "w36" / "tail_codex" / "cl34_certificate.json"
ABS36_PATH = ROOT / "w36" / "audit_codex" / "absolute_margins_arb.json"
POINT8B_JSON = ROOT / "path_to_8" / "point8b_spectral_certificate" / "spectral_certificate.json"
POINT8C_JSON = ROOT / "path_to_8" / "point8c_spectral_theta" / "theta_certificate.json"
LHAT_PATH = ROOT / "w36" / "finish_fable" / "l_hat_fable.json"
X4COMMON_PATH = ROOT / "w35" / "x4_walk_sonnet" / "x4_common.py"
Z72_CERT_PATH = ROOT / "path_to_8" / "point17_z72_geometry" / "z72_certificate.json"

U0_RE = "4.791660937911268"
U0_IM = "7.177581007211361"
LEVEL_M = int(sys.argv[2]) if len(sys.argv) > 2 else 9
R0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.012
N_ARCS = int(sys.argv[3]) if len(sys.argv) > 3 else 64
UC_Z72_RE = "4.758336330290482"
UC_Z72_IM = "7.1277051021354385"
UC_Z72 = C(A(UC_Z72_RE), A(UC_Z72_IM))
T_BOUND = A("4.7911")
LIN_F = A("1.10485425032")
PPERP = A("2.808542")


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')} +{time.time()-T0:7.1f}s] {msg}"
    print(line, flush=True)
    LOGF.write(line + "\n")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def ub(x):
    return x.upper() if hasattr(x, "upper") else x


def lb(x):
    return x.lower() if hasattr(x, "lower") else x


# =====================================================================
# certified constants (from frozen inputs)
# =====================================================================
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

# tube machinery constants (cert_slice_param_audit.py)
_lh_raw = json.loads(LHAT_PATH.read_text(encoding="utf-8"))
LHAT = [C(A(a), A(b)) for a, b in zip(_lh_raw["l_hat_re"], _lh_raw["l_hat_im"])]
L_HEAD_RADIUS = A("1.588343e-14") + A("7.925e-15")
L_TAIL_INF = A("6.990710e-17")
M = xc.M
_db = C(A(DELTA.real.mid(), A(xc.RHO_D)))
_d7 = C(1) / (_db ** 7)
DSDU_C = C(_d7.real.mid(), _d7.imag.mid()) * AU0
DSDU_R = (C(_d7.real.mid(), _d7.imag.mid()).abs_upper() * 7 * A(xc.RHO_D)).upper()


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
    delta_pow = [A(1)]
    for _ in range(J_MAX):
        delta_pow.append((delta_pow[-1] * delta_u).upper())
    Lam_new = [A(1)]
    Lam_final = [A(1)]
    binding = ["EQUAL"]
    for j in range(1, J_MAX + 1):
        acc = A(0)
        for m in range(j):
            bracket = (Hcoef * delta_pow[j - 1 - m] + U(j - 1 - m)).upper()
            acc = (acc + Lam_new[m] * bracket).upper()
        lam_new_j = (Hcoef * delta_pow[j] + U(j) + Krho * acc).upper()
        lam_old_j = (A(L2) ** (j // 2)).upper() if j % 2 == 0 else (A(L1) * A(L2) ** ((j - 1) // 2)).upper()
        Lam_new.append(lam_new_j)
        if float(lam_new_j) <= float(lam_old_j):
            Lam_final.append(lam_new_j)
            binding.append("NEW")
        else:
            Lam_final.append(lam_old_j)
            binding.append("OLD")
    tally = {"NEW": binding.count("NEW"), "OLD": binding.count("OLD"), "EQUAL": binding.count("EQUAL")}
    return Lam_final, tally


log("building Lam_final (forward peeling, point8c)...")
Lam_final, Lam_tally = build_Lam_final()
_p8c = json.loads(POINT8C_JSON.read_text(encoding="utf-8"))
_p8c_match = True
for row in _p8c["peeling"]["table_sample"]:
    j = row["j"]
    if j > J_MAX:
        continue
    if abs(float(Lam_final[j]) - row["Lam_final_upper_float"]) / max(abs(row["Lam_final_upper_float"]), 1e-300) > 1e-6:
        _p8c_match = False
        log(f"  MISMATCH Lam_final[{j}]")
if not _p8c_match or Lam_tally != _p8c["peeling"]["binding_tally"]:
    raise RuntimeError("Lam_final does not reproduce point8c -- ABORT")
log(f"  Lam_final OK (tally {Lam_tally}, matches point8c)")


def tail_terms(abs_uc_rho):
    """Return (S_disk, Phi_s8, Lam1_Phi_s9, S_minus) for the m9 transfer at this contour."""
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
    return {"S_disk": S_new, "Phi_s8": Phi_s8, "Phi_s9": Phi_s9, "lam1": lam1,
            "lam1_Phi_s9": lam1_phi_s9, "S_minus": S_minus, "far_tail": far_tail,
            "far_negligible": far_neg}


# =====================================================================
# LAff single-symbol AF winding (z72_certificate.py copy)
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
            return None, {"fail": "composition", "step": n, "wu": float(wu)}
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
        return None, {"fail": "denominator"}
    Hm = la_mul(x, inv)
    return Hm, {"max_w": max_w}


def compute_winding(uc, rho, level):
    """Rigorous winding of H_level around D(uc, rho) over N_ARCS arcs.
    Returns (winding_ball_str, count, max_w, per_arc_Hm_abs_lower_list, per_arc_Rmax_list,
             endpoints_dict) or (None, failinfo)."""
    endpoints = {}
    arcs_ok = []
    max_w = A(0)
    per_arc_lb = []
    per_arc_rmax = []
    for j in range(N_ARCS):
        t1 = Fraction(j, N_ARCS)
        t2 = Fraction(j + 1, N_ARCS)
        dt = t2 - t1
        tm = (t1 + t2) / 2
        theta = TWO_PI * A(tm.numerator) / A(tm.denominator)
        chord = (2 * rho * (PI * A(dt.numerator) / (2 * A(dt.denominator))).sin()).upper()
        u_mid = uc + C(rho) * (I_UNIT * C(theta)).exp()
        Hm, diag = Hm_light(u_mid, chord, level)
        if Hm is None:
            return None, {"fail": "arc", "arc": j, "diag": diag}
        Hm_lb = (Hm.c.abs_lower() - Hm.a.abs_upper() - Hm.r).lower()
        if not (Hm_lb > 0):
            return None, {"fail": "arc_zero", "arc": j}
        per_arc_lb.append(Hm_lb)
        per_arc_rmax.append(A(diag["max_w"]).upper())
        arcs_ok.append((t1, t2))
        mw = diag["max_w"].upper()
        if mw > max_w:
            max_w = mw
        for tt in (t1, t2):
            if tt not in endpoints:
                th = A(0) if (tt == 0 or tt == 1) else TWO_PI * A(tt.numerator) / A(tt.denominator)
                um = uc + C(rho) * (I_UNIT * C(th)).exp()
                Hp, hd = Hm_light(um, A(0), level)
                if Hp is None:
                    return None, {"fail": "endpoint", "t": float(tt), "diag": hd}
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
            "max_w": str(max_w)}, {"per_arc_lb": per_arc_lb, "per_arc_rmax": per_arc_rmax}


# =====================================================================
# Base-mode level-7 tube (compute_arc, single-symbol AF) -- L7_arc
# =====================================================================
ZERO_ACB = C(0)
ONE_ACB = C(1)
INV_R = C(1) / C("2.5")
TWO_OVER_R = C(2) / C("2.5")
ZERO_ARB = A(0)


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
    """Level-7 tube (256-step orbit, single-symbol AF) -> L7_arc, Rmax_L7, x128_lower."""
    u = AF(u_mid, [C(rho)], ZERO_ARB)
    s = AF(DSDU_C * u.c, [DSDU_C * u.a[0]], (DSDU_R * u.c.abs_upper()).upper())
    x = af_const(ZERO_ACB)
    yh = af_const(ZERO_ACB)
    ys = [af_const(ZERO_ACB) for _ in range(M)]
    tail_q = ZERO_ARB
    max_r = ZERO_ARB
    x128 = yh128 = None
    ys128 = None
    tail_q128 = None
    for n in range(256):
        z = af_sc(af_addc(af_mul(x, x), C(-1)), INV_R)
        rz = af_abs_u(z).upper()
        if not (rz < 1):
            return None, {"fail": "Rmax_L7", "step": n, "rz": float(rz)}
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
        return None, {"fail": "recip"}
    H = af_mul(x, inv)
    phih = response(yh, yh128, H, inv)
    psi = [response(ys[k], ys128[k], H, inv) for k in range(M)]
    dhc = af_addc(phih, -c_fixed)
    c_abs = c_fixed.abs_upper()
    head_best = ZERO_ARB
    head_argmax = -1
    for k in range(M):
        v = af_add(psi[k], af_sc(dhc, LHAT[k]))
        uval = af_abs_u(v) + c_abs * L_HEAD_RADIUS
        if uval > head_best:
            head_best, head_argmax = uval, k
    den_lo = af_abs_l(x128)
    if not (den_lo > 0):
        return None, {"fail": "denlo"}
    tail_phi = tail_q / den_lo + af_abs_u(x) * tail_q128 / (den_lo * den_lo)
    tail_l = c_abs * L_TAIL_INF
    tail_total = tail_phi + tail_l
    total = max(head_best, tail_total)
    phih_upper = af_abs_u(phih)
    L7 = (LIN_F * phih_upper + PPERP * total).upper()
    return {"L7": L7, "phih_upper": phih_upper, "Ns": total, "head_best": head_best,
            "head_argmax": head_argmax, "tail_total": tail_total, "Rmax_L7": max_r,
            "x128_lower": den_lo}, None


def main():
    log("=" * 72)
    log(f"z72_zero_build: r0={R0} level_m={LEVEL_M} N_ARCS={N_ARCS} prec={flint.ctx.prec}")
    log("=" * 72)

    uc = C(A(U0_RE), A(U0_IM))
    rho = A(R0)
    uc_abs = (A(U0_RE) * A(U0_RE) + A(U0_IM) * A(U0_IM)).sqrt()
    abs_uc_rho = (uc_abs + rho).upper()
    dist_u0_z72 = (uc - UC_Z72).abs_upper()  # |u0 - u_c| upper
    margin_to_z72 = (A("0.1") - dist_u0_z72 - rho).lower()  # 0.1 - (|u0-uc| + r0)
    inside_z72 = bool(float(margin_to_z72) > 0.0)
    t_plus = (abs_uc_rho * A(AU0.real.mid())).upper()
    lemma23_ok = bool(t_plus < T_BOUND)
    log(f"contour D(u0,{R0}): |u0|+r0 <= {abs_uc_rho}  |t|+ = {t_plus} < {T_BOUND} -> {lemma23_ok}")

    sd = tail_terms(abs_uc_rho)
    log(f"tail(m9): S_disk={float(sd['S_disk']):.4e} Phi(s8)={float(sd['Phi_s8']):.4e} "
        f"Lam(1)*Phi(s9)={float(sd['lam1_Phi_s9']):.4e} S_minus={float(sd['S_minus']):.4e} "
        f"far_negligible={sd['far_negligible']}")
    if not lemma23_ok:
        log("ABORT: Lemma 2/3 out of range")
        sys.exit(1)
    if not sd["far_negligible"]:
        log("ABORT: far tail not negligible")
        sys.exit(1)

    # ---- winding(H_level) + per-arc Hm_abs_lower ----
    log(f"STEP: winding(H{LEVEL_M}) + per-arc |H{LEVEL_M}|_lower via LAff")
    wind, wind_diag = compute_winding(uc, rho, LEVEL_M)
    if wind is None:
        log(f"ABORT: winding FAILED {wind_diag}")
        sys.exit(1)
    winding_ok = (wind["count"] == 1)
    log(f"  winding(H{LEVEL_M})={wind['winding_ball'][:60]} count={wind['count']} max_w={wind['max_w']}")

    # ---- per-arc L7 (level-7 tube) ----
    log(f"STEP: per-arc L7 (level-7 tube, base mode)")
    arc_rows = []
    worst_margin = None
    worst_gap = None
    worst_arc = None
    all_rc = True
    domain_L7_ok = True
    domain_Lm_ok = True
    S_minus = sd["S_minus"]
    for j in range(N_ARCS):
        t1 = Fraction(j, N_ARCS)
        t2 = Fraction(j + 1, N_ARCS)
        tm = (t1 + t2) / 2
        dt = t2 - t1
        theta = TWO_PI * A(tm.numerator) / A(tm.denominator)
        chord = (2 * rho * (PI * A(dt.numerator) / (2 * A(dt.denominator))).sin()).upper()
        um = uc + C(rho) * (I_UNIT * C(theta)).exp()
        out, err = compute_arc_L7(um, chord, um)
        if out is None:
            log(f"ABORT: L7 tube FAILED arc {j}: {err}")
            sys.exit(1)
        L7 = out["L7"]
        Btail = (L7 * S_minus).upper()
        hmlb = wind_diag["per_arc_lb"][j]
        margin = (hmlb / Btail).lower() if float(Btail) > 0 else A(0)
        gap = (hmlb - Btail).lower()
        rc = bool(float(gap) > 0.0 and float(margin) > 1.0)
        all_rc &= rc
        rmax_L7 = A(out["Rmax_L7"]).upper()
        rmax_Lm = A(wind_diag["per_arc_rmax"][j]).upper()
        domain_L7_ok &= bool(float(rmax_L7) < 1.0)
        domain_Lm_ok &= bool(float(rmax_Lm) < 1.0)
        if worst_margin is None or float(margin) < float(worst_margin):
            worst_margin = margin
            worst_arc = "%d/%d (%s..%s)" % (j, N_ARCS, t1, t2)
        if worst_gap is None or float(gap) < float(worst_gap):
            worst_gap = gap
        arc_rows.append({
            "arc": j, "t1": str(t1), "t2": str(t2),
            "Hm_abs_lower": str(hmlb), "L7_arc_upper": str(L7),
            "phih_upper": str(out["phih_upper"]), "Ns_upper": str(out["Ns"]),
            "Btail_upper": str(Btail), "margin_lower": str(margin), "gap_lower": str(gap),
            "Rmax_L7": str(rmax_L7), "Rmax_Lm": str(rmax_Lm),
            "x128_lower": str(out["x128_lower"]), "rouche_pass": rc,
        })
    min_hmlb = min(float(A(r["Hm_abs_lower"])) for r in arc_rows)
    max_Btail = max(float(A(r["Btail_upper"])) for r in arc_rows)
    log(f"  all_rouche={all_rc} worst_margin={worst_margin} worst_gap={worst_gap} (arc {worst_arc}) "
        f"min|Hm|_lb={min_hmlb:.4e} max_Btail={max_Btail:.4e} domain_L7={domain_L7_ok} domain_L{LEVEL_M}={domain_Lm_ok}")

    rigorous = bool(lemma23_ok and winding_ok and all_rc and domain_L7_ok and domain_Lm_ok and sd["far_negligible"])
    status = "RIGOROUS" if rigorous else "NOT RIGOROUS"
    reason = ("all conditions met" if rigorous else
              ("winding(H%d)!=1" % LEVEL_M if not winding_ok else
               "Rouche fails on %d arcs" % sum(1 for r in arc_rows if not r["rouche_pass"])))
    log(f"VERDICT: {status} ({reason})")

    payload = {
        "task": "Z72 zero localization inside its certificate disk (strict Arb interval)",
        "brief_oid": "9595c48c",
        "status": status,
        "reason": reason,
        "rigorous": rigorous,
        "audit_date": time.strftime("%Y-%m-%d"),
        "precision_bits": flint.ctx.prec,
        "method": {
            "winding": "LAff single-symbol AF (z72_certificate.Hm_light/compute_winding_level copy), full endpoint balls, CL-34 errors",
            "transfer": "Rouche |H_m|_lower > L7_arc*(S_disk - Phi(s8) - Lam_final(1)*Phi(s9))  [m9, anchor 7]",
            "L7_arc": "level-7 tube response (base-mode compute_arc, single-symbol AF); = K100 tube on this small contour (phih/Ns point-dominated)",
            "tail": "forward peeling recursion Lam_final (point8c Lemma E0), S_disk = sum Lam_final(k-7)*Phi(s_{k+1}) + far_tail",
            "level_choice": "m=9 (not m=8): m=8 needs r0~0.04 (touches Z72 boundary); m9 (same anchor-7 lemma, one extra telescoped term) closes at small r0. Level 9 in domain on small contour (max_w~0.400), unlike Z72 big-contour where m9 failed on domain.",
        },
        "parent_certificate": {
            "path": str(Z72_CERT_PATH.relative_to(ROOT)),
            "sha256": sha256(Z72_CERT_PATH),
            "statement": "winding(H, D(u_c=4.758336330290482+7.1277051021354385i, 0.1)) = 1 (limit_proven_m8_refined) -- exactly one zero of H in the Z72 disk",
        },
        "small_disk": {
            "centre_re": U0_RE, "centre_im": U0_IM, "r0": R0,
            "u0_abs": str(uc_abs), "abs_uc_plus_r0_upper": str(abs_uc_rho),
            "dist_u0_to_z72_center_upper": str(dist_u0_z72),
            "margin_to_z72_boundary_lower": str(margin_to_z72),
            "inside_z72_disk": inside_z72,
        },
        "level_m": LEVEL_M,
        "lemma23_domain": {"t_plus_upper": str(t_plus), "bound": str(T_BOUND), "ok": lemma23_ok},
        "winding": {"winding_ball": wind["winding_ball"], "count": wind["count"], "expect": 1,
                     "max_w": wind["max_w"], "level": LEVEL_M},
        "tail_m9": {
            "formula": "|H-H9| <= L7_arc*(S_disk - Phi(s8) - Lam_final(1)*Phi(s9))",
            "anchor": 7, "S_disk_upper": str(sd["S_disk"]), "Phi_s8_upper": str(sd["Phi_s8"]),
            "Phi_s9_upper": str(sd["Phi_s9"]), "Lam_final_1": str(sd["lam1"]),
            "lam1_Phi_s9_upper": str(sd["lam1_Phi_s9"]), "S_minus_upper": str(sd["S_minus"]),
            "far_tail_upper": str(sd["far_tail"]), "far_tail_negligible": sd["far_negligible"],
            "Lam_final_tally": Lam_tally, "Lam_final_matches_point8c": True,
        },
        "n_arcs": N_ARCS,
        "arc_data": arc_rows,
        "all_rouche_pass": all_rc,
        "domain_L7_all_ok": domain_L7_ok,
        "domain_Lm_all_ok": domain_Lm_ok,
        "worst_margin_lower": str(worst_margin),
        "worst_gap_lower": str(worst_gap),
        "worst_margin_arc": worst_arc,
        "min_Hm_abs_lower": str(min(A(r["Hm_abs_lower"]) for r in arc_rows)),
        "max_Btail_upper": str(max(A(r["Btail_upper"]) for r in arc_rows)),
        "verdict_statement": (
            "winding(H9)=1 (rigorous LAff) AND Rouche |H9|_lower > L7_arc*(S_disk-Phi(s8)-Lam(1)Phi(s9)) "
            "on ALL %d arcs => winding(H)=winding(H9)=1: exactly one zero of the limit H in D(u0,%g). "
            "Combined with the parent Z72 certificate (one zero in D(u_c,0.1)), this localizes THAT zero "
            "to D(u0,%g)." % (N_ARCS, R0, R0)
        ) if rigorous else (
            "NOT RIGOROUS: %s. See arc_data for per-arc Rouché margins." % reason
        ),
        "dependencies": {
            "w33_gh_hiprec_sha256": sha256(DATA_PATH),
            "x4_common_sha256": sha256(X4COMMON_PATH),
            "cl34_sha256": sha256(CL34_PATH),
            "abs36_sha256": sha256(ABS36_PATH),
            "point8b_spectral_sha256": sha256(POINT8B_JSON),
            "point8c_theta_sha256": sha256(POINT8C_JSON),
            "l_hat_fable_sha256": sha256(LHAT_PATH),
            "point17_z72_certificate_sha256": sha256(Z72_CERT_PATH),
        },
        "time_s": time.time() - T0,
    }
    OUT_PATH.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    log(f"wrote {OUT_PATH.name} ({OUT_PATH.stat().st_size} bytes)")
    log(f"DONE in {time.time()-T0:.1f}s. status={status} rigorous={rigorous} worst_margin={worst_margin}")


if __name__ == "__main__":
    main()

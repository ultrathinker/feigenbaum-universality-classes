# -*- coding: utf-8 -*-
# [oid:kraw6H] _probe_paths.py -- one-shot measurement (NOT RIGOROUS, dev aid):
# ratios between the generator's square-path magnitudes and the verifier's
# hand-rolled tuple/disk-path magnitudes on selected roots.  Justifies the
# frozen comparison bands in krawczyk6_verify.py (observed worst factors are
# quoted in brackets there).  Kept for reproducibility; run after a
# certificate regeneration:
#     python -X utf8 _probe_paths.py
import json
import sys
from pathlib import Path

import flint

HERE = Path(__file__).resolve().parent
CERT = json.loads((HERE / "krawczyk6_certificates.json").read_text(encoding="utf-8"))


def recompute(spec):
    """Tuple (re, im, r) disk-path recompute -- the verifier core (deg 6)."""
    flint.ctx.prec = spec["prec_bits"]
    A = flint.arb

    def rabs_ub(x):
        a, b = abs(x.lower()), abs(x.upper())
        return b if bool(b > a) else a

    def vconst(re_, im_=None):
        return (A(re_), A(0) if im_ is None else A(im_), A(0))

    def vadd(x, y):
        re_, im_ = x[0] + y[0], x[1] + y[1]
        r = x[2] + y[2] + re_.rad() + im_.rad()
        return (re_, im_, r)

    def vaddc(x, c_re, c_im):
        re2, im2 = x[0] + A(c_re), x[1] + A(c_im)
        r = x[2] + re2.rad() + im2.rad()
        return (re2, im2, r)

    def vneg(x):
        return (-x[0], -x[1], x[2])

    def cmul_mid(a, b):
        return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])

    def mabs_ub(m):
        return (m[0] * m[0] + m[1] * m[1]).sqrt().upper()

    def vmul(x, y):
        m = cmul_mid((x[0], x[1]), (y[0], y[1]))
        r = (mabs_ub((x[0], x[1])) * y[2] + mabs_ub((y[0], y[1])) * x[2]
             + x[2] * y[2] + m[0].rad() + m[1].rad())
        return (m[0], m[1], r)

    def vscal6(x):
        re_, im_ = 6 * x[0], 6 * x[1]
        r = 6 * x[2] + re_.rad() + im_.rad()
        return (re_, im_, r)

    def vabs_ub(x):
        return ((rabs_ub(x[0]) + x[2]) ** 2
                + (rabs_ub(x[1]) + x[2]) ** 2).sqrt().upper()

    def vabs_lb(x):
        return (x[0] * x[0] + x[1] * x[1]).sqrt().lower() - x[2]

    n = spec["n"]
    m_lvl = spec["m"]
    r = A(spec["box_radius"])

    c0 = vconst(spec["box_center_re"], spec["box_center_im"])
    Abox = (c0[0], c0[1], c0[2] + r)              # DISK B(c0, r_reported)

    def step(z, dz, a):
        """x_{k+1} = 1 - a*x^6 ; dx_{k+1} = -x^6 - 6*a*x^5*dx (dx first)."""
        z2 = vmul(z, z)
        z3 = vmul(z2, z)
        z4 = vmul(z2, z2)                          # NOT z3*z: other association
        z5 = vmul(z4, z)
        z6 = vmul(z3, z3)                           # NOT z5*z
        az5dx = vmul(vmul(a, z5), dz)
        dz2 = vneg(vadd(z6, vscal6(az5dx)))
        az6 = vmul(a, z6)
        z2n = vaddc(vneg(az6), 1, 0)
        return z2n, dz2

    z = vconst(0)
    dz = vconst(0)
    for _ in range(n):
        z, dz = step(z, dz, c0)
    F0, Df0 = z, dz

    zb = vconst(0)
    dzb = vconst(0)
    lowers = []
    j = 0
    half = 1 << m_lvl
    for k in range(1, n + 1):
        zb, dzb = step(zb, dzb, Abox)
        if k <= half and k == (1 << j):
            lowers.append(vabs_lb(zb))
            j += 1
    D = dzb

    df_lb = vabs_lb(D)
    df_nonzero = bool(df_lb > 0)
    f0_ub = vabs_ub(F0)

    n_dist = (f0_ub / df_lb).upper()
    newton_inside = bool(n_dist < r)

    d0re, d0im = Df0[0], Df0[1]
    den = (d0re * d0re + d0im * d0im)
    if not bool(den.lower() > 0):
        raise RuntimeError("verifier: |f'(c0)| interval contains 0")
    yre_int = d0re / den
    yim_int = (-d0im) / den
    yre, yim = yre_int.mid(), yim_int.mid()
    yrad = yre_int.rad() + yim_int.rad() + yre.rad() + yim.rad()
    Yball = (yre, yim, yrad)
    Y_abs_ub = vabs_ub(Yball)
    E = vadd(vmul(Yball, D), vconst(-1))          # Y*D - 1
    e_abs_ub = vabs_ub(E)
    sqrt2 = A(2).sqrt().upper()
    k_dist = (Y_abs_ub * f0_ub + e_abs_ub * sqrt2 * r).upper()
    k_inside = bool(k_dist < r)
    kappa = k_dist / r

    return {
        "f0_abs_ub": f0_ub, "df_abs_lb": df_lb, "newton_dist_ub": n_dist,
        "krawczyk_dist_ub": k_dist, "krawczyk_kappa_ub": kappa.upper(),
        "e_abs_ub": e_abs_ub, "k_inside": k_inside, "newton_inside": newton_inside,
        "df_nonzero": df_nonzero,
        "period_abs_lb": lowers,
    }


def mid(x):
    return flint.arb(str(x) if isinstance(x, str) else x).mid()


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sel = ["S001", "S013", "S067", "S076", "S145", "S172", "S180", "S188",
           "S196", "S220", "S078", "S090", "S287"]
    roots = {r["label"]: r for r in CERT["roots"]}
    print("%-6s %-9s | rep/rec ratios: f0, df_lb, n_dist, k_dist, kappa, min_period"
          % ("label", "|E|"))
    for lbl in sel:
        rec = roots[lbl]
        spec = {"label": lbl, "n": rec["n"], "m": rec["m"],
                "box_center_re": rec["box_center_re"],
                "box_center_im": rec["box_center_im"],
                "box_radius": rec["box_radius"], "prec_bits": rec["prec_bits"]}
        v = recompute(spec)
        f0r = mid(rec["f0_abs_ub"]) / mid(v["f0_abs_ub"])
        dfr = mid(rec["df_abs_lb"]) / mid(v["df_abs_lb"])
        ndr = mid(rec["newton_dist_ub"]) / mid(v["newton_dist_ub"])
        kdr = mid(rec["krawczyk_dist_ub"]) / mid(v["krawczyk_dist_ub"])
        kapr = mid(rec["krawczyk_kappa_ub"]) / mid(v["krawczyk_kappa_ub"])
        mp_rep = min(mid(pc["abs_lb"]) for pc in rec["period_checks"][1:])
        mp_rec = min(mid(x) for x in v["period_abs_lb"][1:])
        perr = mp_rep / mp_rec
        print("%-6s | E=%.1e | %.4f %.4f %.4f %.4g %.4g %.4f | vK_in=%s vpass=%s"
              % (lbl, mid(v["e_abs_ub"]), f0r, dfr, ndr, kdr, kapr, perr,
                 v["k_inside"], v["k_inside"] and v["newton_inside"]
                 and v["df_nonzero"] and all(bool(x > 0) for x in v["period_abs_lb"])))


if __name__ == "__main__":
    main()

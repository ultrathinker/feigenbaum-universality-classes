# -*- coding: utf-8 -*-
"""F7: комбинаторные инварианты орбит суперстабильных параметров.

Для каждого суперстабильного a (и комплексной динамики x_{n+1}=1-a*x^d, d∈{4,6})
считаем:
  (a) глубина центрального возврата j*(a);
  (б) знаковый маршрут первых 16 точек (Re x_n, Im x_n);
  (в) масштабный код log|x_{2^j}|/log α, j=0..m.

mpmath dps=30, x_0=0. Период P=2^(m+1). Проверка |x_P|<1e-20.
"""
import json
import math
import os
import sys
import time
from pathlib import Path

import mpmath as mp


HERE = Path(__file__).parent
PROG = HERE / "progress.log"


def log_progress(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
    with open(PROG, "a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="")


# ---------- orbit helpers ----------
def orbit_points(a, d, P):
    """x_0=0, x_n=1-a*x_{n-1}^d, return list of mpc of length P+1."""
    mp.mp.dps = 30
    x = mp.mpf(0)
    pts = [x]
    d_mp = mp.mpf(d)
    for _ in range(P):
        # use pow for complex; but x is mpc here
        # x^d via x**d using mpmath power
        x = mp.mpf(1) - a * (x ** d_mp)
        pts.append(x)
    return pts


def orbit_dps30(a_re, a_im, d, P, a_re_str=None, a_im_str=None):
    """High-precision orbit at dps=30. Returns list of (re, im) Python float pairs.

    a_re_str / a_im_str are optional string precisions (mp.mpf-parseable) used if supplied.
    """
    mp.mp.dps = 30
    if a_re_str is not None:
        a_re_mp = mp.mpf(a_re_str)
    else:
        a_re_mp = mp.mpf(a_re)
    if a_im_str is not None:
        a_im_mp = mp.mpf(a_im_str) if a_im_str else mp.mpf(0)
    else:
        a_im_mp = mp.mpf(a_im) if a_im else mp.mpf(0)
    a = mp.mpc(a_re_mp, a_im_mp)
    d_mp = mp.mpf(d)
    x = mp.mpc(0, 0)
    out = []
    for _ in range(P + 1):
        out.append((float(mp.re(x)), float(mp.im(x))))
        x = mp.mpc(1, 0) - a * (x ** d_mp)
    return out


def check_period(a_re, a_im, d, P, a_re_str=None, a_im_str=None):
    """Return |x_P| at dps=30 (start from x_0=0)."""
    mp.mp.dps = 30
    if a_re_str is not None:
        a_re_mp = mp.mpf(a_re_str)
    else:
        a_re_mp = mp.mpf(a_re)
    if a_im_str is not None:
        a_im_mp = mp.mpf(a_im_str) if a_im_str else mp.mpf(0)
    else:
        a_im_mp = mp.mpf(a_im) if a_im else mp.mpf(0)
    a = mp.mpc(a_re_mp, a_im_mp)
    d_mp = mp.mpf(d)
    x = mp.mpc(0, 0)
    for _ in range(P):
        x = mp.mpc(1, 0) - a * (x ** d_mp)
    return abs(x)


def j_star(a_re, a_im, d, m, c=1.5, a_re_str=None, a_im_str=None):
    """Depth of central return (single-pass; returns max j with all suborbit points < c·r_j)."""
    P = 1 << (m + 1)
    orbit = orbit_dps30(a_re, a_im, d, P, a_re_str=a_re_str, a_im_str=a_im_str)
    last = 0
    for j in range(0, m + 1):
        stride = 1 << j
        if stride > P:
            break
        n_pts = P // stride
        if n_pts < 2:
            continue
        re0, im0 = orbit[stride]
        r_j = (re0 * re0 + im0 * im0) ** 0.5
        if r_j <= 0:
            r_j = 1e-300
        ok = True
        for k in range(1, n_pts):
            rr, ii = orbit[k * stride]
            mod2 = rr * rr + ii * ii
            mod = mod2 ** 0.5
            if mod > c * r_j:
                ok = False
                break
        if ok:
            last = j
    return last


def j_star_robust(a_re, a_im, d, m, c=1.5, a_re_str=None, a_im_str=None):
    """Same: max j such that all suborbit points < c · |x_{2^j}|."""
    P = 1 << (m + 1)
    orbit = orbit_dps30(a_re, a_im, d, P, a_re_str=a_re_str, a_im_str=a_im_str)
    last_ok = 0
    for j in range(0, m + 1):
        stride = 1 << j
        n_pts = P // stride
        if n_pts < 2:
            continue
        re0, im0 = orbit[stride]
        m0 = (re0 * re0 + im0 * im0) ** 0.5
        if m0 <= 0:
            m0 = 1e-300
        ok = True
        for k in range(1, n_pts):
            rr, ii = orbit[k * stride]
            mod2 = rr * rr + ii * ii
            mod = mod2 ** 0.5
            if mod > c * m0:
                ok = False
                break
        if ok:
            last_ok = j
    return last_ok


def j_star_strict(a_re, a_im, d, m, a_re_str=None, a_im_str=None):
    """Strict depth: max j such that |x_{2^{j+1}}| < |x_{2^j}| (monotonic decrease)."""
    P = 1 << (m + 1)
    orbit = orbit_dps30(a_re, a_im, d, P, a_re_str=a_re_str, a_im_str=a_im_str)
    mags = []
    for (r, i) in orbit:
        mags.append((r * r + i * i) ** 0.5)
    last = 0
    for j in range(0, m):
        a0 = mags[1 << j]
        a1 = mags[1 << (j + 1)]
        if a0 > 0 and a1 < a0:
            last = j + 1
        else:
            break
    return last


def orbit_first_mag(a_re, a_im, d, m, a_re_str=None, a_im_str=None):
    """Return list |x_{2^j}| for j=0..m."""
    P = 1 << (m + 1)
    orbit = orbit_dps30(a_re, a_im, d, P, a_re_str=a_re_str, a_im_str=a_im_str)
    out = []
    for j in range(0, m + 1):
        stride = 1 << j
        r, i = orbit[stride]
        out.append((r * r + i * i) ** 0.5)
    return out


def sign_route(a_re, a_im, d, P, a_re_str=None, a_im_str=None):
    """Sign of Re(x_n) and Im(x_n) for n=1..min(P,16). Returns (re_signs, im_signs)."""
    orbit = orbit_dps30(a_re, a_im, d, min(P, 16), a_re_str=a_re_str, a_im_str=a_im_str)
    re_signs = []
    im_signs = []
    for k in range(1, len(orbit)):
        r, i = orbit[k]
        re_signs.append("+" if r > 0 else ("-" if r < 0 else "0"))
        im_signs.append("+" if i > 0 else ("-" if i < 0 else "0"))
    return re_signs, im_signs


def scale_code(a_re, a_im, d, m, alpha, a_re_str=None, a_im_str=None):
    """log|x_{2^j}| / log(alpha) for j=0..m."""
    P = 1 << (m + 1)
    orbit = orbit_dps30(a_re, a_im, d, P, a_re_str=a_re_str, a_im_str=a_im_str)
    code = {}
    log_alpha = math.log(float(alpha))
    for j in range(0, m + 1):
        stride = 1 << j
        if stride > P:
            break
        r, i = orbit[stride]
        m2 = r * r + i * i
        m_abs = math.sqrt(m2) if m2 > 0 else 0.0
        if m_abs <= 0:
            ratio = None
        else:
            ratio = math.log(m_abs) / log_alpha
        code[j] = ratio
    return code


# ---------- drivers ----------
def _research_root():
    """D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип — но ищем автоматически."""
    p = Path(__file__).resolve()
    cur = p
    for _ in range(10):
        cur = cur.parent
        # match if this dir is named "Исследование-3-Принцип"
        if cur.name.startswith("Исследование-3-Принцип"):
            return cur
    # fallback: assume cwd is research root
    return Path("D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип")


def process_d6():
    """Process all 287 d=6 tier-1 (non-late) and 9 late."""
    log_progress("d=6: загрузка match_results.json + krawczyk_late")
    base = _research_root()
    mr_path = base / "agentB_deg6" / "match_results.json"
    late_path = base / "wave3_lateboxes" / "krawczyk6_late_certificates.json"
    targets_late_path = base / "wave3_lateboxes" / "targets_late.json"

    mr = json.load(open(mr_path, "r", encoding="utf-8"))
    late_cert = json.load(open(late_path, "r", encoding="utf-8"))
    targets_late = json.load(open(targets_late_path, "r", encoding="utf-8"))

    A = mp.mpf(mr["A6"])
    delta = mp.mpf(mr["delta6"])
    a_F = mp.mpf(mr["a_F"])
    alpha = (delta / A) ** mp.mpf(2)  # alpha_cascade (used in scale_code ratio)
    log_progress(f"  a_F={float(a_F):.12f}, A={float(A):.6f}, delta={float(delta):.6f}, alpha={float(alpha):.6f}")

    # Map id -> m for late
    late_ids = {tid: info["m"] for tid, info in targets_late["targets"].items()}
    log_progress(f"  late IDs: {sorted(late_ids.keys())} (m={sorted(set(late_ids.values()))})")

    # Build entry list — accepted only (have a_ss_re/im)
    rows = []
    seen_ids = set()
    for entry in mr["catalog_results"]:
        sid = entry["id"]
        if sid.startswith("Q"):
            continue
        if not (entry.get("accepted") and entry.get("exact_period_ok")):
            continue
        try:
            a_re_str = str(entry["a_ss_re"])
            a_im_str = str(entry.get("a_ss_im", "0.0") or "0.0")
            a_re = float(a_re_str)
            a_im = float(a_im_str) if a_im_str else 0.0
        except (KeyError, TypeError, ValueError):
            continue
        u_cat_re = float(entry["u_cat_re"])
        u_cat_im = float(entry["u_cat_im"])
        u_act_re = float(entry["u_actual_re"])
        u_act_im = float(entry["u_actual_im"])
        rows.append({
            "id": sid, "a_re": a_re, "a_im": a_im,
            "a_re_str": a_re_str, "a_im_str": a_im_str,
            "u_cat_re": u_cat_re, "u_cat_im": u_cat_im,
            "u_act_re": u_act_re, "u_act_im": u_act_im,
            "is_late": sid in late_ids, "m_cert": late_ids.get(sid, 9), "d": 6,
        })
        seen_ids.add(sid)

    # Late (not in match_results accepted set, but in late certificates)
    for root in late_cert["roots"]:
        sid = root["label"]
        a_re_str = str(root["box_center_re"])
        a_im_str = str(root["box_center_im"])
        a_re = float(a_re_str)
        a_im = float(a_im_str) if a_im_str else 0.0
        u_cat_re = float(root["u_cat_re"])
        u_cat_im = float(root["u_cat_im"])
        if sid in seen_ids:
            continue  # already covered (not the case in our data, but safe)
        rows.append({
            "id": sid, "a_re": a_re, "a_im": a_im,
            "a_re_str": a_re_str, "a_im_str": a_im_str,
            "u_cat_re": u_cat_re, "u_cat_im": u_cat_im,
            "u_act_re": u_cat_re, "u_act_im": u_cat_im,
            "is_late": True, "m_cert": root["m"], "d": 6,
        })

    log_progress(f"  catalog+late entries: {len(rows)} (d=6)")
    return rows, alpha


def _parse_complex_from_repr(s):
    """Parse '(re + im*j)' or '(re + 0.0j)' style string. Returns (re, im)."""
    if "+" in s[1:] and "j" in s:
        # find last '+' or '-' before 'j'
        idx_j = s.index("j")
        # split sign
        # split by '+'
        body = s.replace("(", "").replace(")", "")
        # last '+' or '-'
        for op in ("+", "-")[-1::-1]:
            pass
        # Easier: find the operator that separates re and im
        # im part starts with '+' or '-' after the real number
        # try regex
        import re
        m = re.match(r"\(?\s*([-+0-9.eE]+)\s*([+-])\s*([-+0-9.eE]+)\s*j\s*\)?", body)
        if m:
            re_p = float(m.group(1))
            sign = 1.0 if m.group(2) == "+" else -1.0
            im_p = float(m.group(3)) * sign
            return re_p, im_p
    # Fallback: try float
    return float(s), 0.0


def process_d4():
    base = _research_root()
    mr_path = base / "quartic_pilot" / "match_results.json"
    mr = json.load(open(mr_path, "r", encoding="utf-8"))
    A = mp.mpf(mr["A4"])
    delta = mp.mpf(mr["delta4"])
    a_F = mp.mpf(mr["a_F"])
    alpha = (delta / A) ** mp.mpf(2)

    rows = []
    for entry in mr["catalog_results"]:
        sid = entry["id"]
        # d=4 entries: id like Q01..
        if not sid.startswith("Q"):
            continue
        a_ss = entry.get("a_ss")
        if not a_ss:
            continue
        try:
            a_re, a_im = _parse_complex_from_repr(a_ss)
            a_re_str = str(a_re)
            a_im_str = str(a_im) if a_im else "0.0"
        except Exception:
            continue
        u_cat_re = float(entry["u_cat_re"])
        u_cat_im = float(entry["u_cat_im"])
        u_act_re = float(entry["u_actual_re"])
        u_act_im = float(entry["u_actual_im"])
        rows.append({
            "id": sid,
            "a_re": a_re, "a_im": a_im,
            "a_re_str": a_re_str, "a_im_str": a_im_str,
            "u_cat_re": u_cat_re, "u_cat_im": u_cat_im,
            "u_act_re": u_act_re, "u_act_im": u_act_im,
            "is_late": False,
            "m_cert": 9,
            "d": 4,
        })
    log_progress(f"  catalog entries with a_ss: {len(rows)} (d=4)")
    return rows, alpha


def lattice_distance(u_re, u_im, q):
    """Per-coordinate distance from u to nearest lattice node k/q."""
    f_re = u_re * q
    f_im = u_im * q
    # periodic distance to nearest integer
    n_re = round(f_re)
    n_im = round(f_im)
    p_re = (f_re - n_re) / q
    p_im = (f_im - n_im) / q
    d_per = (min(abs(f_re - n_re), 1.0 - abs(f_re - n_re)) / q,
             min(abs(f_im - n_im), 1.0 - abs(f_im - n_im)) / q)
    d2 = (p_re ** 2 + p_im ** 2) ** 0.5
    return d_per, d2


def classify(u_re, u_im, q=16, thr=0.004):
    """Classify: node (both <thr), semi (one), off."""
    d_per, d2 = lattice_distance(u_re, u_im, q)
    if d_per[0] < thr and d_per[1] < thr:
        return "node"
    elif d_per[0] < thr or d_per[1] < thr:
        return "semi"
    return "off"


def main():
    log_progress("=== f7_invariants.py ===")
    rows6, alpha6 = process_d6()
    rows4, alpha4 = process_d4()

    out = []
    t0 = time.time()
    total = len(rows6) + len(rows4)
    done = 0
    log_progress(f"Всего нулей к обработке: {total}")

    for r in rows6 + rows4:
        d = r["d"]
        m = r["m_cert"]
        P = 1 << (m + 1)
        a_re = r["a_re"]
        a_im = r["a_im"]
        try:
            xP = check_period(a_re, a_im, d, P, a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str"))
        except Exception as e:
            log_progress(f"  ! exception in check_period {r['id']}: {e}")
            xP = float("nan")
        try:
            js = j_star_robust(a_re, a_im, d, m, c=1.5, a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str"))
        except Exception as e:
            log_progress(f"  ! exception in j_star {r['id']}: {e}")
            js = -1
        try:
            sr_re, sr_im = sign_route(a_re, a_im, d, P, a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str"))
        except Exception as e:
            sr_re, sr_im = [], []
        alpha = alpha6 if d == 6 else alpha4
        try:
            sc = scale_code(a_re, a_im, d, m, alpha, a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str"))
        except Exception as e:
            sc = {}

        # lattice distance to u_actual
        u_re, u_im = r["u_act_re"], r["u_act_im"]
        d1_per16, d2_16 = lattice_distance(u_re, u_im, 16)
        d1_per32, d2_32 = lattice_distance(u_re, u_im, 32)
        cls16 = classify(u_re, u_im, q=16, thr=0.004)

        # lattice distance at u_cat (catalog prediction)
        uc_re, uc_im = r["u_cat_re"], r["u_cat_im"]
        cls_cat = classify(uc_re, uc_im, q=16, thr=0.004)
        d2_16_cat = lattice_distance(uc_re, uc_im, 16)[1]
        d2_32_cat = lattice_distance(uc_re, uc_im, 32)[1]

        out.append({
            "id": r["id"],
            "d": d,
            "m_cert": m,
            "is_late": r["is_late"],
            "a_re": a_re,
            "a_im": a_im,
            "u_act_re": u_re,
            "u_act_im": u_im,
            "u_cat_re": uc_re,
            "u_cat_im": uc_im,
            "xP_abs": float(xP),
            "j_star_c15": js,
            "j_star_c12": j_star_robust(a_re, a_im, d, m, c=1.2,
                                        a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str")),
            "j_star_c20": j_star_robust(a_re, a_im, d, m, c=2.0,
                                        a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str")),
            "j_star_dec": j_star_strict(a_re, a_im, d, m,
                                         a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str")),
            "first_mags": orbit_first_mag(a_re, a_im, d, m,
                                          a_re_str=r.get("a_re_str"), a_im_str=r.get("a_im_str")),
            "sign_re": sr_re,
            "sign_im": sr_im,
            "scale_code": sc,
            "d1_16_re_im": d1_per16,
            "d2_to_node16": d2_16,
            "d2_to_node32": d2_32,
            "class16": cls16,
            "class16_cat": cls_cat,
            "d2_to_node16_cat": d2_16_cat,
            "d2_to_node32_cat": d2_32_cat,
        })
        done += 1
        if done % 30 == 0:
            dt = time.time() - t0
            log_progress(f"  обработано {done}/{total}, прошло {dt:.1f} с")

    log_progress(f"Запись f7_invariants.json ({len(out)} записей)")
    with open(HERE / "f7_invariants.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    log_progress("=== f7_invariants.py DONE ===")


if __name__ == "__main__":
    main()

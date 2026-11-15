# -*- coding: utf-8 -*-
"""F7 — вещественные нули H при u > 0 (хаотическая сторона).

Одномерный скан H_d(u) для вещественных u>0 на уровне m=8:
  - d=6: u∈(0, 48],  шаг 1/512
  - d=4: u∈(0, 20],  шаг 1/512

Кандидаты:
  - знаковые смены H в скане (грубые)
  - каталожные нули H6_zeros.json / H4_zeros.json с Im=0 и Re>0 (точный начальный)

Каждый кандидат полируется Ньютоном (mpmath dps=40) либо оставляется каталожным.
Сохранение расстояния до ближайшего k/16, k/32, k/8.
"""
import json
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).parent
PROG = HERE / "progress.log"


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
    with open(PROG, "a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="")


def H_at_u_real(u_val, a_F, A, delta, d, m, mp_dps=40):
    """H(u) = x_{2^m} for a = a_F + u/(A*delta^m), x_0=0."""
    mp.mp.dps = mp_dps
    a_mp = mp.mpf(a_F) + mp.mpf(u_val) / (mp.mpf(A) * mp.mpf(delta) ** m)
    d_mp = mp.mpf(d)
    x = mp.mpf(0)
    n = 1 << m
    for _ in range(n):
        x = mp.mpf(1) - a_mp * (x ** d_mp)
    return x


def H_scan_array(u_lo, u_hi, step, a_F, A, delta, d, m):
    mp.mp.dps = 25
    ns = int(round((u_hi - u_lo) / step)) + 1
    us = np.linspace(u_lo, u_hi, ns)
    d_mp = mp.mpf(d)
    a_F_mp = mp.mpf(a_F)
    A_mp = mp.mpf(A)
    delta_mp = mp.mpf(delta)
    Hs = []
    for u_val in us:
        a_mp = a_F_mp + mp.mpf(float(u_val)) / (A_mp * delta_mp ** m)
        x = mp.mpf(0)
        n = 1 << m
        for _ in range(n):
            x = mp.mpf(1) - a_mp * (x ** d_mp)
        Hs.append(float(x))
    return us, Hs


def newton_polish_safe(u_init, a_F, A, delta, d, m, max_iter=20, dps=40, bracket_window=0.05):
    """Newton polish only if start is near a simple zero.

    Step capping: if |step| > bracket_window, abort.
    Stop early if not improving.
    """
    mp.mp.dps = dps
    u_mp = mp.mpf(u_init)
    delta_mp = mp.mpf(delta)
    a_F_mp = mp.mpf(a_F)
    A_mp = mp.mpf(A)
    d_mp = mp.mpf(d)
    n_iter = 1 << m
    last_h_abs = None
    for it in range(max_iter):
        a_mp = a_F_mp + u_mp / (A_mp * delta_mp ** m)
        # H
        x = mp.mpf(0)
        for _ in range(n_iter):
            x = mp.mpf(1) - a_mp * (x ** d_mp)
        # dH/da via dy/da
        y = mp.mpf(0)
        cur_x = mp.mpf(0)
        for _ in range(n_iter):
            x_new = mp.mpf(1) - a_mp * (cur_x ** d_mp)
            y_new = -(cur_x ** d_mp) - a_mp * d_mp * (cur_x ** (d_mp - 1)) * y
            cur_x = x_new
            y = y_new
        dH_du = y / (A_mp * delta_mp ** m)
        if abs(dH_du) < mp.mpf("1e-300"):
            return float(u_mp), float(cur_x), "dH/du tiny", it
        step_mp = cur_x / dH_du
        if abs(step_mp) > bracket_window:
            return float(u_mp), float(cur_x), f"step too large at iter {it}", it
        u_new = u_mp - step_mp
        if abs(u_new - u_mp) < mp.mpf("1e-25"):
            return float(u_new), float(cur_x), f"converged iter={it+1}", it
        u_mp = u_new
        h_abs = abs(cur_x)
        if last_h_abs is not None and h_abs >= last_h_abs:
            return float(u_mp), float(cur_x), f"non-improving iter={it+1}", it
        last_h_abs = h_abs
    return float(u_mp), float(cur_x), "max_iter", max_iter


def lattice_dist_real(u, q):
    f = u * q
    n = round(f)
    p = min(abs(f - n), 1.0 - abs(f - n)) / q
    return abs(p), p


def main():
    log("=== f7_realaxis.py ===")
    mr6 = json.load(open(_root() / "agentB_deg6" / "match_results.json", "r", encoding="utf-8"))
    mr4 = json.load(open(_root() / "quartic_pilot" / "match_results.json", "r", encoding="utf-8"))
    h6 = json.load(open(_root() / "agentB_deg6" / "H6_zeros.json", "r", encoding="utf-8"))
    h4 = json.load(open(_root() / "quartic_pilot" / "H4_zeros.json", "r", encoding="utf-8"))

    # === d=6 ===
    log("[d=6] Скан H(u) на u∈(0, 48] при m=8, шаг 1/512")
    a_F = mr6["a_F"]; A = mr6["A6"]; delta = mr6["delta6"]
    u_lo, u_hi, step, d, m = 0.0, 48.0, 1.0/512, 6, 8
    t0 = time.time()
    us, Hs = H_scan_array(u_lo, u_hi, step, a_F, A, delta, d, m)
    log(f"  H рассчитан на {len(us)} точках за {time.time()-t0:.1f} с")
    zeros_brackets = []
    for i in range(len(us) - 1):
        u0 = float(us[i]); u1 = float(us[i + 1])
        if Hs[i] * Hs[i + 1] < 0:
            zeros_brackets.append(0.5 * (u0 + u1))
        elif Hs[i] == 0.0:
            zeros_brackets.append(u0)
    log(f"  знаковых смен: {len(zeros_brackets)}")

    catalog_reals = [z["re"] for z in h6["zeros_level8"] if abs(z["im"]) < 1e-9 and z["re"] > 0]
    log(f"  реалы u>0 из каталога (m=8 в смысле полировки на этом уровне): {len(catalog_reals)}")
    log(f"  ВНИМАНИЕ: каталожные значения — предельные, они НЕ обязаны быть нулями H при m=8.")
    zeros_init = sorted(set([round(z, 6) for z in zeros_brackets]))
    log(f"  Стартовых точек (только знаковые смены): {len(zeros_init)}")

    # Полируем Ньютоном только стартовые точки от знаковых смен
    log(f"  Newton polish (safety: bracket_window=0.1)")
    polished = []
    for z_init in zeros_init:
        u_pol, h_pol, msg, it = newton_polish_safe(z_init, a_F, A, delta, d, m, max_iter=50, bracket_window=0.1)
        src = f"newton_iter={it}"
        polished.append({
            "u_init": z_init, "u_pol": u_pol, "h_pol": h_pol,
            "msg": msg, "src": src
        })

    # Де-дуп по u_pol
    seen = set()
    pol_dedup = []
    for p in sorted(polished, key=lambda x: x["u_pol"]):
        k = round(p["u_pol"], 5)
        if k in seen:
            continue
        seen.add(k)
        pol_dedup.append(p)
    polished = pol_dedup

    log(f"  Полированных (уникальных): {len(polished)}")
    for p in polished:
        u = p["u_pol"]
        p16 = lattice_dist_real(u, 16)
        p32 = lattice_dist_real(u, 32)
        p8 = lattice_dist_real(u, 8)
        log(f"    u_pol={u:.10f} H={p['h_pol']:.2e} dist16={p16[1]:.5f} dist32={p32[1]:.5f} dist8={p8[1]:.5f} src={p['src']}")

    near_16_d6 = sum(1 for p in polished if lattice_dist_real(p["u_pol"], 16)[1] < 0.004)
    near_32_d6 = sum(1 for p in polished if lattice_dist_real(p["u_pol"], 32)[1] < 0.004)
    near_8_d6 = sum(1 for p in polished if lattice_dist_real(p["u_pol"], 8)[1] < 0.004)
    log(f"[d=6] Полированных нулей: {len(polished)}, < 1/16: {near_16_d6}, < 1/32: {near_32_d6}, < 1/8: {near_8_d6}")

    # === d=4 ===
    log("[d=4] Скан H(u) на u∈(0, 20] при m=8, шаг 1/512")
    a_F = mr4["a_F"]; A = mr4["A4"]; delta = mr4["delta4"]
    u_lo, u_hi, step, d, m = 0.0, 20.0, 1.0/512, 4, 8
    t0 = time.time()
    us, Hs = H_scan_array(u_lo, u_hi, step, a_F, A, delta, d, m)
    log(f"  H рассчитан на {len(us)} точках за {time.time()-t0:.1f} с")
    zeros_brackets4 = []
    for i in range(len(us) - 1):
        if Hs[i] * Hs[i + 1] < 0:
            zeros_brackets4.append(0.5 * (float(us[i]) + float(us[i + 1])))
        elif Hs[i] == 0.0:
            zeros_brackets4.append(float(us[i]))
    log(f"  знаковых смен: {len(zeros_brackets4)}")

    catalog_reals4 = [z["re"] for z in h4["zeros_level8"] if abs(z["im"]) < 1e-9 and z["re"] > 0]
    log(f"  реалы u>0 из каталога d=4 (предельные, не обязаны быть m=8 нулями): {len(catalog_reals4)}: {catalog_reals4}")
    zeros_init4 = sorted(set([round(z, 6) for z in zeros_brackets4]))
    log(f"  Стартовых точек (только знаковые смены): {len(zeros_init4)}")

    polished4 = []
    for z_init in zeros_init4:
        u_pol, h_pol, msg, it = newton_polish_safe(z_init, a_F, A, delta, d, m, max_iter=50, bracket_window=0.1)
        src = f"newton_iter={it}"
        polished4.append({
            "u_init": z_init, "u_pol": u_pol, "h_pol": h_pol, "msg": msg, "src": src
        })

    seen = set()
    pol_dedup4 = []
    for p in sorted(polished4, key=lambda x: x["u_pol"]):
        k = round(p["u_pol"], 5)
        if k in seen:
            continue
        seen.add(k)
        pol_dedup4.append(p)
    polished4 = pol_dedup4

    log(f"  Полированных (уникальных): {len(polished4)}")
    for p in polished4:
        u = p["u_pol"]
        p16 = lattice_dist_real(u, 16)
        p32 = lattice_dist_real(u, 32)
        p8 = lattice_dist_real(u, 8)
        log(f"    u_pol={u:.10f} H={p['h_pol']:.2e} dist16={p16[1]:.5f} dist32={p32[1]:.5f} dist8={p8[1]:.5f} src={p['src']}")

    near_16_d4 = sum(1 for p in polished4 if lattice_dist_real(p["u_pol"], 16)[1] < 0.004)
    near_32_d4 = sum(1 for p in polished4 if lattice_dist_real(p["u_pol"], 32)[1] < 0.004)
    near_8_d4 = sum(1 for p in polished4 if lattice_dist_real(p["u_pol"], 8)[1] < 0.004)
    log(f"[d=4] Полированных нулей: {len(polished4)}, < 1/16: {near_16_d4}, < 1/32: {near_32_d4}, < 1/8: {near_8_d4}")

    out = {
        "d6_realaxis": {
            "a_F": a_F, "A": A, "delta": delta, "m": 8,
            "u_range": [0.0, 48.0],
            "step": 1.0/512,
            "polished": polished,
            "n_near_node16": near_16_d6,
            "n_near_node32": near_32_d6,
            "n_near_node8": near_8_d6,
        },
        "d4_realaxis": {
            "a_F": a_F, "A": A, "delta": delta, "m": 8,
            "u_range": [0.0, 20.0],
            "step": 1.0/512,
            "polished": polished4,
            "n_near_node16": near_16_d4,
            "n_near_node32": near_32_d4,
            "n_near_node8": near_8_d4,
        }
    }
    with open(HERE / "f7_realaxis.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    log("=== f7_realaxis.py DONE ===")


def _root():
    p = Path(__file__).resolve()
    cur = p
    for _ in range(10):
        cur = cur.parent
        if cur.name.startswith("Исследование-3-Принцип"):
            return cur
    return Path("D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип")


if __name__ == "__main__":
    main()

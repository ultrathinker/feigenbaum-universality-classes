# -*- coding: utf-8 -*-
"""F7 — сопоставление инвариантов с решёткой 1/16.

Для каждого нуля из f7_invariants.json:
  - класс по отношению к 1/16 (per-coordinate < 0.004), 1/32
  - кросс-таблица j_star_dec × class16
  - уникальность знакового кода в каждом классе
  - масштабный код: медиана log|x_{2^j}|/log α, медиана slope
  - соседи на расстоянии 1/16·(p,q), |p|,|q|<=2: есть ли?
"""
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).parent
PROG = HERE / "progress.log"


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
    with open(PROG, "a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="")


THR_NODE = 0.004  # per-coordinate ≤ 0.004 → узел


def lattice_dist(u_re, u_im, q):
    f_re = u_re * q
    f_im = u_im * q
    n_re = round(f_re)
    n_im = round(f_im)
    p_re = abs(f_re - n_re)
    p_im = abs(f_im - n_im)
    p_re = min(p_re, 1.0 - p_re) / q
    p_im = min(p_im, 1.0 - p_im) / q
    return (p_re, p_im, math.sqrt(p_re * p_re + p_im * p_im))


def classify(u_re, u_im, q, thr=THR_NODE):
    p_re, p_im, _ = lattice_dist(u_re, u_im, q)
    if p_re < thr and p_im < thr:
        return "node"
    elif p_re < thr or p_im < thr:
        return "semi"
    return "off"


def nearest_neighbor_dist(u_re, u_im, q):
    """Distance to nearest k/q lattice node."""
    return lattice_dist(u_re, u_im, q)[2]


def main():
    log("=== f7_lattice.py ===")
    data = json.load(open(HERE / "f7_invariants.json", "r", encoding="utf-8"))
    log(f"loaded {len(data)} entries")

    # ----- 1. classification per zero (1/16 and 1/32) -----
    cls16 = Counter()
    cls32 = Counter()
    d2_16_dist = []
    d2_32_dist = []
    for e in data:
        ur, ui = e["u_act_re"], e["u_act_im"]
        c16 = classify(ur, ui, 16)
        c32 = classify(ur, ui, 32)
        e["class16_v2"] = c16
        e["class32"] = c32
        e["d2_node16_v2"] = nearest_neighbor_dist(ur, ui, 16)
        e["d2_node32"] = nearest_neighbor_dist(ur, ui, 32)
        cls16[c16] += 1
        cls32[c32] += 1
        d2_16_dist.append(e["d2_node16_v2"])
        d2_32_dist.append(e["d2_node32"])

    log(f"class16 overall: {dict(cls16)}")
    log(f"class32 overall: {dict(cls32)}")

    # separate late / non-late d=6 / d=4
    groups = {
        "d6_late": [e for e in data if e["d"] == 6 and e["is_late"]],
        "d6_nonlate": [e for e in data if e["d"] == 6 and not e["is_late"]],
        "d6_all": [e for e in data if e["d"] == 6],
        "d4_all": [e for e in data if e["d"] == 4],
    }

    out_summary = {}
    for gname, grp in groups.items():
        log(f"--- {gname} (n={len(grp)}) ---")
        c16 = Counter(e["class16_v2"] for e in grp)
        c32 = Counter(e["class32"] for e in grp)
        log(f"  class16: {dict(c16)}; class32: {dict(c32)}")

        # median distances
        d16 = sorted(e["d2_node16_v2"] for e in grp)
        d32 = sorted(e["d2_node32"] for e in grp)
        log(f"  d2_node16 quantiles 5/25/50/75/95%: "
            f"{d16[len(d16)//20]:.2e}, {d16[len(d16)//4]:.2e}, "
            f"{d16[len(d16)//2]:.2e}, {d16[3*len(d16)//4]:.2e}, {d16[19*len(d16)//20]:.2e}")
        log(f"  d2_node32 quantiles: "
            f"{d32[len(d32)//20]:.2e}, {d32[len(d32)//4]:.2e}, "
            f"{d32[len(d32)//2]:.2e}, {d32[3*len(d32)//4]:.2e}, {d32[19*len(d32)//20]:.2e}")

        out_summary[gname] = {
            "n": len(grp),
            "class16": dict(c16),
            "class32": dict(c32),
            "d2_node16_quantiles": [d16[len(d16)//20], d16[len(d16)//4], d16[len(d16)//2], d16[3*len(d16)//4], d16[19*len(d16)//20]],
            "d2_node32_quantiles": [d32[len(d32)//20], d32[len(d32)//4], d32[len(d32)//2], d32[3*len(d32)//4], d32[19*len(d32)//20]],
        }

        # ----- 2. j_star_dec x class16 cross-tab -----
        log(f"  j_star_dec × class16:")
        cell = defaultdict(int)
        for e in grp:
            cell[(e["j_star_dec"], e["class16_v2"])] += 1
        js_vals = sorted(set(e["j_star_dec"] for e in grp))
        cls_vals = ["off", "semi", "node"]
        header = "        " + "  ".join(f"{c:>6s}" for c in cls_vals) + "   n"
        log(header)
        for j in js_vals:
            row = [cell[(j, c)] for c in cls_vals]
            log(f"   j*={j:2d} " + "  ".join(f"{v:>6d}" for v in row) + f"  {sum(row):4d}")
        out_summary[gname]["jstar_dec_class16"] = {
            f"{j},{c}": cell[(j, c)] for j in js_vals for c in cls_vals
        }

        # ----- 3. sign-route uniqueness -----
        log(f"  sign_route uniqueness (n={len(grp)}):")
        sr_strings = ["".join(e["sign_re"]) for e in grp]
        si_strings = ["".join(e["sign_im"]) for e in grp]
        sr_full = ["".join(e["sign_re"]) + "|" + "".join(e["sign_im"]) for e in grp]
        log(f"    Re-strings: {len(set(sr_strings))} unique / {len(sr_strings)} (покрытие {len(set(sr_strings))/len(sr_strings):.3%})")
        log(f"    Im-strings: {len(set(si_strings))} unique / {len(si_strings)}")
        log(f"    joint strings (Re+Im): {len(set(sr_full))} unique / {len(sr_full)}")
        # entropy of Re-sign string
        from collections import Counter as _C
        cc = _C(sr_full)
        ent = -sum((v / len(sr_full)) * math.log2(v / len(sr_full)) for v in cc.values())
        log(f"    joint sign-string entropy = {ent:.3f} бит; max log2(n) = {math.log2(len(sr_full)):.3f}")
        out_summary[gname]["sign_route_unique"] = {
            "re_unique": len(set(sr_strings)),
            "im_unique": len(set(si_strings)),
            "joint_unique": len(set(sr_full)),
            "joint_entropy_bits": ent,
        }

        # ----- 3b. sign codes per class16 -----
        log(f"  sign_route uniqueness по class16:")
        for c in cls_vals:
            sub = [e for e in grp if e["class16_v2"] == c]
            if not sub:
                continue
            sf = ["".join(e["sign_re"]) + "|" + "".join(e["sign_im"]) for e in sub]
            n_uniq = len(set(sf))
            top = _C(sf).most_common(3)
            log(f"    class={c:>5s} n={len(sub)} unique={n_uniq} "
                f"({n_uniq/len(sub):.3%}) top3={top}")

        # ----- 4. j_star_dec per |u| bin -----
        log(f"  j_star_dec × |u| bins:")
        bins = [(0, 8), (8, 16), (16, 24), (24, 32), (32, 50), (50, 1e9)]
        for lo, hi in bins:
            sub = [e for e in grp if lo <= math.hypot(e["u_act_re"], e["u_act_im"]) < hi]
            if not sub:
                continue
            cc = Counter(e["j_star_dec"] for e in sub)
            log(f"    |u|∈[{lo},{hi}) n={len(sub)} j_dec={dict(cc)}")

    # ----- 5. lattice neighbors at 1/16·(p,q) -----
    # For each zero, find zeros within lattice steps of 1/16·(p,q) with |p|,|q|<=2.
    # Test whether paired zeros (at lattice shifts) have correlated orbits.
    log("=== lattice neighbors test ===")
    for gname in ["d6_all", "d4_all"]:
        grp = groups[gname]
        pairs = defaultdict(list)
        for i, e in enumerate(grp):
            for j, ee in enumerate(grp):
                if j <= i:
                    continue
                # distance between u's along lattice steps
                du_re = (e["u_act_re"] - ee["u_act_re"]) * 16
                du_im = (e["u_act_im"] - ee["u_act_im"]) * 16
                p, q = round(du_re), round(du_im)
                if abs(p) <= 2 and abs(q) <= 2 and abs(p) + abs(q) > 0:
                    # exactly at 1/16 lattice step
                    err_re = abs(du_re - p)
                    err_im = abs(du_im - q)
                    if err_re < 0.1 and err_im < 0.1:  # tolerance
                        key = (p, q)
                        pairs[key].append((i, j))
        # collect counts
        ne = sum(1 for v in pairs.values() if v)
        log(f"  {gname}: {sum(len(v) for v in pairs.values())} pairs at lattice steps")
        for key in sorted(pairs.keys())[:8]:
            ps = pairs[key]
            log(f"    step (p,q)={key}: {len(ps)} pairs; example: {ps[0]}")

    # ----- 6. save summary -----
    with open(HERE / "f7_lattice.json", "w", encoding="utf-8") as f:
        # strip non-serializable bit (it was added in-place; we serialize)
        serializable = []
        for e in data:
            se = {k: v for k, v in e.items() if not isinstance(v, dict) or k != "scale_code"}
            serializable.append(se)
        json.dump({"summary": out_summary, "data_with_class": serializable},
                  f, ensure_ascii=False, indent=1)
    log("=== f7_lattice.py DONE ===")


if __name__ == "__main__":
    main()

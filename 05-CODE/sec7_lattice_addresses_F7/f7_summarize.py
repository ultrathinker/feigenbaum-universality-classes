# -*- coding: utf-8 -*-
"""F7 — финальный анализ и сборка таблиц/выводов.

Использует f7_invariants.json, f7_lattice.json, f7_realaxis.json.
"""
import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

HERE = Path(__file__).parent
PROG = HERE / "progress.log"


def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
    with open(PROG, "a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="")


THR_NODE = 0.004  # per-coord threshold for "node"


def nearest_node_dist(u_re, u_im, q):
    f_re = u_re * q; f_im = u_im * q
    n_re = round(f_re); n_im = round(f_im)
    p_re = abs(f_re - n_re); p_im = abs(f_im - n_im)
    p_re = min(p_re, 1.0 - p_re) / q
    p_im = min(p_im, 1.0 - p_im) / q
    return math.sqrt(p_re * p_re + p_im * p_im)


def main():
    log("=== f7_summarize.py ===")
    data = json.load(open(HERE / "f7_invariants.json", "r", encoding="utf-8"))
    rl = json.load(open(HERE / "f7_realaxis.json", "r", encoding="utf-8"))

    # ----- Section A: j_star_dec — descriptive stats -----
    log("--- j_star_dec distribution ---")
    d6_late = [d for d in data if d["d"] == 6 and d["is_late"]]
    d6_non = [d for d in data if d["d"] == 6 and not d["is_late"]]
    d4 = [d for d in data if d["d"] == 4]
    log(f"  d6 late: j_star_dec={Counter(e['j_star_dec'] for e in d6_late)}")
    log(f"  d6 non-late: j_star_dec={Counter(e['j_star_dec'] for e in d6_non)}")
    log(f"  d4 all: j_star_dec={Counter(e['j_star_dec'] for e in d4)}")

    # ----- Section B: j_star_dec × hierarchical lattice class -----
    log("--- j_star_dec × hierarchical lattice (1/16 близость vs 1/32 близость) ---")
    def hclass(e):
        d16 = e["d2_to_node16"]
        d32 = e["d2_to_node32"]
        if d16 < 0.004:
            return "1/16"
        elif d32 < 0.004:
            return "1/32-not-1/16"
        else:
            return "off"
    for grp_name, grp in [("d6_late", d6_late), ("d6_non", d6_non), ("d4", d4)]:
        cell = defaultdict(int)
        for e in grp:
            cell[(e["j_star_dec"], hclass(e))] += 1
        js_vals = sorted(set(e["j_star_dec"] for e in grp))
        cls_vals = ["1/16", "1/32-not-1/16", "off"]
        log(f"  {grp_name} (n={len(grp)}):")
        header = "        " + "  ".join(f"{c:>18s}" for c in cls_vals) + "  sum"
        log(header)
        for j in js_vals:
            row = [cell[(j, c)] for c in cls_vals]
            log(f"   j*={j:2d} " + "  ".join(f"{v:>18d}" for v in row) + f"  {sum(row):4d}")

    # ----- Section C: late 9 vs non-late 287 — feature contrasts -----
    log("--- Сравнение 9 поздних vs 333 не-поздних (d=6) ---")
    log(f"  Метрика                     late_mean         non-late_mean   ratio")
    for key in ["d2_to_node16", "d2_to_node32"]:
        v_late = [e[key] for e in d6_late]
        v_non = [e[key] for e in d6_non]
        log(f"  {key:25s} {mean(v_late):.4e}       {mean(v_non):.4e}        {mean(v_late)/mean(v_non):.3f}")
    # Лог-расстояние
    log_late = [math.log10(e["d2_to_node16"]) for e in d6_late]
    log_non = [math.log10(e["d2_to_node16"]) for e in d6_non]
    log(f"  log10(d2_to_node16) {mean(log_late):.3f}            {mean(log_non):.3f}          "
        f"{10**(mean(log_late)-mean(log_non)):.2f}x")

    # Коды
    late_codes = [("".join(e["sign_re"]) + "|" + "".join(e["sign_im"])) for e in d6_late]
    non_codes = [("".join(e["sign_re"]) + "|" + "".join(e["sign_im"])) for e in d6_non]
    log(f"  late знаковых кодов уникальных: {len(set(late_codes))} / {len(late_codes)}")
    log(f"  non-late знаковых кодов уникальных: {len(set(non_codes))} / {len(non_codes)}")

    # j_star_dec per late
    log(f"  Поздние (9):")
    for e in d6_late:
        d16 = e["d2_to_node16"]
        s8 = "".join(e["sign_re"][:8])
        log(f"    {e['id']:>5s} m={e['m_cert']:2d} j*={e['j_star_dec']:2d} "
            f"d2_node16={d16:.3e} sign_re[0:8]={s8}")

    # ----- Section D: distance to 1/16 vs 1/32 — are they different populations? -----
    log("--- d2_to_node16 vs d2_to_node32 ---")
    for grp_name, grp in [("d6_late", d6_late), ("d6_non", d6_non), ("d4", d4)]:
        ds16 = [e["d2_to_node16"] for e in grp]
        ds32 = [e["d2_to_node32"] for e in grp]
        log(f"  {grp_name}: median(d_16)={median(ds16):.3e}, median(d_32)={median(ds32):.3e}")
        # сколько из тех, у кого d_32 < thr не имеют d_16 < thr (т.е. узлы 1/32, не 1/16)
        n_pure_32 = sum(1 for e in grp if e["d2_to_node16"] >= THR_NODE and e["d2_to_node32"] < THR_NODE)
        n_pure_16 = sum(1 for e in grp if e["d2_to_node16"] < THR_NODE)
        n_both = sum(1 for e in grp if e["d2_to_node16"] < THR_NODE and e["d2_to_node32"] < THR_NODE)
        n_off = sum(1 for e in grp if e["d2_to_node16"] >= THR_NODE and e["d2_to_node32"] >= THR_NODE)
        log(f"    1/16-only: {n_pure_16}, 1/32-not-1/16: {n_pure_32}, оба: {n_both}, off: {n_off}")

    # ----- Section E: знаковая уникальность среди разных классов -----
    log("--- Знаковая уникальность: per lattice class ---")
    for grp_name, grp in [("d6_non", d6_non), ("d6_late", d6_late), ("d4", d4)]:
        codes_by_cls = defaultdict(list)
        for e in grp:
            codes_by_cls[hclass(e)].append("".join(e["sign_re"]) + "|" + "".join(e["sign_im"]))
        log(f"  {grp_name}:")
        for cls, lst in sorted(codes_by_cls.items(), key=lambda kv: -len(kv[1])):
            top = Counter(lst).most_common(1)[0]
            tok_short = (top[0][:32] + "...") if len(top[0]) > 32 else top[0]
            log(f"    {cls:18s} n={len(lst):3d}  unique={len(set(lst)):3d}  "
                f"top1={tok_short} count={top[1]}")

    # ----- Section F: real-axis scan summary -----
    log("--- Вещественная ось u>0 (m=8) ---")
    d6_real = rl["d6_realaxis"]
    d4_real = rl["d4_realaxis"]
    log(f"  d=6: {len(d6_real['polished'])} нуль(ей) в (0, 48]:")
    for p in d6_real["polished"]:
        log(f"    u={p['u_pol']:.6f}, dist16={nearest_node_dist(p['u_pol'], 0, 16):.4e}, "
            f"dist32={nearest_node_dist(p['u_pol'], 0, 32):.4e}")
    log(f"  d=4: {len(d4_real['polished'])} нуль(ей) в (0, 20]")
    for p in d4_real["polished"]:
        log(f"    u={p['u_pol']:.6f}, dist16={nearest_node_dist(p['u_pol'], 0, 16):.4e}, "
            f"dist32={nearest_node_dist(p['u_pol'], 0, 32):.4e}")

    # ----- Section G: hierarchical relationship j* ↔ 1/16 -----
    log("--- Иерархия j* ↔ расстояние к 1/16 ---")
    # Кубируем: распределение d_16 в зависимости от j_star_dec
    cell_dist = defaultdict(list)
    for e in data:
        cell_dist[(e["d"], e["j_star_dec"])].append(e["d2_to_node16"])
    for d_val in [4, 6]:
        log(f"  d={d_val}:")
        for j_key in sorted(k for (d, k) in cell_dist if d == d_val):
            ds = cell_dist[(d_val, j_key)]
            log(f"    j*={j_key:2d}  n={len(ds):3d}  median(d_16)={median(ds):.3e}  "
                f"Q25/Q75={sorted(ds)[len(ds)//4]:.3e}/{sorted(ds)[3*len(ds)//4]:.3e}")

    log("=== f7_summarize.py DONE ===")


if __name__ == "__main__":
    main()

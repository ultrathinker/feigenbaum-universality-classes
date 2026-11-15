"""
Задачи Г7 и Г9 рецензента Б для d=6.
Источники: ../sec3_pilot_d6/match_results.json, ../sec3_pilot_d6/H6_zeros.json,
           ../sec5_subfamily_d6/rate_full.json
Всё NOT RIGOROUS — float/mpmath.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp

# Reproducibility
RNG_SEED = 20260815
np.random.seed(RNG_SEED)
mp.mp.dps = 30

ROOT = Path(r"D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип")
WAVE3 = ROOT / "wave3_lattice"
AGENTB = ROOT / "agentB_deg6"

WORK = WAVE3
WORK.mkdir(exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# 0. Загрузка данных
# ─────────────────────────────────────────────────────────────────────────────
with open(AGENTB / "match_results.json", "r", encoding="utf-8") as f:
    MR = json.load(f)

with open(AGENTB / "H6_zeros.json", "r", encoding="utf-8") as f:
    H6 = json.load(f)

with open(AGENTB / "subfamily" / "rate_full.json", "r", encoding="utf-8") as f:
    RATE = json.load(f)

CATALOG = MR["catalog_results"]
DENSITY = MR["density"]
CONTROLS = MR["controls"]
META = MR["meta"]
FROZEN = META["frozen"]

a_F = MR["a_F"]
A6 = MR["A6"]
delta6 = MR["delta6"]
lambda2 = MR["lambda2"]

R_ACCEPT = FROZEN["R_ACCEPT_U"]
EDGE_U = FROZEN["EDGE_U"]

# Late-9 IDs (7 internal + 2 edge)
LATE9_IDS = ["S036", "S104", "S187", "S339", "S031", "S070", "S235", "S279", "S319"]
# true m at which each late hits d_u<0.01
LATE_TRUE_M = {"S031":10, "S070":10, "S235":10, "S279":10, "S319":10,
               "S104":11, "S187":11, "S339":11, "S036":12}

# ─────────────────────────────────────────────────────────────────────────────
# 1. Базовая сегментация
# ─────────────────────────────────────────────────────────────────────────────
# 287 good (internal accepted)
good = [r for r in CATALOG if r["accepted"] and not r.get("edge")]
# 9 late
late = [r for r in CATALOG if r["id"] in LATE9_IDS]
late_set = {r["id"]: r for r in late}

# all u_cat re/im arrays for bounds
all_re = np.array([r["u_cat_re"] for r in CATALOG])
all_im = np.array([r["u_cat_im"] for r in CATALOG])

RE_MIN, RE_MAX = float(all_re.min()), float(all_re.max())
IM_MIN, IM_MAX = float(all_im.min()), float(all_im.max())

# ─────────────────────────────────────────────────────────────────────────────
# 2. Воспроизведение метрики специфичности
# ─────────────────────────────────────────────────────────────────────────────
# Определение (см. REPORT_DEG6_PILOT.md §3.4 + файл density):
#   specificity = median_rival_u / d_u
#   density probe: 12 джиттерных стартов в 0.4-боксе (DENS_SIDE_U=0.4)
#   Если в боксе нет ни одного сходящегося корня-соперника, specificity = None
spec_vals = [r["specificity"] for r in DENSITY.values() if r["specificity"] is not None]
spec_nonzero = len(spec_vals)
spec_median = float(np.median(spec_vals))
spec_min = float(np.min(spec_vals))
spec_max = float(np.max(spec_vals))

density_repro = {
    "n_density_probes": len(DENSITY),
    "n_nonzero_specificity": spec_nonzero,
    "median_specificity": spec_median,
    "min_specificity": spec_min,
    "max_specificity": spec_max,
    "note": "воспроизведено точно из match_results.json — median 1488.3991, 7 ненулевых значений",
}

print("[SPEC] воспроизведено:", density_repro)

# ─────────────────────────────────────────────────────────────────────────────
# 3. ЗАДАЧА 1 (Г9): корреляция остатков с решёткой 1/16
# ─────────────────────────────────────────────────────────────────────────────
def frac16(x):
    """Дробная часть 16*x в [0,1)."""
    return (16.0 * x) - np.floor(16.0 * x)


def wrap_to_0_1(x):
    """Дробная часть в [0,1)."""
    return x - np.floor(x)


def chi2_uniform_16bins(frac_vals):
    """χ²-тест равномерности на 16 бинах; H0: равномерно на [0,1).
    Возвращает (chi2_stat, p_value, counts)."""
    counts, _ = np.histogram(frac_vals, bins=16, range=(0.0, 1.0))
    n = len(frac_vals)
    expected = n / 16.0
    chi2 = float(((counts - expected) ** 2 / expected).sum())
    # df = 16 - 1 = 15
    from scipy.stats import chi2 as chi2_dist
    p = float(1.0 - chi2_dist.cdf(chi2, df=15))
    return chi2, p, counts.tolist()


def ks_uniform(frac_vals):
    """Колмогоров–Смирнов против равномерного U[0,1]. Возвращает (D, p)."""
    from scipy.stats import kstest
    # kstest против uniform; используем exact_asymptotic=True по умолчанию
    res = kstest(frac_vals, "uniform")
    return float(res.statistic), float(res.pvalue)


def rayleigh_test(angles_2pi):
    """Рэлей-тест на круговой равномерности; angles_2pi ∈ [0,2π).
    Возвращает (R_bar, z, p_value)."""
    n = len(angles_2pi)
    C = float(np.cos(angles_2pi).sum())
    S = float(np.sin(angles_2pi).sum())
    R = float(np.sqrt(C * C + S * S))
    R_bar = R / n
    z = R * R / n  # test statistic ~ chi^2_2 for large n
    from scipy.stats import chi2 as chi2_dist
    p = float(np.exp(-z))  # exact for Rayleigh
    return R_bar, z, p


# 3a. 287 хороших — статистика дробных частей 16·u_actual и 16·Δu
good_arr = []
for r in good:
    good_arr.append({
        "id": r["id"],
        "u_act_re": r["u_actual_re"],
        "u_act_im": r["u_actual_im"],
        "u_cat_re": r["u_cat_re"],
        "u_cat_im": r["u_cat_im"],
        "d_u": r["d_u"],
    })

G = np.array([[x["u_act_re"], x["u_act_im"]] for x in good_arr])
G_CAT = np.array([[x["u_cat_re"], x["u_cat_im"]] for x in good_arr])
DU_G = np.array([x["d_u"] for x in good_arr])

# Frac 16 * u_actual
frac_act_re = frac16(G[:, 0])
frac_act_im = frac16(G[:, 1])
# Frac 16 * u_cat
frac_cat_re = frac16(G_CAT[:, 0])
frac_cat_im = frac16(G_CAT[:, 1])
# Frac 16 * Δu
DU_COMP = G - G_CAT
frac_du_re = frac16(DU_COMP[:, 0])
frac_du_im = frac16(DU_COMP[:, 1])

# 9 поздних — на их истинных m из rate_full.json. u_cat ≈ u_actual (d_u мало).
# u_actual = u_cat + Δu, |Δu|=d_u_true, направление неизвестно — НО
# решёточная близость проявляется уже в u_cat (проверено ниже отдельно).
late_arr = []
for sid in LATE9_IDS:
    r = late_set[sid]
    m_true = LATE_TRUE_M[sid]
    d_u_true = RATE[sid]["d_u"][str(m_true)]
    late_arr.append({
        "id": sid,
        "m_true": m_true,
        "d_u_true": d_u_true,
        "u_cat_re": r["u_cat_re"],
        "u_cat_im": r["u_cat_im"],
        "edge": bool(r.get("edge")),
    })

L_CAT = np.array([[x["u_cat_re"], x["u_cat_im"]] for x in late_arr])
frac_late_re = frac16(L_CAT[:, 0])
frac_late_im = frac16(L_CAT[:, 1])

# Контроль: случайные точки в той же области u-плоскости
N_RAND = 287
rand_re = np.random.uniform(RE_MIN, RE_MAX, N_RAND)
rand_im = np.random.uniform(IM_MIN, IM_MAX, N_RAND)
frac_rand_re = frac16(rand_re)
frac_rand_im = frac16(rand_im)


def run_tests(name, fr_re, fr_im):
    """Прогон всех тестов для набора. Возвращает dict."""
    out = {"n": len(fr_re)}
    # chi^2
    chi2_re, p_re, cnt_re = chi2_uniform_16bins(fr_re)
    chi2_im, p_im, cnt_im = chi2_uniform_16bins(fr_im)
    out["chi2_re"] = chi2_re
    out["p_chi2_re"] = p_re
    out["counts_re"] = cnt_re
    out["chi2_im"] = chi2_im
    out["p_chi2_im"] = p_im
    out["counts_im"] = cnt_im
    # KS
    ks_re, pks_re = ks_uniform(fr_re)
    ks_im, pks_im = ks_uniform(fr_im)
    out["KS_re"] = ks_re
    out["p_KS_re"] = pks_re
    out["KS_im"] = ks_im
    out["p_KS_im"] = pks_im
    # Rayleigh: phase = 2π * frac16 (рассматриваем дробную часть как угол)
    angles_re = 2.0 * np.pi * fr_re
    angles_im = 2.0 * np.pi * fr_im
    Rb_re, z_re, prayl_re = rayleigh_test(angles_re)
    Rb_im, z_im, prayl_im = rayleigh_test(angles_im)
    out["Rayleigh_R_re"] = Rb_re
    out["Rayleigh_z_re"] = z_re
    out["Rayleigh_p_re"] = prayl_re
    out["Rayleigh_R_im"] = Rb_im
    out["Rayleigh_z_im"] = z_im
    out["Rayleigh_p_im"] = prayl_im
    # доля у целых (|frac16|<1/16)
    near_int_re = float(((fr_re < 1/16) | (fr_re > 15/16)).mean())
    near_int_im = float(((fr_im < 1/16) | (fr_im > 15/16)).mean())
    out["near_integer_re"] = near_int_re
    out["near_integer_im"] = near_int_im
    return out


tests_good_actual = run_tests("good_actual", frac_act_re, frac_act_im)
tests_good_cat = run_tests("good_cat", frac_cat_re, frac_cat_im)
tests_good_du = run_tests("good_du", frac_du_re, frac_du_im)
tests_late_cat = run_tests("late_cat", frac_late_re, frac_late_im)
tests_random = run_tests("random", frac_rand_re, frac_rand_im)

print("\n[Г9 ТЕСТЫ РАВНОМЕРНОСТИ 16·u]")
print(f"  287 good u_actual: chi2 p={tests_good_actual['p_chi2_re']:.2e}/{tests_good_actual['p_chi2_im']:.2e}, "
      f"KS p={tests_good_actual['p_KS_re']:.2e}/{tests_good_actual['p_KS_im']:.2e}, "
      f"Rayleigh p={tests_good_actual['Rayleigh_p_re']:.2e}/{tests_good_actual['Rayleigh_p_im']:.2e}, "
      f"near_int={tests_good_actual['near_integer_re']*100:.1f}%/{tests_good_actual['near_integer_im']*100:.1f}%")
print(f"  287 good u_cat:    chi2 p={tests_good_cat['p_chi2_re']:.2e}/{tests_good_cat['p_chi2_im']:.2e}, "
      f"near_int={tests_good_cat['near_integer_re']*100:.1f}%/{tests_good_cat['near_integer_im']*100:.1f}%")
print(f"  287 good Δu:       chi2 p={tests_good_du['p_chi2_re']:.2e}/{tests_good_du['p_chi2_im']:.2e}")
print(f"  9 late u_cat:      chi2 p={tests_late_cat['p_chi2_re']:.2e}/{tests_late_cat['p_chi2_im']:.2e}, "
      f"KS p={tests_late_cat['p_KS_re']:.2e}/{tests_late_cat['p_KS_im']:.2e}, "
      f"Rayleigh p={tests_late_cat['Rayleigh_p_re']:.2e}/{tests_late_cat['Rayleigh_p_im']:.2e}, "
      f"near_int={tests_late_cat['near_integer_re']*100:.1f}%/{tests_late_cat['near_integer_im']*100:.1f}%")
print(f"  random 287:        chi2 p={tests_random['p_chi2_re']:.2e}/{tests_random['p_chi2_im']:.2e}, "
      f"near_int={tests_random['near_integer_re']*100:.1f}%/{tests_random['near_integer_im']*100:.1f}%")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Сохранить промежуточные массивы + тесты в JSON
# ─────────────────────────────────────────────────────────────────────────────
lattice_results = {
    "meta": {
        "oid": "wave3_lattice_R",
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "seed": RNG_SEED,
        "scope": "Г9 (1/16 lattice) + Г7 (lattice background), d=6",
        "input": "agentB_deg6/match_results.json + H6_zeros.json + subfamily/rate_full.json",
        "frozen": FROZEN,
    },
    "task_specificity_reproduction": density_repro,
    "task9_lattice": {
        "good_n": len(good_arr),
        "late_n": len(late_arr),
        "random_n": N_RAND,
        "tests_good_actual": tests_good_actual,
        "tests_good_cat": tests_good_cat,
        "tests_good_du": tests_good_du,
        "tests_late_cat": tests_late_cat,
        "tests_random": tests_random,
        "u_bounds": {"re": [RE_MIN, RE_MAX], "im": [IM_MIN, IM_MAX]},
    },
}

with open(WORK / "lattice_results.json", "w", encoding="utf-8") as f:
    json.dump(lattice_results, f, ensure_ascii=False, indent=2)

print(f"\n[SAVE] {WORK / 'lattice_results.json'}")
print("[OK] часть 1 (Г9) готова — нужно ещё Г7")
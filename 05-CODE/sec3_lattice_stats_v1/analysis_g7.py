"""
Часть 2: ЗАДАЧА 2 (Г7) — решёточный фон для специфичности.

(a) Воспроизведено в analysis.py: median 1488.4×, 7 ненулевых.
(b) Фон A: K=1000 случайных наборов из 287 точек, равномерных в той же области
    u-плоскости. Для каждого считаем медианную невязку d_u (расстояние до
    ближайшего каталожного нуля H₆). P(median ≤ 1.13e-3).
(c) Фон B: точки на решётках ℤ+iℤ и (ℤ/2+iℤ/2) в той же области — та же
    статистика.
(d) Итог: во сколько раз реальные совпадения лучше фона, с p-значением.

NOT RIGOROUS — float.
"""
import json
import time
from pathlib import Path

import numpy as np
from scipy.stats import median_abs_deviation

WORK = Path(r".")
AGENTB = Path(r"../sec3_pilot_d6")

# Загрузить результаты части 1
with open(WORK / "lattice_results.json", "r", encoding="utf-8") as f:
    lattice_results = json.load(f)

with open(AGENTB / "match_results.json", "r", encoding="utf-8") as f:
    MR = json.load(f)
CATALOG = MR["catalog_results"]

# Реальные u_actual 287 good
good = [r for r in CATALOG if r["accepted"] and not r.get("edge")]
real_u = np.array([[r["u_actual_re"], r["u_actual_im"]] for r in good])
real_du = np.array([r["d_u"] for r in good])
REAL_MEDIAN_DU = float(np.median(real_du))
REAL_MEAN_DU = float(np.mean(real_du))
print(f"[REAL] 287 good: median d_u = {REAL_MEDIAN_DU:.4e}, mean = {REAL_MEAN_DU:.4e}")

# Каталожные нули (m=9 limit extrapolated) — куда смотрим расстояния
cat_zeros = np.array([[r["u_cat_re"], r["u_cat_im"]] for r in CATALOG])
RE_MIN = float(cat_zeros[:, 0].min())
RE_MAX = float(cat_zeros[:, 0].max())
IM_MIN = float(cat_zeros[:, 1].min())
IM_MAX = float(cat_zeros[:, 1].max())
print(f"[BOUNDS] re∈[{RE_MIN:.2f},{RE_MAX:.2f}], im∈[{IM_MIN:.2f},{IM_MAX:.2f}]")

# ─────────────────────────────────────────────────────────────────────────────
# Ближайший сосед по KDTree (быстро)
# ─────────────────────────────────────────────────────────────────────────────
from scipy.spatial import cKDTree

CAT_TREE = cKDTree(cat_zeros)


def nearest_distance(points):
    """Для каждой точки массива (N,2) вернуть расстояние до ближайшего нуля H₆."""
    d, _ = CAT_TREE.query(points, k=1)
    return d


# ─────────────────────────────────────────────────────────────────────────────
# (b) Случайный фон: K=1000 наборов по 287 точек
# ─────────────────────────────────────────────────────────────────────────────
K_RAND = 1000
N = 287
np.random.seed(20260815)

t0 = time.time()
rand_medians = np.empty(K_RAND)
for k in range(K_RAND):
    pts = np.column_stack([
        np.random.uniform(RE_MIN, RE_MAX, N),
        np.random.uniform(IM_MIN, IM_MAX, N),
    ])
    d = nearest_distance(pts)
    rand_medians[k] = float(np.median(d))
t_rand = time.time() - t0

P_RAND_LE_REAL = float((rand_medians <= REAL_MEDIAN_DU).mean())
RATIO_RAND = float(rand_medians.mean() / REAL_MEDIAN_DU)
print(f"\n[BG RANDOM] K={K_RAND}, N={N}, {t_rand:.1f}s")
print(f"  mean median = {rand_medians.mean():.4e}, std = {rand_medians.std():.4e}")
print(f"  P(median ≤ {REAL_MEDIAN_DU:.2e}) = {P_RAND_LE_REAL:.4e}")
print(f"  ratio (random/real) = {RATIO_RAND:.1f}×")
print(f"  min median = {rand_medians.min():.4e}")

# ─────────────────────────────────────────────────────────────────────────────
# (c) Решёточный фон: ℤ+iℤ и (ℤ/2+iℤ/2) в той же области
# ─────────────────────────────────────────────────────────────────────────────
def lattice_points(step, re_min, re_max, im_min, im_max):
    re_first = np.ceil(re_min / step) * step
    re_last = np.floor(re_max / step) * step
    re = np.arange(re_first, re_last + 0.5 * step, step)
    im_first = np.ceil(im_min / step) * step
    im_last = np.floor(im_max / step) * step
    im = np.arange(im_first, im_last + 0.5 * step, step)
    R, I = np.meshgrid(re, im, indexing="ij")
    return np.column_stack([R.ravel(), I.ravel()])


lat_Z = lattice_points(1.0, RE_MIN, RE_MAX, IM_MIN, IM_MAX)
lat_half = lattice_points(0.5, RE_MIN, RE_MAX, IM_MIN, IM_MAX)
print(f"\n[LATTICE] ℤ+iℤ in area: {len(lat_Z)} points")
print(f"[LATTICE] ℤ/2+iℤ/2 in area: {len(lat_half)} points")


def sample_with_replacement(pop, n, rng):
    """Равномерная выборка с возвратом n точек из pop."""
    idx = rng.integers(0, len(pop), n)
    return pop[idx]


def lattice_medians(pop, n=287, K=1000, seed=20260816):
    rng = np.random.default_rng(seed)
    meds = np.empty(K)
    for k in range(K):
        pts = sample_with_replacement(pop, n, rng)
        meds[k] = float(np.median(nearest_distance(pts)))
    return meds


t0 = time.time()
lat_Z_meds = lattice_medians(lat_Z, n=287, K=1000, seed=20260817)
t_Z = time.time() - t0
P_Z_LE_REAL = float((lat_Z_meds <= REAL_MEDIAN_DU).mean())
RATIO_Z = float(lat_Z_meds.mean() / REAL_MEDIAN_DU)
print(f"\n[BG LATTICE Z] K=1000, {t_Z:.1f}s")
print(f"  mean median = {lat_Z_meds.mean():.4e}")
print(f"  P(median ≤ {REAL_MEDIAN_DU:.2e}) = {P_Z_LE_REAL:.4e}")
print(f"  ratio (Z-lattice/real) = {RATIO_Z:.1f}×")

t0 = time.time()
lat_half_meds = lattice_medians(lat_half, n=287, K=1000, seed=20260818)
t_half = time.time() - t0
P_HALF_LE_REAL = float((lat_half_meds <= REAL_MEDIAN_DU).mean())
RATIO_HALF = float(lat_half_meds.mean() / REAL_MEDIAN_DU)
print(f"\n[BG LATTICE Z/2] K=1000, {t_half:.1f}s")
print(f"  mean median = {lat_half_meds.mean():.4e}")
print(f"  P(median ≤ {REAL_MEDIAN_DU:.2e}) = {P_HALF_LE_REAL:.4e}")
print(f"  ratio (Z/2-lattice/real) = {RATIO_HALF:.1f}×")

# ─────────────────────────────────────────────────────────────────────────────
# (d) Дополнительно: моделирование «density probe» для фона — медианная
# «специфичность» 287 точек, размещённых случайно/на решётке.
# Берём для каждой точки 12 rival-стартов в 0.4-боксе, считаем median dist.
# Если rivals пересекают каталожный нуль — specificity = median_rival/d_u.
# Если ни один не ложится на нуль — specificity = None (как в данных).
# ─────────────────────────────────────────────────────────────────────────────
DENS_SIDE = 0.4
DENS_GRID = 12


def density_probe_for_set(points, seed=20260819):
    """Для каждой точки в points делаем density probe (12 стартов в 0.4-боксе),
    считаем median_rival_u — медианное расстояние от rivals до ближайшего
    нуля H₆. Возвращаем массив median_rival_u (None для проб без нулей)."""
    rng = np.random.default_rng(seed)
    out = np.empty(len(points))
    out[:] = np.nan
    # сетка 12 стартов в 0.4-боксе, та же что в пилоте (DENS_GRID=12)
    # Пилот использовал 12 «джиттерных стартов» — это набор из 12 точек в
    # квадрате 0.4×0.4. Простейшая дискретизация — равномерная 4×3 или 3×4.
    side = DENS_SIDE
    # 4 × 3 = 12 точек внутри бокса (со смещением от центра)
    offsets = []
    n_per_side = int(np.sqrt(DENS_GRID))
    # 4×3 = 12
    for i in range(4):
        for j in range(3):
            x = -side/2 + (i + 0.5) * side / 4
            y = -side/2 + (j + 0.5) * side / 3
            offsets.append((x, y))
    offsets = np.array(offsets[:DENS_GRID])

    for i, p in enumerate(points):
        rivals = p + offsets
        # расстояние каждого rival до ближайшего каталожного нуля
        d = nearest_distance(rivals)
        # если ВСЕ rivals очень далеко (>0.5) — эквивалент «n_rival=0»
        if d.min() > 0.5:
            out[i] = np.nan
        else:
            out[i] = float(np.median(d))
    return out


# Это для 287 точек, не 1000 итераций — слишком медленно. Сделаем K=200.
K_DENS = 200
t0 = time.time()
rand_spec_medians = np.empty(K_DENS)
rand_spec_nonnan_count = np.zeros(K_DENS, dtype=int)
for k in range(K_DENS):
    pts = np.column_stack([
        np.random.uniform(RE_MIN, RE_MAX, 287),
        np.random.uniform(IM_MIN, IM_MAX, 287),
    ])
    med_riv = density_probe_for_set(pts, seed=20260819 + k)
    valid = med_riv[~np.isnan(med_riv)]
    rand_spec_medians[k] = float(np.median(valid)) if len(valid) > 0 else np.nan
    rand_spec_nonnan_count[k] = int(len(valid))
t_dens_rand = time.time() - t0

valid_rand_spec = rand_spec_medians[~np.isnan(rand_spec_medians)]
print(f"\n[DENS PROBE RANDOM BG] K={K_DENS}, {t_dens_rand:.1f}s")
print(f"  median(spec) over iters: median={np.nanmedian(rand_spec_medians):.4e}, mean={np.nanmean(rand_spec_medians):.4e}")
print(f"  iters with ≥1 valid: {(~np.isnan(rand_spec_medians)).sum()}/{K_DENS}")

# Реальная «специфичность» из данных
real_spec_vals = np.array([v for v in [r["specificity"] for r in MR["density"].values()] if v is not None])
REAL_MEDIAN_SPEC = float(np.median(real_spec_vals))
print(f"  real median specificity = {REAL_MEDIAN_SPEC:.1f}×")

# ─────────────────────────────────────────────────────────────────────────────
# Сохранить результаты
# ─────────────────────────────────────────────────────────────────────────────
g7_results = {
    "task7_background": {
        "real_median_d_u": REAL_MEDIAN_DU,
        "real_mean_d_u": REAL_MEAN_DU,
        "n_real": 287,
        "bounds": {"re": [RE_MIN, RE_MAX], "im": [IM_MIN, IM_MAX]},
        "background_random": {
            "K": K_RAND,
            "n_per_set": N,
            "mean_median_d_u": float(rand_medians.mean()),
            "std_median_d_u": float(rand_medians.std()),
            "min_median_d_u": float(rand_medians.min()),
            "max_median_d_u": float(rand_medians.max()),
            "P_median_le_real": P_RAND_LE_REAL,
            "ratio_to_real": RATIO_RAND,
            "elapsed_sec": t_rand,
        },
        "background_lattice_Z": {
            "step": 1.0,
            "n_in_area": int(len(lat_Z)),
            "K": K_RAND,
            "mean_median_d_u": float(lat_Z_meds.mean()),
            "P_median_le_real": P_Z_LE_REAL,
            "ratio_to_real": RATIO_Z,
        },
        "background_lattice_Z2": {
            "step": 0.5,
            "n_in_area": int(len(lat_half)),
            "K": K_RAND,
            "mean_median_d_u": float(lat_half_meds.mean()),
            "P_median_le_real": P_HALF_LE_REAL,
            "ratio_to_real": RATIO_HALF,
        },
        "density_probe_simulation": {
            "K": K_DENS,
            "real_median_specificity": REAL_MEDIAN_SPEC,
            "random_bg_median_spec": float(np.nanmedian(rand_spec_medians)),
            "random_bg_mean_spec": float(np.nanmean(rand_spec_medians)),
        },
    }
}

lattice_results["task7_background"] = g7_results["task7_background"]

with open(WORK / "lattice_results.json", "w", encoding="utf-8") as f:
    json.dump(lattice_results, f, ensure_ascii=False, indent=2)

print(f"\n[SAVE] {WORK / 'lattice_results.json'} обновлён")
"""
Картинки для отчёта: распределения дробных частей 16·u (Г9) и распределения
фоновой медианной невязки (Г7).
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WORK = Path(r".")
AGENTB = Path(r"../sec3_pilot_d6")

with open(WORK / "lattice_results.json", "r", encoding="utf-8") as f:
    RES = json.load(f)

with open(AGENTB / "match_results.json", "r", encoding="utf-8") as f:
    MR = json.load(f)
CATALOG = MR["catalog_results"]

# Reproduce frac16 etc.
def frac16(x):
    return (16.0 * x) - np.floor(16.0 * x)


np.random.seed(20260815)

good = [r for r in CATALOG if r["accepted"] and not r.get("edge")]
late_ids = ["S036", "S104", "S187", "S339", "S031", "S070", "S235", "S279", "S319"]
late = [r for r in CATALOG if r["id"] in late_ids]

G_act = np.array([[r["u_actual_re"], r["u_actual_im"]] for r in good])
G_cat = np.array([[r["u_cat_re"], r["u_cat_im"]] for r in good])
L_cat = np.array([[r["u_cat_re"], r["u_cat_im"]] for r in late])

cat_zeros = np.array([[r["u_cat_re"], r["u_cat_im"]] for r in CATALOG])
RE_MIN, RE_MAX = float(cat_zeros[:, 0].min()), float(cat_zeros[:, 0].max())
IM_MIN, IM_MAX = float(cat_zeros[:, 1].min()), float(cat_zeros[:, 1].max())

rand_re = np.random.uniform(RE_MIN, RE_MAX, 287)
rand_im = np.random.uniform(IM_MIN, IM_MAX, 287)


# ─── Рис. 1: дробные части 16·u для четырёх выборок ──────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)

datasets = [
    ("287 good (u_actual)", frac16(G_act[:, 0]), frac16(G_act[:, 1]), axes[0, 0]),
    ("287 good (u_catalog)", frac16(G_cat[:, 0]), frac16(G_cat[:, 1]), axes[0, 1]),
    ("9 late (u_catalog)", frac16(L_cat[:, 0]), frac16(L_cat[:, 1]), axes[1, 0]),
    ("287 random (uniform)", frac16(rand_re), frac16(rand_im), axes[1, 1]),
]

bins = np.linspace(0, 1, 17)
for title, fr_re, fr_im, ax in datasets:
    counts_re, _, _ = ax.hist(fr_re, bins=bins, alpha=0.55, label="Re", color="C0", edgecolor="k")
    counts_im, _, _ = ax.hist(fr_im, bins=bins, alpha=0.55, label="Im", color="C3", edgecolor="k")
    # ожидание равномерного
    ax.axhline(len(fr_re) / 16, ls="--", color="k", lw=1.0, label="uniform exp.")
    ax.set_title(title)
    ax.set_xlabel("дробная часть 16·u")
    ax.set_ylabel("counts")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.3)

fig.suptitle("Г9: проверка решётки 1/16 в u-плоскости (n=16 bins, χ² на равномерность)",
             fontsize=12)
fig.tight_layout()
fig.savefig(WORK / "fig_g9_lattice16.png", dpi=120)
plt.close(fig)
print("[SAVE] fig_g9_lattice16.png")

# ─── Рис. 2: scatter u-плоскости — все 342 каталожных + 9 late выделены ──
fig, ax = plt.subplots(figsize=(11, 8))
ax.scatter(cat_zeros[:, 0], cat_zeros[:, 1], s=8, c="lightgray", label="каталог H₆ (342)")
ax.scatter(G_cat[:, 0], G_cat[:, 1], s=14, c="C0", label="287 good (matched)", alpha=0.7)
ax.scatter(L_cat[:, 0], L_cat[:, 1], s=70, c="red", marker="*",
           label="9 late (near 1/16 lattice)", edgecolor="k", linewidths=0.5, zorder=10)
# нарисовать решётку 1/16 в области
step = 1/16
re_lat = np.arange(np.ceil(RE_MIN/step)*step, np.floor(RE_MAX/step)*step + step/2, step)
im_lat = np.arange(np.ceil(IM_MIN/step)*step, np.floor(IM_MAX/step)*step + step/2, step)
for r in re_lat:
    ax.axvline(r, color="gray", alpha=0.15, lw=0.4)
for i in im_lat:
    ax.axhline(i, color="gray", alpha=0.15, lw=0.4)
ax.set_xlabel("Re(u)")
ax.set_ylabel("Im(u)")
ax.set_title("Каталог H₆: 287 хороших (синие) и 9 поздних (красные звёзды) на решётке 1/16")
ax.set_aspect("equal")
ax.legend(loc="upper left", fontsize=9)
ax.grid(alpha=0.2)
fig.tight_layout()
fig.savefig(WORK / "fig_g9_uplane.png", dpi=120)
plt.close(fig)
print("[SAVE] fig_g9_uplane.png")

# ─── Рис. 3: распределение расстояний до ближайшего каталожного нуля ────
# Реальные 287, случайный фон (1000 итераций), решётка Z (1000), решётка Z/2 (1000)
from scipy.spatial import cKDTree
TREE = cKDTree(cat_zeros)

# уже в JSON есть агрегаты; для гистограмм воспроизведём
K = 1000
N = 287
real_d = np.linalg.norm(G_act - G_cat, axis=1)  # = d_u
# random
rand_meds = []
for k in range(K):
    pts = np.column_stack([
        np.random.uniform(RE_MIN, RE_MAX, N),
        np.random.uniform(IM_MIN, IM_MAX, N),
    ])
    d, _ = TREE.query(pts, k=1)
    rand_meds.append(np.median(d))
rand_meds = np.array(rand_meds)

# lattice Z
def lattice_points(step, re_min, re_max, im_min, im_max):
    re = np.arange(np.ceil(re_min/step)*step, np.floor(re_max/step)*step + step/2, step)
    im = np.arange(np.ceil(im_min/step)*step, np.floor(im_max/step)*step + step/2, step)
    R, I = np.meshgrid(re, im, indexing="ij")
    return np.column_stack([R.ravel(), I.ravel()])

lat_Z = lattice_points(1.0, RE_MIN, RE_MAX, IM_MIN, IM_MAX)
rng = np.random.default_rng(20260817)
Z_meds = []
for k in range(K):
    idx = rng.integers(0, len(lat_Z), N)
    pts = lat_Z[idx]
    d, _ = TREE.query(pts, k=1)
    Z_meds.append(np.median(d))
Z_meds = np.array(Z_meds)

lat_half = lattice_points(0.5, RE_MIN, RE_MAX, IM_MIN, IM_MAX)
rng = np.random.default_rng(20260818)
half_meds = []
for k in range(K):
    idx = rng.integers(0, len(lat_half), N)
    pts = lat_half[idx]
    d, _ = TREE.query(pts, k=1)
    half_meds.append(np.median(d))
half_meds = np.array(half_meds)

real_med = float(np.median(real_d))

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
# Panel A: распределение d_u для реальных 287
axes[0].hist(real_d, bins=40, color="C0", edgecolor="k", alpha=0.8)
axes[0].axvline(real_med, color="red", ls="--", lw=1.5, label=f"median = {real_med:.2e}")
axes[0].set_xlabel("d_u (distance to nearest catalog zero)")
axes[0].set_ylabel("counts (out of 287)")
axes[0].set_title("Реальные 287: d_u")
axes[0].legend()
axes[0].grid(alpha=0.3)

# Panel B: медианные d_u для K=1000 фонов
bins_b = np.linspace(0, 6, 60)
axes[1].hist(rand_meds, bins=bins_b, color="C1", edgecolor="k", alpha=0.55, label="random uniform")
axes[1].hist(Z_meds, bins=bins_b, color="C2", edgecolor="k", alpha=0.55, label="lattice ℤ+iℤ")
axes[1].hist(half_meds, bins=bins_b, color="C4", edgecolor="k", alpha=0.55, label="lattice ℤ/2+iℤ/2")
axes[1].axvline(real_med, color="red", ls="--", lw=2.0,
                label=f"real median = {real_med:.2e}")
axes[1].set_xlabel("median d_u по 287 точкам")
axes[1].set_ylabel("counts (out of K=1000)")
axes[1].set_title("Фон: медиана d_u (лог-шкала по x)")
axes[1].set_xlim(0, 6)
axes[1].legend(fontsize=8)
axes[1].grid(alpha=0.3)

# Panel C: same but log y
axes[2].hist(rand_meds, bins=bins_b, color="C1", edgecolor="k", alpha=0.55, label="random uniform")
axes[2].hist(Z_meds, bins=bins_b, color="C2", edgecolor="k", alpha=0.55, label="lattice ℤ+iℤ")
axes[2].hist(half_meds, bins=bins_b, color="C4", edgecolor="k", alpha=0.55, label="lattice ℤ/2+iℤ/2")
axes[2].axvline(real_med, color="red", ls="--", lw=2.0,
                label=f"real median = {real_med:.2e}")
axes[2].set_xlabel("median d_u по 287 точкам")
axes[2].set_ylabel("counts (out of K=1000)")
axes[2].set_title("Фон: медиана d_u (log y)")
axes[2].set_yscale("log")
axes[2].set_xlim(0, 6)
axes[2].legend(fontsize=8)
axes[2].grid(alpha=0.3)

fig.suptitle("Г7: реальные совпадения (287 d_u << 10⁻³) против случайного и решёточного фона",
             fontsize=12)
fig.tight_layout()
fig.savefig(WORK / "fig_g7_background.png", dpi=120)
plt.close(fig)
print("[SAVE] fig_g7_background.png")

print("[OK] all 3 figures saved")
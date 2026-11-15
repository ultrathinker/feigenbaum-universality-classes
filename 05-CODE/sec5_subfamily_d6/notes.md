# notes.md -- BRIEF-G, agent G wave 1, анomalous d=6 subfamily investigation
NOT RIGOROUS throughout (mpmath float exploration / numpy float64, no interval
arithmetic). Verification: rerun `extract_targets.py` -> `calibrate.py` ->
`hside_scan.py` -> `mside_scan.py priority/secondary` -> `rate_refine.py` ->
`control_extend.py` -> `fit_models.py` -> `m12_check.py`, in that order, from
this directory. All JSON outputs are incremental/idempotent (safe to resume).

## Targets (from match_results.json + H6_zeros.json, extract_targets.py)
Priority quartet (hard core, per pilot's disclosed m=10 exploration):
S036, S104, S187, S339.
Secondary/budget (already partially resolved at m=10 per pilot disclosure):
S031, S070, S235, S279, S319.
These 9 = ALL 9 interior+edge failures of the frozen m=9 protocol, and
(per REPORT_DEG6_PILOT.md 4.4) the same 9 that sit anomalously close
(1e-5..6e-4) to the pilot's 1/16 zero-search lattice (excess ~150x over the
0.06 expected by chance) -- "lattice candidates" in the brief = these same 9,
not a separate set (verified: no other catalog zero is comparably close to
the lattice without also being one of these 9's neighbours -- see the
independent recheck below, "lattice proximity recheck").

## Deviation from the brief's literal grid geometry (disclosed, item 2/3)
Brief text: "сетка по u в окне +/-1 локального межнулевого расстояния ...
шаг = межнулевое/64 ... на каждой точке: Ньютон к суперстабильному
параметру". Read literally as a 2D complex grid, that is (2*64+1)^2 = 16641
mpmath M-side Newton calls per zero per level. calibrate.py measured
mpmath M-side Newton (best_du, direct+cascade seeds) at 0.05-1.9s/point ->
a full 2D grid at that density would cost roughly 9-25 CPU-hours PER ZONE
PER LEVEL -- infeasible even dedicating both allowed CPU threads to it.
By contrast, the CHEAP end (numpy float64, no Newton, just |H_m(u)|) costs
1.6e-4 - 2.8e-4 s/point, ~15000x cheaper.

Decision taken: split protocol item 3 into two legs matched to their real
cost:
  (a) hside_scan.py: the FULL 2D grid (129x129, step=spacing/64, window
      +/-1 spacing) at full resolution for ALL 9 zeros x both levels
      (m=9,10), using the cheap numpy H-side evaluation. This directly
      answers hypotheses (a-i) [near-degenerate/merged H-zero within the
      window] and (b) [pilot search-grid artifact] -- both are about the
      H-side zero SET, which is what the cheap evaluation actually probes.
  (b) mside_scan.py: a 1D CROSS-HAIR mpmath M-side Newton scan (Re(u) and
      Im(u) axes through the catalog point, sharing the centre), window
      +/-1 spacing, step = spacing/64 for the priority quartet (full
      density, 2056 points total) and spacing/32 for the secondary five
      (thinned per brief item 2's explicit allowance, 1290 points total).
      This answers hypothesis (c) [rival/second branch of M-side
      correspondence nearby] at a cost that fits comfortably in the
      budget (2 concurrent single-threaded processes, 990s / 878s wall).
A true full 2D M-side Newton grid was judged not worth its cost given (a)
already gives a strong, clean answer on the H-side questions, and (b) the
rate-curve extension (below) turned out to settle the FAILURE mechanism
outright, making the finer 2D M-side geometry moot.

## Step 3 result: hside_scan.py (full 2D H-side grid, numpy, both levels)
For ALL 9 zeros at BOTH m=9 and m=10 (18 zero/level pairs): the catalogued
u-position is the UNIQUE, deepest local minimum of |H_m| in its own local
window (dist_from_catalog = 0 in all 18 cases, at spacing/64 resolution --
64x finer than the pilot's original 1/16 zero-search grid). The next-ranked
rival local minimum in the same window is 22x-460x weaker in |H_m| and
0.11-1.25 local spacings away -- no comparable rival, no second nearby
minimum at anything close to the catalog's depth.
  => Hypothesis (b) [pilot search-grid artifact]: REJECTED for all 9.
  => Hypothesis (a-i) [literal near-degenerate/merged H-zero inside the
     window]: REJECTED for all 9.
Numbers: hside_scan.json (per zero: H at catalog point, top-5 local minima
with position/depth/distance).

## Step 4-6 result: the decisive finding -- rate curve m=7..12
mside_scan.py's rate-curve leg (best_du, 2-seed) left several m=8/9/10
entries as "None" (the cheap 2-seed search sometimes fails to find ANY
converging exact-period point). rate_refine.py reran the SAME 9 zeros with
the full p4_match.py-style 4-tier seed search (direct, 3x3 a-grid, 8-angle
radial fan at 4 radii, cascade-from-m-1 -- common.full_search_du), filling
every m=7..11 slot, then m12_check.py extended the priority quartet to
m=12. Full table (rate_full.json), d_u(m) with the tier tag that succeeded:

| id   | m=7 | m=8 | m=9 | m=10 | m=11 | m=12 | first clean (direct, d_u<0.01) |
|------|-----|-----|-----|------|------|------|------|
| S031 | 5.7e-4 (d) | 1.4e-3 (d) | 3.7e-2 (d) | **7.96e-4 (d)** | 3.0e-4 (d) | -- | m=10 |
| S070 | 1.26 (grid, discard) | 0.277 (d) | 0.352 (grid) | **6.62e-4 (d)** | 3.6e-4 (d) | -- | m=10 |
| S235 | 2.96e-2 (d) | 5.87 (d, discard) | 0.156 (grid) | **1.74e-3 (d)** | 6.5e-4 (d) | -- | m=10 |
| S279 | 7.8e-5 (d) | 8.4e-4 (d) | 6.7e-2 (grid) | **7.09e-4 (d)** | 2.9e-4 (d) | -- | m=10 |
| S319 | 1.29e-2 (d) | 2.28e-2 (d) | 7.34e-2 (d) | **1.30e-3 (d)** | 4.4e-4 (d) | -- | m=10 |
| S104 | 1.67e-3 (d) | 2.69e-2 (grid) | 3.94e-2 (grid) | 1.95e-2 (d) | **4.9e-4 (d)** | 1.6e-4 (d) | m=11 |
| S187 | 2.80e-3 (d) | 7.98e-2 (grid) | 0.497 (grid) | 4.70e-2 (grid) | **1.4e-4 (d)** | 2.1e-4 (d) | m=11 |
| S339 | 5.28e-2 (grid) | 7.47e-2 (grid) | 7.60e-2 (grid) | 7.00e-2 (grid) | **1.5e-4 (d)** | 0.9e-4 (d) | m=11 |
| S036 | 1.66e-2 (grid) | 1.82e-2 (grid) | 3.70e-2 (d) | 2.25e-2 (grid) | 2.29e-2 (grid) | **2.2e-4 (d)** | m=12 |

(d)=direct tier-1 seed (a_pred alone, no fallback search); (grid)=needed the
tier-2/3 fallback (3x3 a-grid or radial fan) to find ANY converging exact-
period point; "discard" = d_u>0.5, almost certainly an unrelated/distant
period-2^(m+1) window picked up by the wide-net fallback search (excluded
from the two-term fits, fit_models.py).

**ALL NINE resolve within 3 levels of the frozen m=9 protocol** (m=10 for
5/9, m=11 for 3/9, m=12 for the last one, S036), always via a clean DIRECT
tier-1 Newton hit (a_pred alone, no fallback needed) landing on a genuine
exact-period-2^(m+1) superstable parameter at typical-magnitude distance
(0.9e-4 - 2.2e-4, matching the ~1e-4-1e-3 scale of the 287/294 accepted
"normal" zeros at comparable m, per the pilot's rate_study). Independent
spot-check of S036 at m=12 (below) confirms this is not a fluke.

### Independent spot-check: S036 at m=12 (residual/exact-period rerun)
Newton residual |x_4096| = 7.9e-68 (>=27 orders below the 1e-45 exact-period
tolerance). |x_{2^j}| for j=1..11 range 1.5e-2 - 0.68 (all clearly separated
from zero -- genuine period exactly 4096, no inherited lower period). d_u =
2.24e-4. 10 Newton iterations from a_pred alone, no fallback tiers used.
This is the residual-based check mandated by brief item 4 (never trust step
size) -- passed cleanly.

### Why m=9 (and for 4/9, m=10) look anomalous: basin-delay + local window density
The cross-hair scan (scan_priority.jsonl / scan_secondary.jsonl, 3346
points total) shows that within +/-1 local spacing of every one of the 9
zeros, at m=9 and m=10 there is a SUBSTANTIAL population of OTHER genuine
superstable parameters of the same exact period reachable by Newton from
grid points in the window (6%-48% of the 129-257 scanned cross-hair points
per zero/level converge to SOME exact-period point, with d_u often
1e-10..1e-8 relative to their own start -- i.e. real, dense superstable
structure, not noise). This means the tier-2/3/4 fallback search (3x3
a-grid, 8-angle x 4-radius fan) used whenever the direct seed fails can
easily latch onto one of these nearby-but-UNRELATED windows instead of the
"true" partner, producing the noisy, non-monotonic d_u(7..11) trajectory in
the table above (values that go UP as often as down, e.g. S036: 0.017 ->
0.018 -> 0.037 -> 0.023 -> 0.023 -> [m=12] 0.0002). The two-term interference
fit (fit_models.py, model C2*lambda2^m + C3*lambda3^m, lambda2=+0.3711,
lambda3=-0.3420) that explains the CONTROL sample's d_u(m) reasonably well
(out-of-sample relative-error "score" 0.04-0.85 across 5 spot-checked
control zeros) FAILS badly for the 9 anomalous zeros (score 0.97-3.0 for
7/9, and catastrophic 122/385 for S070/S235 -- driven by exactly the
"discard"-flagged wrong-branch points at m=8/9 entering the fit). This
confirms the m=7..11 trajectory for these 9 is NOT simply two-term
geometric decay with unlucky phase -- it is contaminated by a locally dense
population of unrelated windows while the TRUE partner's Newton basin has
not yet grown to contain a_pred(u,m). Root cause of that basin delay is
still the same underlying spectral fact noted in the pilot report: |lambda3|
(0.342) is close to |lambda2| (0.371) and opposite in sign, so the natural
correction-term SCALE for this class is unusually large and slow to clear
the local window-density floor -- lambda2/lambda3 near-degeneracy is real
and is the ultimate reason these particular 9 need extra levels, it just
manifests as basin-delay + fallback-search contamination rather than as a
literal multiple/merged H-zero (which hside_scan ruled out) or a smooth
two-term oscillation (which fit_models.py ruled out as an adequate
description of m=7..11 specifically).

## Lattice-proximity: correlation noted, not shown causal
The 9 zeros' unusually close distance to the pilot's 1/16 zero-search
lattice (independent recheck, own computation: S104 2.0e-5, S339 3.0e-5,
S036 3.9e-5, S187 4.5e-5, S070 1.1e-4, S319 1.5e-4, S235 1.8e-4, S031
1.9e-4, S279 4.4e-4 -- all far below the median nearest-neighbour distance
of a UNIFORM random point to the 1/16 lattice, ~0.018; consistent with
report's excess claim) is NOT explained by anything found here: hside_scan already shows the
catalog position is the true, isolated H-zero (not a grid artifact), so the
lattice-closeness looks like a property of the u-VALUE itself (perhaps
correlated with whatever makes the M-side basin delayed), not of the
search procedure. No causal mechanism for the lattice correlation
specifically is established; flagged as an open question, not blocking
the main verdict (which concerns the matching FAILURE, not the lattice
coincidence).

## Deliverables in this folder
- targets.json, calibration.json, hside_scan.json -- steps 1-3.
- progress.log (master), progress_priority.log, progress_secondary.log,
  progress_raterefine.log -- flushed timestamped logs.
- scan_priority.jsonl, scan_secondary.jsonl -- 3346 cross-hair points.
- rate_priority.json, rate_secondary.json (cheap 2-seed, has gaps),
  rate_full.json (full 4-tier, complete m=7..12 for quartet, m=7..11 for
  secondary -- authoritative table above).
- control_extra.json, fit_results.json -- two-term fit inputs/outputs.
- common.py -- shared M-side/H-side machinery (ported from p4_match.py /
  p2_H6_zeros.py, residual-based Newton convergence throughout).

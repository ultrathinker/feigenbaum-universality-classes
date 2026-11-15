# Explicit realisation of the "parameter ↔ unstable manifold" correspondence for period-doubling universality classes

Self-contained package of research materials: the article, two appendices, references, and a replication kit with code and data under a sha256 manifest. All definitions, all numbers entering the conclusions, and all steps of the proofs are contained within these files. External links point either elsewhere within the package or to published literature.

## Contents

| file | contents |
|---|---|
| `01-ARTICLE.md` | Main article (~70 KB) — "Explicit realisation of the 'parameter ↔ unstable manifold' correspondence for period-doubling universality classes" |
| `02-APPENDIX-A-LEMMA.md` | Appendix A. Proposition on the gauge tower for the class \|x\|^d |
| `03-APPENDIX-B-DATA.md` | Appendix B. Data and constants for the three classes (full tabular appendix-dataset, ~189 KB; all tables are direct reads of JSON from `05-CODE/`; B.1.0 — external inputs; B.2.2.1 — layer 2 at d = 6; B.6.1 — the `unique222` row for d = 2) |
| `04-REFERENCES.md` | Reference list (verified bibliography with DOI) |
| `05-SUPPLEMENT-METHODS.md` | Supplement: methodological notes S1–S4, lattice forensics S5, layer 2 S6, Bridge-2 S7, OOS/LOOCV S8, operator grid S9 |
| `05-CODE/` | Replication kit: 35 subfolders with independent implementations, certificates, reports, tamper batteries; index — `README-CODE.md`; path-relink journal — `RELINK-LOG.md` |
| `06-EXTERNAL-INPUTS/` | External inputs for the d = 2 class (copies of artefacts from sibling projects; described in `README-EXTERNAL.md`) |
| `MANIFEST.sha256` | sha256 manifest of the entire package (LF: `<sha256><spaces/tabs><relpath>`, sorted by path) |
| `_tools/verify_manifest.py` | cross-platform manifest verifier |

## What is and is not proved (summary by the rigor scale of §2.3 and the "Dependency map" at the end of `01-ARTICLE.md`)

**RIGOROUS (interval Krawczyk boxes, computer-assisted proofs):**

- **357 Krawczyk boxes:** 61 at d = 4 (m = 9, period 1024), 287 at d = 6 (m = 9, period 1024), 9 late ones at d = 6 (m = 10/11/12, periods 2048/4096/8192). For each box — existence, uniqueness and exact period (Krawczyk + interval Newton, `pass=True` for all 357 boxes).
- **Layer 2 (S078, S090) at d = 6 — re-certification** with radius r = 1·10⁻³⁴ at prec 512: the new certificate passes in both orders of computing powers (`x⁶ = x**6` and `x⁶ = x⁵·x`) with relative width of the F′(X) ball ≤ 1.5·10⁻³ and κ < 2·10⁻³ (see `05-CODE/sec3_krawczyk_d6_layer2_recheck/README.md`). No additional refinement of the centre or mean-value form for F′ is required.

**Numerically verified (certified eigenvalue of the collocation discretization; N → ∞ truncation error estimated by a ladder but not formally certified):**

- **λ₂(2) ± 1.82·10⁻³³ — numerically verified:** certified eigenvalue of the collocation discretization N = 48 (sibling project "Bridge").
- **λ₂(4), λ₂(6) — numerically verified:** certified eigenvalues of the collocation discretization N = 40 (d = 4) and N = 48 (d = 6); Chebyshev basis, ball radii ~10⁻⁴⁰; N → ∞ truncation estimated by a ladder (~10⁻²⁶ at the last step). The three quantities (λ₂(2), λ₂(4), λ₂(6)) share the same status "numerically verified" (§2.3, rigor status scale).

**Proposition on the gauge tower (Appendix A):** the formal identity DT·Φ_k = α⁻ᵈᵏ·Φ_k on the odd real branch when dk ∈ 2ℤ ∩ ℝ₊ (and pointwise on the odd branch for any real dk > 0). The identity is formal for the fixed point g_d; existence at even d — Epstein (1986), uniqueness and hyperbolicity — Lyubich (1999) and Avila–Lyubich (2011); d = 2 — Lanford (1982); d = 4 — CAP bounds of Burbanks–Osbaldestin–Thurlby (2021). Numerical verification of the identity was performed on a collocation approximation of g_d (residual ≤ 2·10⁻²⁰ for integer even dk and ≤ 8·10⁻²⁴ for fractional/odd dk on the odd real branch; A.5.1, A.5.2). Attribution to the classics: see §A.6 of Appendix A.

**NOT RIGOROUS** — numerically supported but not certified (per the "Dependency map" at the end of `01-ARTICLE.md`):

- **d_u-layer (§3):** 22/23 (d = 2), 61/62 (d = 4), 287/294 + 296/296 (d = 6 at m = 9…12) — correspondence of zeros of H_d to catalogued superstable parameters; status — **NOT RIGOROUS (numerical correspondence; only the pair ends are certified)**, consistent with §3.2. The d_u-layer itself depends on uncertified rescaling constants (c_F, A, δ); the boxes themselves are certified (see above).
- **Completeness of catalogues:** not claimed for any of the three classes (146 vs 62 at d = 4; 641 vs 342 at d = 6, see §3.2, §3.3).
- **Lattice of step 1/16 — scanner-selection artefact** (§3.3, §7.2): established for d = 4 and d = 6 by a single protocol (rescan with an incommensurate grid 1/24); the values at the nodes are noise from candidate selection among local minima of |H| on the grid; the nodes themselves are not "near-zeros" of H.
- **Continuation of the B⁺ branch towards the tower mode Φ₁ near d = 4 and d = 6** (previously — "shadow"): symmetric, cosine of eigenfunctions = 1 − O(10⁻⁸) for Δd ∈ [0.0005, 0.01], N = 48/64, two bases (`05-CODE/sec2_tower_shadow_small_dd/`); the cosine of the Φ₁ and B⁺ coefficients in the monomial basis is a basis-dependent quantity, not evidence of asymmetry d = 4 / d = 6 (§2.5, Supplement §S1); NOT RIGOROUS.
- **Rate law rate = min(δ, 1/|λ₂|):** the evidence in our data — NOT RIGOROUS (float/mpmath); as a fact of the literature — conditionally RIGOROUS (Mao–Hu 1987; Reick 1992; Briggs 1994).
- **Spectral grid and crossover d\* = 2.3616 ± 0.0008** — NOT RIGOROUS; resolving the sign change of λ₂ near d = 2 — NOT RIGOROUS, but consistent across three independent trials.
- **Status of the pattern principle as such** — no formal theorem about H_d; sketch in §7.3.
- **OOS "within ×2" model** (§7.3): d = 6 — confirmed (~92 % vs baseline ~90.5 %); d = 4 — no structural superiority over baseline shown.

Full list of limitations — §7.4 of the article.

## Environment

- **Python 3.13.13**, invoked as `python -X utf8 <script>` from the script's folder.
- **python-flint 0.9.0** (Arb/FLINT) — for certificates of section 4 (`05-CODE/sec3_krawczyk_d4`, `05-CODE/sec3_krawczyk_d6`, `05-CODE/sec3_krawczyk_d6_late`, `05-CODE/sec3_krawczyk_d6_layer2_recheck`) and Bridge-2 (`05-CODE/sec4_bridge2_lambda2_cert`).
- **mpmath 1.3.0**, **numpy 2.5.0**, **scipy 1.18.0**.
- OS **Windows 11**. All scripts were run by the authors "in place"; when assembling the package, absolute paths pointing into the included folders were replaced by relative ones (journal `05-CODE/RELINK-LOG.md`); remaining external paths are listed there as well.

## How to reproduce

Checks are performed by running the scripts from `05-CODE/` (the full mapping of subfolders to article sections is in `README-CODE.md`). The list of commands by tables is below in the section "Reproduction commands by tables".

## Open questions

Consolidated list of open questions from §7.4 and §8 of the article (status as of today in parentheses):

1. **Completeness of the H_d zero catalogues** (§7.4, item 1) — not claimed for any of the three classes; quantitative bounds of incompleteness: a rescan with incommensurate step 1/24 finds 146 zeros at level 8 vs 62 (d = 4, ≈2.35×) and 641 vs 342 (d = 6, ≈1.87×). Catalogues are extendable but completeness is not achieved.
2. **Behaviour of the physical branch B⁺(d) near even d — analytically** (§7.4, item 5) — numerically the branch continues the tower mode Φ₁ symmetrically at d = 4 and d = 6 (cosine of eigenfunctions = 1 − O(10⁻⁸) for Δd ∈ [0.0005, 0.01]; NOT RIGOROUS, mpmath dps 40); the Proposition of §2.5 gives only an exact identity at even points, not an expansion of the branch in Δd. The cosine of the Φ₁ and B⁺ coefficients in the monomial basis is a basis-dependent quantity, not evidence of asymmetry d = 4 / d = 6.
3. **Tail N → ∞ for λ₂(4), λ₂(6)** (§7.4, item 7) — currently only a heuristic estimate from a ladder N = 16 → 24 → 32 → 40 (d = 4) and up to 48 (d = 6); truncation error ~10⁻²⁶ at the last step, but no formal certification of the limiting operator.
4. **Formal theorem about H_d** (§8) — the "zero ↔ parameter" correspondence remains an empirical NOT RIGOROUS layer; the sketch "paths to a theorem" in §7.3 contains four open subproblems (interval anchors c_F, A, δ; simplicity of zeros of H; rigorous version of the λ₂^m remainder; formal control of u → 2u). The work contains no formal statement about the zero set of H_d or about its link to parameter space in the nonlinear case.
5. **Mechanism of attraction of H_d zeros to the 1/16 lattice** (§7.4, item 8) — scanner-selection artefact established for d = 4 and d = 6; why the "valley" of small |H| coincides in scale with the step of the catalogue-building grids — open question.
6. **Rigorous identification of our collocation g_d with the exact fixed point of the operator T** (estimate of the N → ∞ truncation error) — not carried out; numerical verification of identity (2) was performed on the collocation approximation (A.5.1, A.5.2). Existence of the fixed point g_d at even d is covered by the classical theorems of Epstein (1986), Lyubich (1999) and Avila–Lyubich (2011).
7. **Novelty of the correspondence principle in the literature** (§7.1) — "probably new" with a flagged residual priority risk (kneading theory / MSS combinatorics, journal sources from the 1980s–1990s — not checked line by line).
8. **The law min(δ, 1/|λ₂|) as a formal theorem** (§8) — current confirmation in the data is NOT RIGOROUS; there is no proof of the min rule as a theorem in the paper independent of Briggs (1994).
9. **d_u-layer and specificity metric** — non-uniform between classes (per-zero for d = 2, per-rival for d = 4, 6); direct comparison of the numbers 82339× and 1488× is not correct.
10. **OOS "within ×2" model** (§7.4, item 6; §7.3) — for d = 4 at n = 61, no structural superiority over baseline by this criterion was shown.
11. **Possible extension to a fourth universality class or to classes with an odd maximum exponent** (§8) — a separate task; the paper does not anticipate it.

## Reproduction commands by tables

Commands are run from the corresponding subfolders or from the package root with absolute/relative paths as indicated. **Part III (new in the current package):** marked as "current-version package". All other commands were carried over from previous versions without changes to the paths; verified by `ls` against the contents of the package folders. Times — known from reports/logs; "n/a" — not in the available journals.

### Table T2.0 (definitions and constants, §1, §2 of the article)
- Script: `python -X utf8 05-CODE/sec5_lambda2_deg2/p1_deg2.py`
- Inputs: none (frozen constants — inside the script)
- Outputs: `05-CODE/sec5_lambda2_deg2/p1_deg2.json`
- Time: n/a

### Table T2.1 (constants λ₂(2), α, δ, §2.3)
- Scripts:
  - `python -X utf8 05-CODE/sec5_lambda2_deg2/p2_build_zeros.py`
  - `python -X utf8 05-CODE/sec5_lambda2_deg2/p2_ratecheck.py`
  - `python -X utf8 05-CODE/sec5_lambda2_deg2/p3_analysis.py`
- Inputs: none
- Outputs: `zeros_deg2.json`, `rate_deg2.json`, `p3_analysis.json`
- Time: n/a

### Table T3.1 (correspondence of zeros of H_d to superstable parameters, §3.2)
- Script: `python -X utf8 05-CODE/sec7_lattice_d2_H2scan/run_scan.py`
- Inputs: fixed points from external packages (read by the script from its folder)
- Outputs: `h2_zeros_scan.json`
- Time: n/a

### Table T4.1 (357 Krawczyk boxes, §4.2, §4.3)
- d = 4 (61 boxes): `python -X utf8 05-CODE/sec3_krawczyk_d4/krawczyk4_certificate.py` → `krawczyk4_certificates.json`
- d = 6 (287 boxes): `python -X utf8 05-CODE/sec3_krawczyk_d6/krawczyk6_certificate.py` → `krawczyk6_certificates.json`
- d = 6, late 9 boxes: `python -X utf8 05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_certificate.py` → `krawczyk6_late_certificates.json`
- Verification: `krawczyk*_verify.py` in the same folders
- Requirements: python-flint

### Table T5.1 (law min(δ, 1/|λ₂|) on the d-grid, §5)
- `python -X utf8 05-CODE/sec5_sweep_d/sweep_rate_law.py`
- `python -X utf8 05-CODE/sec5_sweep_d/sweep_rate_law_v2.py`
- `python -X utf8 05-CODE/sec5_sweep_d/sweep_rate_law_v2b.py`
- `python -X utf8 05-CODE/sec5_sweep_d/sweep_rate_law_v2c.py`
- `python -X utf8 05-CODE/sec5_sweep_d/sweep_rate_law_v3_fractional.py`
- Outputs: `sweep_results.json`, `sweep_v2_results.json`, `sweep_v2b_results.json`, `sweep_v2c_results.json`, `sweep_v3_results.json`

### Table T5.2 (two-term model, §5)
- `python -X utf8 05-CODE/sec5_twoterm/p1_deg2_calib.py`
- `python -X utf8 05-CODE/sec5_twoterm/p1_modelsel.py`
- `python -X utf8 05-CODE/sec5_twoterm/p2_deg2_realzeros.py`
- `python -X utf8 05-CODE/sec5_twoterm/p3_quartic_twoterm.py`
- `python -X utf8 05-CODE/sec5_twoterm/p3b_quartic_hiprec.py`
- `python -X utf8 05-CODE/sec5_twoterm/p3c_quartic_tail.py`

### Table T5.3 (subfamily d = 6 late, §3.3)
- `python -X utf8 05-CODE/sec5_subfamily_d6/extract_targets.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/mside_scan.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/hside_scan.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/fit_models.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/control_extend.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/m12_check.py`
- `python -X utf8 05-CODE/sec5_subfamily_d6/rate_refine.py`

### Table T6.1 (spectral grid NOT RIGOROUS, 45 points d ∈ [1.9, 7], §6, Appendix B.5)
- `python -X utf8 05-CODE/sec2_spectra_grid_d/spectra_grid.py` → `spectra_grid.json` (45 points; spectral grid — NOT RIGOROUS, mpmath float exploration)

### Tables B.1–B.6 (Appendix B)
- B.1 catalogue of H₂ (222 unique zeros): `python -X utf8 05-CODE/sec7_lattice_d2_H2scan/run_scan.py`
- B.2 catalogue of H₄: `python -X utf8 05-CODE/sec3_pilot_d4/p2_H4_zeros.py`
- B.3 catalogue of H₆: `python -X utf8 05-CODE/sec3_pilot_d6/p2_H6_zeros.py`
- B.4 ladders of λ₂(4), λ₂(6): `python -X utf8 05-CODE/sec4_bridge2_lambda2_cert/build_certificates.py` (python-flint)
- B.5 spectral grid: see T6.1
- B.6 lattice statistics: `python -X utf8 05-CODE/sec7_lattice_stats_F10/f10_stats.py`

### Table A.5 (Appendix A, summary numerical check)
- current package: `python -X utf8 05-CODE/appA_lemma_check_v05/check_lemma_v05.py` → `check_lemma_v05.out`

## External inputs (d = 2)

Contents of `06-EXTERNAL-INPUTS/` (in full — `README-EXTERNAL.md`): copies of artefacts from two sibling projects, from which the numbers for the d = 2 row of Table T3.1 and the λ₂(2) certificate are taken. **They were not recomputed in this package** — their rigor status is preserved as in the source projects: "NOT RIGOROUS" (experiment, frozen protocol) for the catalogue of 23 zeros of H₂; "RIGOROUS" (Arb balls) for the 23 superstable parameters of the family z² + c and for the Arb localization of the Z72 zero; "RIGOROUS on the discretization N = 48" for λ₂(2) ± 1.82·10⁻³³ (same status as for λ₂(4), λ₂(6) in §4.4).

## Integrity check

```bash
# Cross-platform manifest verifier (format LF: <sha256><spaces><relpath>)
python _tools/verify_manifest.py

# Or (Linux/macOS/git-bash only):
sha256sum -c MANIFEST.sha256
```

The manifest covers 100 % of the package files. Verification instructions and the list of exclusions are in `README-CODE.md`, §4, §5.
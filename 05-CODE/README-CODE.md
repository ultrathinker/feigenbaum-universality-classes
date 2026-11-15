# README-CODE.md — replication package 05-CODE (current package)

Purpose of the package: (a) confirmation of the pattern principle on three universality classes (d = 2, 4, 6); (b) rigorous interval certification of 357 Krawczyk boxes + re-certification of layer 2 (S078, S090) at d = 6; (c) reproduction of Bridge-2 ladders for λ₂(4), λ₂(6); (d) spectral grid and statistics of the 1/16 lattice; (e) rescan of the H₄, H₆ catalogues by an incommensurate grid.

Supplement to the article — `05-SUPPLEMENT-METHODS.md` (S1–S9).

## 1. Composition of the package

Below is a table of **35 subfolders** of `05-CODE/` (31 from previous packages + 3 new in the intermediate package + 1 new in the current package) with mapping to the sections of `01-ARTICLE.md`. Files were verified by `ls` against `CODE-FOLDERS-current.txt`. The package also includes `README-CODE.md` (this file) and `RELINK-LOG.md` (journal of path-relinking during assembly).

| # | Subfolder | Article section | Main scripts | Key data (json) | Report | Status |
|---|---|---|---|---|---|---|
| 1 | `appA_lemma_derivation_A/` | Appendix A, §2.5 | `f4_main.py` → `f4_extras.py` → `f4_test1…f4_test6` | `f4_results.json`, `f4_predcheck.json`, `f4_extras.json` | `LEMMA-F4.md`, `REPORT-F4.md` | NOT RIGOROUS (identity) |
| 2 | `appA_lemma_identity_check/` | Appendix A | `check_phi.py` | — | — | NOT RIGOROUS (mpmath) |
| 3 | `appA_lemma_spectral_C/` | Appendix A, §2.5 | `search_candidate*.py`, `eig_candidate.py`, `search_alpha_candidates.py`, `search_monomial.py`, `search_eigenvalue.py` | `f4_results.json` | `LEMMA-F4.md`, `REPORT-F4.md`, `verify_log.txt` | NOT RIGOROUS |
| 4 | `appA_lemma_verify_B/` | Appendix A, §2.5 | `verify_identity.py`, `verify_frechet.py`, `verify_cosines.py` | `identity_check.json`, `frechet_check.json`, `cosine_check.json` | `LEMMA-F4-VERIFIED.md`, `REPORT-S1-VERIFY.md` | NOT RIGOROUS |
| 5 | `appA_lemma_check_v05/` **(intermediate package)** | Appendix A, A.5 | `check_lemma_v05.py` (integer even dk + odd branch [G]^{dk+1} for fractional dk + x-representation + g′(1) = −α^{d−1}) | — | `check_lemma_v05.out` | NOT RIGOROUS (mpmath) |
| 6 | `sec2_check_d205/` | §2.4, §6.2 | `run_full.py`, `quicktest*.py` | `d205_result.json` | `REPORT_D205_CHECK.md` | NOT RIGOROUS |
| 7 | `sec2_dynamics_verdict/` | §6.2 | `cascade_gen.py`, `analyze_final.py` | `cascade_d275.json`, `cascade_d300.json`, `fits.json` | `REPORT_DYNAMICS_VERDICT.md` | NOT RIGOROUS |
| 8 | `sec2_spectra_grid_d/` | §6, Appendix B.5 | `spectra_grid.py`, `fix_pos_n40.py`, `retrack_pos.py`, `validate_dt.py` | `spectra_grid.json` (45 points) | `REPORT-SPECTRA-GRID.md`, `compare_sweep.json` | NOT RIGOROUS (mpmath float) |
| 9 | `sec2_spectrum_d6_T22/` | §2.3, §4.4 | `build_spectrum_d6.py`, `eig_scan_d6.py`, `cross_check_engine.py`, `refine_and_report.py` | `spectrum_d6.json`, `eig_scan_M32.json`, `eig_scan_M48.json`, `cross_check_engine_N32.json`, `cross_check_engine_N48.json` | `REPORT-T22-D6.md` | NOT RIGOROUS |
| 10 | `sec2_tower_depth_A/` *(former `sec2_tower_depth_F5`)* | §2.4, §7.4 | `f5_d8.py`, `f5_analyze.py`, `f5_consolidate.py`, `check_d8.py` | `f5_d8.json`, `f5_results.json`, `f5_depth_table.json`, `f5_neighborhood.json` | `REPORT-F5.md` | NOT RIGOROUS |
| 11 | `sec2_tower_depth_B/` *(former `sec2_tower_depth_S3`)* | §2.4, §7.4 | `step1_d8_fixedpoint.py`, `step2_d8_spectrum.py`, `step2b_n64_confirm.py`, `step2c_full_confirm.py`, `step3_neighborhood.py`, `step3c_d6_cosine_matched.py`, `step4_consolidate.py` | `s3_f5_results.json`, `step1_d8_fixedpoint.json`, `step2_d8_spectrum.json`, `step3_neighborhood.json`, `step3b_cosine_fix.json`, `step3c_d6_cosine_matched.json` | `REPORT-S3-F5.md`, `ERRATUM.md` (cosines of the report — coefficient metric in monomial basis; see #35) | NOT RIGOROUS |
| 12 | `sec2_vector_homotopy/` | §2.4, §2.5 | `homotopy_track.py`, `pairwise_map.py`, `spectrum_window.py`, `multiplicity_d2.py`, `refine_bridge.py` | `homotopy_track.json`, `pairwise_map.json`, `spectrum_window.json`, `multiplicity_d2.json` | `REPORT_VECTOR_HOMOTOPY.md` | NOT RIGOROUS |
| 13 | `sec3_krawczyk_d4/` | §4.2 | `krawczyk4_certificate.py`, `krawczyk4_verify.py`, `run_tampers4.py` | `krawczyk4_certificates.json` | `REPORT-KRAWCZYK4.md`, `TAMPER_LOG.md` | **RIGOROUS** (Arb) |
| 14 | `sec3_krawczyk_d6/` | §4.3 | `krawczyk6_certificate.py`, `krawczyk6_verify.py`, `calibrate6.py`, `refine6.py`, `run_tampers6.py` | `krawczyk6_certificates.json`, `calibration6.json`, `refined_seeds6.json` | `REPORT-KRAWCZYK6.md`, `TAMPER_LOG.md` | **RIGOROUS** (Arb) |
| 15 | `sec3_krawczyk_d6_late/` | §3.3, §4.3 | `krawczyk6_late_certificate.py`, `krawczyk6_late_verify.py`, `calibrate_late.py`, `refine_late.py`, `extract_targets_late.py`, `run_tampers_late.py` | `krawczyk6_late_certificates.json`, `calibration_late.json`, `refined_seeds_late.json`, `targets_late.json` | `REPORT-LATEBOXES.md`, `TAMPER_LATE.md` | **RIGOROUS** (Arb) |
| 16 | `sec3_krawczyk_d6_layer2_recheck/` **(intermediate package)** | §4.3 (layer 2) | `r2_sweep.py`, `r2_recert.py` (attempts `pow` and `mul`; prec ∈ {256, 512, 1024}; r ∈ {1e-30…1e-34}) | `r2_table.json` (168 rows), `S090_S078_recert.json` | `README.md` | **RIGOROUS** (Arb), re-certification of S078/S090 |
| 17 | `sec3_lattice_stats_v1/` | §7.2 | `analysis.py`, `analysis_g7.py`, `make_figs.py` | `lattice_results.json` | `REPORT-LATTICE.md` | NOT RIGOROUS |
| 18 | `sec3_pilot_d4/` | §3.2, §2.6 | `p1_fixedpoint.py`, `p2_H4_zeros.py`, `p3_mside.py`, `p4_match.py` | `quartic_fixedpoint.json`, `H4_zeros.json`, `mside_cascade.json`, `match_results.json` | `REPORT_QUARTIC_PILOT.md` | NOT RIGOROUS |
| 19 | `sec3_pilot_d6/` | §3.2, §2.6 | `p1_fixedpoint.py`, `p2_H6_zeros.py`, `p3_mside.py`, `p4_match.py`, `diag_*.py`, `extract.py` | `sextic_fixedpoint.json`, `H6_zeros.json`, `mside_cascade.json`, `match_results.json` | `REPORT_DEG6_PILOT.md` | NOT RIGOROUS |
| 20 | `sec3_rescan_patch_d6/` | §3.3 | `rescan.py`, `repolish.py`, `final_validate.py`, `make_fig.py`, `validate.py` | — | — | NOT RIGOROUS |
| 21 | `sec4_bridge1_lambda2_partial/` | §4.4, §5.3 | `run_bridge.py`, `bridge_lib.py`, `corrupt_battery.py`, `test_krawczyk_small.py`, `test_small.py` | `bridge_d4.json`, `bridge_d6.json`, `bridge_lambda2_certificates.json`, `corrupt_battery_results.json` | `REPORT-BRIDGE-LAMBDA2.md` | conditionally RIGOROUS |
| 22 | `sec4_bridge2_lambda2_cert/` | §4.4 | `build_certificates.py`, `verify_most2.py`, `diagn_cond.py`, `proto_cheb.py`, `smoke.py`, `corrupt_battery_most2.py` | `most2_run_d4.json`, `most2_run_d6.json`, `most2_certificates.json`, `verify_most2_results.json`, `corrupt_battery_most2_results.json` | `REPORT-MOST2.md`, `TAMPER_MOST2.md` | **RIGOROUS on the discretization** (Chebyshev basis) |
| 23 | `sec5_lambda2_deg2/` | §2.3, §5.3 | `p1_deg2.py`, `p2_build_zeros.py`, `p2_ratecheck.py`, `p3_analysis.py` | `p1_deg2.json`, `rate_deg2.json`, `zeros_deg2.json`, `p3_analysis.json` | `REPORT_LAMBDA2_DEG2.md`, `audit_results.json` | NOT RIGOROUS |
| 24 | `sec5_subfamily_d6/` | §3.3 | `extract_targets.py`, `mside_scan.py`, `hside_scan.py`, `fit_models.py`, `control_extend.py`, `m12_check.py`, `rate_refine.py` | `targets.json`, `mside_scan.json`, `hside_scan.json`, `fit_results.json`, `rate_*.json`, `calibration.json`, `control_extra.json` | `REPORT_SUBFAMILY.md`, `notes.md` | NOT RIGOROUS |
| 25 | `sec5_sweep_d/` | §5, §6 | `sweep_rate_law.py`, `sweep_rate_law_v2.py`, `sweep_rate_law_v2b.py`, `sweep_rate_law_v2c.py`, `sweep_rate_law_v3_fractional.py` | `sweep_results.json`, `sweep_v2_results.json`, `sweep_v2b_results.json`, `sweep_v2c_results.json`, `sweep_v3_results.json` | — | NOT RIGOROUS |
| 26 | `sec5_twoterm/` | §5 | `p1_deg2_calib.py`, `p1_modelsel.py`, `p2_deg2_realzeros.py`, `p3_quartic_twoterm.py`, `p3b_quartic_hiprec.py`, `p3c_quartic_tail.py` | `p1_deg2_calib.json`, `p1_modelsel.json`, `p2_deg2_realzeros.json`, `p3_quartic_twoterm.json`, `p3b_quartic_hiprec.json`, `p3c_quartic_tail.json` | `REPORT_TWOTERM.md` | NOT RIGOROUS |
| 27 | `sec7_lattice_addresses_F7/` | §7.2 | `f7_lattice.py`, `f7_realaxis.py`, `f7_invariants.py`, `f7_summarize.py` | `f7_lattice.json`, `f7_realaxis.json`, `f7_invariants.json` | `REPORT-F7.md` | NOT RIGOROUS |
| 28 | `sec7_lattice_d2_H2scan/` | §3.2, §7.2 | `run_scan.py`, `run_scan2_*.py`, `diag_*.py`, `run_stats.py` | `h2_zeros_scan.json`, `scan_candidates_v2.json`, `stats_lattice.json` | `REPORT-S2-H2SCAN.md` | NOT RIGOROUS |
| 29 | `sec7_lattice_period_F8/` | §7.2 | `rayleigh_lattice.py`, `quasiperiodicity.py`, `Hcalc.py`, `zeros_vs_m.py`, `local_nodes.py`, `figs.py`, `consolidate.py` | `rayleigh_results.json`, `quasiperiodicity_results.json`, `zeros_vs_m.json`, `f8_results.json`, `hcalc_meta.json`, `local_nodes.json` | `REPORT-F8.md` | NOT RIGOROUS |
| 30 | `sec7_lattice_stats_F10/` | §7.2 | `f10_stats.py` | `f10_results.json` | `REPORT-F10.md` | NOT RIGOROUS |
| 31 | `sec7_lattice_rescan_incommensurable/` **(intermediate package)** | §7.2 (incommensurate grid rescan) | `r1_*`, `r1a_*`, `r1b_*` (1/24 rescan for d = 4 and d = 6) | `H4_step24.json`, `H6_step24.json` (final), `H4_rescan_*.json`, `H6_rescan_*.json`, `r1_stats.json` | `REPORT-R1-LATTICE-RESCAN.md` | NOT RIGOROUS |
| 32 | `sec7_principle_status_F1F2/` | §7.3 | `f1_analysis.py`, `f1_analysis2.py`, `f1_dump.py`, `f2c_experiment.py` | `f1f2_results.json` | `REPORT-F1F2.md`, `THEOREM-SKETCH-F1.md` | NOT RIGOROUS |
| 33 | `sec7_principle_status_check/` | §7.3 | `check_f2c.py`, `check_f2c_fix.py`, `check_oos.py`, `check_rate.py` | `f2c_hits.txt`, `f2c_log_dump.txt`, `rate_full_dump.txt` | `REVIEW-F1F2.md` | NOT RIGOROUS |
| 34 | `sec7_principle_status_v2/` | §7.3 | — | — | `REPORT-F1F2-v2.md`, `THEOREM-SKETCH-F1-v2.md` | NOT RIGOROUS |
| 35 | `sec2_tower_shadow_small_dd/` **(current package)** | §2.5, §7.4, Appendix A §A.4.2 | `r9_main.py` (first run, cosine of coefficients in monomial basis), `check_run.py` + `check_shadow_d6.py` (recheck: own collocation, two bases, cosine of eigenfunctions, mode tracing by continuity) | `r9_results.json`, `check_results.json` | `README-shadow.md` | NOT RIGOROUS (mpmath dps 40) |

**Service files:**

- `05-CODE/_tools/verify_manifest.py` — cross-platform check of `MANIFEST.sha256` (line format: `<sha256><spaces/tabs><relpath>`, LF; see §4).
- `05-CODE/RELINK-LOG.md` — journal of path relinking during package assembly (absolute paths of the authors' machines → relative paths inside the package; remaining external paths are listed there).

## 2. Dependencies

- **Python 3.13** with packages: **mpmath** ≥ 1.3, **numpy**, **scipy**.
- **python-flint 0.9.0** (Arb intervals) — for `sec3_krawczyk_d4`, `sec3_krawczyk_d6`, `sec3_krawczyk_d6_late`, `sec3_krawczyk_d6_layer2_recheck`, `sec4_bridge2_lambda2_cert`.
- All scripts are single-process, except for small worker pools in `sec4_bridge2_lambda2_cert` and `sec3_krawczyk_d6_late`.
- No network access; all inputs are local files of the package.

## 3. Order of reproduction

Dependency graph (intermediate-package steps added):

```
§2.5 (Appendix A):  appA_lemma_derivation_A → appA_lemma_spectral_C
                                    └─→ appA_lemma_verify_B
                                    └─→ appA_lemma_identity_check
                                    └─→ appA_lemma_check_v05        (intermediate package)

§3.3 d=6 anomaly:  sec3_pilot_d6 → sec3_rescan_patch_d6 → sec5_subfamily_d6 → sec3_krawczyk_d6_late
                                                                              └─→ sec3_krawczyk_d6_layer2_recheck  (intermediate package)

§4.2-4.3 Krawczyk:  sec3_krawczyk_d4 → sec3_krawczyk_d6 → sec3_krawczyk_d6_late → sec3_krawczyk_d6_layer2_recheck

§4.4 Bridge-2:  sec3_krawczyk_d4 + sec3_krawczyk_d6 + sec5_lambda2_deg2 → sec4_bridge2_lambda2_cert

§6 grid:  sec2_spectra_grid_d

§7.2 lattice:  sec7_lattice_d2_H2scan + sec3_pilot_d4 + sec3_pilot_d6
                                └─→ sec7_lattice_addresses_F7
                                └─→ sec7_lattice_period_F8
                                └─→ sec7_lattice_stats_F10
                                └─→ sec7_lattice_rescan_incommensurable  (intermediate package: incommensurate grid rescan)

§7.3 status:  sec3_pilot_d4 + sec3_pilot_d6 → sec7_principle_status_F1F2
                                              └─→ sec7_principle_status_check
                                              └─→ sec7_principle_status_v2
```

## 4. Integrity check: sha256

```bash
# Cross-platform check (recommended)
python _tools/verify_manifest.py

# Or git-bash / Linux / macOS
sha256sum -c MANIFEST.sha256
```

The manifest is in LF format: `<sha256><spaces/tabs><relpath>`, sorted by path. Covers 100 % of package files, except the manifest itself.

## 5. What is not included

The publication package does not include: service by-products (`__pycache__/`, `*.pyc`); terminal session logs; regenerable diagnostic files. The full list of exclusions is in the reports of the corresponding subfolders (`TAMPER_LOG.md`, `TAMPER_LATE.md`, `TAMPER_MOST2.md`, `REPORT-*.md`).

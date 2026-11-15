# README-CODE.md — репликационный пакет 05-CODE (текущий пакет)

Назначение пакета: (а) подтверждение принципа узора на трёх классах универсальности (d = 2, 4, 6); (б) строгая интервальная сертификация 357 боксов Кравчика + повторная сертификация слоя 2 (S078, S090) при d = 6; (в) воспроизведение лестниц Мост-2 для λ₂(4), λ₂(6); (г) спектральная сетка и статистика решётки 1/16; (д) перескан каталогов H₄, H₆ несоизмеримой сеткой.

Дополнение к статье — `05-SUPPLEMENT-METHODS.md` (S1–S9).

## 1. Состав пакета

Ниже — таблица **35 подпапок** `05-CODE/` (31 из предыдущих пакетов + 3 новые в промежуточном пакете + 1 новая в текущем пакете) с привязкой к разделам `01-ARTICLE.md`. Файлы проверены `ls` против `CODE-FOLDERS-current.txt`. В пакет также входят `README-CODE.md` (этот файл) и `RELINK-LOG.md` (журнал перелинковки путей при сборке).

| # | Подпапка | Раздел статьи | Главные скрипты | Ключевые данные (json) | Отчёт | Статус |
|---|---|---|---|---|---|---|
| 1 | `appA_lemma_derivation_A/` | Приложение A, §2.5 | `f4_main.py` → `f4_extras.py` → `f4_test1…f4_test6` | `f4_results.json`, `f4_predcheck.json`, `f4_extras.json` | `LEMMA-F4.md`, `REPORT-F4.md` | NOT RIGOROUS (тождество) |
| 2 | `appA_lemma_identity_check/` | Приложение A | `check_phi.py` | — | — | NOT RIGOROUS (mpmath) |
| 3 | `appA_lemma_spectral_C/` | Приложение A, §2.5 | `search_candidate*.py`, `eig_candidate.py`, `search_alpha_candidates.py`, `search_monomial.py`, `search_eigenvalue.py` | `f4_results.json` | `LEMMA-F4.md`, `REPORT-F4.md`, `verify_log.txt` | NOT RIGOROUS |
| 4 | `appA_lemma_verify_B/` | Приложение A, §2.5 | `verify_identity.py`, `verify_frechet.py`, `verify_cosines.py` | `identity_check.json`, `frechet_check.json`, `cosine_check.json` | `LEMMA-F4-VERIFIED.md`, `REPORT-S1-VERIFY.md` | NOT RIGOROUS |
| 5 | `appA_lemma_check_v05/` **(промежуточный пакет)** | Приложение A, A.5 | `check_lemma_v05.py` (целые чётные dk + нечётная ветвь [G]^{dk+1} для дробных dk + x-представление + g′(1) = −α^{d−1}) | — | `check_lemma_v05.out` | NOT RIGOROUS (mpmath) |
| 6 | `sec2_check_d205/` | §2.4, §6.2 | `run_full.py`, `quicktest*.py` | `d205_result.json` | `REPORT_D205_CHECK.md` | NOT RIGOROUS |
| 7 | `sec2_dynamics_verdict/` | §6.2 | `cascade_gen.py`, `analyze_final.py` | `cascade_d275.json`, `cascade_d300.json`, `fits.json` | `REPORT_DYNAMICS_VERDICT.md` | NOT RIGOROUS |
| 8 | `sec2_spectra_grid_d/` | §6, Приложение B.5 | `spectra_grid.py`, `fix_pos_n40.py`, `retrack_pos.py`, `validate_dt.py` | `spectra_grid.json` (45 точек) | `REPORT-SPECTRA-GRID.md`, `compare_sweep.json` | NOT RIGOROUS (mpmath float) |
| 9 | `sec2_spectrum_d6_T22/` | §2.3, §4.4 | `build_spectrum_d6.py`, `eig_scan_d6.py`, `cross_check_engine.py`, `refine_and_report.py` | `spectrum_d6.json`, `eig_scan_M32.json`, `eig_scan_M48.json`, `cross_check_engine_N32.json`, `cross_check_engine_N48.json` | `REPORT-T22-D6.md` | NOT RIGOROUS |
| 10 | `sec2_tower_depth_A/` *(бывшая `sec2_tower_depth_F5`)* | §2.4, §7.4 | `f5_d8.py`, `f5_analyze.py`, `f5_consolidate.py`, `check_d8.py` | `f5_d8.json`, `f5_results.json`, `f5_depth_table.json`, `f5_neighborhood.json` | `REPORT-F5.md` | NOT RIGOROUS |
| 11 | `sec2_tower_depth_B/` *(бывшая `sec2_tower_depth_S3`)* | §2.4, §7.4 | `step1_d8_fixedpoint.py`, `step2_d8_spectrum.py`, `step2b_n64_confirm.py`, `step2c_full_confirm.py`, `step3_neighborhood.py`, `step3c_d6_cosine_matched.py`, `step4_consolidate.py` | `s3_f5_results.json`, `step1_d8_fixedpoint.json`, `step2_d8_spectrum.json`, `step3_neighborhood.json`, `step3b_cosine_fix.json`, `step3c_d6_cosine_matched.json` | `REPORT-S3-F5.md`, `ERRATUM.md` (косинусы отчёта — коэффициентная метрика в мономиальном базисе; см. № 35) | NOT RIGOROUS |
| 12 | `sec2_vector_homotopy/` | §2.4, §2.5 | `homotopy_track.py`, `pairwise_map.py`, `spectrum_window.py`, `multiplicity_d2.py`, `refine_bridge.py` | `homotopy_track.json`, `pairwise_map.json`, `spectrum_window.json`, `multiplicity_d2.json` | `REPORT_VECTOR_HOMOTOPY.md` | NOT RIGOROUS |
| 13 | `sec3_krawczyk_d4/` | §4.2 | `krawczyk4_certificate.py`, `krawczyk4_verify.py`, `run_tampers4.py` | `krawczyk4_certificates.json` | `REPORT-KRAWCZYK4.md`, `TAMPER_LOG.md` | **СТРОГО** (Arb) |
| 14 | `sec3_krawczyk_d6/` | §4.3 | `krawczyk6_certificate.py`, `krawczyk6_verify.py`, `calibrate6.py`, `refine6.py`, `run_tampers6.py` | `krawczyk6_certificates.json`, `calibration6.json`, `refined_seeds6.json` | `REPORT-KRAWCZYK6.md`, `TAMPER_LOG.md` | **СТРОГО** (Arb) |
| 15 | `sec3_krawczyk_d6_late/` | §3.3, §4.3 | `krawczyk6_late_certificate.py`, `krawczyk6_late_verify.py`, `calibrate_late.py`, `refine_late.py`, `extract_targets_late.py`, `run_tampers_late.py` | `krawczyk6_late_certificates.json`, `calibration_late.json`, `refined_seeds_late.json`, `targets_late.json` | `REPORT-LATEBOXES.md`, `TAMPER_LATE.md` | **СТРОГО** (Arb) |
| 16 | `sec3_krawczyk_d6_layer2_recheck/` **(промежуточный пакет)** | §4.3 (слой 2) | `r2_sweep.py`, `r2_recert.py` (попытки `pow` и `mul`; prec ∈ {256, 512, 1024}; r ∈ {1e-30…1e-34}) | `r2_table.json` (168 строк), `S090_S078_recert.json` | `README.md` | **СТРОГО** (Arb), повторная сертификация S078/S090 |
| 17 | `sec3_lattice_stats_v1/` | §7.2 | `analysis.py`, `analysis_g7.py`, `make_figs.py` | `lattice_results.json` | `REPORT-LATTICE.md` | NOT RIGOROUS |
| 18 | `sec3_pilot_d4/` | §3.2, §2.6 | `p1_fixedpoint.py`, `p2_H4_zeros.py`, `p3_mside.py`, `p4_match.py` | `quartic_fixedpoint.json`, `H4_zeros.json`, `mside_cascade.json`, `match_results.json` | `REPORT_QUARTIC_PILOT.md` | NOT RIGOROUS |
| 19 | `sec3_pilot_d6/` | §3.2, §2.6 | `p1_fixedpoint.py`, `p2_H6_zeros.py`, `p3_mside.py`, `p4_match.py`, `diag_*.py`, `extract.py` | `sextic_fixedpoint.json`, `H6_zeros.json`, `mside_cascade.json`, `match_results.json` | `REPORT_DEG6_PILOT.md` | NOT RIGOROUS |
| 20 | `sec3_rescan_patch_d6/` | §3.3 | `rescan.py`, `repolish.py`, `final_validate.py`, `make_fig.py`, `validate.py` | — | — | NOT RIGOROUS |
| 21 | `sec4_bridge1_lambda2_partial/` | §4.4, §5.3 | `run_bridge.py`, `bridge_lib.py`, `corrupt_battery.py`, `test_krawczyk_small.py`, `test_small.py` | `bridge_d4.json`, `bridge_d6.json`, `bridge_lambda2_certificates.json`, `corrupt_battery_results.json` | `REPORT-BRIDGE-LAMBDA2.md` | условно-СТРОГО |
| 22 | `sec4_bridge2_lambda2_cert/` | §4.4 | `build_certificates.py`, `verify_most2.py`, `diagn_cond.py`, `proto_cheb.py`, `smoke.py`, `corrupt_battery_most2.py` | `most2_run_d4.json`, `most2_run_d6.json`, `most2_certificates.json`, `verify_most2_results.json`, `corrupt_battery_most2_results.json` | `REPORT-MOST2.md`, `TAMPER_MOST2.md` | **СТРОГО на дискретизации** (чебышёвский базис) |
| 23 | `sec5_lambda2_deg2/` | §2.3, §5.3 | `p1_deg2.py`, `p2_build_zeros.py`, `p2_ratecheck.py`, `p3_analysis.py` | `p1_deg2.json`, `rate_deg2.json`, `zeros_deg2.json`, `p3_analysis.json` | `REPORT_LAMBDA2_DEG2.md`, `audit_results.json` | NOT RIGOROUS |
| 24 | `sec5_subfamily_d6/` | §3.3 | `extract_targets.py`, `mside_scan.py`, `hside_scan.py`, `fit_models.py`, `control_extend.py`, `m12_check.py`, `rate_refine.py` | `targets.json`, `mside_scan.json`, `hside_scan.json`, `fit_results.json`, `rate_*.json`, `calibration.json`, `control_extra.json` | `REPORT_SUBFAMILY.md`, `notes.md` | NOT RIGOROUS |
| 25 | `sec5_sweep_d/` | §5, §6 | `sweep_rate_law.py`, `sweep_rate_law_v2.py`, `sweep_rate_law_v2b.py`, `sweep_rate_law_v2c.py`, `sweep_rate_law_v3_fractional.py` | `sweep_results.json`, `sweep_v2_results.json`, `sweep_v2b_results.json`, `sweep_v2c_results.json`, `sweep_v3_results.json` | — | NOT RIGOROUS |
| 26 | `sec5_twoterm/` | §5 | `p1_deg2_calib.py`, `p1_modelsel.py`, `p2_deg2_realzeros.py`, `p3_quartic_twoterm.py`, `p3b_quartic_hiprec.py`, `p3c_quartic_tail.py` | `p1_deg2_calib.json`, `p1_modelsel.json`, `p2_deg2_realzeros.json`, `p3_quartic_twoterm.json`, `p3b_quartic_hiprec.json`, `p3c_quartic_tail.json` | `REPORT_TWOTERM.md` | NOT RIGOROUS |
| 27 | `sec7_lattice_addresses_F7/` | §7.2 | `f7_lattice.py`, `f7_realaxis.py`, `f7_invariants.py`, `f7_summarize.py` | `f7_lattice.json`, `f7_realaxis.json`, `f7_invariants.json` | `REPORT-F7.md` | NOT RIGOROUS |
| 28 | `sec7_lattice_d2_H2scan/` | §3.2, §7.2 | `run_scan.py`, `run_scan2_*.py`, `diag_*.py`, `run_stats.py` | `h2_zeros_scan.json`, `scan_candidates_v2.json`, `stats_lattice.json` | `REPORT-S2-H2SCAN.md` | NOT RIGOROUS |
| 29 | `sec7_lattice_period_F8/` | §7.2 | `rayleigh_lattice.py`, `quasiperiodicity.py`, `Hcalc.py`, `zeros_vs_m.py`, `local_nodes.py`, `figs.py`, `consolidate.py` | `rayleigh_results.json`, `quasiperiodicity_results.json`, `zeros_vs_m.json`, `f8_results.json`, `hcalc_meta.json`, `local_nodes.json` | `REPORT-F8.md` | NOT RIGOROUS |
| 30 | `sec7_lattice_stats_F10/` | §7.2 | `f10_stats.py` | `f10_results.json` | `REPORT-F10.md` | NOT RIGOROUS |
| 31 | `sec7_lattice_rescan_incommensurable/` **(промежуточный пакет)** | §7.2 (перескан несоизмеримой сеткой) | `r1_*`, `r1a_*`, `r1b_*` (перескан 1/24 для d = 4 и d = 6) | `H4_step24.json`, `H6_step24.json` (финальные), `H4_rescan_*.json`, `H6_rescan_*.json`, `r1_stats.json` | `REPORT-R1-LATTICE-RESCAN.md` | NOT RIGOROUS |
| 32 | `sec7_principle_status_F1F2/` | §7.3 | `f1_analysis.py`, `f1_analysis2.py`, `f1_dump.py`, `f2c_experiment.py` | `f1f2_results.json` | `REPORT-F1F2.md`, `THEOREM-SKETCH-F1.md` | NOT RIGOROUS |
| 33 | `sec7_principle_status_check/` | §7.3 | `check_f2c.py`, `check_f2c_fix.py`, `check_oos.py`, `check_rate.py` | `f2c_hits.txt`, `f2c_log_dump.txt`, `rate_full_dump.txt` | `REVIEW-F1F2.md` | NOT RIGOROUS |
| 34 | `sec7_principle_status_v2/` | §7.3 | — | — | `REPORT-F1F2-v2.md`, `THEOREM-SKETCH-F1-v2.md` | NOT RIGOROUS |
| 35 | `sec2_tower_shadow_small_dd/` **(текущий пакет)** | §2.5, §7.4, Приложение A §A.4.2 | `r9_main.py` (первый прогон, косинус коэффициентов в мономиальном базисе), `check_run.py` + `check_shadow_d6.py` (перепроверка: своя коллокация, два базиса, косинус собственных функций, трассировка моды по непрерывности) | `r9_results.json`, `check_results.json` | `README-shadow.md` | NOT RIGOROUS (mpmath dps 40) |

**Служебные файлы:**

- `05-CODE/_tools/verify_manifest.py` — кросс-платформенный контроль `MANIFEST.sha256` (формат строк: `<sha256><пробелы/табы><relпуть>`, LF; см. §4).
- `05-CODE/RELINK-LOG.md` — журнал перелинковки путей при сборке пакета (абсолютные пути авторских машин → относительные пути внутри пакета; оставшиеся внешние пути перечислены там же).

## 2. Зависимости

- **Python 3.13** с пакетами: **mpmath** ≥ 1.3, **numpy**, **scipy**.
- **python-flint 0.9.0** (Arb-интервалы) — для `sec3_krawczyk_d4`, `sec3_krawczyk_d6`, `sec3_krawczyk_d6_late`, `sec3_krawczyk_d6_layer2_recheck`, `sec4_bridge2_lambda2_cert`.
- Все скрипты однопроцессные, кроме небольших пулов воркеров в `sec4_bridge2_lambda2_cert` и `sec3_krawczyk_d6_late`.
- Без сетевого доступа; все входы — локальные файлы пакета.

## 3. Порядок воспроизведения

Граф зависимостей (добавлены шаги промежуточного пакета):

```
§2.5 (приложение A):  appA_lemma_derivation_A → appA_lemma_spectral_C
                                       └─→ appA_lemma_verify_B
                                       └─→ appA_lemma_identity_check
                                       └─→ appA_lemma_check_v05        (промежуточный пакет)

§3.3 d=6 аномалия:  sec3_pilot_d6 → sec3_rescan_patch_d6 → sec5_subfamily_d6 → sec3_krawczyk_d6_late
                                                                                └─→ sec3_krawczyk_d6_layer2_recheck  (промежуточный пакет)

§4.2-4.3 Кравчик:  sec3_krawczyk_d4 → sec3_krawczyk_d6 → sec3_krawczyk_d6_late → sec3_krawczyk_d6_layer2_recheck

§4.4 Мост-2:  sec3_krawczyk_d4 + sec3_krawczyk_d6 + sec5_lambda2_deg2 → sec4_bridge2_lambda2_cert

§6 сетка:  sec2_spectra_grid_d

§7.2 решётка:  sec7_lattice_d2_H2scan + sec3_pilot_d4 + sec3_pilot_d6
                                    └─→ sec7_lattice_addresses_F7
                                    └─→ sec7_lattice_period_F8
                                    └─→ sec7_lattice_stats_F10
                                    └─→ sec7_lattice_rescan_incommensurable  (промежуточный пакет: перескан несоизмеримой сеткой)

§7.3 статус:  sec3_pilot_d4 + sec3_pilot_d6 → sec7_principle_status_F1F2
                                              └─→ sec7_principle_status_check
                                              └─→ sec7_principle_status_v2
```

## 4. Проверка целостности: sha256

```bash
# Кросс-платформенный контроль (рекомендуется)
python _tools/verify_manifest.py

# Или git-bash / Linux / macOS
sha256sum -c MANIFEST.sha256
```

Манифест в формате LF: `<sha256><пробелы/табы><relпуть>`, сортировка по пути. Покрывает 100 % файлов пакета, кроме самого манифеста.

## 5. Что не включено

В публикационный комплект не входят: служебные производные (`__pycache__/`, `*.pyc`); терминальные журналы сессий; регенерируемые диагностические файлы. Полный перечень исключений — в отчётах соответствующих подпапок (`TAMPER_LOG.md`, `TAMPER_LATE.md`, `TAMPER_MOST2.md`, `REPORT-*.md`).

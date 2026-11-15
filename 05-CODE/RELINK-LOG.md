# RELINK-LOG — замена абсолютных путей авторской машины на пути внутри пакета

Скрипты запускались авторами «на месте» с абсолютными путями. При сборке пакета каждое вхождение пути,
указывающего в папку, включённую в пакет, заменено на путь относительно расположения скрипта (прямые слэши).
Пути, не имеющие аналога в пакете, оставлены как есть и перечислены ниже как ВНЕШНИЕ.

| файл | было | стало |
|---|---|---|
| 05-CODE/appA_lemma_check_v05/check_lemma_v05.py | `../../wave1_spectra_d/spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/appA_lemma_identity_check/check_phi.py | `../../wave1_spectra_d/spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_dynamics_verdict/analyze.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_dynamics_verdict` | `.` |
| 05-CODE/sec2_dynamics_verdict/analyze_final.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_dynamics_verdict` | `.` |
| 05-CODE/sec2_dynamics_verdict/cascade_gen.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_dynamics_verdict` | `.` |
| 05-CODE/sec2_dynamics_verdict/cascade_gen.py | `../quick_sweep_d/sweep_rate_law_v2c.py` | `../sec5_sweep_d/sweep_rate_law_v2c.py` |
| 05-CODE/sec2_spectra_grid_d/audit_orchestrator.py | `../quartic_pilot/quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec2_spectrum_d6_T22/build_spectrum_d6.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave1_spectra_d\spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_spectrum_d6_T22/cross_check_engine.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_vector_homotopy` | `../sec2_vector_homotopy` |
| 05-CODE/sec2_spectrum_d6_T22/eig_scan_d6.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave1_spectra_d\spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_spectrum_d6_T22/refine_and_report.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave1_spectra_d\spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_tower_shadow_small_dd/check_run.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave5_stage3\r9_shadow_d6_mx` | `.` |
| 05-CODE/sec2_tower_shadow_small_dd/check_shadow_d6.py | `D:\\Исследования\\И-3 — Классы универсальности\\Исследование-3-Принцип\\wave1_spectra_d\\spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_tower_shadow_small_dd/check_shadow_d6.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave1_spectra_d\spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_vector_homotopy/audit_orchestrator.py | `../wave1_spectra_d/spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec2_vector_homotopy/audit_orchestrator.py | `../wave1_spectra_d/spectra_grid.json` | `../sec2_spectra_grid_d/spectra_grid.json` |
| 05-CODE/sec3_krawczyk_d4/audit_orchestrator.py | `../quartic_pilot/match_results.json` | `../sec3_pilot_d4/match_results.json` |
| 05-CODE/sec3_krawczyk_d4/krawczyk4_certificate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/match_results.json` | `../sec3_pilot_d4/match_results.json` |
| 05-CODE/sec3_krawczyk_d4/krawczyk4_verify.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/match_results.json` | `../sec3_pilot_d4/match_results.json` |
| 05-CODE/sec3_krawczyk_d4/run_tampers4.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/match_results.json` | `../sec3_pilot_d4/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/audit_orchestrator.py | `../agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/calibrate6.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/krawczyk6_certificate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/krawczyk6_verify.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/refine6.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6/run_tampers6.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_krawczyk_d6_late/audit_orchestrator.py | `../wave2_krawczyk6/krawczyk6_certificates.json` | `../sec3_krawczyk_d6/krawczyk6_certificates.json` |
| 05-CODE/sec3_krawczyk_d6_late/extract_targets_late.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/subfamily` | `../sec5_subfamily_d6` |
| 05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_certificate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave2_krawczyk6/` | `../sec3_krawczyk_d6` |
| 05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_verify.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave2_krawczyk6/` | `../sec3_krawczyk_d6` |
| 05-CODE/sec3_krawczyk_d6_layer2_recheck/check_pairwise_isolation.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_krawczyk6\krawczyk6_certificates.json` | `../sec3_krawczyk_d6/krawczyk6_certificates.json` |
| 05-CODE/sec3_krawczyk_d6_layer2_recheck/r2_recert.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_krawczyk6\krawczyk6_certificates.json` | `../sec3_krawczyk_d6/krawczyk6_certificates.json` |
| 05-CODE/sec3_krawczyk_d6_layer2_recheck/r2_recert.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave2_krawczyk6\krawczyk6_certificates.json` | `../sec3_krawczyk_d6/krawczyk6_certificates.json` |
| 05-CODE/sec3_lattice_stats_v1/analysis.py | `../agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec3_lattice_stats_v1/analysis.py | `../agentB_deg6/H6_zeros.json` | `../sec3_pilot_d6/H6_zeros.json` |
| 05-CODE/sec3_lattice_stats_v1/analysis.py | `../agentB_deg6/subfamily/rate_full.json` | `../sec5_subfamily_d6/rate_full.json` |
| 05-CODE/sec3_lattice_stats_v1/analysis_g7.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave3_lattice` | `.` |
| 05-CODE/sec3_lattice_stats_v1/analysis_g7.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6` | `../sec3_pilot_d6` |
| 05-CODE/sec3_lattice_stats_v1/make_figs.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave3_lattice` | `.` |
| 05-CODE/sec3_lattice_stats_v1/make_figs.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6` | `../sec3_pilot_d6` |
| 05-CODE/sec3_pilot_d6/audit_orchestrator.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6\match_results.json` | `match_results.json` |
| 05-CODE/sec3_pilot_d6/audit_orchestrator.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6\mside_cascade.json` | `mside_cascade.json` |
| 05-CODE/sec3_pilot_d6/diag_orbit.py | `../quartic_pilot/quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec3_pilot_d6/diag_orbit.py | `../quartic_pilot/H4_zeros.json` | `../sec3_pilot_d4/H4_zeros.json` |
| 05-CODE/sec4_bridge1_lambda2_partial/corrupt_battery.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec4_bridge1_lambda2_partial/corrupt_battery.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6\sextic_fixedpoint.json` | `../sec3_pilot_d6/sextic_fixedpoint.json` |
| 05-CODE/sec4_bridge1_lambda2_partial/test_small.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec4_bridge1_lambda2_partial/test_small.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentB_deg6\sextic_fixedpoint.json` | `../sec3_pilot_d6/sextic_fixedpoint.json` |
| 05-CODE/sec4_bridge2_lambda2_cert/corrupt_battery_most2.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave3_most2\.` | `.` |
| 05-CODE/sec5_lambda2_deg2/audit_orchestrator.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentA_lambda2\audit_progress.log` | `audit_progress.log` |
| 05-CODE/sec5_lambda2_deg2/audit_orchestrator.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_lambda2_deg2/audit_orchestrator.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentA_lambda2\audit_results.json` | `audit_results.json` |
| 05-CODE/sec5_lambda2_deg2/p1_deg2.py | `D:\\Исследования\\И-3 — Классы универсальности\\Исследование-3-Принцип\\quartic_pilot\\p1_fixedpoint.py` | `../sec3_pilot_d4/p1_fixedpoint.py` |
| 05-CODE/sec5_lambda2_deg2/p2_build_zeros.py | `D:\\Исследования\\И-1 — Нули функции H\\Исследование-1\\09-CODE\\section7_mandelbrot\\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_lambda2_deg2/p2_build_zeros.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_lambda2_deg2/p2_ratecheck.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_subfamily_d6/audit_orchestrator.py | `../sextic_fixedpoint.json` | `../sec3_pilot_d6/sextic_fixedpoint.json` |
| 05-CODE/sec5_subfamily_d6/audit_orchestrator.py | `../mside_cascade.json` | `../sec3_pilot_d6/mside_cascade.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_results.json` | `sweep_results.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_progress.log` | `sweep_progress.log` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2_results.json` | `sweep_v2_results.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2_progress.log` | `sweep_v2_progress.log` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2b.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2b_results.json` | `sweep_v2b_results.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2b.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2b_progress.log` | `sweep_v2b_progress.log` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2c.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2c_results.json` | `sweep_v2c_results.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v2c.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v2c_progress.log` | `sweep_v2c_progress.log` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v3_fractional.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v3_results.json` | `sweep_v3_results.json` |
| 05-CODE/sec5_sweep_d/sweep_rate_law_v3_fractional.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quick_sweep_d\sweep_v3_progress.log` | `sweep_v3_progress.log` |
| 05-CODE/sec5_twoterm/audit_a14_check.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\mside_cascade.json` | `../sec3_pilot_d4/mside_cascade.json` |
| 05-CODE/sec5_twoterm/audit_aF_error.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\mside_cascade.json` | `../sec3_pilot_d4/mside_cascade.json` |
| 05-CODE/sec5_twoterm/audit_twoterm.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_twoterm/common.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\progress.log` | `progress.log` |
| 05-CODE/sec5_twoterm/p1_deg2_calib.py | `D:\\Исследования\\И-1 — Нули функции H\\Исследование-1\\09-CODE\\section7_mandelbrot\\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_twoterm/p1_deg2_calib.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p1_deg2_calib.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_twoterm/p1_deg2_calib.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p1_deg2_calib.json` | `p1_deg2_calib.json` |
| 05-CODE/sec5_twoterm/p1_modelsel.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p1_modelsel.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p1_deg2_calib.json` | `p1_deg2_calib.json` |
| 05-CODE/sec5_twoterm/p1_modelsel.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p1_modelsel.json` | `p1_modelsel.json` |
| 05-CODE/sec5_twoterm/p2_deg2_realzeros.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p2_deg2_realzeros.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec5_twoterm/p2_deg2_realzeros.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p2_deg2_realzeros.json` | `p2_deg2_realzeros.json` |
| 05-CODE/sec5_twoterm/p3b_quartic_hiprec.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p3b_quartic_hiprec.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\mside_cascade.json` | `../sec3_pilot_d4/mside_cascade.json` |
| 05-CODE/sec5_twoterm/p3b_quartic_hiprec.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec5_twoterm/p3b_quartic_hiprec.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p3b_quartic_hiprec.json` | `p3b_quartic_hiprec.json` |
| 05-CODE/sec5_twoterm/p3c_quartic_tail.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p3c_quartic_tail.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\mside_cascade.json` | `../sec3_pilot_d4/mside_cascade.json` |
| 05-CODE/sec5_twoterm/p3c_quartic_tail.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec5_twoterm/p3c_quartic_tail.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p3c_quartic_tail.json` | `p3c_quartic_tail.json` |
| 05-CODE/sec5_twoterm/p3_quartic_twoterm.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm` | `.` |
| 05-CODE/sec5_twoterm/p3_quartic_twoterm.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\mside_cascade.json` | `../sec3_pilot_d4/mside_cascade.json` |
| 05-CODE/sec5_twoterm/p3_quartic_twoterm.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentC_twoterm\p3_quartic_twoterm.json` | `p3_quartic_twoterm.json` |
| 05-CODE/sec5_twoterm/p3_quartic_twoterm.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\quartic_pilot\quartic_fixedpoint.json` | `../sec3_pilot_d4/quartic_fixedpoint.json` |
| 05-CODE/sec7_lattice_d2_H2scan/h2_common.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../../06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/results.json` |
| 05-CODE/sec7_lattice_period_F8/consolidate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/figs.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/Hcalc.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/local_nodes.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/quasiperiodicity.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/rayleigh_lattice.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_period_F8/zeros_vs_m.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/wave5_stage3/b1_3_F8_period` | `.` |
| 05-CODE/sec7_lattice_rescan_incommensurable/r1a_common.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/H4_zeros.json` | `../sec3_pilot_d4/H4_zeros.json` |
| 05-CODE/sec7_lattice_rescan_incommensurable/r1a_common.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/H4_zeros.json` | `../sec3_pilot_d4/H4_zeros.json` |
| 05-CODE/sec7_lattice_rescan_incommensurable/r1b_late9.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec7_lattice_rescan_incommensurable/r1_stats.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/quartic_pilot/H4_zeros.json` | `../sec3_pilot_d4/H4_zeros.json` |
| 05-CODE/sec7_lattice_rescan_incommensurable/r1_stats.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentB_deg6/match_results.json` | `../sec3_pilot_d6/match_results.json` |
| 05-CODE/sec7_lattice_stats_F10/f10_stats.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave5_stage3\b1_1_F10_stats` | `.` |
| 05-CODE/sec7_principle_status_check/_extract_log_chunk.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave5_stage3\b3_1_F1F2_deepseek` | `../sec7_principle_status_F1F2` |
| 05-CODE/sec7_principle_status_check/_extract_log_chunk2.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave5_stage3\b3_1_F1F2_deepseek` | `../sec7_principle_status_F1F2` |
| 05-CODE/sec7_principle_status_check/_extract_log_chunk3.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\wave5_stage3\b3_1_F1F2_deepseek` | `../sec7_principle_status_F1F2` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2b_predict.py | `D:/Исследования/И-1 — Нули функции H/Исследование-1/09-CODE/section7_mandelbrot/results.json` | `../study1_d2_frozen_experiment/results.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2b_predict.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentA_lambda2/zeros_deg2.json` | `../../05-CODE/sec5_lambda2_deg2/zeros_deg2.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2_rate.py | `D:/Исследования/И-1 — Нули функции H/Исследование-1/09-CODE/section7_mandelbrot/results.json` | `../study1_d2_frozen_experiment/results.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2_rate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentA_lambda2/zeros_deg2.json` | `../../05-CODE/sec5_lambda2_deg2/zeros_deg2.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2_rate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип/agentA_lambda2/rate_deg2.json` | `../../05-CODE/sec5_lambda2_deg2/rate_deg2.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-4 — Мост (реализация)\Исследование-4-Реализация\kilo_verify` | `.` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mandelbrot\results.json` | `../study1_d2_frozen_experiment/results.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_z72_localization\z72_zero_certificate.json` | `../study1_d2_z72_localization/z72_zero_certificate.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\section7_mside_superstable\mside_certificate.json` | `../study1_d2_superstable_certificates/mside_certificate.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип\agentA_lambda2\rate_deg2.json` | `../../05-CODE/sec5_lambda2_deg2/rate_deg2.json` |

## ВНЕШНИЕ (не в пакете, оставлены как есть)

| файл | путь |
|---|---|
| 05-CODE/appA_lemma_derivation_A/f4_common.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec3_krawczyk_d4/krawczyk4_certificate.py | `D:\Islanded-1` |
| 05-CODE/sec3_lattice_stats_v1/analysis.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec3_rescan_patch_d6/rescan.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec5_lambda2_deg2/audit_orchestrator.py | `D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json` |
| 05-CODE/sec5_lambda2_deg2/p2_build_zeros.py | `D:\\Исследования\\И-1 — Нули функции H\\Feigenbaum-archive\\ideation_sweep1\\codex\\zero_map.json` |
| 05-CODE/sec5_lambda2_deg2/p2_build_zeros.py | `D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json` |
| 05-CODE/sec5_twoterm/p2_deg2_realzeros.py | `D:\\Исследования\\И-1 — Нули функции H\\Feigenbaum-archive\\ideation_sweep1\\codex\\zero_map.json` |
| 05-CODE/sec5_twoterm/p2_deg2_realzeros.py | `D:\Исследования\И-1 — Нули функции H\Feigenbaum-archive\ideation_sweep1\codex\zero_map.json` |
| 05-CODE/sec7_lattice_addresses_F7/f7_invariants.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_addresses_F7/f7_realaxis.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_d2_H2scan/h2_common.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\archive_inputs\w33\t5_cert_sonnet\w33_gh_hiprec.json` |
| 05-CODE/sec7_lattice_d2_H2scan/h2_common.py | `D:\Исследования\И-1 — Нули функции H\Исследование-1\09-CODE\archive_inputs\path_to_8\point15_layer1_addressing\address_certificate.json` |
| 05-CODE/sec7_lattice_period_F8/consolidate.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/figs.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/Hcalc.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/local_nodes.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/quasiperiodicity.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/rayleigh_lattice.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_lattice_period_F8/zeros_vs_m.py | `D:/Исследования/И-3 — Классы универсальности/Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_check/check_f2c.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_check/check_f2c_fix.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_check/check_oos.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_check/check_rate.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_F1F2/f1_analysis.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_F1F2/f1_analysis2.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_F1F2/f1_dump.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 05-CODE/sec7_principle_status_F1F2/f2c_experiment.py | `D:\Исследования\И-3 — Классы универсальности\Исследование-3-Принцип` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v2_rate.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/ideation_sweep1/codex/zero_map.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-4 — Мост (реализация)\Исследование-4-Реализация\sonnet_veha_b\step6c_m11_krawczyk_result.json` |
| 06-EXTERNAL-INPUTS/bridge_lambda2_d2_certificate/v6_vehabcheck.py | `D:\Исследования\И-4 — Мост (реализация)\Исследование-4-Реализация\sonnet_veha_b\step6e_m12_result.json` |
| 06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/frozen_experiment.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/path_to_8/` |
| 06-EXTERNAL-INPUTS/study1_d2_frozen_experiment/make_report.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/path_to_8/` |
| 06-EXTERNAL-INPUTS/study1_d2_superstable_certificates/mside_certificate.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/` |
| 06-EXTERNAL-INPUTS/study1_d2_superstable_certificates/mside_verify.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/` |
| 06-EXTERNAL-INPUTS/study1_d2_superstable_certificates/run_tampers.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive/` |
| 06-EXTERNAL-INPUTS/study1_d2_z72_localization/z72_zero_build.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive` |
| 06-EXTERNAL-INPUTS/study1_d2_z72_localization/z72_zero_verify.py | `D:/Исследования/И-1 — Нули функции H/Feigenbaum-archive` |

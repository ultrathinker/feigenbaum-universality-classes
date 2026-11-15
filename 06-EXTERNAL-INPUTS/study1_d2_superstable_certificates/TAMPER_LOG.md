# TAMPER_LOG -- corruption battery for mside_certificate.json [oid:msideK1]

Date: 2026-08-14 17:28:11

Subject: `mside_certificate.json` (23 Krawczyk/interval-Newton
certificates of superstable parameters) verified by the independent
tuple-path verifier `mside_verify.py`. Principle: every tamper is
internally consistent ("lie plausibly") and must be caught by
recomputation from primary sources or by exact source-pin reconciliation.

Clean control: exit code **0** (expected 0) -- OK

| # | tamper | what was corrupted | exit | caught | first verifier complaint |
|---|--------|--------------------|------|--------|--------------------------|
| T1 | center_last_digit | Z05 box_center_re last digit +1 (shift ~1e-30, still inside the certified box; caught only by the exact source pin) | 1 | YES | FAIL: Z05 box_center_re drift vs source: cert='-1.401155625434136143045155574' source='-1.401155625434136143045155573' |
| T2 | z72_legacy_center | Z72 center replaced by the legacy certified-centre prediction c_pred (rejected provenance), other fields untouched | 1 | YES | FAIL: Z72 box_center_re drift vs source: cert='-1.40115560339843378903373358727' source='-1.40115560629999364024782953707' |
| T3 | newton_dist_scaled | Z17 newton_dist_ub forged 0.9x (plausibly tighter) | 1 | YES | FAIL: Z17 newton_dist_ub drift: reported=[2.27050889452687e-30 +/-  recomputed~[2.52278766058539e-30 +/-  (rel diff > 1e-6) |
| T4 | kappa_scaled | Z40 krawczyk_kappa_ub forged 0.01x (plausibly tighter) | 1 | YES | FAIL: Z40 kappa_ub drift: reported=[4.12408228992031e-5 +/- 5 recomputed~[0.00412408228991833 +/- 2 (rel diff > 1e-6) |
| T5 | period_check_dropped | Z63: j=10 (d=1024) period check removed | 1 | YES | FAIL: Z63: period_checks length 10 != m+1 = 11 |
| T6 | period_downgrade | Z19: n 2048->1024, m 10->9, period checks trimmed to j<=8 (self-consistent downgrade) | 1 | YES | FAIL: Z19 period n drift: cert=1024 source=2048 |
| T7 | pin_zeroed | sha256 pin of results.json zeroed | 1 | YES | FAIL: sha256 pin mismatch results.json: cert=0000000000000000000000000000000000000000000000000000000000000000 actual=021c378df159524c62cd5a694c258b19a9f0d7eeb22611c22b192876f8cf825d |
| T8 | pass_flag_flip | Z57 pass -> false while all_pass stays true | 1 | YES | FAIL: Z57 pass drift: cert=False recomputed=True |
| T9 | radius_tamper | Z66 box_radius -> "1e-29" (in-ladder value; only the recomputed radius-dependent fields can catch it) | 1 | YES | FAIL: Z66 kappa_ub drift: reported=[0.00439168125506909 +/- 5 recomputed~[0.439168126237115 +/- 2.9 (rel diff > 1e-6) |
| T10 | label_center_swap | Z68 <-> Z69 box centers swapped | 1 | YES | FAIL: Z68 box_center_re drift vs source: cert='-1.40115582181237413940961579485' source='-1.40115580341004315212881751609' |

Follow-up complaints (2nd FAIL line), for the record:
- T1: FAIL: Z05 f0_abs_ub drift: reported=[7.66113747142414e-24 +/-  recomputed~[1.99196557110729e-21 +/-  (rel diff > 1e-6)
- T2: FAIL: Z72 box_center_im drift vs source: cert='-0.00000062060634564680869303119869225' source='-0.0000006249490201460275312368120738'
- T6: FAIL: Z19 m drift: cert=9 source=10
- T8: FAIL: all_pass drift: cert=True recomputed=False
- T10: FAIL: Z68 box_center_im drift vs source: cert='-0.000000558427797053044061367699459913' source='-0.000000551874621760202881948411638124'

Summary: **10/10 tampers caught**; clean control exit 0. Battery PASSED.

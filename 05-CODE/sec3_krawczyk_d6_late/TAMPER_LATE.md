# TAMPER_LATE -- corruption battery for krawczyk6_late_certificates.json [oid:kraw6L]

Date: 2026-08-15 13:20:49

Subject: 9 late-subfamily Krawczyk certificates (sextic f_a(x)=1-a*x^6, periods 2048/4096/8192; refined seeds from refine_late.py) produced by `krawczyk6_late_certificate.py` (Arb acb ball path) and verified by the independent `krawczyk6_late_verify.py` (hand-rolled (re,im,r) tuple disk path). Principle: every tamper lies plausibly -- internally consistent formats and magnitudes -- so only recomputation or exact source-pin reconciliation can catch it.

Engines: **A** = rebuild via the generator with corrupted input (caught = no ladder configuration certifies); **B** = mutated certificate JSON (caught = verifier exit 1).

Controls (must pass): OK
- C1 clean_rebuild: OK -- certified at r=1e-26 prec=256 (!)
- C2 clean_verify: OK -- exit=0

| # | tamper | engine | what was corrupted | caught | detail / first verifier complaint |
|---|--------|--------|--------------------|--------|----------------------------------|
| T1a | center_shift(rebuild) | A | refined center Re + 1e-20 (>> box radius 1e-26): true root leaves the box, K(X) not inside X at all radii | YES | N(X) not inside X; K(X) not inside X (r=1e-26) |
| T2 | foreign_seed(rebuild) | A | S104 certified with S031's center (distance ~2e-9 vs radius 1e-26): root outside every box | YES | exact-period lower bound = 0 (r=1e-26) |
| T3a | period_1024(rebuild) | A | claim period 1024 for the S031 period-2048 root: Newton/Krawczyk distance >> r | YES | N(X) not inside X; K(X) not inside X (r=1e-26) |
| T3b | period_4096(rebuild) | A | claim period 4096 for the S036 period-8192 root: the box contains the period-8192 root, exact-period lower bound at d=4096 collapses to 0 | YES | N(X) not inside X; K(X) not inside X (r=1e-26) |
| T4a | radius_1e-8(rebuild) | A | nearest new pair S104/S187 (sep ~2e-10) inflated to r=1e-8: boxes overlap and/or K(X) leaves X | YES | isolation(S104,S187)@1e-8=False, K-inside=False/False |
| T5 | mantissa_flip(pin) | B | S235 box_center_im last decimal digit +1 (~1e-70 shift, below box radius -- the EXACT string pin on refined_seeds_late.json catches it) | YES | FAIL: S235 box center drift vs refined seed |
| T1b | center_shift(pin) | B | mutated certificate JSON | YES | FAIL: S036 box center drift vs refined seed |
| T3c | period_checks_trimmed | B | mutated certificate JSON | YES | FAIL: S187: period_checks length 11 != m+1 = 12 |
| T4b | radius_inflate | B | mutated certificate JSON | YES | FAIL: S031: box_radius outside sanity band (0, 1e-20]: [1.00000000000000e-8 |
| T6 | pin_zeroed | B | mutated certificate JSON | YES | FAIL: sha256 pin mismatch refined_seeds_late.json: cert=0000000000000000 actual=a035a35ae455abcc |
| T7 | period_mark_forged | B | mutated certificate JSON | YES | FAIL: S031 period_checks[j=3].abs_lb drift: reported/recomputed ratio [0.0010000000000 outside [[0.90909, 1.1] |
| T8 | pass_flag_flip | B | mutated certificate JSON | YES | FAIL: S339 pass drift: cert=False recomputed=True |
| T9 | period_check_dropped | B | mutated certificate JSON | YES | FAIL: S279: period_checks length 10 != m+1 = 11 |
| T10 | center_swap_S104_S187 | B | mutated certificate JSON | YES | FAIL: S104 box center drift vs refined seed |
| T11 | u_cat_drift | B | mutated certificate JSON | YES | FAIL: S031 u_cat drift vs source |
| T12 | ladder_forged | B | mutated certificate JSON | YES | FAIL: frozen radii_ladder drift: ['1e-30', '1e-29', '1e-28'] |

Summary: **16/16 tampers caught**; controls OK. Battery PASSED.

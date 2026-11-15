# TAMPER_LOG -- corruption battery for krawczyk6_certificates.json [oid:kraw6H]

Date: 2026-08-15 10:52:55

Subject: 287 Krawczyk certificates (sextic f_a(x)=1-a*x^6, period 1024; 285 source-string seeds + 2 refined-layer seeds S078/S090) produced by `krawczyk6_certificate.py` (Arb acb ball path) and verified by the independent `krawczyk6_verify.py` (hand-rolled (re,im,r) tuple disk path). Principle: every tamper lies plausibly -- internally consistent formats and magnitudes -- so only recomputation or exact source-pin reconciliation can catch it.

Engines: **A** = rebuild via the generator with corrupted input (caught = no ladder configuration certifies); **B** = mutated certificate JSON (caught = verifier exit 1).

Controls (must pass): OK
- C1 clean_rebuild: OK -- certified at r=1e-24 prec=256 (!)
- C2 clean_verify: OK -- exit=0

| # | tamper | engine | what was corrupted | caught | detail / first verifier complaint |
|---|--------|--------|--------------------|--------|----------------------------------|
| T1a | seed_shift(rebuild) | A | seed a_ss Re + 1e-20: true root leaves every box, K(X) not inside X at all radii | YES | N(X) not inside X; K(X) not inside X (r=1e-24) |
| T2 | family_coeff(rebuild) | A | orbit of g_a(x)=1-(1+1e-9)a*x^6 at 4 seeds: no ladder configuration certifies | YES | S001:fail; S100:fail; S200:fail; S287:fail |
| T3a | period_512(rebuild) | A | claim period 512 for a period-1024 root: Newton/Krawczyk distance >> r | YES | N(X) not inside X; K(X) not inside X (r=1e-24) |
| T3b | period_2048(rebuild) | A | claim period 2048: the box contains the period-1024 root, exact-period lower bound at d=1024 collapses to 0 | YES | exact-period lower bound = 0 (r=1e-24) |
| T4a | radius_1e-8(rebuild) | A | nearest pair S143/S145 (sep 2.56e-10) inflated to r=1e-8: boxes overlap (pairwise isolation violated) and/or K(X) leaves X | YES | isolation(S143,S145)=False, K-inside=False/False |
| T5 | fake_parameters(rebuild) | A | rebuild at predicted positions of the float-rejected zeros S036/S104 (anomalous subfamily), S070 (secondary lattice) and control CTRL1: none certifies | YES | REJ:S036:fail; REJ:S104:fail; REJ:S070:fail; CTRL1:fail |
| T1b | seed_shift(pin) | B | mutated certificate JSON | YES | FAIL: S005 box center drift vs source seed: cert=('1.68326020501997666794689551783','0.0') source=('1.68326020 |
| T3c | period_downgrade | B | mutated certificate JSON | YES | FAIL: S019 period/m drift: cert=(512,8) source=(1024,9) |
| T4b | radius_inflate | B | mutated certificate JSON | YES | FAIL: S001: box_radius outside sanity band (0, 1e-20]: [1.00000000000000e-8 |
| T6 | box_swap | B | mutated certificate JSON | YES | FAIL: S084 box center drift vs source seed: cert=('1.68326021644426793214326653246','3.1106885264870566e-08') |
| T7 | pin_zeroed | B | mutated certificate JSON | YES | FAIL: sha256 pin mismatch match_results.json: cert=0000000000000000 actual=6c367ac1176a193c |
| T8 | kappa_forged | B | mutated certificate JSON | YES | FAIL: S006 kappa_ub drift: reported/recomputed ratio [0.0100000000024 outside [[0.33333, 3] |
| T9 | period_check_dropped | B | mutated certificate JSON | YES | FAIL: S287: period_checks length 9 != m+1 = 10 |
| T10 | pass_flag_flip | B | mutated certificate JSON | YES | FAIL: S030 pass drift: cert=False recomputed=True |
| T11 | refined_center_shift | B | mutated certificate JSON | YES | FAIL: S078 refined center drift vs refined_seeds6.json |
| T12 | seed_dist_forged | B | mutated certificate JSON | YES | FAIL: S090 seed_dist_ub drift: reported=[1.85371711575298286588318 recomputed~[1.85371711575298281658239 (rel |

Summary: **16/16 tampers caught**; controls OK. Battery PASSED.

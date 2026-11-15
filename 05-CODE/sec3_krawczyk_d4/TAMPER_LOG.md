# TAMPER_LOG -- corruption battery for krawczyk4_certificates.json [oid:kraw4E]

Date: 2026-08-15 10:19:31

Subject: 61 Krawczyk certificates (quartic f_a(x)=1-a*x^4, period 1024) produced by `krawczyk4_certificate.py` (Arb acb ball path) and verified by the independent `krawczyk4_verify.py` (hand-rolled (re,im,r) tuple path). Principle: every tamper lies plausibly -- internally consistent formats and magnitudes -- so only recomputation or exact source-pin reconciliation can catch it.

Engines: **A** = rebuild via the generator with corrupted input (caught = no ladder configuration certifies); **B** = mutated certificate JSON (caught = verifier exit 1).

Controls (must pass): OK
- C1 clean_rebuild: OK -- certified at r=1e-26 prec=256 (!)
- C2 clean_verify: OK -- exit=0

| # | tamper | engine | what was corrupted | caught | detail / first verifier complaint |
|---|--------|--------|--------------------|--------|----------------------------------|
| T1a | seed_shift(rebuild) | A | seed a_ss Re + 1e-20: true root leaves every box, K(X) not inside X at all radii | YES | N(X) not inside X; K(X) not inside X (r=1e-26) |
| T2 | family_coeff(rebuild) | A | orbit of g_a(x)=1-(1+1e-9)a*x^4 at 4 seeds: no ladder configuration certifies | YES | Q01:fail; Q28:fail; Q39:fail; Q62:fail |
| T3a | period_512(rebuild) | A | claim period 512 for a period-1024 root: x_512(c0)~1e-2, Newton/Krawczyk distance >> r | YES | N(X) not inside X; K(X) not inside X (r=1e-26) |
| T3b | period_2048(rebuild) | A | claim period 2048: the box contains the period-1024 root, exact-period lower bound at d=1024 collapses to 0 | YES | exact-period lower bound = 0 (r=1e-26) |
| T4a | radius_1e-8(rebuild) | A | nearest pair Q51/Q52 inflated to r=1e-8: boxes overlap (pairwise isolation violated) and/or K(X) leaves X | YES | isolation(Q51,Q52)=False, K-inside=False/False |
| T5 | fake_parameters(rebuild) | A | rebuild at the rejected Q57 predicted position and at CTRL1/CTRL2 control points: none certifies | YES | Q57pred:fail; CTRL1:fail; CTRL2:fail |
| T1b | seed_shift(pin) | B | mutated certificate JSON | YES | FAIL: Q05 box center drift vs source seed: cert=('1.59490137873494180679273309559','0.000000066035182006668247 |
| T3c | period_downgrade | B | mutated certificate JSON | YES | FAIL: Q19 period/m drift: cert=(512,8) source=(1024,9) |
| T4b | radius_inflate | B | mutated certificate JSON | YES | FAIL: Q01: box_radius outside sanity band (0, 1e-10]: [1.00000000000000e-8 |
| T6 | box_swap | B | mutated certificate JSON | YES | FAIL: Q02 box center drift vs source seed: cert=('1.59490141222668164936124504894','0.0') source=('1.594901394 |
| T7 | pin_zeroed | B | mutated certificate JSON | YES | FAIL: sha256 pin mismatch match_results.json: cert=00000000000000000000000000000000000000000000000000000000000 |
| T8 | kappa_forged | B | mutated certificate JSON | YES | FAIL: Q06 kappa_ub drift: reported=[4.22415405557193e-6 +/- 3 recomputed~[0.000422415405550884 +/-  (rel diff |
| T9 | period_check_dropped | B | mutated certificate JSON | YES | FAIL: Q61: period_checks length 9 != m+1 = 10 |
| T10 | pass_flag_flip | B | mutated certificate JSON | YES | FAIL: Q30 pass drift: cert=False recomputed=True |

Summary: **14/14 tampers caught**; controls OK. Battery PASSED.

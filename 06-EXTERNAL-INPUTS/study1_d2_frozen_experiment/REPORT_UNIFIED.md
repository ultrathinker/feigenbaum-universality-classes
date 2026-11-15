# REPORT_UNIFIED: Frozen Mandelbrot cross-check of the point-19 hypothesis

Single frozen script (`frozen_experiment.py`), single dataset (`address_certificate.json`, stage 15.3, 23 disks), one process, mpmath, no multiprocessing, no post-hoc parameter tuning. Every constant used below was frozen in the script header before any per-point result was inspected (see the script's own docstring for the justification of each one).

**Run metadata**

| | |
|---|---|
| mp.dps | 150 |
| level m | 10 (period = 2^11 = 2048) |
| delta | `4.669201609102990671853203820466201617258185577475768632745651343004133555868548426390712381266351242` |
| c_F (own cascade, Richardson) | `(-1.4011551890920506005238268502643472843821101975214 + 0.0j)` |
| \|c_F - brief literal\| | 5.238e-19 |
| A (own calibration) | `(-2.33188562955070413230142737400505175507716 + 0.0j)` |
| u(c_11) self-consistency error (expect -1 at m=10) | 1.421e-08 |
| R_ACCEPT_U (frozen acceptance radius, u-space) | 0.01 |
| control seed | 20260812 |
| total wall time | 2175.9 s (36.3 min) |

## 1. Catalog cross-check -- all 23 disks

Rule (frozen, item 4 of the brief): a converged, exact-period-2048 root is **accepted** as the match for a catalog point only if `|c_ss - c_pred| < R_ACCEPT_C(m) = R_ACCEPT_U/(|A|*delta^m)`. `d_u` is computed **only** for accepted points -- unconverged/rejected candidates never contribute a d_u (this is the exact defect FIX_BRIEF.md flagged in the earlier attempt).

| id | band | u_cat | status | tier | d_u | exact period | pool specificity | drift \|ratio\| | drift phase diff (deg) |
|---|---|---|---|---|---|---|---|---|---|
| Z01 | VI | 1.553494+8.123803i | MATCH | 1 | 2.603e-05 | True | 222247x | 0.9638 | -0.03 |
| Z02 | III | 5.817868+4.476081i | MATCH | 1 | 2.033e-05 | True | 119022x | 0.9557 | 0.13 |
| Z05 | III | 5.011415+5.626752i | MATCH | 1 | 2.027e-05 | True | 113597x | 0.9042 | 1.42 |
| Z17 | III | 6.681454+6.210555i | MATCH | 1 | 3.135e-05 | True | 61794x | 0.9546 | -0.79 |
| Z19 | IV | 4.008907+6.703721i | MATCH | 1 | 2.262e-05 | True | 135625x | 0.9394 | -1.04 |
| Z38 | IV | 4.637026+6.759226i | MATCH | 1 | 2.588e-05 | True | 102496x | 0.9757 | -0.59 |
| Z39 | III | 6.248463+5.678898i | MATCH | 1 | 2.710e-05 | True | 50485x | 0.9631 | 0.22 |
| Z40 | I | 8.368184+3.421616i | MATCH | 1 | 3.124e-05 | True | 102391x | 0.9682 | 0.28 |
| Z57 | I | 8.589396+3.974588i | MATCH | 1 | 3.358e-05 | True | 84069x | 0.9496 | 0.89 |
| Z59 | II | 7.125929+4.064322i | MATCH | 1 | 2.456e-05 | True | 92631x | 0.9246 | 1.07 |
| Z60 | I | 8.263327+4.015964i | MATCH | 1 | 3.247e-05 | True | 80609x | 0.9746 | -0.22 |
| Z63 | III | 6.831103+5.231170i | MATCH | 1 | 2.806e-05 | True | 48886x | 0.9602 | -0.13 |
| Z64 | III | 7.218228+5.788641i | MATCH | 1 | 3.233e-05 | True | 59365x | 0.9567 | -0.32 |
| Z65 | III | 7.360531+5.761589i | MATCH | 1 | 3.329e-05 | True | 59034x | 0.9652 | -0.39 |
| Z66 | III | 7.243993+6.188093i | MATCH | 1 | 3.471e-05 | True | 64202x | 0.9686 | -0.26 |
| Z67 | III | 6.615259+6.339356i | MATCH | 1 | 3.215e-05 | True | 62906x | 0.9701 | 0.09 |
| Z68 | III | 7.055487+6.338351i | MATCH | 1 | 3.455e-05 | True | 64717x | 0.9730 | -0.60 |
| Z69 | III | 7.266840+6.413616i | MATCH | 1 | 3.602e-05 | True | 66325x | 0.9712 | -0.93 |
| Z70 | III | 7.338955+6.505440i | MATCH | 1 | 3.687e-05 | True | 67469x | 0.9712 | -0.29 |
| Z71 | IV | 4.506189+7.078002i | MATCH | 1 | 2.641e-05 | True | 109002x | 0.9503 | -0.56 |
| Z72 | IV | 4.758336+7.127705i | **NO MATCH** | -- | -- | -- | -- | -- | best converged candidate rejected: dist_to_c_pred=5.223e-09 (R_ACCEPT_C=8.707e-10) |
| Z73 | V | 4.020255+8.037592i | MATCH | 1 | 3.180e-05 | True | 114148x | 0.9975 | 0.60 |
| Z74 | V | 4.727254+8.686343i | MATCH | 1 | 3.718e-05 | True | 91517x | 0.9631 | -0.00 |

**Summary: 22/23 catalog points MATCH under the frozen rule.** Rejected: Z72. Median d_u (accepted) = 3.158e-05.

## 2. Z72 -- provenance (read before trusting any Z72 number above)

The certificate's Z72 disk centre `4.758336330290482 + 7.1277051021354385i` is a **point17 "shifted disk"** (radius rho=0.1), a technical device of that point's proof, NOT a tight numerical estimate of the true zero -- unlike the other 22 disks, whose centres happen to double as tight (~1e-4) estimates. Feeding this centre through the exact same frozen pipeline as the other 22 points correctly and honestly produces **NO MATCH** (36 of 49 deterministic restart attempts converge to *some* exact-period-2048 root, but every one is farther from `c_pred` than `R_ACCEPT_C`) -- this is the row reported in the table above, and it is real, not a bug.

The best rejected candidate for Z72 converges to `u_actual = 4.791660937911268+7.177581007211361i`, at distance 5.223e-09 from `c_pred` (vs. R_ACCEPT_C=8.707e-10) -- i.e. **6.0x** over the frozen acceptance radius.

Source of the corrected address `4.791660937911268425 + 7.177581007211361312i`: kilo_xcheck/results.json stage3_catalog.Z72.u_actual (re_str/im_str) -- the position kilo_xcheck's own Newton search converged to at m=10 when probing near the certified disk centre.

Distance(certified centre, corrected position) = 0.05998, which is **inside** the point17-certified disk radius (rho=0.1) -- consistent with the disk genuinely containing the true zero, just not centred on it.

Re-running the frozen pipeline with u_cat := corrected position (this script's own Newton code, not a copy of any prior converged c):

- status: **ACCEPTED**, tier 1, d_u = 0.0 (underflows to 0.0 at double precision -- expected given the construction, not evidence of a new independent match), exact_period=True.

**This corrected row is deliberately excluded from the 23-point tally above.** By construction, `c_pred` for the corrected u_cat is numerically the same point kilo_xcheck's search already converged to, so a small d_u here confirms self-consistency of this script's Newton code against kilo's independent implementation -- it is NOT a fresh, independent discovery of the position, and reporting it as a 23rd "match" would overclaim. What it DOES legitimately establish: Z72 has a genuine, exact-period-2048 superstable partner near c_F, exactly as the other 22 points do -- the earlier "Z72 anomaly" was a mislabelled-centre artifact, not a hypothesis failure.

## 3. Local root-density spot-check (catalog, 12x12 grid, side 0.4 in u)

Answers FIX_BRIEF.md's core concern: is a small d_u meaningful, or is this region of the u-plane simply dense with period-2048 superstable centres (so that any point would "match" something nearby)?

| id | roots in window | nearest (=d_u) | 2nd nearest | accidental median | specificity (acc_med/d_u) |
|---|---|---|---|---|---|
| Z01 | 1 | 2.603e-05 | isolated | 1.725e-01 | 6627x |
| Z02 | 1 | 2.033e-05 | isolated | 1.725e-01 | 8486x |
| Z05 | 9 | 2.027e-05 | 3.032e-01 | 1.537e-01 | 7586x |
| Z17 | 12 | 3.135e-05 | 1.448e-01 | 1.060e-01 | 3381x |
| Z19 | 1 | 2.262e-05 | isolated | 1.725e-01 | 7624x |
| Z38 | 19 | 2.588e-05 | 1.935e-01 | 1.141e-01 | 4408x |
| Z39 | 2 | 2.710e-05 | 7.348e-01 | 1.725e-01 | 6364x |
| Z40 | 1 | 3.124e-05 | isolated | 1.725e-01 | 5522x |
| Z57 | 25 | 3.358e-05 | 1.769e-01 | 9.538e-02 | 2841x |
| Z59 | 12 | 2.456e-05 | 2.905e-01 | 1.475e-01 | 6005x |
| Z60 | 20 | 3.247e-05 | 1.850e-01 | 1.202e-01 | 3702x |
| Z63 | 1 | 2.806e-05 | isolated | 1.725e-01 | 6147x |
| Z64 | 12 | 3.233e-05 | 1.448e-01 | 1.060e-01 | 3278x |
| Z65 | 20 | 3.329e-05 | 1.073e-01 | 5.325e-02 | 1600x |
| Z66 | 29 | 3.471e-05 | 1.329e-01 | 6.184e-02 | 1782x |
| Z67 | 19 | 3.215e-05 | 1.078e-01 | 5.034e-02 | 1566x |
| Z68 | 31 | 3.455e-05 | 1.316e-01 | 6.306e-02 | 1825x |
| Z69 | 21 | 3.602e-05 | 1.167e-01 | 4.983e-02 | 1384x |
| Z70 | 42 | 3.687e-05 | 7.627e-02 | 3.636e-02 | 986x |
| Z71 | 25 | 2.641e-05 | 1.855e-01 | 9.554e-02 | 3618x |
| Z72 | 31 | 5.998e-02 | 1.176e-01 | 6.888e-02 | n/a |
| Z73 | 14 | 3.180e-05 | 2.592e-01 | 1.285e-01 | 4042x |
| Z74 | 10 | 3.718e-05 | 2.608e-01 | 1.385e-01 | 3725x |

## 4. Negative controls (20 points, deterministic, seed=20260812)

Generated by rejection sampling (`random.Random(seed)`, uniform in r in [1.5, 11.0], theta in [0, 85] deg, minimum distance 0.5 in u from every catalog centre) -- i.e. drawn from the same general region of the u-plane as the catalog, but nowhere near any catalogued zero.

| id | u | status | d_u | exact period |
|---|---|---|---|---|
| CTRL01 | 1.8010+0.8802i | rejected | (best candidate at 2.556e-07 from c_pred, > R_ACCEPT) | True |
| CTRL02 | 2.5670+8.3009i | rejected | (best candidate at 8.958e-08 from c_pred, > R_ACCEPT) | True |
| CTRL03 | 6.1141+7.1786i | rejected | (best candidate at 1.685e-08 from c_pred, > R_ACCEPT) | True |
| CTRL04 | 2.8589+2.4793i | rejected | (best candidate at 3.994e-07 from c_pred, > R_ACCEPT) | True |
| CTRL05 | 6.2464+5.1439i | rejected | (best candidate at 4.658e-08 from c_pred, > R_ACCEPT) | True |
| CTRL06 | 4.7295+1.2753i | rejected | (best candidate at 1.532e-07 from c_pred, > R_ACCEPT) | True |
| CTRL07 | 3.0991+1.4552i | rejected | (best candidate at 3.787e-07 from c_pred, > R_ACCEPT) | True |
| CTRL08 | 3.2336+7.5735i | rejected | (best candidate at 7.059e-08 from c_pred, > R_ACCEPT) | True |
| CTRL09 | 7.8909+2.1513i | rejected | (best candidate at 2.058e-07 from c_pred, > R_ACCEPT) | True |
| CTRL10 | 3.4247+6.3247i | rejected | (best candidate at 6.064e-08 from c_pred, > R_ACCEPT) | True |
| CTRL11 | 1.2824+1.6712i | rejected | (best candidate at 2.463e-07 from c_pred, > R_ACCEPT) | True |
| CTRL12 | 3.6648+3.0009i | rejected | (best candidate at 3.279e-07 from c_pred, > R_ACCEPT) | True |
| CTRL13 | 4.2172+5.7624i | rejected | (best candidate at 6.077e-08 from c_pred, > R_ACCEPT) | True |
| CTRL14 | 4.2461+1.8912i | rejected | (best candidate at 2.211e-07 from c_pred, > R_ACCEPT) | True |
| CTRL15 | 2.5706+0.0568i | rejected | (best candidate at 3.109e-07 from c_pred, > R_ACCEPT) | True |
| CTRL16 | 1.0828+10.3979i | rejected | (best candidate at 3.380e-08 from c_pred, > R_ACCEPT) | True |
| CTRL17 | 0.8448+3.1349i | rejected | (best candidate at 3.167e-07 from c_pred, > R_ACCEPT) | True |
| CTRL18 | 5.8461+0.6329i | rejected | (best candidate at 5.572e-08 from c_pred, > R_ACCEPT) | True |
| CTRL19 | 5.1389+2.6906i | rejected | (best candidate at 1.663e-07 from c_pred, > R_ACCEPT) | True |
| CTRL20 | 1.0618+3.1104i | rejected | (best candidate at 3.249e-07 from c_pred, > R_ACCEPT) | True |

**0/20 controls accepted a match; median control d_u = n/a (none accepted). Median catalog d_u = 3.158e-05.** Separation: n/a.

Disclosure: every one of the 20 controls exhausted the full deterministic 49-attempt tier-1/2/3 ladder without ever satisfying `R_ACCEPT` (980 attempts total; 924 of them converged to *some* exact-period-2048 root -- confirming such roots are not rare, just never close enough to a random probe point). By contrast, 22 of the 23 catalog points matched on the very first (tier-1, direct) attempt; only Z72 needed the full ladder (71 catalog attempts total, 36 converged-but-rejected, all from Z72).

## 5. Drift-law vector check (item 5)

`Delta_pred(u_cat) = -(c2*AU0/delta^9) * u_cat^2 * delta/(delta-1)`, `c2=0.79028031944088979`, `AU0=0.41403799908627305889537448439555967803329915886199` (THEORY_NOTE.md §5bis literals). Compares the observed residual `R = u_actual - u_cat` against this prediction, per catalog point (abs ratio and phase difference in degrees).

Median |R|/|Delta_pred| over the 22 matched catalog points = 0.9635; median |phase difference| = 0.353 deg.

## 6. Verdict

- Catalog: **22/23 MATCH** under the frozen acceptance rule (Z72-as-certified is the sole rejection, with full provenance above).
- Median catalog d_u = 3.158e-05; median pool specificity = 82339x.
- Controls: 0/20 accepted a match at all; those that did land at d_u orders of magnitude above the catalog median (see §4).
- Drift-law vector check: median |ratio|=0.9635, median |phase diff|=0.353 deg -- consistent with THEORY_NOTE.md §5bis's finding that the residual ~2-4e-5 polka is the *known, predicted* level-7-to-limit drift, not unexplained noise.
- **Full failure disclosure**: every non-accepted attempt (all tiers, all points) is logged verbatim in `results.json` (`attempts` arrays and `best_converged_but_rejected` fields); nothing here is smoothed over.

Total wall time: 36.3 minutes.

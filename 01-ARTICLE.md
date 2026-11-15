# Explicit realisation of the "parameter ↔ unstable manifold" correspondence for period-doubling universality classes

**Author:** Yevhen Borzenkov

**Affiliation:** independent researcher; correspondence — author's personal email (to be supplied at submission).

**Target journal:** Experimental Mathematics

**Keywords:** Feigenbaum universality; renormalization group; period doubling; unstable manifold; scaling function; superstable parameters; corrections to scaling; interval arithmetic; Krawczyk operator.

## Abstract

For the period-doubling classes f_a(x) = 1 − a·|x|^d, d ∈ {2, 4, 6}, the paper gives an explicit computable realisation of the "parameter space ↔ unstable manifold W^u(g_d)" correspondence: the zero set of the scaling function H_d(u) is a catalogue of superstable parameters of period 2^{m+1} along W^u. The correspondence is a consequence of Feigenbaum universality (Feigenbaum 1978, 1979; Coullet–Tresser 1978; for d = 2 — Milnor 1989, Lyubich 1999); the contribution of the paper is the explicit H_d, its verification at d = 4, 6 with interval-certified parameters and the rate of convergence, including the λ₂^m layer.

The correspondence is confirmed on 22/23 zeros (d = 2; external input from a companion paper, not rechecked here), 61/62 (d = 4) and 287/294 (d = 6, m = 9), with nine parameters at m = 10–12. The rigorous part — Krawczyk certificates for 357 parameters and numerically verified λ₂(4) (N = 40), λ₂(6) (N = 48); λ₂(2) — separately, ±1.8·10⁻³³. The law rate = min(δ, 1/|λ₂|) is a re-discovery of the corrections-to-scaling structure (Mao–Hu 1987; Reick 1992; Briggs 1994); a quantitative blind confirmation — d = 6, for d = 4 at the level of order and sign (no point prediction because of the discrepancy in a_F(4)). Out-of-sample verification of the remainder confirmed at d = 6 (92.0 % vs baseline 90.5 %); for d = 4 (n = 61) the two-term model does not show superiority over baseline in the out-of-sample test. The attraction of zeros of H to the 1/16 lattice is a scanner-selection artefact, shown for both classes by a single protocol.

**Genre and rigor.** The "zero ↔ parameter" correspondence is verified numerically (NOT RIGOROUS — for the values of the d_u metric, controls and rate of convergence); the rigorous part is only the theorem on Krawczyk boxes of §4.1 and the Proposition on the gauge tower of §2.5. There is no formal theorem about H_d — the correspondence and the rate of convergence remain numerically verified rather than proved.

## 1. Introduction and contribution

Feigenbaum universality — the infinite sequence of period-doubling bifurcations in one-parameter families of one-dimensional maps with an even maximum converges to a limit point c_F with geometric rate δ, and near it the orbits of different families look the same after suitable scaling. The formal carrier of this sameness is the universal scaling function H(u): the fixed point of the doubling operator T, to which the sequence of renormalised iterates g_m(x) of the map f_a converges as a → c_F. Existence and hyperbolicity of the fixed point for d = 2 — Lanford (1982); for even d ≥ 4 — Epstein (1986) (existence), Lyubich (1999) and Avila–Lyubich (2011) (uniqueness and hyperbolicity); CAP bounds at d = 4 — Burbanks–Osbaldestin–Thurlby (2021).

A standard fact of the theory is the correspondence between parameter space of the class and the unstable manifold W^u(g_d) in the linearising coordinate (Feigenbaum 1978, 1979; Coullet, Tresser 1978); for d = 2 — a standard result (Milnor 1989; Lyubich 1999).

The classical theory of period-doubling renormalisation is also developed in Sullivan (1992), McMullen (1996), Collet–Eckmann (1980, book), de Melo–van Strien (1993) and Vul–Sinai–Khanin (1984) — the general background of the present work; the detailed content of those works is neither used nor paraphrased here.

**The contribution of this work** is not a re-discovery of the correspondence itself, but its **explicit computable realisation**: the universal scaling function H_d(u), whose zero set serves as a catalogue of points of W^u(g_d) in the linearising coordinate, identified with superstable parameters of period 2^{m+1}. The catalogue is built by an explicit recipe (§2.2, §2.6) and verified empirically on the classes d = 4, 6 with specificity controls (§3); it is backed by interval Krawczyk certificates for the superstable parameters themselves (§4). Separately verified is the rate of convergence: it is given not by δ alone but by a two-term law

    e(m) ~ −(C·δ^(−m) + D·λ₂^m),

in which the slower of the two terms dominates (§5). The structure of this correction has been known since 1987–1994 (Mao–Hu 1987; Reick 1992 — "crossover"; Briggs 1994) — the paper independently re-discovered it and presents it here as **a re-discovery, modern certification and extension**, not as a new result.

In addition to the two main lines: interval Krawczyk certificates for 357 superstable parameters (§4); numerically verified λ₂(2), λ₂(4), λ₂(6) ("Bridge"/"Bridge-2", §4.4); an extended but uncertified λ₂(d), δ(d) curve over d ∈ [1.3, 20] (for λ₂ skipping d = 8; see §6.3) with a resolution of the sign-change question of λ₂(d) (§6); the Proposition on the gauge tower (§2.5). Sections S1–S9 are in the Supplement; source codes — README-CODE.md.

---

## 2. Definitions and constants

### 2.1 Universal scaling function H(u) and its zeros

Let f_a be a one-parameter family of one-dimensional maps with a single maximum, undergoing a Feigenbaum period-doubling cascade as a → c_F. The doubling operator T: (Tg)(x) = −α·g(g(x/α)), with α chosen so that (Tg)(0) = 1 when g(0) = 1. The sequence of renormalised iterates of the cascade g_m(x) converges as m → ∞ along the unstable direction of T to the fixed point H:

    g_m(x) → H(x),    TH = H,    H(0) = 1.

H is the universal scaling function of the class, determined by the order of the maximum d up to scale. Existence and hyperbolicity for d = 2 — Lanford (1982); for general even d — Epstein (1986), Lyubich (1999), Avila–Lyubich (2011); CAP bounds at d = 4 — Burbanks–Osbaldestin–Thurlby (2021).

The central object of the paper is the zero set of H(u), where u is the coordinate on the unstable manifold W^u(g_d) in the linearising parametrisation (§2.2). Besides the substantive zeros, H(u) carries a coordinate series — the geometric tower κ_k = α^(−d·k), k = 1, 2, 3, … Only the substantive zeros enter the matching count (§3); the tower ones are excluded (justification — §7.1).

### 2.2 Classes f_a(x) = 1 − a·|x|^d, superstable parameters and H_d(u)

The paper considers the family of universality classes generated by the maps

    f_a(x) = 1 − a·|x|^d,    d = 2, 4, 6, …

(for even d the modulus can be dropped: |x|^d = x^d). Each family at a suitable a undergoes a period-doubling cascade to a limit point c_F(d).

A parameter a is called **superstable of level m** if the orbit of the point x = 0 is periodic with minimal period P = 2^(m+1) (the critical orbit does not return to 0 at the intermediate steps 2, 4, …, P/2, §4.1(iv)); m = 0 — period 2, m = 1 — period 4.

Scaled coordinate: for a superstable a_i of level m

    u_m(a_i) = A·(a_i − c_F)·δ^m,

where A is an amplitude constant of the class, calibrated so that the zero u = −1 of H_d corresponds to the calibration cascade (values of A, c_F, δ for the three classes — Table T2.0, §2.6).

**Formal construction of H_d(u).** The fixed point g_d is perturbed along its eigenvector h_d, the family g_d + s·h_d with s = u·AU0_d/δ^m; the orbit x_{n+1} = g_d(x_n) + s·h_d(x_n), x_0 = 0; at level m define

    H_{d,m}(u) = x_{2^(m+1)} / x_{2^m},

and H_d(u) := lim_{m→∞} H_{d,m}(u) (a frozen finite m in practice, §2.6). The catalogue of zeros of H_d is built by a rough scan over u followed by polishing with Newton's method in mpmath; construction — §2.6.

**Complex extension.** **Holomorphic continuation** is used: the fixed point and its perturbation are realised as power series in t = z^d for complex z with integer d (in the code — `t = z**4`, `t = z**6`), i.e. the continuation is analytic, not |z|^d.

The correspondence (a hypothesis, verified in §3): the zero set {u : H_d(u) = 0} and the set of scaled superstable parameters {u_m(a_i)} converge pointwise as m → ∞.

### 2.3 Spectrum of the linearisation of the doubling operator: δ, λ₂, λ₃, …

The linearisation DT_H of T at the fixed point H is a linear operator on the space of perturbations; its spectrum consists of the leading eigenvalue δ (unique, real, > 1, corresponding to the unstable direction along which the cascade is built — the classical Feigenbaum constant of the class) and a subleading spectrum λ₂, λ₃, … (in decreasing modulus), corresponding to decaying directions.

**Table T2.1 — spectra of the three universality classes.**

| class d | α | δ | λ₂ | λ₃ | coordinate tower (k = 1) | external anchor |
|---|---|---|---|---|---|---|
| 2 | 2.502907875095892822284 | 4.66920160910299067185320382… | **B⁻** = −0.123652712552687024335 (numerically verified, ±1.8×10⁻³³, N = 48) | −0.0573070210666818468843650322586 | α₂^(−2) = 0.159628… | Christiansen–Cvitanović–Rugh 1990; Briggs 1994 |
| 4 | 1.6903029714052449 | 7.28468621707334336431 | **B⁺** = +0.29183840898354510453 (1/λ₂ = 3.4266; numerically verified, N = 40) | **B⁻** = −0.255664558474297406325347479489 | α₄^(−4) = 0.122501605427 | Burbanks–Osbaldestin–Thurlby 2021: δ₄ ∈ [7.28468621706, 7.28468621709], conditionally RIGOROUS; Briggs 1994; a preprint for degree 2 by the same group of authors separately — arXiv:2103.05991, not used here |
| 6 | 1.4677424503199009 | 9.29624683277137008283 | **B⁺** = +0.37107146768193013247 (1/λ₂ = 2.6949; numerically verified, N = 48) | **B⁻** = −0.341968322547802105099859066308 | α₆^(−6) = 0.10002322865… | Briggs 1994 |

The notations B⁺ and B⁻ are introduced because, as shown in §2.4, for d ≳ 1.9 the linearisation spectrum contains TWO real subleading branches — a positive B⁺ and a negative B⁻ — and which of them is called "λ₂" depends on d in a non-trivial way.

_Footnote to T2.1._ In the file `spectra_grid.json` the field `lam2` stores the B⁻ branch (the paper's λ₃), the field `lam3` the next mode (λ₄); the field correspondence is a table in App. B.1.1.

**Sign of λ₃(4).** An early pilot computation with a deflated power iteration (a method that breaks down on this class, Supplement §S2) gave the wrong sign of λ₃(4) while the modulus was correct; the sign was fixed by an independent full eig-diagonalisation.

**Unified scale of rigor statuses**, used throughout the paper:

1. **RIGOROUS** — interval certificate for an exactly posed problem without discretisation of the limiting object; example — existence, uniqueness and exact period of the 357 superstable parameters (§4).
2. **conditionally RIGOROUS** — interval or classical rigorous result under explicitly stated assumptions of the method, or an external interval anchor; example — δ(2) (lineage of Lanford 1982), δ(4) (Burbanks–Osbaldestin–Thurlby 2021).
3. **numerically verified** — certified eigenvalue of a collocation discretisation at a fixed N; the N → ∞ truncation error is estimated by a ladder of neighbouring N but is not certified interval-wise. This status applies to **all three** λ₂(2), λ₂(4), λ₂(6): λ₂(2) — N = 48 (±1.82·10⁻³³), λ₂(4) — N = 40, λ₂(6) — N = 48 ("Bridge-2", §4.4).
4. **NOT RIGOROUS** — exploratory computations (mpmath/float) without interval guarantees; example — the curve λ₂(d), δ(d) in §6 beyond the anchors of T2.1, the whole d_u layer in §3.

**Status of the rows of T2.1:** λ₂(2), λ₂(4), λ₂(6) — **numerically verified** (level 3). δ(2), λ₃(2), δ(4) — conditionally RIGOROUS (δ(4) — external anchor Burbanks–Osbaldestin–Thurlby 2021). δ(6), λ₃(4), λ₃(6) and the rest of the subleading spectrum of d = 6 — NOT RIGOROUS, validated by agreement with Briggs (1994). Working route — Supplement §S2.

### 2.4 Two subleading branches B⁺/B⁻ and pointwise gauge behavior at even d

**The non-gauge part of the spectrum** is the spectrum of the linearisation DT_H minus all eigenvalues generated by changes of coordinates. At every d ≳ 1.9 we investigated, the linearisation spectrum contains two real subleading branches: a positive B⁺(d) > 0 and a negative B⁻(d) < 0; over the whole interval d ∈ [1.9, 7] both are continuous, do not cross, and keep their sign (grid — §6.2). At fractional and odd d both branches are substantive as ordinary eigenvalues without any special role: λ₂(d) is the one larger in modulus, λ₃(d) the smaller; over the whole range except for the single point d = 2 it is B⁺.

**One of the branches coincides in modulus with the gauge value α^(−d) exactly at even d** — at d = 2 it is the second-in-modulus branch, B⁺; at d = 4, 6, 8 it is the d-th-in-modulus branch without counting δ. At these points that branch is a gauge mode, not a stand-alone substantive eigenvalue.

At d = 2 the role of the coordinate mode is taken by B⁺: B⁺(2) = α(2)^(−2) = 0.15962844038269977 coincides with α⁻² at the 10⁻³⁰ level: relative tracking residual 1.63·10⁻³⁰ (the `tower1.match_rel` field in `spectra_grid.json`); the difference of the 30-digit tower1 and α⁻ᵈ entries in B.5 is 4.9·10⁻³¹, i.e. at the level of the last stored digit. An independent check of multiplicity (Supplement §S2) showed: α^(−2) is a simple eigenvalue, B⁺ coincides with α⁻² as a (simple) eigenvalue — the coordinate mode Φ₁ takes the place of the subleading one. Hence λ₂(2) = −0.123652712552687024335 is B⁻, not B⁺. At any fractional d near 2 both branches are substantive and λ₂ = B⁺ > 0 — the sign jumps exactly at d = 2 (§6.2). At d = 4, 6, 8 the gauge value passes through a deeper branch (T2.2) — the labels λ₂/λ₃ are unaffected.

**Table T2.2 — gauge events at even d.**

| even d | gauge value α^(−d) | position without δ / with δ | n_>(d) | effect on the λ₂/λ₃ labels |
|---|---|---|---|---|
| 2 | α₂^(−2) = 0.15962844038269977 | 1 / 2 (boundary case — B⁺ coincides with α⁻² as a simple eigenvalue; the coordinate mode Φ₁ takes the place of the subleading one) | 0 | λ₂(2) = B⁻; the sign of λ₂ jumps exactly here |
| 4 | α₄^(−4) = 0.122501605427 | 4 / 5 | 3 | unchanged: λ₂(4) = B⁺, λ₃(4) = B⁻ |
| 6 | α₆^(−6) = 0.10002322865… | 6 / 7 | 5 | unchanged: λ₂(6) = B⁺, λ₃(6) = B⁻ |
| 8 | α₈^(−8) = 0.08644853405611… | 8 / 9 | 7 | unchanged (Note M6) |

**Observation (NOT RIGOROUS, at d = 4, 6, 8):** n_>(d) = d − 1 for even d ≥ 4; at d = 2 — a boundary case, n_>(2) = 0. Observation at three points, not an established pattern. The mechanism is the Proposition of §2.5.

### 2.5 Proposition on the gauge tower

Consider the class f(x) = F(|x|^d) = G(t), t = |x|^d, with normalisation G(0) = 1; the doubling operator T(f)(x) = −α_f·f(f(x/α_f)), α_f = −1/f(1); in the t-representation G(t) = −α·G(W(t)). The linearisation DT on perturbations Φ(t) with Φ(0) = 0 splits into gauge + J + S; the explicit formula and its finite-difference check against the nonlinear operator (residual ≤ 3·10⁻²⁷ for d ∈ {2, 3, 4, 6}) — source Appendix A, §A.1.

**Proposition (pointwise identity + analyticity criterion).** Let d > 1, α > 0, g be real on [−1, 1], k ≥ 1. Define w_k(x) := x·|x|^{dk} − x (in the x-representation; equivalently: Φ_k(t) at t = |x|^d) and

    Φ_k(t) := G(t)^{dk+1} − d·t^{k+1}·G′(t) − K(t),     Φ_k(0) = 0,      (1)

where K(t) := G(t) − d·t·G′(t). Then for **any** real dk > 0 for Φ_k in the form (1), understood on the odd real branch [G]^{dk+1} := sign(G)·|G|^{dk+1} when dk ∉ 2ℤ (the case of odd integer dk — separately, see the Corollaries below and the numerical check in A.5), the **pointwise identity**

    DTΦ_k(t) = α^{−dk}·Φ_k(t)     for all t ∈ [0, 1].      (2)

holds. When dk ∈ 2ℤ, Φ_k ∈ C^ω([0,1]) and is an analytic eigenvector of DT with eigenvalue α^{−dk}; the converse (Φ_k ∈ C^ω ⇒ dk ∈ 2ℤ) holds under simple zeros of G on (0, 1). On the entire branch G^{dk+1} with dk ∉ 2ℤ the identity (2) fails; on the odd branch it holds, but Φ_k ∉ C^ω.

**Rigor status (conditionally RIGOROUS).** Two sub-steps: (a) the **algebraic identity** DTΦ_k = α^{−dk}·Φ_k on the odd real branch is formally rigorous for a given fixed point G (proof — Appendix A §A.3, step (a)); (b) **existence and uniqueness of the fixed point g_d** at even d — from the literature (Epstein 1986; Lyubich 1999; Avila–Lyubich 2011). Identity (2) is a formal identity for the fixed point g_d; its numerical check was performed on a collocation approximation of g_d (residual ≤ 2·10⁻²⁰ for integer even dk and ≤ 8·10⁻²⁴ for fractional/odd dk). In total: **conditionally RIGOROUS** (formal identity for the fixed point g_d, whose existence at even d is covered by Epstein (1986), Lyubich (1999) and Avila–Lyubich (2011)); a rigorous estimate of the N → ∞ truncation error was not carried out.

**Corollaries.** d even integer ⇒ dk even for all k ⇒ the whole {α^{−dk}}_{k≥1} are eigenvalues of DT in the analytic class. d odd integer ⇒ the tower is sparse, the first level is k = 2. d irrational ⇒ there are no analytic members of the tower.

#### Exact result vs empirical continuation of the branch (d = 4 and d = 6)

The exact result is identity (2), proved above; this is a separate question from the empirical continuation of the physical branch B⁺(d) in d towards Φ_k near an even point.

- **d = 2 — exact coincidence**, not a limit. At d = 2, k = 1 the numerical check of identity (2) (App. A.5.1) simultaneously confirms that B⁺ coincides with Φ₁ (§2.4).
- **d = 4 and d = 6 — the continuation of the branch B⁺ coincides with Φ₁ (k = 1, dk = d) symmetrically.** Cosine of eigenfunctions = 1 − O(10⁻⁸).

**Remark (NOT RIGOROUS, mpmath dps 40).** The cosine of the Φ₁ and B⁺ coefficients in the monomial basis (§7.4, item 5) is a basis-dependent quantity, not evidence of asymmetry d = 4 / d = 6.

**Attribution to the classics.** Eigenvalues of coordinate changes are a classical fact; the explicit form and the generalisation of the classical coordinate mode to the class |x|^d via the construction χ_ε = id + ε·(x|x|^{dk} − x) (App. A.6.1) is carried out in this work.

### 2.6 Class constants and construction of catalogues

_Remark._ Here and below H_d(u) denotes the universal scaling function of the period-doubling universality class of order d; H in §2.1 is the special case d = 2 (H = H_2); for d = 2, 4, 6 the function H_d is realised as a power series in t = |z|^d, d ∈ ℕ.

**Table T2.0 — class constants needed for the rescaling u_m(a_i) = A·(a_i − c_F)·δ^m** (values NOT RIGOROUS, mpmath, except where separately noted):

| class | c_F | A | δ | AU0 | provenance |
|---|---|---|---|---|---|
| d = 2 | −1.40115518909205060052382685026… | −2.33188562955070413230142737400505175507716… | 4.66920160910299067185320382… | 0.41403799908627305889… | external input, `06-EXTERNAL-INPUTS/study1_d2_frozen_experiment` (README-EXTERNAL); dps 150 |
| d = 4 | 1.594901356228820564497830482717… | 1.49043254635406998722066572645949… | 7.28468621707334336430893056799… | (calibration constant of the script; see below) | Note M6, sec3_pilot_d4; dps 60 |
| d = 6 | 1.6832601982046430005119730057818… | 1.32131970544285165944683814471791… | 9.29624683277137008283447656636… | (calibration constant of the script; see below) | Note M6, sec5_subfamily_d6; dps 60 |

A is a calibration constant, chosen so that u = −1 corresponds to the calibration cascade. For d = 4, 6: **AU0_d** := −L_s, where L_s is the limit of s*_m·δ_d^m as m → ∞.

_Remark._ For d = 2 the normalisation differs from AU0_d for d = 4, 6: AU0(2) := 0.41403… ≠ L(2) := 0.428… — two different normalisations (AU0_d := −L_s for d = 4, 6; for d = 2 AU0 is the calibration constant of the script). See Appendix B.1.0, remark on naming. For d = 2 the normalising constant AU₀(2) = 0.41403799908627305889… is a calibration constant of the script and does not coincide with L(2) = 0.428837498429404068…; substituting one for the other is a normalisation error (see App. B.1.0).

**"Match" rule:** for a superstable a_i of level m we compute u_m(a_i); the closest zero of the H_d catalogue (Euclidean metric) is chosen; if d_u < 0.01 (threshold R_accept), the match is recorded. The exact-period condition (min_{j≤m}|x_{2^j}| ≥ 10⁻⁴⁵) is a separate mandatory check.

**Obtained catalogues:** for d = 2 — **23 certified zeros of H₂ in the annulus 7 < |u| < 10** (external input, catalogue from Study 1; §2.6, App. B.1.0, B.3.3); the extended uncertified scan of H₂ (B.3.3) yields 222 unique zeros (441 with conjugate pairs). For d = 4 — 62 zeros (catalogue of H₄), for d = 6 — 342 zeros. **Completeness of the catalogues is not claimed** for any class.

External inputs of the d = 2 class (23 zeros of H₂, 23 superstable parameters of the family z² + c, the certificate of λ₂(2), the Z72 localisation) are taken from the Study-1 manuscript and were not independently re-checked in the present work (README-EXTERNAL).

---

## 3. Correspondence results

### 3.1 Statement

For the class d consider the zero set {u : H_d(u) = 0} and the set of scaled superstable parameters of level m, {u_m(a_i)}. The correspondence claims that as m → ∞ these sets converge pointwise:

    d_u(a_i, m) = |u_m(a_i) − (closest zero of H_d)|  →  0.

This is verified empirically: at a fixed m we present a catalogue of zeros of H_d, a catalogue of superstable parameters of the same level and a measure of match d_u, with specificity controls (§3.4).

### 3.2 Results on the three classes

For d = 6: (i) **287 out of 294** at the frozen level m = 9; (ii) **296 out of 296** when unfrozen (m = 10…12, §3.3).

**Table T3.1 — summary table of the correspondence on the three classes** ("61/62" — match of zeros of H₄ to parameters; the criterion differs from "61/61 pass" of §4.2).

| class d | level m (period P = 2^(m+1)) | match (frozen protocol) | continuation when unfrozen | median d_u | controls (matches/total) | specificity |
|---|---|---|---|---|---|---|
| 2 (z² + c) | m = 10 | 22/23 | — (for d = 2 continuation when unfrozen was not performed by design) | 3.16·10⁻⁵ | 0/20 (for d = 2 — 20 control parameters of the same range, not superstable, compared with the extended catalogue of 441 zeros of H₂ with conjugate pairs; §3.4) | median 82339× (per-zero, 22 accepted; range 48886…222247×; external input) |
| 4 (1 − a·x⁴) | m = 9 | 61/62 | — (for d = 4 continuation when unfrozen was not performed by design) | 7.3·10⁻⁵ | 0/6 | median 2164× (n = 8 nonzero; upper median/Q10 = 3487×; range 212…38478×) |
| 6 (1 − a·x⁶) | m = 9 | 287/294 | 296/296 at m = 10…12 (§3.3) | 1.13·10⁻³ | 0/6 | median 1488× (n = 7 nonzero; range 22…3639900×) |

**Non-uniformity of metrics between classes.** Comparing "82339× (d = 2)" and "1488× (d = 6)" directly is not correct; by T3.1 the medians are — 82339× (d = 2, n = 22), 2164× (d = 4, n = 8), 1488× (d = 6, n = 7), with a lower tail down to 212× (d = 4) and 22× (d = 6, n = 7).

_Footnote to T3.1 (definition of "n")._ In the "specificity" column, n is the number of catalogue zeros for which at least one competing candidate root was found nearby (per-rival; for d = 4: n = 8 of 61 accepted, for d = 6: n = 7 of 287 accepted); for d = 2 the specificity is per-zero (by neighbours on the common scan grid). In the "controls" column, n is the number of tested control parameter candidates of the same range, fixed by the protocol (n = 6 for d = 4, 6; n = 20 for d = 2) and unrelated to the rival count. The denominators in the two columns refer to different sets and are not directly comparable.

_Remark._ For d = 2 the numbers are from the external input (folder `06-EXTERNAL-INPUTS`, family z² + c; the link with 1 − a·|x|² is a conformal reparametrisation, see Dudko–Lyubich (2023)).

The three classes give one qualitative result: the overwhelming majority of parameters have a partner with vanishingly small residual; the controls — none. The cases of incomplete match (1/23, 1/62, 7/294 interior) are not counter-examples; for d = 6 an analysis (§3.3) reveals the mechanism of delay.

**OOS caveat d = 4.** Out-of-sample verification of the remainder "within a factor 2" (§7.3) for d = 4 at n = 61 **does not show superiority of the two-term model over baseline**: median 65.1 %, spread 9.3–97.7 %; baseline — 76.7 %. This is a consequence of the small sample size, not a failure of the model; for d = 6 (n = 287) the confirmation was obtained (92.0 % vs 90.5 %).

**Rigor status.** The d_u layer in all three classes — NOT RIGOROUS: the rescaling constants (c_F, A, δ) are uncertified float/mpmath values (§2.6). The rigorous part — existence of the parameter as a root in an interval box (§4).

### 3.3 Anomalous d = 6 subfamily and the 1/16-step lattice

**Mechanism of delay of the nine "holdouts".** At the frozen level m = 9, seven superstable parameters of the d = 6 class did not find a partner among the zeros of H₆ (interior failures), plus two more at the edge. All nine are genuine zeros of the correspondence; their non-match at m = 9 is explained by **finite-m delay**: the seed had not yet fallen into the basin of attraction of the true partner. When continuing the same direct method (Newton, 1–3 levels of m deeper) all nine converge to a partner. The spectral cause is the notable closeness of the moduli |λ₃| and |λ₂| at d = 6 (|λ₃/λ₂| = 0.922, §2.3), making the pair prone to interference: one compares the moduli of the eigenvalues λ₂ (B⁺ = +0.371) and λ₃ (B⁻ = −0.342) on the same branch of the spectrum at d = 6. All nine parameters are certified interval-wise by Krawczyk (§4.3).

**Table T3.2 — late convergence of the 9 "holdouts".**

| identifiers | first clean direct hit (d_u < 0.01) | d_u at that level |
|---|---|---|
| S031 | m = 10 | 8.0·10⁻⁴ |
| S070 | m = 10 | 6.6·10⁻⁴ |
| S235 | m = 10 | 1.7·10⁻³ |
| S279 | m = 10 | 7.1·10⁻⁴ |
| S319 | m = 10 | 1.3·10⁻³ |
| S104, S187, S339 | m = 11 | 1.4·10⁻⁴ … 4.9·10⁻⁴ |
| S036 | m = 12 | 2.2·10⁻⁴ (exact period certified, 7.9·10⁻⁶⁸) |

**The 1/16-step lattice is a scanner-selection artefact (d = 4 and d = 6).** The catalogue positions of the nine holdouts lie close to the nodes of the 1/16 lattice; on the whole catalogue of d = 4, d = 6 — per-coordinate 36.1–52.9 %, 2D joint 22.9–39.2 % against a background of 1.56 % (p ≪ 10⁻¹⁵); in ORIG-ONLY — up to 69.0 % joint. The decisive test is a rescan with step 1/24; forensics — Supplement §S5. For d = 2 the lattice is weaker (joint ≈ 12.9 % on 441 zeros with conjugate pairs; 13.51 % on 222 zeros without conjugates, Appendix B.6.1); the known23 control does NOT show the lattice (p = 0.30).

### 3.4 Controls, specificity and incompleteness of catalogues

In addition to comparing superstable parameters with zeros of H, for each class a control protocol was run — the same comparison for parameters of the same range that are not superstable. Result — 0/20 (d = 2), 0/6 (d = 4), 0/6 (d = 6).

---

## 4. Rigorous certificates

### 4.1 Method: Krawczyk interval box (common to d = 2, 4, 6)

The rigorous part answers a narrower question than the correspondence as a whole: "does there exist, and is it unique, a superstable parameter of exactly the stated period in a neighbourhood of the approximation" — solved by interval arithmetic without dependence on c_F, A, δ.

**Setup.** Complex a-plane; for a parameter of level m, period P = 2^(m+1), with seed c₀ build a square box X = c₀ + [±r]². F(a) := x_P(a), F′(a) := dx_P/da — both recursively along the orbit. Arithmetic — Arb balls via python-flint.

On the box X compute F₀ = F(c₀), Df₀ = F′(c₀) and their interval versions D = F′(X); introduce the interval Newton operator N(X) := c₀ − F₀/D and the Krawczyk operator K(X) := c₀ − Y·F₀ + (1 − Y·D)·(X − c₀), Y := 1/Df₀; both act pointwise on every point of the box X and yield an interval image N(X), K(X); lower bounds |x_{2^j}(X)| for j = 0, …, m.

**PASS criterion** (all comparisons are strict comparisons of balls): (i) 0 ∉ D; (ii) |N(X) − c₀|_upper < r; (iii) |K(X) − c₀|_upper < r; (iv) inf_X |x_{2^j}| > 0 for all j = 0, …, m.

**Theorem.** Let X = c₀ + [±r]² ⊂ ℂ and let (i)–(iv) hold. There exists a unique root a*_i of F(a) = 0 in X, the critical orbit does not return to 0 at steps 2, 4, …, P/2 — the exact period of a*_i equals P. The statement is standard: the Krawczyk operator and the interval Newton method — see Krawczyk (1969), Neumaier (1990); it is used here without proof.

The theorem about each box **does not depend** on the uncertified constants of the pilot (c_F, A, δ) — the seeds are exact decimal strings, the theorem is formulated in terms of the orbit of f_a. The correspondence of boxes to zeros of H is a separate NOT RIGOROUS layer on top of the rigorous result.

### 4.2 Degree 4: 61 of 61 certified

**Table T4.1 — certification summary for d = 4.**

| metric | min | median | max |
|---|---|---|---|
| box radius r (nominal r = 1·10⁻²⁶; exact radius — in B.2.1) | 1.0×10⁻²⁶ | 1.0×10⁻²⁶ | 1.0×10⁻²⁶ |
| κ (Krawczyk contraction, \|K−c₀\|/r) | 1.2×10⁻⁵ | 2.9×10⁻⁴ | 5.0×10⁻⁴ |
| inf \|x_{2^j}\| on the box, j = 1…9 | 0.00734 (Q01, j = 9, period 1024) | — | 0.0214 |
| d_u (float, NOT RIGOROUS layer, 8 rivals) | 5.8×10⁻⁷ | 7.3×10⁻⁵ | 3.5×10⁻³ |

All 61 certificates passed (the criterion differs from "61/62" in §3.2); the worst contraction κ = 5.0×10⁻⁴ leaves a margin of 2000×; the worst isolation — Q51/Q52, ≈ 1.57×10⁻⁹ in max-norm (|ΔRe|, |ΔIm|) (Euclidean ≈ 2.07×10⁻⁹), margin ≈ 7.8×10¹⁶ relative to the sum of radii; all 1830 pairs do not intersect. An independent second code path confirmed the result; a battery of 14 tamper types — all caught (Supplement §S4).

**Limitations.** Completeness of the catalogue of zeros of H₄ is not claimed (Q57 was rejected by the pilot). The "box ↔ zero of H₄" correspondence uses uncertified float constants — a NOT RIGOROUS layer.

### 4.3 Degree 6: 287 + 9 = 296 boxes certified

Direct continuation of the protocol of §4.2, adapted to F(a) = x_1024(a). The verdict set is 287 "interior" superstable parameters of the blind protocol (|u| ≤ 44; the catalogue as a whole contains 342 = 294 interior plus 48 edge entries, the edges intentionally not covered by the certificate).

**Theorem (for each of the 287).** There exists a unique root of x_1024(a) = 0 in the disk B(c_i, r_i), r_i ∈ {1×10⁻²², 1×10⁻²³, 1×10⁻²⁴} for the main array ("layer 1", 285 parameters) or {1×10⁻³⁴} for the two special cases ("layer 2"); the period is exactly 1024; all 287 boxes are pairwise disjoint. Isolation of boxes (certificate, pairwise_isolation field, distance of centres in max-norm (|ΔRe|, |ΔIm|) relative to the sum of radii): the worst pair by margin is S172/S180, distance ≈ 3.7×10⁻¹⁰ with r_i + r_j = 1.1×10⁻²² (margin ≈ 3.3×10¹²); the smallest distance between centres — S143/S145, ≈ 1.8×10⁻¹⁰ (Euclidean ≈ 2.6×10⁻¹⁰); all 41 041 pairs of boxes are disjoint.

**Layer 2 — S078 and S090.** A box around the stored seed for S078, S090 — points with extreme |F′(a*)| ≈ 8.1×10¹⁴ and 1.1×10¹⁸ — does not certify already at radius 10⁻²⁵ (S078) / 10⁻²⁹ (S090). Solution: eight Newton iterations at dps = 120, truncation of the centre to 60 significant digits; |c_refined − seed|_upper = 1.74×10⁻²⁴ (S078) / 1.85×10⁻²⁴ (S090), both < 5×10⁻²⁴. Reducing the radius to r = 10⁻³⁴ removes the fragility. → §S6 for details on pow/mul fragility.

The rigor of the certificate does not depend on the algebraic order of exponentiation: the interval estimate under any order is a legal enclosure; the order only affects the width of the enclosure and hence the passability of the criterion, not its rigor. All 357 boxes are rigorous without caveats; the pow/mul cross-check is an additional stability check, not a condition of rigor.

**Table T4.1-b — certification summary for d = 6.**

| metric | min | median | max |
|---|---|---|---|
| box radius r | 1.0×10⁻³⁴ (S078, S090, layer 2) | 1.0×10⁻²³ | 1.0×10⁻²² (12 parameters) |
| κ (Krawczyk contraction) | 4.2×10⁻⁷ | 1.5×10⁻¹ | 5.0×10⁻¹ (S196; boundary case: the rule κ ≤ 0.5 admits exactly the boundary value) |
| inf \|x_{2^j}\| on the box, j = 1…9 | 0.0182 (S203) | 0.049 | 0.0686 |
| d_u (float, NOT RIGOROUS layer, 287 interior_accepted) | 1.4248×10⁻¹⁰ | 1.1321×10⁻³ | 7.9790×10⁻³ |

Distribution of layer-1 radii: 98 at r = 1×10⁻²⁴, 175 at r = 1×10⁻²³, 12 at r = 1×10⁻²². The worst pairwise isolation by margin — S172/S180 (distance of centres ≈ 3.7×10⁻¹⁰ with sum of radii 1.1×10⁻²²; margin ≈ 3.3×10¹²). An independent second code path confirmed the result; a battery of 16 tamper types — all caught (Supplement §S4).

_Remark._ The worst indicators (κ = 5.0×10⁻¹ for S196, isolation S172/S180) never coincide on the same box: S196 has good isolation, S172/S180 — a typical Krawczyk contraction.

**Summary.** In total the rigorous boxes of the work are **357** (61 at d = 4, 296 = 287 + 9 at d = 6): 348 — layer m = 9, period 1024; nine "late" at d = 6 (§3.3) — periods 2048/4096/8192, radius 1×10⁻²⁶.

**Important distinction.** "296 out of 296" of the d_u layer (§3.3) is a NOT RIGOROUS correspondence of zeros of H to parameters; "296 rigorous boxes" of this section is the existence of the parameters themselves as roots of exact period. The coincidence of the numbers is a consequence of the fact that all 9 holdouts of level m = 9 find partners at m = 10–12 (finite-m delay), and not independent evidence for the correspondence. The coincidence of the two "296" stems from the fact that every parameter presented has been proved. Completeness of the catalogue of H₆ is not claimed; of the 48 edge entries (|u| > 44), only S319 and S339 (m = 10/11) are certified, the other 46 are not; in total 287 (m = 9, interior) + 7 (late, interior) + 2 (late, edge) = 296.

### 4.4 Numerically verified λ₂(2), λ₂(4), λ₂(6) — "Bridge" and "Bridge-2"

**λ₂(2).** An interval certificate of an eigenpair of the linearisation of the doubling operator for the discretised (N = 48) approximation of the true operator. Result: λ₂(2) = −0.1236527125526870243353936525855 ± 1.82×10⁻³³ — **numerically verified** (§2.3).

**λ₂(4), λ₂(6) — "Bridge-2" (Chebyshev basis).** The first attempt (monomial basis) gave rigorous certificates for N = 16, 24, but with a gap of 18–30 orders of magnitude between the certified centre and the float anchor of T2.1. Switching to the Chebyshev basis reduces the growth of the condition number to polynomial (roughly linear growth of log(cond) with N: ~7·10⁵ at N = 16 to ~1.5·10⁸ at N = 48, factor ≈ 21 over N = 16 → 48); K(X) ⊂ int(X) passes on the ladder N = 16, 24, 32, 40 (d = 4) and up to 48 (d = 6), ball radii of λ₂ ~10⁻⁴⁰…10⁻⁴¹. Detailed walk-through — Supplement §S7.

Certified centres: λ₂(4) = 0.291838408983545104530996961… (N = 40), λ₂(6) = 0.371071467681930132470983381… (N = 48). Both float anchors of T2.1 are reproduced in all 20 significant digits.

---

## 5. Rate of convergence: the law min(δ, 1/|λ₂|)

**Framing caveat.** The material of this section is an independent re-discovery, modern certification and extension of the well-known corrections-to-scaling structure, not a new result of this work.

### 5.1 History: the law has been known since 1987–1994

The first work to fix the corrections-to-scaling structure was Mao, Hu (1987, *J. Stat. Phys.* 46, 111). Reick (1992, *Phys. Rev. A* 45, 777, "Universal corrections… Multiple scaling and crossover") developed this line. The most complete exposition is Briggs (1994, arXiv:chao-dyn/9407019): his tables reproduce the law min(δ, 1/|λ₂|) exactly. The external anchor at d = 2 is Christiansen, Cvitanović, Rugh (1990).

**Table T5.1 — comparison with Briggs (1994).**

| d | leading exponent (Briggs) | second | third | fourth | our match |
|---|---|---|---|---|---|
| 2 | 4.6692 = δ | −8.0872 = 1/λ₂ | −17.4499 = 1/λ₃ | — | all listed signs agree |
| 4 | **3.4266 = 1/λ₂** | −3.9114 (= 1/λ₃, our λ₃ = −0.2557) | −6.4283 | δ = 7.2846 (only fourth) | 1/λ₂ and the sign of 1/λ₃ agree |
| 6 | **2.6949 = 1/λ₂** | −2.9242 (= 1/λ₃) | −4.7499 (our λ₄ = −0.2105) | 4.8100 (our 1/λ₅ = 4.8112) | three exponents match exactly; 1/λ₅ differs from Briggs by 2.5·10⁻⁴ relative (4th significant digit) |

Conclusion of the section: a rate of convergence equal to δ is correct ONLY for d = 2; at d = 4, 6 δ is only the fourth largest exponent in Briggs, and the much slower 1/λ₂ dominates.

**Status of the law min(δ, 1/|λ₂|):** conditionally RIGOROUS as a fact of the literature (Mao–Hu 1987; Reick 1992; Briggs 1994) — with internal numerical re-verification (NOT RIGOROUS).

### 5.2 Evidence (why a modern verification of the re-discovery is needed)

The residual model of the correspondence is two-term: e(m) ~ −(C·δ^(−m) + D·λ₂^m); the observed behaviour is that the slower of the two terms dominates. Detailed LOOCV/OOS procedures — Supplement §S8.

- **Direct detection of the λ₂-mode under the δ-layer at d = 2.** After subtracting the δ-layer the remainder is sign-alternating (5 sign alternations out of 6; ratios 7.88–7.89 against 1/|λ₂(2)| = 8.087). LOOCV shows an improvement of the two-term model over the one-term one by a factor of 58 (calibration) and 300–530 (real zeros of H).
- **Blind confirmation at d = 6.** The computation was carried out under a frozen procedure with 8 candidates. The protocol itself chose 1/λ₂(6): geometric-mean deviation −0.2 %, closest rival +8.7 %.
- **Interference at d = 6 explained.** The spread of ratios at d = 6 (0.6…7.5) is a consequence of the nearly degenerate pair λ₂ = +0.371, λ₃ = −0.342 (|λ₃/λ₂| = 0.922).

**Table T5.2 — prediction vs observation.**

| class | δ | status | 1/λ₂ | status of λ₂ | slow term | prediction | observation |
|---|---|---|---|---|---|---|---|
| 2 | 4.669 | conditionally RIGOROUS | 8.087 (λ₂ < 0) | numerically verified | δ-term | 4.669 | 4.666–4.671 (by u = −1); 4.31–4.91 (by all real zeros) |
| 4 | 7.285 | conditionally RIGOROUS | 3.4266 | numerically verified | λ₂-term | 3.4266 | 3.36–3.46 after ε-correction of a_F (T5.3) |
| 6 | 9.296 | NOT RIGOROUS | 2.6949 | numerically verified | λ₂-term | 2.6949 | medians of d_u ratios on a pair of levels (m = 12 and m = 13): 2.556 and 2.830 on different subsamples (see Supplement §S8), geometric mean 2.6893 (−0.2 %) |

Remark on accuracy: ball radii of λ₂(4), λ₂(6) ~10⁻⁴⁰.

**Table T5.3 — confirmation of d = 4 after correcting the reference constant.**

| stage | observed ratio | remark |
|---|---|---|
| Raw (before a_F correction) | ρ_m = 3.58492658… (m = 12, k = m + 1 = 13; from ratio_table["12"]) | discrepancy with prediction 3.4266 is 4.6 % |
| After ε-correction of a_F (ε ≈ 8.6·10⁻²² vs the stated 4·10⁻²³) | 3.36–3.46 (bracket of correction methods) | the boundaries of the bracket are the spread between correction methods, not measurement noise |
| Prediction of the law (1/λ₂) | **3.4266** (λ₂(4) = 0.29183840898354510453, N = 40) | matches the bracket (3.36–3.46 ≈ 3.4266 ± 1.9 %) |

m = 14…17 carries an artefact of the Richardson recursion (a_F was built from the same levels) and is not included in the fit. The number 3.585 is not the median over 13 zeros, but the ratio ρ_m := e_{m+1}/e_m (ratio of adjacent correspondence residuals at levels m and m+1, not the Feigenbaum constant δ) at m = 12.

The discrepancy between the claimed accuracy of a_F(4) (4·10⁻²³) and the actual ε-correction (≈8.6·10⁻²²) — about a factor of 20 — is not explained in the present work, nor is what causes it. **After the ε-correction the law gives no point prediction for d = 4; the bracket 3.36–3.46 is a consequence of the uncertainty of the correction method**, verified here at the level of order and sign. The same uncertainty in a_F(4) likely also explains the weak result of the out-of-sample "within ×2" verification for d = 4 (§7.3): without an accurate a_F(4) the λ₂-term of the two-term model is calibrated for d = 4 more roughly than for d = 6.

**External estimate of δ₄.** Burbanks, Osbaldestin, Thurlby (2021) give δ₄ ∈ [7.28468621706, 7.28468621709] (conditionally RIGOROUS anchor).

### 5.3 Role of λ₂(2) in the law

The numerically verified value of λ₂(2) (±1.8·10⁻³³, §4.4) is the most accurate among the three numerically verified values of λ₂ of this work; a rigorous estimate of N → ∞ remains an open problem (§2.3). It serves as a self-consistent anchor for verifying the chain of computations leading to the d = 2 row of Table T5.2.

### 5.4 Artefactual traps

Circularity of the cascade tail (a_F is built by Richardson from the same levels, using it as a reference for the fit on these levels) and the growing parasitic term from the ε error of the reference constant are the two traps of the pipeline behind the numbers of T5.3; the mechanism — Supplement §S1, items 2–3. The accuracy of reference constants must be checked with a "bracket" of independent estimates before any substantive conclusion from the shape of the residual curve, not after the fact.

---

## 6. Spectrum over d: three zones

**Framing caveat.** All material of this section, including T6.1, is NOT RIGOROUS (mpmath/float, without interval certificates); the only interval anchor among the constants of T2.1 is δ(4) (external, Burbanks–Osbaldestin–Thurlby 2021); λ₂(2/4/6) — certificates on a discretisation N; everything else — NOT RIGOROUS. The very course of the curve between anchors and beyond them is uncertified.

### 6.1 δ-zone, d ≤ d* ≈ 2.36

At the nodes d = 1.3, 1.5, 1.7 the observed convergence ratio agrees with δ(d) to 3–4 significant digits (δ ≈ 3.3787, 3.8006, 4.1720); these δ(d) are **NOT RIGOROUS** (exploratory sweep). The δ-zone extends to d* ≈ 2.36.

### 6.2 Transition window 2.4–2.9: crossover d*, sign jump of λ₂ at d = 2

A special operator grid of 45 points over d ∈ [1.9, 7] (double truncation N = 40/48, residuals ≤ 10⁻²⁷). **There is no smooth zero of λ₂(d):** both subleading branches B⁺, B⁻ keep their sign throughout [1.9, 7]. The sign change of the classical λ₂ is a jump exactly at d = 2 (§2.4). Details — Supplement §S9.

**Crossover d*.** The point where the law changes its dominant term: δ(d)·λ₂(d) = 1 on the B⁺ branch; bisection gives **d* = 2.3616 ± 0.0008** (δ = 5.2155, λ₂ = 0.19175). Both components are NOT RIGOROUS.

**Dynamical confirmation.** Direct numerical continuation of the cascade up to m = 16 at d = 2.75 and d = 3.0: rates 4.09–4.40 and 4.14–4.49 agree with the prediction of the B⁺ branch (4.53 and 4.22; error 2–10 %).

### 6.3 λ₂-zone, d ≥ 3: monotonic growth, saturation

**Table T6.1 — λ₂(d), surface level (exploratory sweeps).**

| d | 5 | 7 | 9 | 10 | 13 | 14 | 15 | 16 | 17 | 19 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| λ₂(d), implied | 0.330 | 0.388 | 0.425 | 0.439 | 0.478 | 0.501 | 0.508 | 0.513 | 0.518 | 0.525 | 0.528 |

The node d = 8 ended with a bootstrap-method error (Supplement §S3). Linear interpolation between d = 7 and d = 9 gives λ₂(8) ≈ 0.407 — **not computed, a guide, neither a certificate nor a measured value**.

Outside the transition zone (d ≥ 3.5) the operator curve B⁺ and the implied λ₂ agree to 0.7–3.5 %; δ(d) agrees with the independent cascade route to 10⁻⁸–10⁻¹⁴. The range d ∈ [3.5, 5] in this paragraph is confirmed by the data of the operator grid of Appendix B.5 (mpmath, NOT RIGOROUS); T6.1 (exploratory sweep) covers d = 5 and above.

**Conclusion of the section (saturation in the range 0.527–0.535 is observed by d = 20; the hypothesis of a limit 0.53–0.55 as d → ∞ is extrapolation, not a proved limit).** λ₂(d) grows monotonically with d and saturates to 0.53–0.55 as d → ∞ — a refinement of the Briggs (1994) estimate of |λ₂| → ~1/2. δ(d) ranges from 3.38 to 17.0.

**Sweep vs operator grid in the transition zone — the actual table (replacing the qualitative "sweep underestimates/overestimates by 12–22 %"):**

| d | sweep (implied λ₂) | operator value (grid of §6.2, B⁺ branch) | discrepancy |
|---|---|---|---|
| 2.8 | 0.19367 | 0.22409 | −13.6 % |
| 2.9 | 0.19024 | 0.23073 | −17.6 % |
| 3.0 | no node in sweep | 0.23714611165 | — |

_Remark._ In the transition zone 2.4–2.9 the sweep systematically underestimates λ₂ relative to the operator value (negative discrepancy); outside the transition zone (d ≥ 3.5) the operator curve B⁺ and the sweep agree to 0.7–3.5 %.

---

## 7. Discussion

### 7.1 Relation with corrections-to-scaling theory and priority

The material of §5 fits into the perturbation theory of the fixed point of the doubling operator: δ determines the orbital convergence of the cascade always, and the subleading spectrum corresponds to decaying directions, whose contribution to the parametric convergence shows up as a correction to the δ-layer — a perturbative derivation in §3 of Briggs (1994). The present work has no formal proof of the min-rule as a theorem.

**Novelty of the correspondence principle.** Verdict: **probably new**, with an incomplete priority search (kneading/MSS literature and Perrier–Girault 2022 were not compared in detail) — the construction "zero set of H_d(u) corresponds pointwise to superstable parameters via u = A·(a − c_F)·δ^m" was not found in any of the works we checked (targeted search across Isaeva–Kuznetsov (2005); Levin–Świątek (2018/2021); Gorbovickis–Yampolsky (2017/2018); van der Weele–Capel–Kluiving (1987); Perrier–Girault (2022, scaling and fine structure of superstable orbits of the logistic map) — the closest thematically among the works we know; a detailed comparison with the H_d construction in the present work was not carried out; see also Varin (2012) — mentioned only as a caveat); the literature search is limited to the listed works, a systematic review of the combinatorial (kneading/MSS) literature was not carried out; a recent review of renormalisation theory and the MLC conjecture (Dudko 2025, arXiv:2512.24171 — preprint submitted 30.12.2025, after the submission of the manuscript; "Proceedings of the ICM 2026" status not independently confirmed, mentioned only as a negative priority control) does not contain such a construction. Briggs (1994), Mao–Hu (1987), Reick (1992) establish the scaling rate via the spectrum but not the pointwise correspondence of zeros of H to parameters.

_Remark._ New in this work — the recipe for constructing H_d and the catalogue of zeros; everything else (the law rate = min(δ, 1/|λ₂|)) is a re-discovery/extension of known structure.

### 7.2 The 1/16-step lattice in the u-plane

Full statistics, the decisive control and the established mechanism — §3.3: the attraction of zeros of H to the nodes of the 1/16-step lattice in the u-plane — **established for d = 4; for d = 6 — by the same protocol, but the statistics differ (joint 39.2 % vs 22.9 % at d = 4; see S5)** — an identified and explained scanner-selection artefact, proved for both classes (d = 4, d = 6) by a single protocol (rescan with an incommensurate shifted grid); d = 2 is weaker (joint 13.51 % on 222 zeros — Appendix B.6.1). Here — a summary; details — Supplement §S5.

The "zero ↔ parameter" correspondence is not undermined: the median d_u (1.13·10⁻³ for d = 6) is 7× smaller than the median distance to the nearest node.

### 7.3 Status of the correspondence principle

**Perturbative derivation.** Within the two-term model of §5 the scaled position has the form u = A·(a − a_F)·δ^m + R_m, R_m scales as **λ₂^m**.

**Out-of-sample check: "within a factor 2".** The model e_9 = c₀ + c₁·u + c₂·u²·δ^(−9) is trained on 30 % of the catalogue and tested on 70 %. The λ₂-term is not explicitly included in the OOS model; the test checks the δ⁻ᵐ correction.

- **d = 6** (n = 287): "within ×2" — 92.0 % (median over 500 splits); baseline — 90.5 %. An advantage of ~1.5 pp is given without an estimate of statistical significance — **the model is confirmed as working**.
- **d = 4** (n = 61): "within ×2" — median 65.1 %, spread 9.3–97.7 %; baseline — 76.7 %. At this size **the superiority of the model over a constant is not shown**.

For d = 4 we give one median without a battery of splits (as opposed to the 500 splits for d = 6); the comparison between classes in this sense is methodologically asymmetric, and the sample of d = 4 (n = 61) itself is small.

**Proposition (on the linearised catalogue — exactly one zero; a trivial special case).** If the operator T is linearised in a neighbourhood of g_d (family g_d + s·h_d by its linear approximation in s), the resulting family has exactly one zero. *Proof.* The numerator is affine in u (A·u + B, A ≠ 0); an affine function of one complex variable with nonzero leading coefficient has exactly one root. ∎ The pointwise but non-unique structure of zeros of H_d is a manifestation of the nonlinear geometry of W^u. _Remark._ The result is trivial and is not used for the nonlinear case; it is given to explain how the scan-and-polish works (after linearisation Newton's iteration for finding the closest zero converges — which explains why Newton works at the d_u level, §3.2).

**Summary.** The principle is confirmed as universality plus an explicit λ₂^m remainder; the content is in the rate of convergence (min(δ, 1/|λ₂|)) and in the structure of the residual (verified out-of-sample at d = 6: an excess of 1.5 pp over the baseline; significance was not estimated).

**Genre and rigor of the work.** By genre the work is experimental mathematics. The only formal results are the theorem on Krawczyk boxes (§4.1), the Proposition on the gauge tower (§2.5) and the Proposition on the linearised catalogue (§7.3); the first two are rigorous under the assumption of the existence of g_d; for even d ≥ 4 this assumption is covered by Epstein (1986), Lyubich (1999) and Avila–Lyubich (2011); a rigorous estimate of the N → ∞ truncation error was not carried out. There is no formal theorem about H_d; the "zero ↔ parameter" correspondence remains an empirical NOT RIGOROUS layer.

### 7.4 Limitations

1. **Completeness of the catalogues of zeros of H is not claimed** for any of the three classes; rescan with incommensurate step 1/24: d = 4 — 146 zeros at level 8 vs 62 (≈2.35×); d = 6 — 641 vs 342 (≈1.87×), §3.2, §3.3. The influence of catalogue incompleteness on the specificity numbers (§3.2, §3.4) was not estimated.
2. The "zero ↔ parameter" correspondence (d_u layer, §3) in all three classes — NOT RIGOROUS: uncertified float constants (c_F, A, δ). The rigorous part — existence and uniqueness of the parameter as a root in a box (§4).
3. The subleading (λ₂-) term of the two-term model for d = 4 is not measured quantitatively with sufficient accuracy: |C| ≲ 1.6 with D ≈ 0.038.
4. Transition window d ≈ 2.4–2.9 (§6.2): crossover d* = 2.3616 ± 0.0008 found; the operator grid is NOT RIGOROUS, certified only at the anchors of T2.1.
5. **Continuation of the branch B⁺ towards Φ₁** near even d is established only numerically (cos = 1 − O(10⁻⁸); mpmath, NOT RIGOROUS, §2.5); using the monomial basis cos ≈ 0.21–0.99999982 (range over Table A.4.2) — a basis-dependent quantity, not evidence of asymmetry d = 4 / d = 6 (§2.5, S1).
6. The OOS "within ×2" model (§7.3) is convincing for d = 6, but not for d = 4 at n = 61.
7. A rigorous estimate of the N → ∞ truncation error for λ₂(4), λ₂(6) was not obtained (only a heuristic ladder N = 40/48, §4.4).
8. The mechanism of the 1/16-step lattice is not established even when the artefact status has been proved for d = 4.
9. **The discrepancy between the "claimed accuracy" of a_F(4) ≈ 4·10⁻²³ and the actual ε-correction ≈ 8.6·10⁻²² (§5.2) — a factor ≈ 20.** In the present work this discrepancy is not traced to a source; as a consequence, for d = 4 the two-term model does not give a point prediction; the bracket 3.36–3.46 around 3.4266 reflects the uncertainty of the correction method, not data noise. The same uncertainty in a_F(4) likely also explains the weak result of the OOS "within ×2" check for d = 4 (item 6 below, §7.3).

---

## 8. Conclusion

The paper gives an explicit computable realisation of the "parameter space ↔ unstable manifold W^u(g_d)" correspondence of period doubling: the universal scaling function H_d(u), whose zero set serves as a catalogue of superstable parameters along W^u in the linearising coordinate. The correspondence has been verified on three classes f_a(x) = 1 − a|x|^d — d = 2 (from external inputs), d = 4 and d = 6 — with specificity of match against controls: medians 82339× (d = 2, n = 22), 2164× (d = 4, n = 8), 1488× (d = 6, n = 7), with a lower tail down to 212× (d = 4) and 22× (d = 6, n = 7) (T3.1); for d = 6 the anomalous nine non-matches at the frozen level have received an explained mechanism (§3.3).

The rigorous part of the work: interval Krawczyk certificates for the existence, uniqueness and exact period of 357 superstable parameters (61 at d = 4, 296 at d = 6); numerically verified λ₂(2) ±1.8×10⁻³³ and λ₂(4), λ₂(6) at N = 40/48 ("Bridge-2", §4.4). A by-product is the law rate = min(δ, 1/|λ₂|), re-discovered, certified and extended to d up to 20, with a quantitative blind confirmation at d = 6 and confirmation at the level of order and sign at d = 4 (no point prediction because of the discrepancy of a_F(4)), with a resolution of the behaviour of λ₂(d) near d = 2 (sign jump; crossover d* = 2.3616 ± 0.0008, §6). The mechanism of the jump is the Proposition on the gauge tower (§2.5) — a formal identity for the fixed point g_d, whose existence and uniqueness at even d are established by Epstein (1986), Lyubich (1999) and Avila–Lyubich (2011); the identity is proved for all d > 1 on the odd real branch when dk ∈ 2ℤ ∩ ℝ₊.

The attraction of zeros of H to the 1/16-step lattice is established as a scanner-selection artefact for d = 4 and d = 6 (§3.3, §7.2). Further work: a theoretical proof of the law min(δ, 1/|λ₂|); a rigorous estimate of the N → ∞ truncation error for λ₂(4), λ₂(6); an analytic description of the behaviour of B⁺(d) near even d; the mechanism of the attraction of zeros of H to the 1/16 lattice.

---

## Dependency map: claim → status → where verified

A summary table for the reader: for each key claim of the paper — status on the four-level scale of §2.3 and the section where it is verified.

| claim | status | where |
|---|---|---|
| δ(2), λ₃(2) | conditionally RIGOROUS (lineage of Lanford 1982, CCR 1990) | §2.3 |
| λ₂(2) = −0.123652712552687024335 | numerically verified, ±1.8×10⁻³³, N = 48 | §2.3, §4.4 |
| δ(4) ∈ [7.28468621706, 7.28468621709] | conditionally RIGOROUS (external anchor Burbanks et al. 2021) | §2.3 |
| λ₂(4), λ₂(6) | numerically verified, N = 40 / N = 48 | §2.3, §4.4 |
| δ(6), λ₃(4), λ₃(6) and the rest of the subleading spectrum of d = 6 | NOT RIGOROUS (validated by agreement with Briggs 1994) | §2.3, §5.1 |
| Correspondence d = 2: 22/23, median d_u = 3.16×10⁻⁵ | d_u layer NOT RIGOROUS; 23 boxes RIGOROUS (external input) | §3.2 |
| Correspondence d = 4: 61/62, median d_u = 7.3×10⁻⁵ | d_u layer NOT RIGOROUS; 61/61 boxes RIGOROUS | §3.2, §4.2 |
| Correspondence d = 6: 287/294 (m = 9) and 296/296 (m = 10…12) | d_u layer NOT RIGOROUS; 296 boxes RIGOROUS | §3.2, §3.3, §4.3 |
| Mechanism of the anomalous d = 6 subfamily (finite-m delay) | NOT RIGOROUS (qualitative explanation; the orbit of S036 has been checked exactly, 7.9×10⁻⁶⁸) | §3.3 |
| Construction of the catalogues of H_d, completeness not claimed | NOT RIGOROUS | §2.6, §3.4 |
| Specificity 1488× / 2164× / 82339× across the classes | NOT RIGOROUS; metrics are non-uniform between classes | §3.2, §3.4 |
| 1/16-step lattice in the u-plane | scanner-selection artefact (proved d = 4; d = 6 — by the same protocol) | §3.3, §7.2 |
| Law min(δ, 1/\|λ₂\|) — statement | conditionally RIGOROUS as a fact of the literature (Mao–Hu 1987; Reick 1992; Briggs 1994) — with internal numerical re-verification (NOT RIGOROUS) | §5.1 |
| Law min(δ, 1/\|λ₂\|) — confirmation in this paper's data | NOT RIGOROUS, with internal re-verification (LOOCV 58×/300–530×, blind protocol d = 6) | §5.2 |
| Curve λ₂(d), δ(d), three zones, saturation to ~0.53–0.55 | NOT RIGOROUS | §6 |
| 61 boxes d = 4: existence, uniqueness, period 1024 | RIGOROUS | §4.2 |
| 287 boxes d = 6: existence, uniqueness, period 1024 | RIGOROUS | §4.3 |
| 9 late boxes d = 6: period 2048/4096/8192 | RIGOROUS | §3.3, §4.3 |
| Crossover d* = 2.3616 ± 0.0008; sign of λ₂ jumps at d = 2 | NOT RIGOROUS (three independent trials consistent) | §6.2 |
| Simplicity of the eigenvalue α^(−2) at d = 2 | NOT RIGOROUS (two independent implementations) | §2.4 |
| B⁺/B⁻ branches, pointwise calibrability at even d | NOT RIGOROUS for the position of the branches | §2.4 |
| Cosine of Φ₁ and B⁺ in the monomial basis (0.21–0.99999982, A.4.2) | basis-dependent quantity, not evidence of asymmetry d = 4 / d = 6; in Chebyshev basis 1 − O(10⁻⁸) | §2.4, A.4.2 |
| λ₂-mode under the δ-layer at d = 2 (direct detection, §5.2) | NOT RIGOROUS — an observation, not a verified result; the specific levels m have not been re-confirmed by an independent reproduction | §5.2, Supplement §S8 |
| Proposition on the gauge tower: DTΦ_k = α^(−dk)Φ_k for dk ∈ 2ℤ (and pointwise on the odd branch for any real dk > 0) | conditionally RIGOROUS: (a) the algebraic identity — proved, (b) existence of the fixed point g_d at even d — Epstein 1986; uniqueness and hyperbolicity — Lyubich 1999, Avila–Lyubich 2011; d = 2 — Lanford 1982; d = 4 — Burbanks–Osbaldestin–Thurlby 2021 | §2.5 |
| "Zero ↔ parameter" correspondence — probably a new construction | conditional (priority risks not fully closed) | §7.1 |
| OOS "within ×2" model: d = 6 confirmed (~92 % vs baseline ~90.5 %), d = 4 does not show superiority over baseline | NOT RIGOROUS | §7.3 |
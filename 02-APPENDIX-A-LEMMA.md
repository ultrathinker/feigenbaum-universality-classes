# Appendix A. Proposition on the gauge tower for the class |x|^d: statement, proof, numerical check

*Appendix to the article "Explicit realisation of the 'parameter ↔ unstable manifold' correspondence for period-doubling universality classes". The algebraic identity is formally rigorous at a given fixed point g; the numerical part — NOT RIGOROUS (mpmath float); confirmation of the pointwise identity on the odd real branch [G]^{dk+1} for all dk>0, and numerical cross-check with the analytic tower for dk∈2ℤ.*

---

## A.1 Setup and operator

**Class of objects.** Let d > 1 be a parameter of the Cvitanović–Feigenbaum universality class. We consider functions of x ∈ ℝ of the form f(x) = F(|x|^d), where F is a smooth function of t = |x|^d, with the normalisation F(0) = 1 (i.e. f(0) = 1). The function is even in x; in the x-representation it is an even function of |x|^d. This is the working function space: **real-analytic even functions of x of class C^ω** on [−1, 1] (equivalent to real-analytic functions of t ∈ [0, 1]).

**Functional equation and fixed point.** The Cvitanović doubling operator in the project's realisation:

  T(f)(x) = −α_f · f(f(x/α_f)),   α_f = −1/f(1).

In the t-representation F(t) = G(t) — a fixed point, G(0) = 1, and

  G(t) = −α·G(W(t)),   τ = t/α^d,   U = G(τ),   W = |U|^d,

where α = α_G = −1/G(1) is the Feigenbaum constant of the given degree d. From the definition:

  K(t) := G(t) − d·t·G′(t),    a₀ := G(1) = −1/α.

**Linearisation.** On perturbations Φ(t) with Φ(0) = 0 the Fréchet derivative of T at G has the form

  (DTΦ)(t) = −(Φ(1)/G(1))·K(t) − α·Φ(W(t)) − α·d·sign(U)·|U|^{d−1}·G′(W)·Φ(t/α^d).   (∗)

This is equivalent to splitting DT = gauge + J + S: the first term is the gauge contribution of the variation of α_f; J(Φ) = −α·Φ(W) — inner composition; S(Φ) = m(t)·Φ(τ) with m(t) = α^d·G′(t)/G′(τ) — shift of argument. Formula (∗) has been checked by finite-difference cross-check against the true nonlinear operator T in the package `appA_lemma_verify_B/`.

**Transition to the x-representation.** If we set t = |x|^d, W = |U|^d = |G(τ)|^d, τ = t/α^d = |x|^d/α^d = |x/α|^d, then t/α^d = |x/α|^d, and Φ(W(t)) = Φ(|G(|x/α|^d)|^d). Writing g(x) = G(|x|^d) and perturbation φ(x) := Φ(|x|^d) we obtain the doubling operator in the x-representation; this convention is used in A.3 for a short derivation.

---

## A.2 Proposition (gauge tower for the class |x|^d)

**Proposition (pointwise identity + analyticity criterion).**
Let d > 1, α > 0, g be real on [−1, 1], k ≥ 1, and let w_k(x) := x·|x|^{dk} − x. Define

  **Φ_k(t) := G(t)^{dk+1} − d·t^{k+1}·G′(t) − K(t)**,     Φ_k(0) = 0.        (1)

Then for **any** real dk > 0 for Φ_k in the form (1), understood on the odd real branch [G]^{dk+1} := sign(G)·|G|^{dk+1} when dk ∉ 2ℤ (see A.3.2 on the conditions of differentiability and uniqueness of the choice of branch), the **pointwise identity**

  **DTΦ_k(t) = α^{−dk}·Φ_k(t)**    for all t ∈ [0, 1].                    (2)

holds. The function Φ_k in the form (1) **belongs to the real-analytic class** C^ω([0,1]) (and, equivalently, the analytic class of even functions of x) when dk ∈ 2ℤ (an even integer; this is the requirement that the exponent of G(t)^{dk+1} be even, so as to avoid a kink at the zeros of G and at x = 0); the converse (Φ_k ∈ C^ω ⇒ dk ∈ 2ℤ) is established under simple zeros of G on (0, 1). When dk ∈ 2ℤ, Φ_k is an analytic eigenvector of DT with eigenvalue α^{−dk}. On the ordinary integer branch G^{dk+1} with dk ∉ 2ℤ the identity (2) **does not** hold (numerical check: d=3, k=1, dk=3, ordinary branch — rel ≈ 0.465; A.5); on the odd branch [G]^{dk+1} it holds (rel ≈ 5.4·10⁻³⁰ for the same (d, k); A.5.2), but Φ_k ∉ C^ω.

**Corollary 1 (tower at even d).** If d is an even integer, then for all k = 1, 2, 3, … the product dk is automatically even ⇒ **the whole sequence {α^{−dk}}_{k ≥ 1} are eigenvalues of DT** in the analytic class.

**Corollary 2 (sparse tower at odd integer d).** If d is an odd integer, then dk is even ⇔ k is even. Analytic eigenvectors Φ_k in class C^ω exist only at even k; the first level is k = 2 with eigenvalue α^{−2d}.

**Corollary 3 (K → α^d·K).** The trivial coordinate mode K(t) = G(t) − d·t·G′(t) is an eigenvector of the extended operator (without the condition Φ(0) = 0) with eigenvalue α^d. A two-line computation: (1) compare the coefficients of |x|^d in g(x) = −α·g(g(x/α)) ⇒ g′(1) = −α^{d−1} (A.4.4, numerical check up to machine dps 40 in A.5.4); (2) then DT_full(K)/K = α^d (numerical check: α² = 6.264547831… for d = 2; α_6^6 = 9.997677674… for d = 6; A.5.4). This corollary is a direct computation from the structure of the fixed point; it does not require the conjugation χ_ε.

**Corollary 4 (rational and irrational d).** For d rational p/q (in lowest terms) analytic members of the tower appear at dk = p·k/q ∈ 2ℤ, i.e. at k divisible by 2q (if p is odd) or by q (if p is even). The smallest k₀ ∈ ℕ with p·k₀/q ∈ 2ℤ defines the first level. For d irrational no dk is an integer ⇒ Φ_k for all k is not analytic ⇒ no analytic members of the tower.

---

## A.3 Proof

_Navigation._ §A.3.1–A.3.2 use the conjugation construction χ_ε and work for dk ∈ 2ℤ (give "Corollary 1/2" about the tower). §A.3.3 — an independent general proof for all dk > 0 on the odd real branch; logically independent of A.3.1–A.3.2. The numerical check A.5.2 (4 rows, fractional/odd dk) refers precisely to A.3.3.

### A.3.1. Conjugation χ_ε = id + ε·w (≤ 1 page)

Let χ_ε = id + ε·w, w(x) = x^n − x, n = dk + 1. The conditions w(0) = 0 and w(1) = 0 preserve the normalisation f(0) = 1. The family f_ε = χ_ε ∘ g ∘ χ_ε⁻¹ gives φ := d/dε f_ε|_0 = w ∘ g − g′·w. In the t-representation with t = |x|^d, g(x) = G(t) we obtain (1).

### A.3.2. Generator w̃ for T (≤ 1 page)

For f_ε in the class T ∘ χ_ε = χ̃_ε ∘ T, where χ̃_ε = id + ε·w̃ with w̃(x) = w(a·x)/a − w(a)·x/a, a := f(1). With a₀ = G(1) = −1/α and dk ∈ 2ℤ we have a₀^{dk} = (−1/α)^{dk} = α^{−dk}, whence w̃(x) = a₀^{dk}·w(x) = α^{−dk}·w(x). T(f_ε) ≈ f_{ε·α^{−dk}} ⇒ DT maps Φ_k to itself with weight α^{−dk}. ∎

*Applicability conditions.* (i) g is real of class C^1 on [−1, 1]; (ii) when dk ∉ 2ℤ we choose the continuous real branch G^{dk+1} := sign(G)·|G|^{dk+1} — uniquely defined wherever G ≠ 0; at the zeros of G the branch [G]^{dk+1} has a kink exactly at those points, which form the analytic obstruction (A.2). (iii) Complex exponentiation (a₀)^{dk} is used only when dk ∈ 2ℤ: for dk = 2m we have a₀^{dk} = (G(1))^{2m} = (G(1)^2)^m ≥ 0, and no choice of a root branch is required.

### A.3.3. Proof in the x-representation (≤ 1 page)

We use the notation φ(x) := Φ_k(t = |x|^d) — the same eigenfunction in the x-representation, equivalent to (∗) via the chain rule t = |x|^d. In the x-representation:

  DTφ(x) = α·φ(1)·K(x) − α·φ(g(x/α)) − α·g′(g(x/α))·φ(x/α).     (∗∗)

Substitute φ = w ∘ g − g′·w with w(x) = x·|x|^{dk} − x. The gauge term: φ(1) = w(g(1)) = w(−1/α) = −α^{−(dk+1)} + α^{−1}, in total α·φ(1) = 1 − α^{−dk}, so the gauge term = (1 − α^{−dk})·K(x). In the two middle terms the summands ±α·g′(g(x/α))·w(g(x/α)) cancel mutually; the remaining −α·w(g(g(x/α))) and +α·g′(g(x/α))·g′(x/α)·w(x/α) after expanding w by −g(x)/α (from the fixed-point equation g(g(x/α)) = −g(x)/α) give α^{−dk}·φ(x) − (1 − α^{−dk})·K(x). In total:

  DTφ(x) = (1 − α^{−dk})·K(x) + [α^{−dk}·φ(x) − (1 − α^{−dk})·K(x)] = α^{−dk}·φ(x).   (2)

The K-summands cancel between the gauge term and the remainder in one line.

**Linking element (g-chain rule).** Differentiating the fixed-point identity g(g(x/α)) = −g(x)/α with respect to x (g of class C^1, see A.3.2), we get

  g′(g(x/α))·g′(x/α)·(1/α) = −g′(x)/α,   ⇒   g′(g(x/α))·g′(x/α) = −g′(x).

This equality is substituted into (∗∗) together with the differentiation of g(g(x/α)) and underpins the grouping "−α·w(g(g(x/α))) + α·g′(g(x/α))·g′(x/α)·w(x/α) = α^{−dk}·φ(x) − (1−α^{−dk})·K(x)". ∎

*Remark.* A pointwise numerical check on a grid t ∈ (0, 1] for d ∈ {2, 4, 6} at dk ∈ 2ℤ gives a residual ≤ 10⁻²⁵ (N = 48, mpmath dps 40; A.5). Numerical check at higher N and truncation estimate — A.5.1.

---

## A.4 Position of the tower among subleading modes and the shadow as d → 2

### A.4.1. Numbering for d = 2/4/6/8

A summary numbering of the tower mode α_d^{−d} among the subleading eigenvalues of DT at even d (from `sec2_tower_depth_B/`, `sec2_tower_depth_A/`):

| d  | δ       | n_>(d) | pos. without δ | pos. with δ |
|---:|--------:|-------:|-----------:|---------:|
| 2  | 4.6692  | 0      | **1**      | 2        |
| 4  | 7.2847  | 3      | **4**      | 5        |
| 6  | 9.2962  | 5      | **6**      | 7        |
| 8  | 10.9486 | 7      | **8**      | 9        |

**Pattern**: n_>(d) = d − 1 for even d ≥ 4 (checked at 4, 6, 8). At d = 2 — a boundary case (see below), n_>(2) = 0 (`sec2_vector_homotopy/multiplicity_d2.py`).

**Consistent formulation.** For even d ≥ 4 the eigenvalue α_d^{−d} is the **(d+1)-th subleading including δ, d-th without δ** (4/5, 6/7, 8/9 at d = 4, 6, 8). At d = 2 — α_2^{−2} — 1st without δ / 2nd with δ; the B⁺ branch degenerates into the coordinate mode.

**Remark on d = 8.** δ(8) = 10.9486; see `sec2_tower_depth_B/step3_neighborhood.json`. n_>(8) = 7 is confirmed at N = 48/56/64 with a rel-spread of 3·10⁻¹⁷…2·10⁻¹¹; the WINDOW = 0.42 artefact that cut off a real stable mode |λ| = 0.429465 has been fixed.

### A.4.2. Odd d and the shadow as d → 2

For odd integer d the eigenvector Φ_k in the analytic class exists only at even k (dk even ⇔ k even); the first level is k = 2, λ = α⁻²ᵈ. Numerical values from `appA_lemma_derivation_A/` (rel N40/N48): α_3^{−6} = 0.01948837063, 4.4·10⁻²¹/8.0·10⁻²⁰; α_5^{−10} = 0.01203767415, 1.8·10⁻¹⁴/4.5·10⁻¹⁷; α_7^{−14} = 0.008551691871, 7.1·10⁻⁷/4.1·10⁻¹¹; α_{5/2}^{−10} = 5.037362847·10⁻⁴, 9.7·10⁻¹⁶/8.1·10⁻¹⁵.

At d = 2 the coordinate mode α⁻² and the B⁺ branch degenerate; at d = 2.01, 2.05 the B⁺ branch by continuity continues the coordinate mode vector: cos(Φ₁(2), B⁺(2.01)) = 0.999999, cos(Φ₁(2), B⁺(2.05)) = 0.999980 (`sec2_vector_homotopy/`).

**Continuation of the branch at d = 4 and d = 6 — symmetrically** (`sec2_tower_shadow_small_dd/`, NOT RIGOROUS, mpmath dps 40). The cosine is taken between eigen*functions* (by values on a common grid of 120 points t ∈ (0,1) or by coefficients in a shifted Chebyshev basis e_j(t) = T_j(2t−1) − (−1)^j), the mode is traced by continuity from Φ₁(d₀):

| d | N | λ of the branch | cos(Φ₁), values | cos(Φ₂), values | cos(Φ₁), Chebyshev coeffs. | cos(Φ₁), monomial coeffs. |
|---:|---:|---:|---:|---:|---:|---:|
| 4.001 | 48 / 64 | 0.122549065496717 | 1.0 / 1.0 | 0.986 | 1.0 / 1.0 | 0.99999982 / 0.66 |
| 4.01 | 48 / 64 | 0.122976009744813 | 0.9999999 / 0.9999999 | 0.986 | 0.9999998 / 0.9999998 | 0.99998 / 0.47 |
| 6.0005 | 48 / 64 | 0.10003794816778 | 1.0 / 1.0 | 0.986 | 1.0 / 1.0 | 0.24 / 0.98 |
| 6.001 | 48 / 64 | 0.100052667662639 | 1.0 / 1.0 | 0.986 | 1.0 / 1.0 | 0.46 / 0.21 |
| 6.002 | 48 / 64 | 0.100082106593178 | 1.0 / 1.0 | 0.986 | 1.0 / 1.0 | 0.44 / 0.26 |
| 6.005 | 48 / 64 | 0.100170422901463 | 0.99999999 / 0.99999999 | 0.986 | 0.99999998 / 0.99999998 | 0.96 / 0.23 |
| 6.01 | 48 / 64 | 0.100317615071611 | 0.99999996 / 0.99999996 | 0.986 | 0.99999991 / 0.99999991 | 0.50 / 0.26 |

The B⁺ branch near both even points is the same eigenfunction as the tower mode Φ₁ (cos = 1 − O(10⁻⁸)), and is clearly distinct from Φ₂ (0.986). The last column shows that the cosine of the *vectors of coefficients in the monomial basis* for the same mode (the same λ across all signs) at N = 48–64 does not reflect the closeness of the functions (poor conditioning of the basis) — both for d = 6 and d = 4 (N = 64). The cosine of Φ₁ and B⁺ in the monomial basis depends on the choice of basis and is poorly conditioned at N = 48–64; in the functional/Chebyshev basis (the columns above in this table) cos = 1 − O(10⁻⁸).

### A.4.3. Rational d = p/q — sparse analytic members

The smallest k₀ with p·k₀/q ∈ 2ℤ defines the first tower member. **d = 5/2**: k₀ = 4 ⇒ α⁻¹⁰ = 5.037·10⁻⁴. **d = 2.05 = 41/20**: dk = 41k/20 is integer only when 20 | k, even — when 40 | k; the first level is k = 40, λ = α⁻⁸². Not checked numerically at N = 48 (G(t)⁸³ completely loses accuracy; the check would be an artefact of truncation, not a confirmation/refutation of the identity). _Remark (n/a)._ In the cell "numerical value for d = 2.05 at N = 48" we put the mark **n/a** — no numerical value has been obtained in this regime (see above).

### A.4.4. K → α^d·K (Corollary of A.2)

K(t) = G(t) − d·t·G′(t) is an eigenvector of the extended DT_full (without the condition Φ(0) = 0) with α^d: DT_full(K)/K = α^d (d = 2: 6.264547831… = α²; d = 6: 9.997677674… = α_6⁶; spread < 10⁻⁴⁵). A two-line computation: g(x) = −1/α + c_d·|x|^d + …, from g(x) = −α·g(g(x/α)) ⇒ c_d = −α^{d−1}·c_d, g′(1) = d·c_d = −α^{d−1}. Cross-check: A.5.4.

---

## A.5 Numerical check

Metric: `max|DTΦ_k − λΦ_k| / max|Φ_k|` over 240 points t ∈ (0, 1.2] (step 1/200), mpmath dps 40, N = 48 (from `05-CODE/sec2_spectra_grid_d/spectra_grid.json`). Numbers — from `check_lemma_v05.out`. In A.5.1 + A.5.2 (check of (2)) — **13 rows** of real (d, k): 9 in A.5.1 + 4 in A.5.2; across A.5.1–A.5.4 in total — **22 rows** (9+4+6+3); these 22 rows do not include the separately discussed n/a case d = 2.05, k = 40 (see A.4.3): a numerical check at N = 48 is impossible, the cell is n/a.

### A.5.1. Entire branch G^{dk+1}, integer even dk

| d   | k | dk | max\|Φ\|  | max\|DTΦ−λΦ\| | rel       |
|----:|--:|---:|----------:|--------------:|----------:|
| 2.0 | 1 |  2 | 1.030e+00 | 1.265e-30     | 1.228e-30 |
| 2.0 | 2 |  4 | 1.783e+00 | 1.106e-30     | 6.204e-31 |
| 2.0 | 3 |  6 | 2.666e+00 | 9.978e-31     | 3.742e-31 |
| 4.0 | 1 |  4 | 1.933e+00 | 2.096e-29     | 1.084e-29 |
| 4.0 | 2 |  8 | 2.835e+00 | 2.043e-29     | 7.204e-30 |
| 4.0 | 3 | 12 | 4.351e+00 | 1.976e-29     | 4.541e-30 |
| 6.0 | 1 |  6 | 2.909e+00 | 5.293e-20     | 1.819e-20 |
| 6.0 | 2 | 12 | 4.036e+00 | 5.926e-20     | 1.468e-20 |
| 6.0 | 3 | 18 | 5.492e+00 | 6.002e-20     | 1.093e-20 |

All residuals ≤ 2·10⁻²⁰; at N = 72 we obtain ≤ 10⁻⁴⁹ (for d = 6).

*Truncation estimate.* The residual |DTΦ_k − λΦ_k| is determined by the accuracy of computing DT in a collocation approximation of g at order N. At d = 6, N = 48 for the subbranch λ₅ (lam_M48_f = 0.2078466910229993, cf. A.5.1 row 7) the dense residual of the operator T at M = 48 is ≈1.35·10⁻³² (`05-CODE/sec2_spectrum_d6_T22/spectrum_d6.json` field `lambda5_new.dense_resid_M48`); the numerical rel ≈ 1.8·10⁻²⁰ in the table corresponds to truncating the power series at the tail terms and the convergence of g at N = 48 ≈ mpmath dps 20 for exponent 6 (g^{(40)}(ξ)/40! in the worst case). Since for d = 2, 4 the residual at N = 48 is already ≤ 10⁻²⁹, an additional check at higher N (56/64/72) does not change the qualitative conclusion of A.5.5; measurements at N = 72 give a residual ≤ 10⁻⁴⁹ for d = 6, consistent with the convergence interpolation "doubling N → ~10 orders" (power law by analyticity radius of g — internal `05-CODE/sec2_spectra_grid_d/spectra_grid.py` field `fp_dense_resid_N48`).

### A.5.2. Odd real branch [G]^{dk+1}, fractional/odd dk

| d     | k | dk   | max\|Φ\|  | max\|DTΦ−λΦ\| | rel       |
|------:|--:|-----:|----------:|--------------:|----------:|
| 2.5   | 1 | 2.5  | 1.259e+00 | 2.120e-30     | 1.684e-30 |
| 2.05  | 1 | 2.05 | 1.054e+00 | 6.981e-31     | 6.626e-31 |
| 3.0   | 1 | 3    | 1.482e+00 | 7.970e-30     | 5.376e-30 |
| 5.0   | 1 | 5    | 2.406e+00 | 7.447e-24     | 3.095e-24 |

**Key observation:** on the odd real branch [G]^{dk+1} the residual for fractional and odd dk is ~10⁻³⁰…10⁻²⁴, of the same order as for integer even dk; see the table above (homogeneity of rows 2, 5 and row 4 of the table). On the **main complex branch** G^{dk+1} (in contrast to the odd real [G]^{dk+1}) the residual for the same rows 2, 5 is of order 0.05–0.45. There is **no algebraic** obstruction on the odd branch; the only obstruction to Φ_k ∈ C^ω is **analytic** (a kink at the zeros of G on (0, 1) when dk ∉ 2ℤ; see the analyticity criterion in A.2).

### A.5.3. x-representation: Φ_k(t) ≡ w ∘ g − g′·w

Checking the agreement of Φ_k(t = |x|^d) with φ_x(x) = w(g(x)) − g′(x)·w(x), w(x) = x|x|^{dk} − x, over 100 points x ∈ (0, 1] (domain t = |x|^d ∈ (0, 1], step 1/100):

| d   | k | max\|φ_x\| | max\|φ_x − Φ_t\| | rel       |
|----:|--:|-----------:|-----------------:|----------:|
| 2.0 | 1 | 1.030e+00  | 4.592e-41        | 4.456e-41 |
| 2.0 | 2 | 1.359e+00  | 6.888e-41        | 5.070e-41 |
| 4.0 | 1 | 1.933e+00  | 9.184e-41        | 4.751e-41 |
| 4.0 | 2 | 2.587e+00  | 1.378e-40        | 5.324e-41 |
| 6.0 | 1 | 2.908e+00  | 9.184e-41        | 3.158e-41 |
| 6.0 | 2 | 4.036e+00  | 1.378e-40        | 3.413e-41 |

For all six combinations (d, k) ∈ {2, 4, 6} × {1, 2} the residuals are ≤ 6·10⁻⁴¹ — the identity Φ_k(t = |x|^d) ≡ w(g(x)) − g′(x)·w(x) holds exactly (domain x ∈ (0, 1]; `check_lemma_v05.py`, all six points give ~3–5·10⁻⁴¹).

### A.5.4. g′(1) = −α^{d−1}

| d   | g′(1)            | −α^{d−1}         | difference |
|----:|-----------------:|-----------------:|---------:|
| 2.0 | −2.502907875096  | −2.502907875096  | 0        |
| 4.0 | −4.829405415303  | −4.829405415303  | 0        |
| 6.0 | −6.811602179884  | −6.811602179884  | 0        |

(g(x) = G(|x|^d) ⇒ g′(1) = d·G′(1).) Exact equality in all digits up to machine dps 40. Corollary of A.2 (K → α^d·K).

### A.5.5. Summary picture

- Entire branch, even dk (9 combinations): rel ≤ 2·10⁻²⁰, the identity is confirmed.
- Odd real branch, fractional/odd dk (4 combinations): rel ≤ 8·10⁻²⁴, the identity holds pointwise.
- Analytic spectrum of DT: Φ_k ∈ C^ω ⇔ dk ∈ 2ℤ; Φ_k is an eigenvector in the analytic class ⇔ dk ∈ 2ℤ (the analyticity criterion and the spectral criterion coincide — A.2).

All numbers — from `check_lemma_v05.out` (mpmath dps 40, N = 48, 240 points). Auxiliary numbers from `appA_lemma_derivation_A/` are in A.4.2; no discrepancies.

---

## A.6 Attribution, rigor status and open questions

### A.6.1. Attribution to the classics

Eigenvalues of coordinate changes x → x + εx^n with value α^{1−n} are a classical fact: Feigenbaum (1979), Collet–Eckmann–Lanford (1980), Cvitanović "Universality in Chaos" (1989), Christiansen–Cvitanović–Rugh (1990, *J. Phys. A* 23, L713), Briggs (1994, arXiv:chao-dyn/9407019) explicitly speaks of "gauge". The condition that n be odd is folklore of doubling theory. The generalisation to the class |x|^d is a consequence of the construction χ_ε = id + ε·(x|x|^{dk} − x), carried out in this work.

**What is new in this work:** (1) the explicit t-form Φ_k (1) for direct numerical check; (2) the analyticity criterion Φ_k ∈ C^ω ⇔ dk ∈ 2ℤ (a unified criterion of analyticity and belonging to the spectrum of DT); (3) numerical cross-check of the tower at d = 2, 4, 6, 8 (A.4.1) and at odd d (A.4.2).

### A.6.2. Rigor status

**The algebraic identity** (2) DTΦ_k = α^{−dk}·Φ_k on the odd real branch [G]^{dk+1} is **rigorous** at a given fixed point G: it suffices that G be real and dk > 0, so that substituting Φ_k into (∗) gives a cancellation. Formally proved in A.3 for all d > 1 and dk > 0 on the branch [G]^{dk+1}; in the analytic class — for dk ∈ 2ℤ.

**Rigor overall** (conditionally RIGOROUS) = rigor of the existence of the fixed point g_d of the operator T(f) = −α_f·f(f(x/α_f)). Two sub-steps: (a) the **algebraic identity** (2) is formally rigorous at a given fixed point G; (b) the **existence and uniqueness of g_d** at even d — from the literature. Identity (2) is a formal identity for this fixed point; its numerical check was performed on a collocation approximation of g_d (tables A.5.1/A.5.2; residual ≤ 2·10⁻²⁰ for integer even dk and ≤ 8·10⁻²⁴ for fractional/odd dk on the odd real branch, as before). In the literature: existence at even d ≥ 4 — Epstein (1986, *Comm. Math. Phys.* 106, 395–426; analytic proof of existence of the functional equation in the class of even functions of general degree at the critical point); uniqueness and hyperbolicity of the fixed point for an arbitrary even integer degree — Lyubich (1999, *Ann. of Math.* 149, 319–420) and Avila–Lyubich (2011, *Publ. Math. IHÉS* 114, 171–223; renormalisation horseshoe, exponential contraction along hybrid classes); d = 2 — Lanford (1982, *Bull. Amer. Math. Soc.* 6, 427; computer-assisted proof); d = 4 — Burbanks–Osbaldestin–Thurlby (2021, *J. Math. Phys.* 62, 112701; CAP bounds on the fixed point, α, δ and eigenfunctions at degree 4, more than 320 significant digits). For d = 8 numerical "certificates" (`sec2_tower_depth_A/`, `sec2_tower_depth_B/`): dense residual 6.96·10⁻³⁶ at N = 48/56; agreement with an independent solution to 17 significant digits — numerical, not interval verification.

**Summary without overstating:** the identity (2) is proved pointwise on the odd branch for dk > 0; its validity for the concrete numerical g_d is confirmed to the accuracy of the collocation input data; the existence of the fixed point g_d at even d is covered by the classical theorems of Epstein (1986), Lyubich (1999) and Avila–Lyubich (2011).

### A.6.3. Open questions

1. **Behaviour of the physical branch B⁺(d) near even d — analytically.** Numerically the branch continues the tower mode Φ₁ symmetrically at d = 4 and d = 6 (A.4.2; cosine of eigenfunctions = 1 − O(10⁻⁸) for Δd ∈ [0.0005, 0.01]); Proposition A.2 gives only the identity at the even points themselves and says nothing about the branch at Δd ≠ 0. An analytic description (e.g. an expansion of B⁺(d₀ + Δd) in Δd) has not been obtained. The cosine of the coefficients of Φ₁ and B⁺ in the monomial basis is a basis-dependent quantity and does not testify to an asymmetry between d = 4 and d = 6; in the Chebyshev basis cos = 1 − O(10⁻⁸) (A.4.2).

2. **Rigorous identification of our collocation g_d with the exact fixed point of the operator T** (estimate of the N → ∞ truncation error between the collocation approximation and the limiting object) — not carried out; a numerical check of identity (2) was performed on the collocation approximation (A.5.1, A.5.2).

3. **d = 2.05, k = 40** — the first member of the "sparse" tower (α⁻⁸²); a numerical check is impossible at N = 48 (G(t)⁸³ loses accuracy). Requires N ≥ 100 or a special construction.

---

## References

References — common list `04-REFERENCES.md`.
# Explicit realisation of the "parameter ↔ unstable manifold" correspondence for period-doubling universality classes

**Yevhen Borzenkov** · ORCID [0009-0002-7946-589X](https://orcid.org/0009-0002-7946-589X)
Preprint, version 0.7.4 (frozen 16.08.2026), not peer-reviewed. Public edition 15.11.2026.
Русская версия — [README.ru.md](README.ru.md).

Companion to [feigenbaum-h-zeros](https://github.com/ultrathinker/feigenbaum-h-zeros) (quadratic class d = 2).
Here the same zero pattern of the universal function H_d is followed through the classes |x|^d with d = 2, 4, 6:
zeros of H_d encode the geometry of the parameter space of each class, not only of the quadratic one.

## Main results

Rigorous (computer-assisted, interval Krawczyk operator and interval Newton, python-flint / Arb):

- **357 certified boxes** of superstable parameters: 61 at d = 4 (period 1024), 287 at d = 6 (period 1024),
  9 late ones at d = 6 (periods 2048, 4096, 8192); existence, uniqueness and exact period for each.

Numerically verified (certified eigenvalue of the collocation discretization; the N → ∞ truncation is estimated, not certified):

- λ₂(4) = 0.2918384089835451…, λ₂(6) = 0.3710714676819301… (and λ₂(2) from the companion "Bridge" computation).

Analytical: a proposition on the gauge tower — Φ_k is an exact eigenvector of DT with eigenvalue α^(−dk)
when dk ∈ 2ℤ (Appendix A).

Not rigorous, stated as such: the zero ↔ parameter correspondence itself (22/23 at d = 2, 61/62 at d = 4,
287/294 at d = 6), completeness of the catalogues, the spectral grid d ∈ [1.9, 7] with crossover d* = 2.3616 ± 0.0008,
and the rate law rate = min(δ, 1/|λ₂|) in our data. The rate law itself is known (Mao–Hu 1987, Reick 1992, Briggs 1994);
existence and hyperbolicity of the fixed points at even d — Epstein, Lyubich, Avila–Lyubich, Burbanks–Osbaldestin–Thurlby.
Full list of limitations — §7.4 of the article.

## Repository layout

| Path | Contents |
|---|---|
| [`00-README.md`](00-README.md) | Package overview, what is and is not proved, reproduction commands |
| [`01-ARTICLE.md`](01-ARTICLE.md) | Main article |
| `02-APPENDIX-A-LEMMA.md`, `03-APPENDIX-B-DATA.md` | Appendix A (gauge tower), Appendix B (data and constants for the three classes) |
| [`04-REFERENCES.md`](04-REFERENCES.md), `05-SUPPLEMENT-METHODS.md` | References, methodological supplement |
| [`05-CODE/`](05-CODE/) | Replication kit: 35 subfolders with code, certificates, reports, tamper batteries; index [`05-CODE/README-CODE.md`](05-CODE/README-CODE.md) |
| `06-EXTERNAL-INPUTS/` | Inputs for d = 2 copied from the companion projects |
| [`ru/`](ru/) | Russian version of the texts (code is shared) |
| `MANIFEST.sha256`, `_tools/verify_manifest.py` | sha256 of every file and its verifier |
| [`PUBLICATION-NOTES.md`](PUBLICATION-NOTES.md) | What was changed for this public edition |

## Reproducing

Python 3.13, `python-flint` 0.9.0, `mpmath` 1.3.0, `numpy`, `scipy`. Run scripts from their own folder
(`python -X utf8 <script>`). Quick checks:

```
python _tools/verify_manifest.py
cd 05-CODE/sec3_krawczyk_d4      && python krawczyk4_verify.py        # 61 boxes, ~1 s
cd ../sec3_krawczyk_d6           && python krawczyk6_verify.py        # 287 boxes, ~5 s
cd ../sec3_krawczyk_d6_late      && python krawczyk6_late_verify.py   # 9 boxes, ~1 s
```

All four pass on this edition (15.11.2026).

## Use of AI tools

The computations, code and text were prepared with the help of AI tools; every claim marked rigorous rests
on an interval-arithmetic certificate that is reproduced by the code in this repository. Some reports in
`05-CODE/` name the tool that produced them.

## Related work by the author

- Zeros of the universal scaling function, class d = 2 — [feigenbaum-h-zeros](https://github.com/ultrathinker/feigenbaum-h-zeros).
- Hyperbolicity of period-doubling renormalization at non-integer critical exponent ℓ ∈ [2, 4] —
  [feigenbaum-hyperbolicity](https://github.com/ultrathinker/feigenbaum-hyperbolicity), Zenodo DOI
  [10.5281/zenodo.22306911](https://doi.org/10.5281/zenodo.22306911).

## Citation and license

Citation metadata — [`CITATION.cff`](CITATION.cff). License — MIT ([`LICENSE`](LICENSE)), for code, data and text.

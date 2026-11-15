# Public edition: what was changed

This repository is the package of version 0.7.4 (frozen 16.08.2026), prepared for publication on 15.11.2026.
The frozen manifest is kept as `_tools/MANIFEST.original.sha256`; the current `MANIFEST.sha256` covers every
file of the repository. Of the 486 entries of the frozen manifest, **483 are byte-identical** here,
2 were modified and 1 removed.

## 1. Removed

An internal review note of an intermediate round; the corrections it led to are incorporated in the v2 reports
of `05-CODE/sec7_principle_status_v2/`, which keep the citation of the note as a historical reference. The working folder of the translation and editorial review (never part of the frozen
manifest) was not published either.

- `05-CODE/sec7_principle_status_check/REVIEW-F1F2.md`

## 2. Modified — a path defect of the frozen package, fixed

- `05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_certificate.py`
- `05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_verify.py`

When absolute paths were relinked to package-relative ones (journal `05-CODE/RELINK-LOG.md`), a trailing `/`
was lost in one implicitly concatenated string literal, so both scripts looked for
`../sec3_krawczyk_d6krawczyk6_certificates.json` and the verifier of the 9 late boxes could not run. The fix adds
that one character in each file (`r"../sec3_krawczyk_d6/"`). Neither script is pinned by any certificate.

## 3. Added

- `.gitattributes`
- `CITATION.cff`
- `LICENSE`
- `PUBLICATION-NOTES.md`
- `README.md`
- `README.ru.md`
- `_tools/MANIFEST.original.sha256`
- `ru/00-README.md`
- `ru/01-ARTICLE.md`
- `ru/02-APPENDIX-A-LEMMA.md`
- `ru/03-APPENDIX-B-DATA.md`
- `ru/04-REFERENCES.md`
- `ru/05-CODE/README-CODE.md`
- `ru/05-SUPPLEMENT-METHODS.md`
- `ru/06-EXTERNAL-INPUTS/README-EXTERNAL.md`

`.gitattributes` stores every file byte-for-byte, so the sha256 checks hold after cloning on any platform.
`ru/` holds the Russian text of the article, appendices and the two README files of the code; the code is shared.

## 4. Checks run on this edition (15.11.2026)

On a copy of this repository (Windows, Python 3.13.13, python-flint 0.9.0):

| Check | Result |
|---|---|
| `_tools/verify_manifest.py` | 0 problems |
| `05-CODE/sec3_krawczyk_d4/krawczyk4_verify.py` | `VERDICT: PASS` — 61 roots recomputed, isolation over 1830 pairs |
| `05-CODE/sec3_krawczyk_d6/krawczyk6_verify.py` | `VERDICT: PASS` — 287 roots recomputed, isolation over 41041 pairs |
| `05-CODE/sec3_krawczyk_d6_late/krawczyk6_late_verify.py` | `VERDICT: PASS` — 9 roots recomputed (after the fix of section 2) |

Other components (spectral grid, Bridge-2, tamper batteries) were not rerun for this edition.

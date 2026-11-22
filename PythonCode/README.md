# Birthday Polynomials and Newton's Method: code and data

Code and data for the article

> Americo Cunha Jr., *Birthday Polynomials and the Boundaries of Newton's Method*,
> submitted to The Mathematical Intelligencer (2026).

Laboratório Nacional de Computação Científica, Petrópolis, Brazil;
Universidade do Estado do Rio de Janeiro, Rio de Janeiro, Brazil.
Contact: americo@lncc.br · ORCID 0000-0002-8342-0363

The eight digits of a date, DD/MM/YYYY, are read as the coefficients of a polynomial,
and each point of the complex plane is colored by the root that Newton's method reaches
from it. This folder reproduces every figure, table, and number of the article and of its
Online Resource 1. All computations are deterministic.

The interactive website https://americocunha.org/BirthdayFractal/ is a separate,
simplified JavaScript implementation, kept in the root of this repository; this folder,
`PythonCode/`, holds the code of the article.

## Contents

| Path | Purpose | Output |
|---|---|---|
| `birthday/dates.py` | The 38 dates of the article; date → polynomial | |
| `birthday/newton.py` | Newton's method on a grid, classification, colors, pictures | |
| `birthday/critical.py` | Critical points of the Newton map: exact classification, numerical test | |
| `render.py` | Picture of any date | `figures/newton_DD-MM-YYYY.png` |
| `make_figures.py` | Figs. 1–5 and 7–12 of the article | `figures/` |
| `make_portraits.py` | Fig. 6 (needs the source images) | `figures/` |
| `certify_trapping.py` | Exact certificate for Proposition 3 and Corollary 1 | `results/trapping_certificate.json` |
| `certify_orbit_intervals.py` | Critical-orbit intervals of 20/11/1924 (Table S5) | `results/orbit_intervals_certificate.json` |
| `lyapunov.py` | Lyapunov indicator for 20/11/1924 (Table S4) | `results/lyapunov.json` |
| `classify_pictures.py` | Classification of the points of each picture (Table S6) | `results/picture_classification.json` |
| `make_supplement.py` | Tables S1–S6 and the per-date data | `supplement/ESM_1.tex`, `results/dates.csv` |
| `audit.py` | Recomputes the 87 numbers stated in the article and Online Resource 1 | `results/audit.txt` |

`results/` holds the outputs used in the article, and `supplement/ESM_1.tex` is the source
of Online Resource 1.

## Requirements

Python 3.13 and the packages in `requirements.txt` (`pip install -r requirements.txt`).
`make_figures.py` also needs LaTeX with the packages `helvet`, `sansmath`, and `cm-super`.

## Quick start

```
python render.py 30/08/1984                                   # your birthday here
python render.py 30/08/1984 --order YMD                       # digits in YYYYMMDD order
python render.py 30/08/1984 --relax 0.9+0.2j                  # relaxed step z - a p/p'
python render.py 03/07/1986 --center -0.5 --half-width 0.05   # magnify a region
```

## Reproducing the article

Run from this folder, in this order. Times are for one core.

| Step | Command | Time |
|---|---|---|
| 1 | `python certify_trapping.py` | 5 s |
| 2 | `python certify_orbit_intervals.py` | 1 s |
| 3 | `python lyapunov.py` | 1 min |
| 4 | `python classify_pictures.py` | 7 min |
| 5 | `python make_supplement.py`, then `pdflatex ESM_1` twice in `supplement/` | 5 s |
| 6 | `python audit.py` | 15 s |
| 7 | `python make_figures.py` | 8 min |
| 8 | `python make_portraits.py SOURCE_FOLDER` | 5 s |

Step 4 can run in parallel: `python classify_pictures.py 0 2` and
`python classify_pictures.py 1 2` in two terminals, then `python classify_pictures.py merge 2`.
Step 7 accepts figure names, for example `python make_figures.py Fig7f Fig12`.
Step 8 expects the source images under the file names listed in `make_portraits.py`;
they are not included (see the caption of Fig. 6 for their origin and licenses).

## Notes

- The convergence of the free critical points in Table S3 is a finite numerical criterion,
  not a proof. The nonconvergence for 20/11/1924 is proved, by `certify_trapping.py`.
- Critical poles are identified in exact arithmetic and are not iterated numerically.
- In the certificates, no inequality is decided in ordinary floating point: they use exact
  rational arithmetic, Sturm sequences, and interval arithmetic with outward rounding.
- The results were produced with NumPy 2.5.3, SymPy 1.14.0, mpmath 1.3.0, Matplotlib 3.11.2,
  and Pillow 12.3.0, under Python 3.13.16 and again under 3.13.18, with identical results.

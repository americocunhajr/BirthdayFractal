# Birthday Polynomials and Newton's Method: code and data

Online Resource 2 of the article

> Americo Cunha Jr., *Birthday Polynomials and the Boundaries of Newton's Method*,
> submitted to The Mathematical Intelligencer.

Laboratório Nacional de Computação Científica (LNCC), Petrópolis, Brazil;
Universidade do Estado do Rio de Janeiro (UERJ), Rio de Janeiro, Brazil.
Contact: americo.cunhajr@gmail.com

The eight digits of a date, written as DDMMYYYY, are used as the coefficients of a
polynomial, and each point of the complex plane is colored by the root to which
Newton's method converges. This package contains the programs and data needed to
reproduce the figures of the article, Tables S1–S6 of its Online Resource 1, and the
numbers quoted in the text. The computations are deterministic and use no random numbers.

An interactive browser version of the construction, written separately in JavaScript, is at
https://americocunhajr.github.io/BirthdayFractal/ (source code:
https://github.com/americocunhajr/BirthdayFractal). The figures and numbers of the article
come from the programs in this package.

The programs were first written for a Portuguese version of the article, so most file
names, identifiers, comments, and console messages are in Portuguese. The glossary at
the end translates the recurring terms. The scripts that produce English output are
`supplement_en.py`, `figuras_vetoriais_en.py`, and `make_submission_figures.py`.

## Quick start

Fractal of one date (PNG, 1600 × 1600, in `figs/`):

```
python newton_arte.py 30/08/1984 --saida figs
```

Variations mentioned in Appendix A of the article:

```
python newton_aniversario.py 30/08/1984 --ordem AMD --saida figs     # digits in YYYYMMDD order
python newton_aniversario.py 30/08/1984 --a 0.9+0.2j --saida figs    # relaxed step z - a p/p'
python -c "from newton_aniversario import coeficientes; from newton_arte import renderizar_coef; \
renderizar_coef(coeficientes('03/07/1986'), centro=complex(-0.5, 0), R=0.05)[0].save('figs/zoom.png')"
```

The last command magnifies a stretch of the real axis for a date without real roots
(`centro` and `R` are the center and half-width of the window).

## Dependencies

Versions used for the submitted version:

| Component | Version |
|---|---|
| Python | 3.13.16 |
| NumPy | 2.5.3 |
| Matplotlib | 3.11.2 |
| Pillow | 12.3.0 |
| mpmath | 1.3.0 |
| SymPy | 1.14.0 |
| pdfTeX | TeX Live 2023 |

`figuras_extras.py` and `figuras_vetoriais_en.py` call LaTeX through Matplotlib and need
the packages `lmodern`, `cm-super`, `helvet`, and `sansmath`.

## Reproduction, in order

All commands run in this `research/ESM_2/` folder, unless otherwise indicated.

From a fresh clone of the BirthdayFractal repository:

```bash
git clone https://github.com/americocunhajr/BirthdayFractal.git
cd BirthdayFractal/research/ESM_2
python -m pip install -r requirements.txt
```

1. **Lyapunov indicator for 20/11/1924 (Table S4)**, about 1 min:
   `python analise_mandelbrot.py` → `mandelbrot_orbita.json`

2. **Certificates for 20/11/1924**, under 1 s each:
   `python certifica_aprisionamento.py` → `certificado_aprisionamento.json`
   Exact certificate for Proposition 3 (Table A1): Sturm counts of $p$, $p'$, and
   $q_t(x)=(x-t)p'(x)-p(x)$ on the rational intervals $K_\pm$, exact values of $N_p$ at
   two endpoints, and the distance bounds (root counts on $K_\pm$ enlarged by 0.55 and
   0.21). Also certifies Corollary 1: the degree-49 numerator of $N_p^2(x)-x$ has exactly
   one zero in $K_-$, and the multiplier of the resulting 2-cycle, enclosed by interval
   arithmetic, is less than $-1$.

   **Finer verification (Table S5)**:
   `python certifica_mandelbrot.py` → `certificado_mandelbrot.json`
   Exact rational arithmetic and Sturm sequences (SymPy), plus interval arithmetic with
   80 digits (`mpmath.iv`). The field `exato` holds the exact rational endpoints of all
   intervals, the isolating intervals of $c$, of the real root of $p$, and of the real
   zeros of $p'$, the Sturm counts on $J_-$ and $J_+$, the rational bounds for each
   inequality, and the degree of $\gcd(p',p'')$.

3. **Classification of the points of each fractal picture (Table S6)**, about 50 min on
   one core:
   - sequentially: `python auditoria_figuras.py`
   - or in parallel: `python auditoria_figuras.py 0 2` and `python auditoria_figuras.py 1 2`,
     then `python auditoria_figuras.py juntar 2`

   Output: `auditoria_figuras.json`. The submitted version was produced in parallel.

4. **Tables S1–S6 and the per-date CSV**, a few seconds:
   `python supplement_en.py` → `supplementary_dates.csv` and the LaTeX code of
   Tables S1–S6, written into `esm/ESM_1.tex` between the lines `% BEGIN TABLES S1-S6`
   and `% END TABLES S1-S6`. Then `cd esm; pdflatex ESM_1.tex` (twice) gives Online
   Resource 1.
   Tables S1–S3 are computed here (exact classification of the critical points of $N_p$
   with SymPy, then numerical iteration of the regular free critical points, using the
   functions of `analise_suplementar.py`). Tables S4–S6 are read from the JSON files of
   steps 1–3.

5. **Audit of the numbers in the text**, a few seconds:
   `python auditoria_datas.py` → `auditoria_datas.txt`
   Recomputes 87 numerical statements of the article, including the propositions, and marks
   each one OK or FALHA (failure). All 87 are OK in the submitted version.

6. **Figures**
   - Date fractals, one call per date:
     ```
     for d in 25/12/1642 04/01/1643 30/04/1777 15/04/1707 16/08/1821 29/04/1854 \
              28/02/1878 03/02/1893 20/11/1924 05/07/1687 23/09/1846 30/11/1889 \
              29/05/1919 29/12/1972 16/05/1980 08/08/1900 30/08/1984 03/07/1986 \
              27/04/1952 15/08/1950; do python newton_arte.py $d --saida figs; done
     python newton_arte.py 08/08/1900 --sem-normalizar --saida figs
     ```
   - Cayley examples (Fig. 3), magnifications of the picture for Gauss's date of birth
     (Fig. 4), and the failure example $z^3-2z+2$ (Fig. 11): `python figuras_extras.py`
     (it also writes Portuguese-labelled versions of Figs. 2 and 5, not used in the
     English article).
   - Figs. 2 and 5 with English Helvetica lettering: `python figuras_vetoriais_en.py`
     → `submission/Fig2.pdf`, `submission/Fig5.pdf`
   - Fig. 12 (trapped critical orbit): `python figura_aprisionamento_en.py` → `submission/Fig12.pdf`
   - Portraits of Fig. 6: `python retratos.py mono`. This needs the original images,
     which are not included; their path is set at the top of the program.
   - Files with the journal's names (`Fig1a.jpg`, …, `Fig11.jpg`):
     `python make_submission_figures.py` → `submission/`

## Data files

| File | Content |
|---|---|
| `supplementary_dates.csv` | One row per date (38 dates): polynomial, degree, distinct roots, multiplicity of $z=0$, window, Girard heuristic, classification of the critical points of $N_p$ (total $2k-2$, regular free critical points, critical poles), and convergence of the regular free critical points. |
| `mandelbrot_orbita.json` | Runs of the indicator $\lambda_n$ for 20/11/1924 (precision, starting point, $n$, $\lambda_n$, visited interval). |
| `certificado_aprisionamento.json` | Exact certificate for Proposition 3: intervals, checks, and exact values of $N_p$ at two endpoints. |
| `certificado_mandelbrot.json` | Finer verification for 20/11/1924: summary of each check and the exact record (field `exato`). |
| `auditoria_figuras.json` | Number of points in each class, (a) to (d), for each fractal picture. |
| `auditoria_datas.txt` | Result of the audit of step 5. |

## Programs

| File | Role |
|---|---|
| `newton_aniversario.py` | Date → coefficients (`coeficientes`), window rule with isolated roots (`raizes_da_janela`). |
| `newton_arte.py` | Newton iteration on the grid, classification, smoothed count, OKLCH colors; command line for one date. |
| `figuras_extras.py` | Figs. 3, 4, and 11 (and Portuguese versions of Figs. 2 and 5). |
| `figuras_vetoriais_en.py` | Figs. 2 and 5 in English. |
| `retratos.py` | Cropping and monochrome toning of the portraits of Fig. 6. |
| `make_submission_figures.py` | Copies the figures under the journal's file names. |
| `analise_suplementar.py` | Exact classification of critical points (`pontos_criticos`) and numerical test of the regular free ones (`teste_criticos`); used by `supplement_en.py` and `auditoria_datas.py`. Run on its own, it writes the Portuguese tables. |
| `analise_mandelbrot.py` | Lyapunov indicator for 20/11/1924. |
| `certifica_aprisionamento.py` | Exact certificate for Proposition 3 and Corollary 1. |
| `certifica_mandelbrot.py` | Finer verification of the minimal invariant intervals $I_\pm$. |
| `figura_aprisionamento_en.py` | Fig. 12. |
| `auditoria_figuras.py` | Classification counts of every fractal picture. |
| `auditoria_datas.py` | Audit of the numbers quoted in the article. |
| `supplement_en.py` | English Tables S1–S6 and `supplementary_dates.csv`. |
| `esm/ESM_1.tex` | Source of Online Resource 1. |

## Remarks

- The convergence of the regular free critical points (Table S3) is a finite numerical
  criterion, not a proof.
- Critical poles, of which the only case is $z=-1$ for 11/02/2016, are identified in exact
  arithmetic and are not iterated numerically.
- The confinement of the critical orbit of 20/11/1924 is a computer-assisted verification.
  No inequality was decided in ordinary floating point.

## Glossary of Portuguese terms

| Portuguese | English |
|---|---|
| auditoria | audit |
| bacia | basin |
| certificado, certifica | certificate, certifies |
| datas | dates |
| falha | failure |
| figuras | figures |
| janela, centro, meia largura | window, center, half-width |
| livres regulares | regular free (critical points) |
| polos críticos | critical poles |
| raízes, raiz | roots, root |
| retratos | portraits |
| sem normalizar | without normalization |
| suplemento, suplementar | supplement, supplementary |
| OK / FALHA | passed / failed |

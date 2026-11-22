# Scientific code and reproducibility

BirthdayFractal contains two distinct implementations:

| Component | Location | Purpose |
| --- | --- | --- |
| Browser application | [`app.js`](app.js), [`worker.js`](worker.js) | Interactive, approximate Newton-fractal visualization |
| Scientific Python package | [`PythonCode/`](PythonCode/) | Figures, tables, audits, exact certificates, and data for the companion paper |

The original Portuguese-language research package has been **replaced** by the author's English-language package supplied as `BirthdayFractal_paper.zip`. Its contents are published as individual source files in [`PythonCode/`](PythonCode/). The package's own [README](PythonCode/README.md) documents all scripts and the complete reproduction sequence.

## Installation

Use **Python 3.13** with the [pinned requirements](PythonCode/requirements.txt). From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate             # Windows: .venv\Scripts\activate
python -m pip install -r PythonCode/requirements.txt
cd PythonCode
```

## Key programs

| Command | Purpose | Output |
| --- | --- | --- |
| `python render.py 30/08/1984` | Render a Newton fractal | `figures/newton_30-08-1984.png` |
| `python certify_trapping.py` | Exact proof of the trapping and repelling 2-cycle | `results/trapping_certificate.json` |
| `python certify_orbit_intervals.py` | Certified critical-orbit intervals | `results/orbit_intervals_certificate.json` |
| `python lyapunov.py` | Finite-time Lyapunov indicators | `results/lyapunov.json` |
| `python classify_pictures.py` | Classify image samples for Table S6 | `results/picture_classification.json` |
| `python make_supplement.py` | Supplementary tables and dates | `supplement/ESM_1.tex`, `results/dates.csv` |
| `python audit.py` | Audit 87 numerical statements | `results/audit.txt` |
| `python make_figures.py` | Generate article figures | `figures/` |
| `python make_portraits.py SOURCE_FOLDER` | Prepare Figure 6 portraits | `figures/` |

The reusable [`birthday/`](PythonCode/birthday/) Python package provides date coefficients, Newton rendering, and critical-point classification. The [`results/`](PythonCode/results/) directory contains the provided numerical audits and certificate records. Online Resource 1 has the LaTeX source [`supplement/ESM_1.tex`](PythonCode/supplement/ESM_1.tex).

## Verification

The automated [research-code workflow](.github/workflows/research-check.yml) checks all **13 Python source files**, the **four JSON results**, and **38 dates** in the CSV. It then recomputes both exact certificates and requires their `certified` fields to be `true`.

These certificate calculations use exact rational/Sturm and interval arithmetic. Finite-iteration image classifications and Lyapunov averages are numerical evidence rather than proofs. Historical portraits used in Figure 6 are not included, and some figure scripts require a LaTeX installation.

**Full instructions:** [PythonCode/README.md](PythonCode/README.md).

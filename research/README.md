# Research code and reproducibility

The [Online Resource 2 source package](ESM_2/) accompanies the manuscript *Birthday Polynomials and the Boundaries of Newton’s Method* by Americo Cunha Jr., submitted to *The Mathematical Intelligencer* (under review).

This directory contains the **original Python research programs**, the supplementary LaTeX source, the 38-date CSV, audit results, and machine-readable exact certificates supplied by the author. It is distinct from the approximate **JavaScript** fractal generator in the repository root.

## First steps

```bash
git clone https://github.com/americocunhajr/BirthdayFractal.git
cd BirthdayFractal/research/ESM_2
python -m pip install -r requirements.txt
python certifica_aprisionamento.py
python certifica_mandelbrot.py
```

See **[ESM_2/README.md](ESM_2/README.md)** for the full command sequence (figures, exact certificates, Tables S1–S6, numerical audits, and per-date rendering), dependencies, output-file map, and Portuguese-to-English glossary.

## Scientific scope and validation

- `newton_aniversario.py` constructs birthday polynomials; `newton_arte.py` computes research-grade images.
- `certifica_aprisionamento.py` and `certifica_mandelbrot.py` use exact rational arithmetic, Sturm sequences, and interval arithmetic to certify the trapping properties and period-two orbit for Mandelbrot’s date.
- The programs `supplement_en.py`, `auditoria_datas.py`, and `auditoria_figuras.py` generate supplementary tables and audits.
- Some figure scripts need LaTeX and external portrait images; those original images are **not included** in the supplied archive. The precise requirements are documented in the package README.
- Floating-point fractal images and finite-iteration diagnostics are **numerical illustrations**, not proofs. Do not substitute them for the exact certificates.

**Checks performed for this repository import:** all 14 Python source files parse; all 4 supplied JSON files parse; the supplied CSV contains 38 records; both certification scripts returned `certificado: True` in the available Python environment. These tests do not constitute an independent mathematical verification of every result or figure.

The source and data originate from the author's supplied `ESM_2.zip` package (21 November 2026).

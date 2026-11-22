# Scientific code and reproducibility

BirthdayFractal has two separate computational implementations:

| Component | Purpose | Location |
| --- | --- | --- |
| Interactive JavaScript | Fast, approximate rendering of Newton attraction basins in the browser | [app.js](app.js), [worker.js](worker.js) |
| Python research package (Online Resource 2) | Article figures, exact certificates, tables, date audits, and supporting data | [research/ESM_2](research/ESM_2/) |

## Python research package

The author supplied an archive named `ESM_2.zip`, whose contents include 14 Python programs, a LaTeX supplementary document, a 38-date CSV, four JSON records, and an audit text file. This package is being organized under `research/ESM_2/` with its original detailed README, a requirements file, and a reproduction guide.

Key programs:

- `newton_aniversario.py`: maps dates to polynomial coefficients and defines root/window selection.
- `newton_arte.py`: computes research-grade Newton basins and PNG images.
- `certifica_aprisionamento.py`: exact certificate for the trapping intervals and a repelling period-two orbit of the 20/11/1924 polynomial.
- `certifica_mandelbrot.py`: refined interval/Sturm verification for Mandelbrot’s birthday.
- `supplement_en.py`: tables S1–S6 and `supplementary_dates.csv`.
- `auditoria_datas.py` and `auditoria_figuras.py`: article-value and figure-classification audits.
- `esm/ESM_1.tex`: Online Resource 1 manuscript source.

The supplied README documents the complete execution order, dependencies, and a glossary of Portuguese identifiers.

## Scientific validation

The 14 original Python files passed syntax parsing. The supplied CSV contains 38 date records. Both certified trapping/period-two verification scripts returned `certificado: True` in the available environment. Some figures require external portrait files or a LaTeX installation and cannot be fully reproduced with the supplied archive alone.

The JavaScript viewer is an educational visualization with numerical tolerances; it is not a substitute for the exact certificates. See [the publication companion](https://americocunhajr.github.io/BirthdayFractal/).


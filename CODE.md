# Scientific code and reproducibility

BirthdayFractal has two separate computational implementations:

| Component | Purpose | Location |
| --- | --- | --- |
| Interactive JavaScript | Fast, approximate rendering of Newton attraction basins in the browser | [app.js](app.js), [worker.js](worker.js) |
| Python research package (Online Resource 2) | Article figures, exact certificates, tables, date audits, and supporting data | [PythonCode](PythonCode/) |

## Python research package

The author supplied an archive named `ESM_2.zip`, whose contents include 14 Python programs, a LaTeX supplementary document, a 38-date CSV, four JSON records, and an audit text file. The complete author-supplied source is published under [`PythonCode/`](PythonCode/), including its detailed [README](PythonCode/README.md), [requirements](PythonCode/requirements.txt), audit data, and supplementary material.

Key programs:

- `newton_aniversario.py`: maps dates to polynomial coefficients and defines root/window selection.
- `newton_arte.py`: computes research-grade Newton basins and PNG images.
- `certifica_aprisionamento.py`: exact certificate for the trapping intervals and a repelling period-two orbit of the 20/11/1924 polynomial.
- `certifica_mandelbrot.py`: refined interval/Sturm verification for Mandelbrot’s birthday.
- `supplement_en.py`: tables S1–S6 and `supplementary_dates.csv`.
- `auditoria_datas.py` and `auditoria_figuras.py`: article-value and figure-classification audits.
- `esm/ESM_1.tex`: Online Resource 1 manuscript source.

The supplied README documents the complete execution order, dependencies, and a glossary of Portuguese identifiers.

## Reproduce the certified results

Use Python 3.13 (the original environment used Python 3.13.16). From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate   # On Windows: .venv\\Scripts\\activate
python -m pip install -r PythonCode/requirements.txt
cd PythonCode
python certifica_aprisionamento.py
python certifica_mandelbrot.py
```

Both scripts produce JSON files that contain the field `"certificado": true` when the stated exact/interval verification conditions succeed. Other scripts produce figures, numerical audits, and supplementary tables. **For the complete order of execution and output inventory**, consult the [original research README](PythonCode/README.md).

## Scientific validation

The 14 original Python files passed syntax parsing. The supplied CSV contains 38 date records. Both certified trapping/period-two verification scripts returned `certificado: True` in the available environment. Some figures require external portrait files or a LaTeX installation and cannot be fully reproduced with the supplied archive alone.

The JavaScript viewer is an educational visualization with numerical tolerances; it is not a substitute for the exact certificates. See [the publication companion](https://americocunhajr.github.io/BirthdayFractal/).


# BirthdayFractal

The public companion website for *Birthday Polynomials and the Boundaries of Newton's Method* by Americo Cunha Jr.

**Website:** https://americocunhajr.github.io/BirthdayFractal/

## Concept

The eight digits of a date (DDMMYYYY) are read as polynomial coefficients, highest degree first. For each complex initial condition, Newton's method follows the rational map N_p(z)=z-p(z)/p'(z). Color points by the root approached to reveal their attraction basins and the shared Julia-set boundary.

This repository publishes a lightweight, dependency-free **interactive browser edition**. The numerical gallery is computed in JavaScript using floating-point arithmetic, not the scientific Python renderer or the exact certification code from the associated paper.

## Features

- Enter any valid date, including leading-zero days.
- Automatic polynomial construction and numerical root approximations.
- Newton-fractal rendering in a Web Worker.
- Pan by clicking to zoom; zoom controls and reset.
- Download a PNG and copy a shareable date URL.
- Four color palettes and three resolutions.
- Browse the 38 historical and personal dates considered in the paper.
- Responsive dark interface; no server or JavaScript framework required.

## Local preview

Serve the files over HTTP (required for the Web Worker):

```bash
python -m http.server 8000
```

Visit http://localhost:8000/ .

## GitHub Pages

A GitHub Actions workflow in `.github/workflows/pages.yml` deploys the site on changes to `main`. If GitHub Pages is not already enabled, choose **GitHub Actions** as the build/deployment source under **Settings → Pages**.

## Reproducibility

The interactive generator ([`app.js`](app.js), [`worker.js`](worker.js)) runs entirely in the browser. The associated scientific Python package includes the source of the article's figures, data tables, and exact certificates. The full reproducibility material is in [`research/ESM_2/`](research/ESM_2/), with the [author's documented reproduction commands](research/ESM_2/README.md), input data, exact certificates, and Python source. See [Scientific code and reproducibility](CODE.md) for an overview and the distinction between numerical visualizations and proofs.

## Citation

Americo Cunha Jr. “Birthday Polynomials and the Boundaries of Newton's Method.” Manuscript under review at *The Mathematical Intelligencer*, 2026.

```bibtex
@unpublished{cunha2026birthdaypolynomials,
  author = {Cunha Jr, Americo},
  title  = {Birthday Polynomials and the Boundaries of Newton's Method},
  note   = {Manuscript under review at The Mathematical Intelligencer},
  year   = {2026}
}
```

## Author

**Americo Cunha Jr.** is Research Faculty at the National Laboratory for Scientific Computing (LNCC) and Professor of Applied and Computational Mathematics at Rio de Janeiro State University (UERJ), Brazil. See [his academic website](https://americocunha.org/) or the [About the author](https://americocunhajr.github.io/BirthdayFractal/#author) section.

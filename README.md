# BirthdayFractal

The public companion website for *Birthday Polynomials and the Boundaries of Newton's Method* by Americo Cunha Jr.

**Website:** https://americocunha.org/BirthdayFractal/

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

## Discoverability and publication

Five static articles about Newton, Euler, Gauss, Julia, and Mandelbrot complement the interactive gallery. Metadata, an XML sitemap, and structured data are maintained in the repository. See [SEO.md](SEO.md) for Search Console verification and indexing instructions.

**GitHub Pages:** In **Settings → Pages**, select **Deploy from a branch**, using `main` and `/ (root)`. The built-in GitHub **pages build and deployment** workflow publishes the website. The separate `.github/workflows/pages.yml` workflow validates scripts and search assets without deploying a second Pages artifact.

**Local preview:**

```bash
python -m http.server 8000
```

Open http://localhost:8000/ (the fractal renderer uses a Web Worker and requires HTTP).

## Reproducibility

The JavaScript viewer ([`app.js`](app.js), [`worker.js`](worker.js)) runs entirely in the browser. The **scientific Python package** is in [`PythonCode/`](PythonCode/). It includes research-grade figure generators, exact certificates, numerical audits, the 38-date dataset, and supplementary LaTeX source.

Consult the [PythonCode README](PythonCode/README.md) for the complete reproduction sequence and the [scientific-code guide](CODE.md) for an overview. The browser images are approximate illustrations, not substitutes for the mathematical certificates.

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

## License

**BirthdayFractal** is released under the MIT license, consistent with several of the author's software repositories. See [LICENSE](LICENSE) for details.

<img src="logo/mit_license_red.png" width="10%"> 

### Institutional support

<img src="logo/logo_lncc.png" width="25%"> &nbsp; &nbsp; &nbsp; <img src="logo/logo_uerj.png" width="10%"> 

### Funding

<img src="logo/faperj.png" width="25%"> &nbsp; &nbsp; <img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">


## Author

**Americo Cunha Jr.** is Research Faculty at the National Laboratory for Scientific Computing (LNCC) and Professor of Applied and Computational Mathematics at Rio de Janeiro State University (UERJ), Brazil. See [his academic website](https://americocunha.org/) or the [About the author](https://americocunhajr.github.io/BirthdayFractal/#author) section.

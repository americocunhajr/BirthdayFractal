"""All figures of the article except the portraits of Fig. 6, under the journal's file
names, in figures/: fractal pictures as JPEG (1600 x 1600, quality 95, 600 dpi) and the
vector figures Fig2.pdf, Fig5.pdf, Fig12.pdf (Helvetica lettering, through LaTeX).

Usage:
  python make_figures.py               # all figures (about 8 minutes on one core)
  python make_figures.py Fig7f Fig12   # only the figures named
"""
import functools
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image, ImageDraw

from birthday.dates import coefficients
from birthday.newton import (render, standard_window, hues, oklch_to_oklab, oklab_to_srgb,
                             INK_RGB)

OUT = "figures"
N, SS = 1600, 2

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage[scaled]{helvet}"
                           r"\renewcommand{\familydefault}{\sfdefault}"
                           r"\usepackage{sansmath}\sansmath\usepackage{amsmath}"})
INK = np.array(INK_RGB) / 255
BLUE = oklab_to_srgb(oklch_to_oklab(np.array(0.55), np.array(0.10), np.array(250.0)))
CORAL = oklab_to_srgb(oklch_to_oklab(np.array(0.62), np.array(0.14), np.array(30.0)))
BAND = oklab_to_srgb(oklch_to_oklab(np.array(0.93), np.array(0.04), np.array(80.0)))

# fractal pictures of dates
DATE_FIGURES = {
    "Fig1a": "25/12/1642", "Fig1b": "04/01/1643",
    "Fig7a": "15/04/1707", "Fig7b": "16/08/1821", "Fig7c": "29/04/1854",
    "Fig7d": "28/02/1878", "Fig7e": "03/02/1893", "Fig7f": "20/11/1924",
    "Fig8a": "05/07/1687", "Fig8b": "23/09/1846", "Fig8c": "30/11/1889",
    "Fig8d": "29/05/1919", "Fig8e": "29/12/1972", "Fig8f": "16/05/1980",
    "Fig9a": "08/08/1900",
    "Fig10a": "30/08/1984", "Fig10b": "03/07/1986", "Fig10c": "27/04/1952", "Fig10d": "15/08/1950",
}

# Fig. 4: Gauss's date, a pole of N_p (real zero of p') and a preimage of another pole
GAUSS = coefficients("30/04/1777")
POLE = complex(-0.45128, 0.0)
PREIMAGE = complex(-0.6158, 0.0919)


def save_jpeg(img, name):
    img.convert("RGB").save(os.path.join(OUT, name + ".jpg"), quality=95, subsampling=0,
                            dpi=(600, 600))


@functools.lru_cache(maxsize=None)
def gauss_views():
    """The three views of Fig. 4 (without rectangles) and their windows."""
    center, R = standard_window(GAUSS)
    levels = [(center, R), (POLE, R / 6), (PREIMAGE, R / 48)]
    imgs = [render(GAUSS, n=N, ss=SS, center=z0, R=r, Rref=R, color_center=center)[0]
            for z0, r in levels]
    return imgs, levels


def mark(img, window, sub):
    """Draw on img (showing window) the rectangle of the sub-window magnified next."""
    (z0, r), (z1, r1) = window, sub
    draw = ImageDraw.Draw(img)
    s = (N - 1) / (2 * r)
    x0 = (z1.real - r1 - (z0.real - r)) * s
    x1 = (z1.real + r1 - (z0.real - r)) * s
    y0 = ((z0.imag + r) - (z1.imag + r1)) * s
    y1 = ((z0.imag + r) - (z1.imag - r1)) * s
    for k in range(-1, 6):  # thin white frame outside, ink inside
        draw.rectangle([x0 - k, y0 - k, x1 + k, y1 + k],
                       outline=(255, 255, 255) if k in (-1, 5) else INK_RGB)
    return img


def fig4(panel):
    imgs, levels = gauss_views()
    img = imgs[panel].copy()
    if panel < 2:
        mark(img, levels[panel], levels[panel + 1])
    return img


def fig11():
    """z^3 - 2z + 2 around the cycle {0, 1}; colors and tolerance of the standard window."""
    c = [1, 0, -2, 2]
    center, R = standard_window(c)
    return render(c, n=N, ss=SS, center=0.5 + 0j, R=1.0, Rref=R, color_center=center)[0]


FRACTALS = {name: (lambda d=d: render(coefficients(d), n=N, ss=SS)[0])
            for name, d in DATE_FIGURES.items()}
FRACTALS.update({
    "Fig3a": lambda: render([1, 0, -1], n=N, ss=SS)[0],
    "Fig3b": lambda: render([1, 0, 0, -1], n=N, ss=SS)[0],
    "Fig4a": lambda: fig4(0), "Fig4b": lambda: fig4(1), "Fig4c": lambda: fig4(2),
    "Fig9b": lambda: render(coefficients("08/08/1900"), n=N, ss=SS, normalize=False)[0],
    "Fig11": fig11,
})


def fig2():
    """Tangent construction of Newton's method."""
    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    x = np.linspace(0.95, 2.25, 300)
    xk, xk1, xs = 2.0, 1.5, np.sqrt(2)
    ax.axhline(0, color=INK, lw=0.7)
    ax.plot(x, x ** 2 - 2, color=BLUE, lw=1.8)
    xt = np.linspace(1.35, 2.2, 10)
    ax.plot(xt, 2 + 4 * (xt - 2), color=CORAL, lw=1.2)
    ax.plot([xk, xk], [0, 2], ls=(0, (3, 2)), color=INK, lw=0.7)
    ax.plot([xk], [2], "o", ms=4, color=INK)
    for xv, label, off, va in [(xk, r"$x_k$", (0, -7), "top"), (xk1, r"$x_{k+1}$", (4, -7), "top"),
                               (xs, r"$x^*$", (-6, 7), "bottom")]:
        ax.plot([xv], [0], "o", ms=3.5, color=INK if xv != xs else BLUE, zorder=5)
        ax.annotate(label, (xv, 0), xytext=off, textcoords="offset points",
                    ha="center", va=va, color=INK if xv != xs else BLUE)
    ax.text(2.2, 2.85, r"$y=f(x)$", color=BLUE, ha="right")
    ax.text(1.80, 0.88, "tangent", color=CORAL, fontsize=9, ha="left", va="center",
            rotation=41, rotation_mode="anchor")
    ax.set_xlim(0.95, 2.28)
    ax.set_ylim(-1.4, 3.2)
    ax.axis("off")
    fig.savefig(os.path.join(OUT, "Fig2.pdf"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def fig5():
    """From a date to a picture: digits, polynomial, roots, basins (30/04/1777)."""
    center, _ = standard_window(GAUSS)
    roots = np.roots(GAUSS)
    colors = oklab_to_srgb(oklch_to_oklab(np.full(len(roots), 0.66), np.full(len(roots), 0.13),
                                          hues(roots, center)))
    fig = plt.figure(figsize=(6.4, 2.25))
    a1 = fig.add_axes([0.0, 0.0, 0.40, 1.0]); a1.axis("off")
    a1.text(0.5, 0.86, "30/04/1777", ha="center", va="center", fontsize=12, color=INK)
    a1.text(0.5, 0.64, r"$3\;\;0\;\;0\;\;4\;\;1\;\;7\;\;7\;\;7$", ha="center", va="center",
            fontsize=12, color=INK)
    a1.annotate("", xy=(0.5, 0.42), xytext=(0.5, 0.55),
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8))
    a1.text(0.5, 0.30, r"$p(z)=3z^7+4z^4+z^3$", ha="center", va="center", fontsize=10, color=INK)
    a1.text(0.5, 0.16, r"$+\,7z^2+7z+7$", ha="center", va="center", fontsize=10, color=INK)
    a2 = fig.add_axes([0.45, 0.10, 0.25, 0.80])
    t = np.linspace(0, 2 * np.pi, 200)
    a2.plot(np.cos(t), np.sin(t), color=INK, lw=0.4, ls=(0, (2, 2)))
    a2.axhline(0, color=INK, lw=0.4); a2.axvline(0, color=INK, lw=0.4)
    a2.scatter(roots.real, roots.imag, s=34, c=np.clip(colors, 0, 1), edgecolors=INK,
               linewidths=0.6, zorder=3)
    a2.set_xlim(-1.8, 1.8); a2.set_ylim(-1.8, 1.8); a2.set_aspect("equal")
    a2.set_xticks([]); a2.set_yticks([])
    for s in a2.spines.values():
        s.set_visible(False)
    a2.text(1.75, 0.08, r"$\mathrm{Re}\,z$", ha="right", va="bottom", fontsize=8, color=INK)
    a2.text(0.08, 1.75, r"$\mathrm{Im}\,z$", ha="left", va="top", fontsize=8, color=INK)
    a3 = fig.add_axes([0.745, 0.06, 0.255, 0.88]); a3.axis("off")
    a3.imshow(gauss_views()[0][0].resize((900, 900), Image.LANCZOS))
    for xa in (0.405, 0.71):
        fig.add_artist(FancyArrowPatch((xa, 0.5), (xa + 0.035, 0.5), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=10, color=INK, lw=0.8))
    fig.savefig(os.path.join(OUT, "Fig5.pdf"), dpi=600, pad_inches=0.02)
    plt.close(fig)


def fig12():
    """The trapped critical orbit of 20/11/1924: (a) N_p on the real line with K- and K+
    and the first iterates of c; (b) N_p^2 on K-, with the period-two point x*."""
    plt.rcParams.update({"font.size": 8.5})
    p = np.poly1d(coefficients("20/11/1924"))
    dp = p.deriv()
    Np = lambda x: x - p(x) / dp(x)
    alpha, beta, gamma, delta = -0.7626, -0.6240, 0.1580, 0.2012
    pole_roots = np.roots(dp.coeffs)
    poles = np.sort(pole_roots[np.abs(pole_roots.imag) < 1e-12].real)
    crit = np.roots(dp.deriv().coeffs)
    c = float(crit[np.abs(crit.imag) < 1e-12].real[0])
    x_star = -0.7269194427359652  # period-two point in K-, certified by certify_trapping.py

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 3.0), gridspec_kw=dict(width_ratios=[1.35, 1]))
    # (a) N_p on the real line
    xmin, xmax, ymin, ymax = -1.05, 0.32, -1.05, 0.45
    for (a, b) in [(alpha, beta), (gamma, delta)]:
        a1.axvspan(a, b, color=BAND, lw=0, zorder=0)
        a1.axhspan(a, b, color=BAND, lw=0, zorder=0)
    for (x0, x1, y0, y1) in [(alpha, beta, gamma, delta), (gamma, delta, alpha, beta)]:
        a1.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec=INK, lw=0.6, zorder=3))
    a1.plot([xmin, xmax], [xmin, xmax], color=INK, lw=0.5, ls=(0, (2, 2)), zorder=1)
    cuts = [xmin] + list(poles) + [xmax]
    for a, b in zip(cuts[:-1], cuts[1:]):
        x = np.linspace(a + 1e-4, b - 1e-4, 4000)
        y = Np(x)
        y[(y < ymin - 1) | (y > ymax + 1)] = np.nan
        a1.plot(x, y, color=BLUE, lw=1.1, zorder=2)
    for xp in poles:
        a1.axvline(xp, color=INK, lw=0.5, ls=(0, (1, 2)))
    x = c  # cobweb c -> x1 -> x2 -> x3 -> x4
    a1.plot([x, x], [x, Np(x)], color=CORAL, lw=0.8, zorder=4)
    for _ in range(3):
        y = Np(x)
        a1.plot([x, y], [y, y], color=CORAL, lw=0.8, zorder=4)
        a1.plot([y, y], [y, Np(y)], color=CORAL, lw=0.8, zorder=4)
        x = y
    a1.plot([c], [Np(c)], "o", ms=3, color=CORAL, zorder=5)
    a1.plot([c], [c], "o", ms=3, color=INK, zorder=5)
    a1.annotate(r"$c$", (c, c), xytext=(-9, 3), textcoords="offset points", color=INK)
    a1.text((alpha + beta) / 2, ymin + 0.03, r"$K_-$", ha="center", va="bottom", color=INK)
    a1.text((gamma + delta) / 2, ymin + 0.03, r"$K_+$", ha="center", va="bottom", color=INK)
    a1.text(xmin + 0.03, (gamma + delta) / 2, r"$K_+$", ha="left", va="center", color=INK)
    a1.text(xmin + 0.03, (alpha + beta) / 2, r"$K_-$", ha="left", va="center", color=INK)
    for xp in poles:
        a1.text(xp + 0.015, ymax - 0.03, "pole", rotation=90, ha="left", va="top", color=INK,
                fontsize=7.5)
    a1.text(0.27, 0.38, r"$y=x$", ha="right", va="bottom", color=INK, fontsize=7.5)
    a1.text(-0.50, 0.33, r"$y=N_p(x)$", ha="center", va="bottom", color=BLUE)
    a1.set_xlim(xmin, xmax); a1.set_ylim(ymin, ymax)
    a1.set_xlabel(r"$x$"); a1.set_ylabel(r"$N_p(x)$", labelpad=1)
    a1.set_xticks([-1, -0.5, 0]); a1.set_yticks([-1, -0.5, 0])
    a1.set_title(r"(a)", loc="left", fontsize=8.5)
    # (b) N_p^2 on K-
    x = np.linspace(alpha, beta, 3000)
    a2.add_patch(Rectangle((alpha, alpha), beta - alpha, beta - alpha, fill=False, ec=INK, lw=0.6))
    a2.plot([alpha, beta], [alpha, beta], color=INK, lw=0.5, ls=(0, (2, 2)))
    a2.plot(x, Np(Np(x)), color=BLUE, lw=1.1)
    x = c
    a2.plot([x, x], [x, Np(Np(x))], color=CORAL, lw=0.8)
    for _ in range(4):
        y = Np(Np(x))
        a2.plot([x, y], [y, y], color=CORAL, lw=0.8)
        a2.plot([y, y], [y, Np(Np(y))], color=CORAL, lw=0.8)
        x = y
    a2.plot([x_star], [x_star], "o", ms=3.6, mfc="white", mec=INK, mew=0.8, zorder=6)
    a2.annotate(r"$x_*$", (x_star, x_star), xytext=(-12, 3), textcoords="offset points", color=INK)
    a2.plot([c], [c], "o", ms=3, color=INK, zorder=5)
    a2.annotate(r"$c$", (c, c), xytext=(4, -9), textcoords="offset points", color=INK)
    a2.text(-0.650, -0.672, r"$y=x$", ha="left", va="top", color=INK, fontsize=7.5)
    a2.text(-0.712, -0.640, r"$y=N_p^2(x)$", ha="center", va="bottom", color=BLUE)
    a2.set_xlim(alpha, beta); a2.set_ylim(alpha, beta)
    a2.set_xticks([-0.75, -0.70, -0.65]); a2.set_yticks([-0.75, -0.70, -0.65])
    a2.set_xlabel(r"$x\in K_-$"); a2.set_ylabel(r"$N_p^2(x)$", labelpad=1)
    a2.set_aspect("equal")
    a2.set_title(r"(b)", loc="left", fontsize=8.5)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(os.path.join(OUT, "Fig12.pdf"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    plt.rcParams.update({"font.size": 9})


VECTOR = {"Fig2": fig2, "Fig5": fig5, "Fig12": fig12}


def main(names):
    os.makedirs(OUT, exist_ok=True)
    unknown = set(names) - set(FRACTALS) - set(VECTOR)
    if unknown:
        sys.exit(f"unknown figure(s): {sorted(unknown)}")
    for name in sorted(FRACTALS, key=lambda s: (int(s[3:].rstrip("abcdefgh")), s)):
        if not names or name in names:
            save_jpeg(FRACTALS[name](), name)
            print(name + ".jpg", flush=True)
    for name, make in VECTOR.items():
        if not names or name in names:
            make()
            print(name + ".pdf", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])

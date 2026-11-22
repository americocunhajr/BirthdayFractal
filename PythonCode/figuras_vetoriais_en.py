"""
English versions of the two vector figures (Fig. 2: tangent construction; Fig. 5:
from a date to a picture), with Helvetica lettering as required by the journal.
Same geometry and colours as figuras_extras.py. Output: submission/Fig2.pdf, Fig5.pdf.
Requires LaTeX with helvet, sansmath and cm-super (used by Matplotlib's usetex).
"""
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from newton_arte import janela_padrao, matizes, oklch_para_oklab, oklab_para_srgb, ESTILOS
from newton_aniversario import coeficientes

OUT = "submission"
import os
os.makedirs(OUT, exist_ok=True)
TINTA = (43, 45, 66)
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.linewidth": 0.6,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage[scaled]{helvet}"
                           r"\renewcommand{\familydefault}{\sfdefault}"
                           r"\usepackage{sansmath}\sansmath\usepackage{amsmath}"})

L, C, h0, *_ = ESTILOS["aquarela"]
azul = oklab_para_srgb(oklch_para_oklab(np.array(0.55), np.array(0.10), np.array(250.0)))
coral = oklab_para_srgb(oklch_para_oklab(np.array(0.62), np.array(0.14), np.array(30.0)))
tinta = np.array(TINTA) / 255

# ------------------------------------------------------------- Fig. 2: tangent
fig, ax = plt.subplots(figsize=(3.6, 2.7))
x = np.linspace(0.95, 2.25, 300)
xk, xk1, xs = 2.0, 1.5, np.sqrt(2)
ax.axhline(0, color=tinta, lw=0.7)
ax.plot(x, x ** 2 - 2, color=azul, lw=1.8)
xt = np.linspace(1.35, 2.2, 10)
ax.plot(xt, 2 + 4 * (xt - 2), color=coral, lw=1.2)
ax.plot([xk, xk], [0, 2], ls=(0, (3, 2)), color=tinta, lw=0.7)
ax.plot([xk], [2], "o", ms=4, color=tinta)
for xv, rot, off, va in [(xk, r"$x_k$", (0, -7), "top"), (xk1, r"$x_{k+1}$", (4, -7), "top"),
                         (xs, r"$x^*$", (-6, 7), "bottom")]:
    ax.plot([xv], [0], "o", ms=3.5, color=tinta if xv != xs else azul, zorder=5)
    ax.annotate(rot, (xv, 0), xytext=off, textcoords="offset points",
                ha="center", va=va, color=tinta if xv != xs else azul)
ax.text(2.2, 2.85, r"$y=f(x)$", color=azul, ha="right")
ax.text(1.80, 0.88, "tangent", color=coral, fontsize=9, ha="left", va="center",
        rotation=41, rotation_mode="anchor")
ax.set_xlim(0.95, 2.28)
ax.set_ylim(-1.4, 3.2)
ax.axis("off")
fig.savefig(f"{OUT}/Fig2.pdf", bbox_inches="tight", pad_inches=0.02)
plt.close(fig)

# ------------------------------------------------------------- Fig. 5: date -> picture
c = np.array(coeficientes("30/04/1777"), float)
cen, R = janela_padrao(c)
raizes = np.roots(c)
hs = matizes(raizes, cen, h0)
cores = oklab_para_srgb(oklch_para_oklab(np.full(len(raizes), 0.66), np.full(len(raizes), 0.13), hs))
fig = plt.figure(figsize=(6.4, 2.25))
a1 = fig.add_axes([0.0, 0.0, 0.40, 1.0]); a1.axis("off")
a1.text(0.5, 0.86, "30/04/1777", ha="center", va="center", fontsize=12, color=tinta)
a1.text(0.5, 0.64, r"$3\;\;0\;\;0\;\;4\;\;1\;\;7\;\;7\;\;7$", ha="center", va="center",
        fontsize=12, color=tinta)
a1.annotate("", xy=(0.5, 0.42), xytext=(0.5, 0.55),
            arrowprops=dict(arrowstyle="-|>", color=tinta, lw=0.8))
a1.text(0.5, 0.30, r"$p(z)=3z^7+4z^4+z^3$", ha="center", va="center", fontsize=10, color=tinta)
a1.text(0.5, 0.16, r"$+\,7z^2+7z+7$", ha="center", va="center", fontsize=10, color=tinta)
a2 = fig.add_axes([0.45, 0.10, 0.25, 0.80])
t = np.linspace(0, 2 * np.pi, 200)
a2.plot(np.cos(t), np.sin(t), color=tinta, lw=0.4, ls=(0, (2, 2)))
a2.axhline(0, color=tinta, lw=0.4); a2.axvline(0, color=tinta, lw=0.4)
a2.scatter(raizes.real, raizes.imag, s=34, c=np.clip(cores, 0, 1), edgecolors=tinta,
           linewidths=0.6, zorder=3)
a2.set_xlim(-1.8, 1.8); a2.set_ylim(-1.8, 1.8); a2.set_aspect("equal")
a2.set_xticks([]); a2.set_yticks([])
for s in a2.spines.values():
    s.set_visible(False)
a2.text(1.75, 0.08, r"$\mathrm{Re}\,z$", ha="right", va="bottom", fontsize=8, color=tinta)
a2.text(0.08, 1.75, r"$\mathrm{Im}\,z$", ha="left", va="top", fontsize=8, color=tinta)
a3 = fig.add_axes([0.745, 0.06, 0.255, 0.88]); a3.axis("off")
a3.imshow(Image.open("figs/newton_30-04-1777.png").resize((900, 900), Image.LANCZOS))
for xa in (0.405, 0.71):
    fig.add_artist(matplotlib.patches.FancyArrowPatch((xa, 0.5), (xa + 0.035, 0.5),
                   transform=fig.transFigure, arrowstyle="-|>", mutation_scale=10,
                   color=tinta, lw=0.8))
fig.savefig(f"{OUT}/Fig5.pdf", dpi=600, pad_inches=0.02)
plt.close(fig)
print("Fig2.pdf, Fig5.pdf")

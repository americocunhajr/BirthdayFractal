"""
Fig. 12 of the article: the trapped critical orbit for 20/11/1924.
(a) Graph of the real Newton map N, the trapping intervals K- and K+ (shaded on both
    axes, with the boxes K- x K+ and K+ x K-), the real poles, and the first iterates of
    the free critical point c (cobweb).
(b) N o N on K-, a unimodal map of K- into itself with minimum at c.
Helvetica lettering, as required by the journal. Output: submission/Fig12.pdf
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from newton_arte import oklch_para_oklab, oklab_para_srgb

OUT = "submission"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8.5, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage[scaled]{helvet}"
                           r"\renewcommand{\familydefault}{\sfdefault}"
                           r"\usepackage{sansmath}\sansmath\usepackage{amsmath}"})

tinta = np.array((43, 45, 66)) / 255
azul = oklab_para_srgb(oklch_para_oklab(np.array(0.55), np.array(0.10), np.array(250.0)))
coral = oklab_para_srgb(oklch_para_oklab(np.array(0.62), np.array(0.14), np.array(30.0)))
faixa = oklab_para_srgb(oklch_para_oklab(np.array(0.93), np.array(0.04), np.array(80.0)))

c = np.array([2, 0, 1, 1, 1, 9, 2, 4], float)
p, dp = np.poly1d(c), np.poly1d(c).deriv()
N = lambda x: x - p(x) / dp(x)
alpha, beta, gamma, delta = -0.7626, -0.6240, 0.1580, 0.2012
polos = np.sort(np.roots(dp.coeffs)[np.abs(np.roots(dp.coeffs).imag) < 1e-12].real)
cc = np.roots(dp.deriv().coeffs)
cc = float(cc[np.abs(cc.imag) < 1e-12].real[0])

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.4, 3.0), gridspec_kw=dict(width_ratios=[1.35, 1]))

# ---------------------------------------------------------------- (a) N on the real line
xmin, xmax, ymin, ymax = -1.05, 0.32, -1.05, 0.45
for (a, b) in [(alpha, beta), (gamma, delta)]:
    a1.axvspan(a, b, color=faixa, lw=0, zorder=0)
    a1.axhspan(a, b, color=faixa, lw=0, zorder=0)
for (x0, x1_, y0, y1_) in [(alpha, beta, gamma, delta), (gamma, delta, alpha, beta)]:
    a1.add_patch(Rectangle((x0, y0), x1_ - x0, y1_ - y0, fill=False, ec=tinta, lw=0.6, zorder=3))
a1.plot([xmin, xmax], [xmin, xmax], color=tinta, lw=0.5, ls=(0, (2, 2)), zorder=1)
cortes = [xmin] + list(polos) + [xmax]
for a, b in zip(cortes[:-1], cortes[1:]):
    x = np.linspace(a + 1e-4, b - 1e-4, 4000)
    y = N(x)
    y[(y < ymin - 1) | (y > ymax + 1)] = np.nan
    a1.plot(x, y, color=azul, lw=1.1, zorder=2)
for xp in polos:
    a1.axvline(xp, color=tinta, lw=0.5, ls=(0, (1, 2)))
# cobweb from c: c -> x1 -> x2 -> x3 -> x4
x = cc
a1.plot([x, x], [x, N(x)], color=coral, lw=0.8, zorder=4)
for k in range(3):
    y = N(x)
    a1.plot([x, y], [y, y], color=coral, lw=0.8, zorder=4)
    a1.plot([y, y], [y, N(y)], color=coral, lw=0.8, zorder=4)
    x = y
a1.plot([cc], [N(cc)], "o", ms=3, color=coral, zorder=5)
a1.plot([cc], [cc], "o", ms=3, color=tinta, zorder=5)
a1.annotate(r"$c$", (cc, cc), xytext=(-9, 3), textcoords="offset points", color=tinta)
a1.text((alpha + beta) / 2, ymin + 0.03, r"$K_-$", ha="center", va="bottom", color=tinta)
a1.text((gamma + delta) / 2, ymin + 0.03, r"$K_+$", ha="center", va="bottom", color=tinta)
a1.text(xmin + 0.03, (gamma + delta) / 2, r"$K_+$", ha="left", va="center", color=tinta)
a1.text(xmin + 0.03, (alpha + beta) / 2, r"$K_-$", ha="left", va="center", color=tinta)
a1.text(polos[0] + 0.015, ymax - 0.03, "pole", rotation=90, ha="left", va="top", color=tinta,
        fontsize=7.5)
a1.text(polos[1] + 0.015, ymax - 0.03, "pole", rotation=90, ha="left", va="top", color=tinta,
        fontsize=7.5)
a1.text(0.27, 0.38, r"$y=x$", ha="right", va="bottom", color=tinta, fontsize=7.5)
a1.text(-0.50, 0.33, r"$y=N_p(x)$", ha="center", va="bottom", color=azul)
a1.set_xlim(xmin, xmax); a1.set_ylim(ymin, ymax)
a1.set_xlabel(r"$x$"); a1.set_ylabel(r"$N_p(x)$", labelpad=1)
a1.set_xticks([-1, -0.5, 0]); a1.set_yticks([-1, -0.5, 0])
a1.set_title(r"(a)", loc="left", fontsize=8.5)

# ---------------------------------------------------------------- (b) N o N on K-
x = np.linspace(alpha, beta, 3000)
y = N(N(x))
a2.add_patch(Rectangle((alpha, alpha), beta - alpha, beta - alpha, fill=False, ec=tinta, lw=0.6))
a2.plot([alpha, beta], [alpha, beta], color=tinta, lw=0.5, ls=(0, (2, 2)))
a2.plot(x, y, color=azul, lw=1.1)
x = cc
a2.plot([x, x], [x, N(N(x))], color=coral, lw=0.8)
for k in range(4):
    y = N(N(x))
    a2.plot([x, y], [y, y], color=coral, lw=0.8)
    a2.plot([y, y], [y, N(N(y))], color=coral, lw=0.8)
    x = y
xe = -0.7269194427359652   # fixed point of N o N in K- (2-cycle), certified in certifica_aprisionamento.py
a2.plot([xe], [xe], "o", ms=3.6, mfc="white", mec=tinta, mew=0.8, zorder=6)
a2.annotate(r"$x_*$", (xe, xe), xytext=(-12, 3), textcoords="offset points", color=tinta)
a2.plot([cc], [cc], "o", ms=3, color=tinta, zorder=5)
a2.annotate(r"$c$", (cc, cc), xytext=(4, -9), textcoords="offset points", color=tinta)
a2.text(-0.650, -0.672, r"$y=x$", ha="left", va="top", color=tinta, fontsize=7.5)
a2.text(-0.712, -0.640, r"$y=N_p^2(x)$", ha="center", va="bottom", color=azul)
a2.set_xlim(alpha, beta); a2.set_ylim(alpha, beta)
a2.set_xticks([-0.75, -0.70, -0.65]); a2.set_yticks([-0.75, -0.70, -0.65])
a2.set_xlabel(r"$x\in K_-$"); a2.set_ylabel(r"$N_p^2(x)$", labelpad=1)
a2.set_aspect("equal")
a2.set_title(r"(b)", loc="left", fontsize=8.5)

fig.tight_layout(w_pad=1.5)
fig.savefig(f"{OUT}/Fig12.pdf", bbox_inches="tight", pad_inches=0.02)
print("Fig12.pdf; c =", cc, "poles =", polos)

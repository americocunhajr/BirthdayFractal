"""
Figuras do artigo que nao sao datas: Cayley (z^2-1, z^3-1), falha (z^3-2z+2),
ampliacoes do fractal de Gauss, tangente (vetorial) e pipeline (vetorial + imagem).
Saida em figs/. Dependencias: numpy, Pillow, matplotlib, newton_arte.py.
"""
import numpy as np
from PIL import Image, ImageDraw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

from newton_arte import renderizar_coef, janela_padrao, matizes, oklch_para_oklab, \
    oklab_para_srgb, ESTILOS
from newton_aniversario import coeficientes

N, SS = 1600, 2
TINTA = (43, 45, 66)  # #2b2d42, a mesma tinta do estilo aquarela

# fontes: Latin Modern, como no texto
for f in font_manager.findSystemFonts():
    if "lmroman" in f.lower():
        font_manager.fontManager.addfont(f)
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 0.6,
                     "text.usetex": True,
                     "text.latex.preamble": r"\usepackage[T1]{fontenc}\usepackage{lmodern}\usepackage{amsmath}"})


def salvar(img, nome):
    img.save(f"figs/{nome}.png", optimize=True)
    print("salvo", nome)


# ------------------------------------------------------------- Cayley e falha
for nome, c in [("cayley_z2-1", [1, 0, -1]), ("cayley_z3-1", [1, 0, 0, -1])]:
    img, nc, _, _ = renderizar_coef(c, n=N, ss=SS)
    print(nome, "sem convergencia:", round(nc, 4))
    salvar(img, nome)

# falha: janela em torno do ciclo {0, 1}; cores e tolerancia da janela padrao
cf = [1, 0, -2, 2]
cenf, Rf = janela_padrao(cf)
img, nc, _, _ = renderizar_coef(cf, n=N, ss=SS, centro=0.5 + 0j, R=1.0, Rref=Rf, centro_cor=cenf)
print("falha_z3-2z+2 sem convergencia:", round(nc, 4))
salvar(img, "falha_z3-2z+2")

# ------------------------------------------------------------- ampliacoes (Gauss)
c = np.array(coeficientes("30/04/1777"), float)
cen, R = janela_padrao(c)
polo = complex(-0.45128, 0.0)            # zero real de p' (polo de N_p)
alvo3 = complex(-0.6158, 0.0919)          # pre-imagem de segunda ordem de um polo
niveis = [(cen, R), (polo, R / 6), (alvo3, R / 48)]
imgs = []
for (z0, r) in niveis:
    img, nc, _, _ = renderizar_coef(c, n=N, ss=SS, centro=z0, R=r, Rref=R, centro_cor=cen)
    imgs.append(img)


def retangulo(img, janela, sub):
    (z0, r), (z1, r1) = janela, sub
    dr = ImageDraw.Draw(img)
    esc = (N - 1) / (2 * r)
    x0 = (z1.real - r1 - (z0.real - r)) * esc
    x1 = (z1.real + r1 - (z0.real - r)) * esc
    y0 = ((z0.imag + r) - (z1.imag + r1)) * esc
    y1 = ((z0.imag + r) - (z1.imag - r1)) * esc
    for k in range(-1, 6):  # moldura branca fina por fora, tinta por dentro
        cor = (255, 255, 255) if k in (-1, 5) else TINTA
        dr.rectangle([x0 - k, y0 - k, x1 + k, y1 + k], outline=cor)


retangulo(imgs[0], niveis[0], niveis[1])
retangulo(imgs[1], niveis[1], niveis[2])
for i, img in enumerate(imgs, 1):
    salvar(img, f"zoom_{i}")

# ------------------------------------------------------------- tangente
L, C, h0, *_ = ESTILOS["aquarela"]
azul = oklab_para_srgb(oklch_para_oklab(np.array(0.55), np.array(0.10), np.array(250.0)))
coral = oklab_para_srgb(oklch_para_oklab(np.array(0.62), np.array(0.14), np.array(30.0)))
tinta = np.array(TINTA) / 255

fig, ax = plt.subplots(figsize=(3.6, 2.7))
x = np.linspace(0.95, 2.25, 300)
f = x ** 2 - 2
xk, xk1, xs = 2.0, 1.5, np.sqrt(2)
ax.axhline(0, color=tinta, lw=0.7)
ax.plot(x, f, color=azul, lw=1.8)
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
ax.text(1.80, 0.88, "tangente", color=coral, fontsize=9, ha="left", va="center", rotation=41, rotation_mode="anchor")
ax.set_xlim(0.95, 2.28)
ax.set_ylim(-1.4, 3.2)
ax.axis("off")
fig.savefig("figs/tangente.pdf", bbox_inches="tight", pad_inches=0.02)
plt.close(fig)
print("salvo tangente.pdf")

# ------------------------------------------------------------- pipeline
raizes = np.roots(c)
hs = matizes(raizes, cen, h0)
cores = oklab_para_srgb(oklch_para_oklab(np.full(len(raizes), 0.66), np.full(len(raizes), 0.13), hs))
fig = plt.figure(figsize=(6.4, 2.25))
# painel 1: data e polinomio
a1 = fig.add_axes([0.0, 0.0, 0.40, 1.0]); a1.axis("off")
a1.text(0.5, 0.86, "30/04/1777", ha="center", va="center", fontsize=15, color=tinta)
a1.text(0.5, 0.64, r"$3\;\;0\;\;0\;\;4\;\;1\;\;7\;\;7\;\;7$", ha="center", va="center",
        fontsize=14, color=tinta)
a1.annotate("", xy=(0.5, 0.42), xytext=(0.5, 0.55),
            arrowprops=dict(arrowstyle="-|>", color=tinta, lw=0.8))
a1.text(0.5, 0.30, r"$p(z)=3z^7+4z^4+z^3$", ha="center", va="center", fontsize=11, color=tinta)
a1.text(0.5, 0.16, r"$+\,7z^2+7z+7$", ha="center", va="center", fontsize=11, color=tinta)
# painel 2: raizes
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
# painel 3: bacias
a3 = fig.add_axes([0.745, 0.06, 0.255, 0.88]); a3.axis("off")
a3.imshow(Image.open("figs/newton_30-04-1777.png").resize((900, 900), Image.LANCZOS))
for xa in (0.405, 0.71):
    fig.add_artist(matplotlib.patches.FancyArrowPatch((xa, 0.5), (xa + 0.035, 0.5),
                   transform=fig.transFigure, arrowstyle="-|>", mutation_scale=10,
                   color=tinta, lw=0.8))
fig.savefig("figs/pipeline.pdf", dpi=300, pad_inches=0.02)
plt.close(fig)
print("salvo pipeline.pdf")

"""
Fractal de Newton a partir de uma data de aniversario.

Codificacao: os digitos de DDMMAAAA (ou AAAAMMDD) sao os coeficientes de
p(z) = c_0 z^n + c_1 z^{n-1} + ... + c_n. Cada pixel z_0 do plano complexo e
iterado por z <- z - a p(z)/p'(z); a cor indica a raiz atingida (bacia de
atracao) e o brilho, o numero de iteracoes ate a convergencia.

Dependencias: numpy, matplotlib.
Uso: python newton_aniversario.py 25/12/1990 [--ordem AMD] [--a 1.0] [--n 1600]
"""
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb


def coeficientes(data, ordem="DMA"):
    d, m, a = data.split("/")
    s = {"DMA": d + m + a, "AMD": a + m + d}[ordem]
    c = np.trim_zeros(np.array([int(ch) for ch in s], float), "f")  # remove zeros a esquerda
    if c.size < 3:
        raise ValueError("polinomio de grau < 2")
    return c


def raizes_da_janela(raizes, kappa=3.0):
    """Raizes usadas para enquadrar: exclui raizes isoladas, a mais de kappa vezes a
    distancia mediana ao centro robusto (mediana das partes real e imaginaria),
    calculado sobre as raizes distintas."""
    distintas = []
    for z in raizes:
        if all(abs(z - w) > 1e-6 for w in distintas):
            distintas.append(z)
    distintas = np.array(distintas)
    c = np.median(distintas.real) + 1j * np.median(distintas.imag)
    d = np.abs(distintas - c)
    return distintas[d <= kappa * np.median(d)]


def newton(c, a=1.0, n=1600, maxit=100, tol=1e-10, margem=1.3, foco=False):
    raizes = np.roots(c)
    dc = np.polyder(c)

    # janela centrada no baricentro das raizes (com foco: sem as raizes isoladas)
    rj = raizes_da_janela(raizes) if foco else raizes
    centro = rj.mean()
    R = margem * np.max(np.abs(rj - centro))
    x = np.linspace(centro.real - R, centro.real + R, n)
    y = np.linspace(centro.imag - R, centro.imag + R, n)
    Z = (x[None, :] + 1j * y[:, None]).ravel()

    it = np.full(Z.size, maxit, dtype=float)
    ativo = np.arange(Z.size)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for k in range(maxit):
            za = Z[ativo]
            passo = a * np.polyval(c, za) / np.polyval(dc, za)
            Z[ativo] = za - passo
            fim = ~np.isfinite(passo) | (np.abs(passo) < tol * R)
            it[ativo[fim]] = k
            ativo = ativo[~fim]
            if ativo.size == 0:
                break

    # classificacao: raiz mais proxima, desde que suficientemente proxima
    dist = np.abs(Z[:, None] - raizes[None, :])
    bacia = np.argmin(dist, axis=1)
    bacia[~np.isfinite(Z) | (dist.min(axis=1) > 1e-4 * R)] = -1

    return bacia.reshape(n, n), it.reshape(n, n), raizes, (x[0], x[-1], y[0], y[-1])


def colorir(bacia, it, nraizes):
    h = np.where(bacia >= 0, bacia / nraizes, 0.0)
    s = np.full(bacia.shape, 0.75)
    v = 0.15 + 0.85 * np.exp(-it / 15.0)  # mais escuro onde a convergencia e lenta
    v[bacia < 0] = 0.0                    # sem convergencia: preto
    return hsv_to_rgb(np.dstack([h, s, v]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("data", help="DD/MM/AAAA")
    ap.add_argument("--ordem", default="DMA", choices=["DMA", "AMD"])
    ap.add_argument("--a", type=complex, default=1.0, help="relaxacao (a=1: Newton classico)")
    ap.add_argument("--n", type=int, default=1600, help="resolucao em pixels")
    ap.add_argument("--maxit", type=int, default=100)
    ap.add_argument("--foco", action="store_true",
                    help="enquadrar sem as raizes isoladas (muito distantes das demais)")
    ap.add_argument("--saida", default=".", help="pasta de saida")
    args = ap.parse_args()

    c = coeficientes(args.data, args.ordem)
    bacia, it, raizes, ext = newton(c, a=args.a, n=args.n, maxit=args.maxit, foco=args.foco)
    img = colorir(bacia, it, len(raizes))

    fig, ax = plt.subplots(figsize=(8, 8), dpi=200)
    ax.imshow(img, origin="lower", extent=ext, interpolation="lanczos")
    ax.plot(raizes.real, raizes.imag, "o", ms=3, mfc="white", mec="black", mew=0.5)
    ax.set_xlim(ext[0], ext[1])  # raizes fora da janela (--foco) nao a ampliam
    ax.set_ylim(ext[2], ext[3])
    ax.set_axis_off()
    import os
    nome = os.path.join(args.saida, "newton_" + args.data.replace("/", "-") + ".png")
    fig.savefig(nome, bbox_inches="tight", pad_inches=0)
    print("p(z) coef.:", c.astype(int))
    print("raizes:", np.round(raizes, 4))
    print("fracao sem convergencia:", np.mean(bacia < 0))
    print("salvo em", nome)

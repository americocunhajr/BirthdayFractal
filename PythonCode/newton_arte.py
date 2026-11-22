"""
Fractais de Newton de datas, com esquemas de cor para impressao.

Mesma matematica de newton_aniversario.py (digitos de DD/MM/AAAA como coeficientes,
iteracao de Newton, janela com --foco). Diferencas:
  - contagem de iteracoes suavizada (sem faixas): nu = k - log2(log|passo| / log tol);
  - brilho normalizado por bacia (o "piso" de iteracoes de cada bacia e subtraido),
    para que bacias de raizes multiplas, de convergencia linear, nao saiam pretas;
  - cores definidas em OKLCH (luminosidade perceptual uniforme) e misturadas em OKLab;
  - superamostragem (antisserrilhamento) e escrita direta do PNG, pixel a pixel.

Dependencias: numpy, Pillow.
Uso: python newton_arte.py 30/04/1777 --estilo joias [--n 1600] [--ss 2] [--sem-raizes]
"""
import argparse
import os

import numpy as np
from PIL import Image, ImageDraw

from newton_aniversario import coeficientes, raizes_da_janela


# ---------------------------------------------------------------- calculo
def janela_padrao(c, margem=1.3, foco=True):
    """Centro e meia largura da janela padrao (com a regra das raizes isoladas)."""
    raizes = np.roots(c)
    rj = raizes_da_janela(raizes) if foco else raizes
    centro = rj.mean()
    return centro, margem * np.max(np.abs(rj - centro))


def calcular(c, n=800, ss=2, maxit=300, tol=1e-10, margem=1.3, foco=True,
             centro=None, R=None, Rref=None, com_status=False):
    """centro, R: janela (padrao: janela_padrao). Rref: escala usada na tolerancia
    e na classificacao (padrao: R); em ampliacoes, passar o R da vista completa.
    com_status=True acrescenta ao retorno um dicionario com o estado de cada ponto
    (0 classificado; 1 passo ou iterado nao finito; 2 nao classificado apos maxit
    iteracoes; 3 parou longe de toda raiz) e o ultimo iterado. Nao altera o calculo."""
    raizes = np.roots(c)
    dc = np.polyder(c)
    if centro is None or R is None:
        centro, R = janela_padrao(c, margem, foco)
    Rref = R if Rref is None else Rref
    N = n * ss
    x = np.linspace(centro.real - R, centro.real + R, N)
    y = np.linspace(centro.imag + R, centro.imag - R, N)  # linha 0 = topo
    Z = (x[None, :] + 1j * y[:, None]).ravel()

    kfim = np.full(Z.size, np.nan)    # iteracao de parada
    lpasso = np.full(Z.size, np.nan)  # log do ultimo passo
    ativo = np.arange(Z.size)
    tolR = tol * Rref
    ltol = np.log(tolR)
    with np.errstate(all="ignore"):
        for k in range(maxit):
            za = Z[ativo]
            passo = np.polyval(c, za) / np.polyval(dc, za)
            Z[ativo] = za - passo
            ap = np.abs(passo)
            fim = ~np.isfinite(passo) | (ap < tolR)
            kfim[ativo[fim]] = k
            lpasso[ativo[fim]] = np.log(np.maximum(ap[fim], 1e-300))
            ativo = ativo[~fim]
            if ativo.size == 0:
                break

    # classificacao pela raiz mais proxima (sem matriz N^2 x d)
    melhor = np.full(Z.size, np.inf)
    bacia = np.full(Z.size, -1)
    with np.errstate(invalid="ignore"):
        for j, r in enumerate(raizes):
            d = np.abs(Z - r)
            m = d < melhor
            melhor[m] = d[m]
            bacia[m] = j
    bacia[~np.isfinite(Z) | ~(melhor < 1e-4 * Rref) | np.isnan(kfim)] = -1

    # contagem suavizada: convergencia quadratica (raiz simples) ou linear com
    # taxa (m-1)/m (raiz de multiplicidade m)
    mult = np.array([np.sum(np.abs(raizes - r) < 1e-6) for r in raizes])
    m = np.where(bacia >= 0, mult[np.maximum(bacia, 0)], 1)
    with np.errstate(all="ignore"):
        frac_q = np.log2(np.clip(lpasso / ltol, 1.0, 2.0 ** 1.5))
        rho = np.where(m > 1, (m - 1) / m, 0.5)
        frac_l = np.clip((ltol - lpasso) / -np.log(rho), 0.0, 1.0)
    nu = kfim - np.where(m > 1, frac_l, frac_q)
    if com_status:
        status = np.zeros(Z.size, dtype=np.int8)
        parou = ~np.isnan(kfim)
        with np.errstate(invalid="ignore"):
            status[parou & ~np.isfinite(Z)] = 1
            status[~parou] = 2
            status[parou & np.isfinite(Z) & ~(melhor < 1e-4 * Rref)] = 3
        info = {"status": status.reshape(N, N), "Zfinal": Z.reshape(N, N)}
        return bacia.reshape(N, N), nu.reshape(N, N), raizes, centro, R, info
    return bacia.reshape(N, N), nu.reshape(N, N), raizes, centro, R


# ---------------------------------------------------------------- cor
def oklch_para_oklab(L, C, h):
    h = np.deg2rad(h)
    return np.stack([L, C * np.cos(h), C * np.sin(h)], axis=-1)


def oklab_para_srgb(lab):
    L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]
    l_ = L + 0.3963377774 * a + 0.2158037573 * b
    m_ = L - 0.1055613458 * a - 0.0638541728 * b
    s_ = L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    rgb = np.stack([
        4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s], axis=-1)
    rgb = np.clip(rgb, 0, 1)
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)


def hex_para_oklab(hx):
    v = np.array([int(hx[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    lin = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
    l = 0.4122214708 * lin[0] + 0.5363325363 * lin[1] + 0.0514459929 * lin[2]
    m = 0.2119034982 * lin[0] + 0.6806995451 * lin[1] + 0.1073969566 * lin[2]
    s = 0.0883024619 * lin[0] + 0.2817188376 * lin[1] + 0.6299787005 * lin[2]
    l_, m_, s_ = np.cbrt([l, m, s])
    return np.array([
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_])


def matizes(raizes, centro, h0, espelho=False):
    """Matiz de cada raiz (por indice em `raizes`). Raizes distintas ordenadas pelo
    angulo em torno do centro e espalhadas simetricamente em torno de h0, de modo que
    pares conjugados recebem matizes simetricas (h0 +- delta). Com espelho=True,
    pares conjugados recebem a mesma matiz."""
    dist = []
    for z in raizes:
        if all(abs(z - w) > 1e-6 for w in dist):
            dist.append(z)
    dist = np.array(dist)
    ang = np.angle(dist - centro)
    # chave simetrica: |angulo| e o sinal desempatam; conjugados ficam em posicoes espelhadas
    ordem = np.argsort(ang)
    d = len(dist)
    h = np.empty(d)
    if espelho:
        chave = np.round(np.abs(ang), 6)
        grupos = sorted(set(chave))
        for g, cg in enumerate(grupos):
            h[np.isclose(chave, cg)] = h0 + 300 * (g + 0.5) / len(grupos) - 150
    else:
        for pos, i in enumerate(ordem):
            h[i] = h0 + 360 * (pos + 0.5) / d - 180
    # mapear cada raiz (com repeticoes) para a matiz da raiz distinta correspondente
    return np.array([h[np.argmin(np.abs(dist - z))] for z in raizes])


def suavizar(t):
    return t * t * (3 - 2 * t)


ESTILOS = {
    #            L     C     h0   tinta (borda)  centro claro  modo
    "joias":    (0.70, 0.13, 250, "#151726", 0.25, "claro"),
    "aquarela": (0.88, 0.075, 220, "#2b2d42", 0.40, "claro"),
    "noturno":  (0.78, 0.16, 250, "#07080f", 0.00, "brilho"),
    "espelho":  (0.70, 0.13, 30, "#151726", 0.25, "claro"),
}


def colorir(bacia, nu, raizes, centro, estilo="joias", escala=6.0,
            normalizar=True, pisos=None):
    """centro: ponto de referencia para a ordem angular das matizes.
    normalizar: subtrai o piso de iteracoes de cada bacia (2o percentil), para que
    bacias de convergencia lenta nao saiam escuras. pisos: dicionario {bacia: piso}
    calculado em outra imagem (ampliacoes usam os pisos da vista completa).
    Devolve (rgb, pisos)."""
    L, C, h0, tinta, clarear, modo = ESTILOS[estilo]
    hs = matizes(raizes, centro, h0, espelho=(estilo == "espelho"))
    base = oklch_para_oklab(np.full(len(raizes), L), np.full(len(raizes), C), hs)
    ink = hex_para_oklab(tinta)
    branco = np.array([1.0, 0.0, 0.0])

    t = np.zeros(bacia.shape)
    usados = {}
    validos = bacia >= 0
    piso_global = np.percentile(nu[validos], 2) if validos.any() else 0.0
    for j in np.unique(bacia[validos]):
        m = bacia == j
        if not normalizar:
            piso = piso_global
        elif pisos is not None and j in pisos:
            piso = pisos[j]
        else:
            piso = np.percentile(nu[m], 2)
        usados[int(j)] = piso
        t[m] = 1 - np.exp(-np.maximum(nu[m] - piso, 0) / escala)
    t[~validos] = 1.0

    b = base[np.maximum(bacia, 0)]
    if modo == "claro":
        # perto da raiz: tom claro; meio: cor base; perto da fronteira: tinta
        claro = b + clarear * (branco - b)
        u = suavizar(np.clip(t / 0.35, 0, 1))[..., None]
        v = suavizar(np.clip((t - 0.35) / 0.65, 0, 1))[..., None] ** 1.3
        lab = claro + u * (b - claro)
        lab = lab + v * (ink - lab)
    else:
        # interior escuro, fronteira luminosa na cor da bacia
        u = (t ** 1.6)[..., None]
        lab = ink + u * (b - ink)
        w = suavizar(np.clip((t - 0.85) / 0.15, 0, 1))[..., None] * 0.35
        lab = lab + w * (branco - lab)
    lab[~validos] = ink
    return oklab_para_srgb(lab), usados


# ---------------------------------------------------------------- imagem
def renderizar_coef(c, estilo="aquarela", n=800, ss=2, raizes_visiveis=True, foco=True,
                    centro=None, R=None, Rref=None, centro_cor=None, pisos=None,
                    normalizar=True):
    """Renderiza o metodo de Newton do polinomio de coeficientes c.
    Devolve (imagem, fracao sem convergencia, pisos, (centro, R))."""
    c = np.asarray(c, float)
    bacia, nu, raizes, centro, R = calcular(c, n=n, ss=ss, foco=foco,
                                            centro=centro, R=R, Rref=Rref)
    cc = centro if centro_cor is None else centro_cor
    rgb, pisos = colorir(bacia, nu, raizes, cc, estilo, normalizar=normalizar, pisos=pisos)
    N = n * ss
    img = Image.fromarray((rgb * 255 + 0.5).astype(np.uint8))
    if raizes_visiveis:
        dr = ImageDraw.Draw(img)
        rad = 0.0055 * N
        cor_anel = (21, 23, 38) if estilo != "noturno" else (255, 255, 255)
        for z in raizes:
            px = (z.real - (centro.real - R)) / (2 * R) * (N - 1)
            py = ((centro.imag + R) - z.imag) / (2 * R) * (N - 1)
            if -rad <= px <= N + rad and -rad <= py <= N + rad:
                dr.ellipse([px - rad, py - rad, px + rad, py + rad],
                           fill=(255, 255, 255), outline=cor_anel, width=max(1, int(rad / 3)))
    return img.resize((n, n), Image.LANCZOS), np.mean(bacia < 0), pisos, (centro, R)


def renderizar(data, estilo="aquarela", n=800, ss=2, raizes_visiveis=True, foco=True,
               normalizar=True):
    img, nc, _, _ = renderizar_coef(coeficientes(data), estilo, n, ss, raizes_visiveis,
                                    foco, normalizar=normalizar)
    return img, nc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("data", help="DD/MM/AAAA")
    ap.add_argument("--estilo", default="aquarela", choices=list(ESTILOS))
    ap.add_argument("--n", type=int, default=1600, help="lado da imagem final, em pixels")
    ap.add_argument("--ss", type=int, default=2, help="fator de superamostragem")
    ap.add_argument("--sem-raizes", action="store_true")
    ap.add_argument("--sem-foco", action="store_true")
    ap.add_argument("--sem-normalizar", action="store_true",
                    help="nao normalizar o brilho por bacia (bacias lentas ficam escuras)")
    ap.add_argument("--saida", default=".")
    args = ap.parse_args()
    img, nc = renderizar(args.data, args.estilo, args.n, args.ss,
                         not args.sem_raizes, not args.sem_foco, not args.sem_normalizar)
    sufixo = "_sem-norm" if args.sem_normalizar else ""
    nome = os.path.join(args.saida, "newton_" + args.data.replace("/", "-") + sufixo + ".png")
    img.save(nome, optimize=True)
    print(f"{args.data}: sem convergencia {nc:.2e}; salvo em {nome}")

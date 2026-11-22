"""Recorta os retratos em 4:5 e salva em figs/retrato_<nome>.jpg (800x1000)."""
import sys
from PIL import Image, ImageOps

U = "/root/.claude/uploads/37821492-a744-5e11-bb0f-f9949aecec92/"
# nome: (arquivo, esquerda, topo, largura)  -- altura = 1.25 * largura
CORTES = {
    "cayley":     ("71f71dff-image.jpg", 0, 60, 1700),
    "gauss":      ("f60e4f63-image.jpg", 0, 0, 960),
    "julia":      ("97787b47-image.jpg", 368, 0, 544),
    "poincare":   ("3691523d-image.png", 0, 25, 1692),
    "euler":      ("15398396-image.jpg", 0, 0, 1920),
    "mandelbrot": ("f0ef53df-image.jpg", 0, 10, 276),
    "fatou":      ("720d9b44-image.jpg", 0, 0, 500),
    "newton":     ("7766aa06-image.jpg", 70, 0, 1850),
}


def tom(img, modo):
    if modo == "cor":
        return img
    g = ImageOps.grayscale(img)
    g = ImageOps.autocontrast(g, cutoff=0.5)
    # monocromatico quente: preto -> tinta (#2b2d42), branco -> papel (#f6f1e7)
    return ImageOps.colorize(g, black=(43, 45, 66), white=(246, 241, 231), mid=(150, 140, 128))


def recortar(modo="cor", saida="figs"):
    for nome, (arq, x, y, w) in CORTES.items():
        im = Image.open(U + arq).convert("RGB")
        h = round(1.25 * w)
        assert x + w <= im.width and y + h <= im.height, (nome, im.size, x, y, w, h)
        im = im.crop((x, y, x + w, y + h)).resize((800, 1000), Image.LANCZOS)
        tom(im, modo).save(f"{saida}/retrato_{nome}.jpg", quality=92, subsampling=0)


if __name__ == "__main__":
    recortar(sys.argv[1] if len(sys.argv) > 1 else "cor")

"""
Copy the figures into the flat submission folder under the journal's names
(Fig1a.jpg, ..., Fig11.jpg). Fractal pictures: PNG (1600 x 1600) -> JPEG,
quality 95, no chroma subsampling, 600 dpi tag. Portraits: copied unchanged.
Fig2.pdf and Fig5.pdf (English vector figures) are written by figuras_vetoriais_en.py.
Usage: python make_submission_figures.py   (reads figs/, writes submission/)
"""
import os
import shutil
from PIL import Image

SRC, OUT = "figs", "submission"

FRACTAIS = {
    "Fig1a": "newton_25-12-1642", "Fig1b": "newton_04-01-1643",
    "Fig3a": "cayley_z2-1", "Fig3b": "cayley_z3-1",
    "Fig4a": "zoom_1", "Fig4b": "zoom_2", "Fig4c": "zoom_3",
    "Fig7a": "newton_15-04-1707", "Fig7b": "newton_16-08-1821", "Fig7c": "newton_29-04-1854",
    "Fig7d": "newton_28-02-1878", "Fig7e": "newton_03-02-1893", "Fig7f": "newton_20-11-1924",
    "Fig8a": "newton_05-07-1687", "Fig8b": "newton_23-09-1846", "Fig8c": "newton_30-11-1889",
    "Fig8d": "newton_29-05-1919", "Fig8e": "newton_29-12-1972", "Fig8f": "newton_16-05-1980",
    "Fig9a": "newton_08-08-1900", "Fig9b": "newton_08-08-1900_sem-norm",
    "Fig10a": "newton_30-08-1984", "Fig10b": "newton_03-07-1986",
    "Fig10c": "newton_27-04-1952", "Fig10d": "newton_15-08-1950",
    "Fig11": "falha_z3-2z+2",
}
RETRATOS = {"Fig6a": "newton", "Fig6b": "euler", "Fig6c": "gauss", "Fig6d": "cayley",
            "Fig6e": "poincare", "Fig6f": "fatou", "Fig6g": "julia", "Fig6h": "mandelbrot"}

if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for alvo, nome in FRACTAIS.items():
        im = Image.open(os.path.join(SRC, nome + ".png")).convert("RGB")
        im.save(os.path.join(OUT, alvo + ".jpg"), quality=95, subsampling=0, dpi=(600, 600))
    for alvo, nome in RETRATOS.items():
        shutil.copyfile(os.path.join(SRC, f"retrato_{nome}.jpg"), os.path.join(OUT, alvo + ".jpg"))
    print(len(FRACTAIS) + len(RETRATOS), "figure files written to", OUT)

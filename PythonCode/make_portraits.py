"""Portraits of Fig. 6: crop to 4:5, resize to 800 x 1000 and tone in warm monochrome
(ink #2b2d42 to paper #f6f1e7). The source images (public domain, and the Oberwolfach
photograph of Mandelbrot, CC BY-SA 2.0 DE; see the caption of Fig. 6) are not included.

Usage: python make_portraits.py SOURCE_FOLDER      (writes figures/Fig6a.jpg ... Fig6h.jpg)
"""
import os
import sys

from PIL import Image, ImageOps

# figure: (source file, left, top, width of the crop); height = 1.25 * width
CROPS = {
    "Fig6a": ("SirIsaacNewton.jpg", 70, 0, 1850),
    "Fig6b": ("LeonhardEuler.jpg", 0, 0, 1920),
    "Fig6c": ("CarlFriedrichGauss.jpg", 0, 0, 960),
    "Fig6d": ("ArthurCayley.jpg", 0, 60, 1700),
    "Fig6e": ("HenriPoincare.png", 0, 25, 1692),
    "Fig6f": ("PierreFatou.jpg", 0, 0, 500),
    "Fig6g": ("GastonJulia.jpg", 368, 0, 544),
    "Fig6h": ("Mandelbrot.jpg", 0, 10, 276),
}


def tone(img):
    g = ImageOps.autocontrast(ImageOps.grayscale(img), cutoff=0.5)
    return ImageOps.colorize(g, black=(43, 45, 66), white=(246, 241, 231), mid=(150, 140, 128))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    os.makedirs("figures", exist_ok=True)
    for name, (source, x, y, w) in CROPS.items():
        im = Image.open(os.path.join(sys.argv[1], source)).convert("RGB")
        h = round(1.25 * w)
        assert x + w <= im.width and y + h <= im.height, (name, im.size)
        im = im.crop((x, y, x + w, y + h)).resize((800, 1000), Image.LANCZOS)
        tone(im).save(os.path.join("figures", name + ".jpg"), quality=92, subsampling=0)
        print(name + ".jpg")

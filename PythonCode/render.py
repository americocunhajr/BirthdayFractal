"""Picture of Newton's method for the polynomial of one date.

Usage:
  python render.py 30/08/1984                      # figures/newton_30-08-1984.png
  python render.py 30/08/1984 --order YMD          # digits in YYYYMMDD order
  python render.py 30/08/1984 --relax 0.9+0.2j     # relaxed step z - a p/p'
  python render.py 03/07/1986 --center -0.5 --half-width 0.05   # magnify a region
Other options: --n (side in pixels), --ss (supersampling), --no-normalize, --no-roots,
--out (folder).
"""
import argparse
import os

from birthday.dates import coefficients
from birthday.newton import render


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("date", help="DD/MM/YYYY")
    ap.add_argument("--order", default="DMY", choices=["DMY", "YMD"])
    ap.add_argument("--relax", type=complex, default=1.0, help="step factor a (1: Newton)")
    ap.add_argument("--center", type=complex, help="center of the window (default: roots)")
    ap.add_argument("--half-width", type=float, help="half-width of the window")
    ap.add_argument("--n", type=int, default=1600, help="side of the image in pixels")
    ap.add_argument("--ss", type=int, default=2, help="supersampling factor")
    ap.add_argument("--no-normalize", action="store_true",
                    help="do not normalize brightness per basin (slow basins look dark)")
    ap.add_argument("--no-roots", action="store_true", help="do not mark the roots")
    ap.add_argument("--out", default="figures", help="output folder")
    args = ap.parse_args()
    if (args.center is None) != (args.half_width is None):
        ap.error("--center and --half-width go together")

    img, unclassified, _, _ = render(coefficients(args.date, args.order), n=args.n, ss=args.ss,
                                     show_roots=not args.no_roots, center=args.center,
                                     R=args.half_width, normalize=not args.no_normalize,
                                     relax=args.relax)
    suffix = "".join([f"_{args.order}" if args.order != "DMY" else "",
                      "_relaxed" if args.relax != 1 else "",
                      "_zoom" if args.center is not None else "",
                      "_unnormalized" if args.no_normalize else ""])
    os.makedirs(args.out, exist_ok=True)
    path = os.path.join(args.out, "newton_" + args.date.replace("/", "-") + suffix + ".png")
    img.save(path, optimize=True)
    print(f"{args.date}: unclassified fraction {unclassified:.2e}; saved {path}")


if __name__ == "__main__":
    main()

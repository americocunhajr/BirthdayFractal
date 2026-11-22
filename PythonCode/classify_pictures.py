"""Classification of every point of each fractal picture of the article (Table S6).

Each picture is recomputed with the parameters of make_figures.py (3200 x 3200 grid, at
most 300 iterations), and its points are counted in four classes:
  (a) classified: stopped by the step criterion within 1e-4 R of a root;
  (b) non-finite: the step or the iterate ceased to be a finite floating-point number;
  (c) unclassified after 300 iterations;
  (d) stopped by the step criterion farther than 1e-4 R from every root.
For Fig. 11 (z^3 - 2z + 2), it also counts the unclassified points whose last iterate
lies within 1e-6 of the cycle {0, 1}.

Usage:
  python classify_pictures.py            # all pictures (about 7 minutes on one core)
  python classify_pictures.py i n        # part i of n, to run n processes in parallel
  python classify_pictures.py merge n    # merge the n parts
Output: results/picture_classification.json
"""
import json
import os
import sys

import numpy as np

from birthday.dates import coefficients
from birthday.newton import iterate, standard_window
from make_figures import DATE_FIGURES, GAUSS, POLE, PREIMAGE

OUT = "results/picture_classification.json"
CLASSES = ["classified", "non_finite", "unclassified", "stopped_far"]

PICTURES = [(name, coefficients(d), {}) for name, d in DATE_FIGURES.items()]
PICTURES.insert(2, ("Fig4a", GAUSS, {}))
center, R = standard_window(GAUSS)
_, Rf = standard_window([1, 0, -2, 2])
PICTURES += [
    ("Fig3a", np.array([1., 0, -1]), {}),
    ("Fig3b", np.array([1., 0, 0, -1]), {}),
    ("Fig11", np.array([1., 0, -2, 2]), dict(center=0.5 + 0j, R=1.0, Rref=Rf)),
    ("Fig4b", GAUSS, dict(center=POLE, R=R / 6, Rref=R)),
    ("Fig4c", GAUSS, dict(center=PREIMAGE, R=R / 48, Rref=R)),
]


def classify(pictures):
    result = {}
    for name, c, kw in pictures:
        _, _, _, cen, Rw, info = iterate(np.asarray(c, float), n=1600, ss=2, with_status=True, **kw)
        status = info["status"]
        item = {"points": int(status.size),
                "counts": {CLASSES[k]: int(np.sum(status == k)) for k in range(4)},
                "center": [float(np.real(cen)), float(np.imag(cen))],
                "half_width": float(Rw)}
        if name == "Fig11":
            last = info["Zlast"][status == 2]
            near = np.minimum(np.abs(last - 0), np.abs(last - 1))
            item["unclassified_near_cycle_1e-6"] = int(np.sum(near < 1e-6))
        result[name] = item
        print(name, item["counts"], item.get("unclassified_near_cycle_1e-6", ""), flush=True)
    return result


if __name__ == "__main__":
    args = sys.argv[1:]
    os.makedirs("results", exist_ok=True)
    order = [name for name, _, _ in PICTURES]
    if not args:
        json.dump(classify(PICTURES), open(OUT, "w"), indent=1)
    elif args[0] == "merge":
        result = {}
        for i in range(int(args[1])):
            result.update(json.load(open(OUT.replace(".json", f"_part{i}.json"))))
        assert sorted(result) == sorted(order), "missing parts"
        json.dump({k: result[k] for k in order}, open(OUT, "w"), indent=1)
        print(OUT, "with", len(result), "pictures")
    else:
        i, n = int(args[0]), int(args[1])
        json.dump(classify(PICTURES[i::n]), open(OUT.replace(".json", f"_part{i}.json"), "w"),
                  indent=1)

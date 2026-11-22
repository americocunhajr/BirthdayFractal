"""Finite-time Lyapunov indicator along the orbit of the real free critical point c of
N_p for 20/11/1924 (Table S4):

    lambda_n = (1/n) sum_{j=1}^{n} log |N_p'(x_j)|,   x_j = N_p^j(c),

starting at the first iterate, since N_p'(c) = 0. Sensitivity to n, to the arithmetic
(float64, long double, mpmath with 40 digits) and to the starting point.

Output: results/lyapunov.json   (about 1 minute)
"""
import json
import os

import mpmath as mp
import numpy as np

from birthday.dates import coefficients

c = coefficients("20/11/1924")               # 2 0 1 1 1 9 2 4
dc, d2c = np.polyder(c), np.polyder(c, 2)
roots = np.roots(c)
real_root = float(roots[np.argmin(np.abs(roots.imag))].real)
crit = [z for z in np.roots(d2c) if abs(z.imag) < 1e-12]
c0 = float(crit[0].real)
real_poles = sorted(float(z.real) for z in np.roots(dc) if abs(z.imag) < 1e-12)


def orbit(x0, n, dtype=np.float64):
    """lambda_n, smallest and largest iterate, and smallest distance to the real root."""
    cc, dd, ee = c.astype(dtype), dc.astype(dtype), d2c.astype(dtype)
    x = dtype(x0)
    total = dtype(0)
    xmin = xmax = None
    dmin_root = np.inf
    for _ in range(n):
        x = x - np.polyval(cc, x) / np.polyval(dd, x)      # x_{j+1}
        px, dpx, d2px = np.polyval(cc, x), np.polyval(dd, x), np.polyval(ee, x)
        total += np.log(abs(px * d2px / dpx ** 2))
        xmin = x if xmin is None else min(xmin, x)
        xmax = x if xmax is None else max(xmax, x)
        dmin_root = min(dmin_root, abs(float(x) - real_root))
    return float(total / n), float(xmin), float(xmax), dmin_root


def orbit_mp(x0, n, dps=40):
    mp.mp.dps = dps
    C = [mp.mpf(int(v)) for v in c]
    D = [mp.mpf(int(v)) * (len(c) - 1 - i) for i, v in enumerate(c[:-1])]
    E = [D[i] * (len(D) - 1 - i) for i in range(len(D) - 1)]
    x = mp.mpf(x0)
    total = mp.mpf(0)
    xmin, xmax = mp.inf, -mp.inf
    for _ in range(n):
        x = x - mp.polyval(C, x) / mp.polyval(D, x)
        total += mp.log(abs(mp.polyval(C, x) * mp.polyval(E, x) / mp.polyval(D, x) ** 2))
        xmin, xmax = min(xmin, x), max(xmax, x)
    return float(total / n), float(xmin), float(xmax)


def run(precision, x0_label, n, lam, a, b, dr=None):
    r = {"precision": precision, "x0": x0_label, "n": n, "lambda": lam, "min": a, "max": b}
    if dr is not None:
        r["min_distance_to_root"] = dr
    print(f"{precision:<17} x0={x0_label:<8} n={n:>8}: lambda={lam:.4f}  "
          f"interval=[{a:.4f}, {b:.4f}]")
    return r


if __name__ == "__main__":
    out = {"polynomial": c.tolist(), "critical_point": c0, "real_root": real_root,
           "real_poles": real_poles, "runs": []}
    print("critical point c =", c0, " real root =", real_root, " real poles =", real_poles)
    for n in (10**4, 10**5, 10**6):
        out["runs"].append(run("float64", "c", n, *orbit(c0, n)))
    out["runs"].append(run("long double", "c", 10**5, *orbit(c0, 10**5, np.longdouble)))
    for dx in (1e-8, -1e-8, 1e-4):
        out["runs"].append(run("float64", f"c{dx:+.0e}", 10**5, *orbit(c0 + dx, 10**5)))
    out["runs"].append(run("mpmath 40 digits", "c", 2 * 10**4, *orbit_mp(c0, 2 * 10**4)))
    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/lyapunov.json", "w"), indent=1)

"""Audit of the numbers stated in the article and in Online Resource 1. Each statement is
recomputed and compared with the published value; the report lists OK or FAILED per
statement. Reads results/lyapunov.json, results/orbit_intervals_certificate.json and
results/trapping_certificate.json.

Output: results/audit.txt   (about 30 seconds)
"""
import calendar
import json
import math
from fractions import Fraction

import mpmath as mp
import numpy as np
import sympy as sp

from birthday.critical import distinct_roots, test_critical_points
from birthday.dates import DATES, coefficients
from birthday.newton import framing_roots, standard_window

z = sp.symbols("z")
report = []


def check(statement, condition, detail=""):
    report.append((statement, bool(condition), detail))


def poly(date):
    return sp.Poly([int(v) for v in coefficients(date)], z)


def isolated(date):
    r = np.roots(coefficients(date))
    return len(framing_roots(r)) < len(distinct_roots(r)[0])


def largest_real_root(date):
    r = np.roots(coefficients(date))
    re = r[np.abs(r.imag) < 1e-9].real
    return re[np.argmax(np.abs(re))]


# 1. the dates and their polynomials
dates = [d for d, _ in DATES]
check("38 dates, all distinct", len(dates) == 38 and len(set(dates)) == 38, f"{len(dates)} dates")
check("polynomials = digits of DDMMYYYY without leading zeros",
      all([int(ch) for ch in d.replace("/", "")][8 - len(coefficients(d)):]
          == [int(v) for v in coefficients(d)] for d in dates))
published = {"25/12/1642": 2*z**7+5*z**6+z**5+2*z**4+z**3+6*z**2+4*z+2,
             "04/01/1643": 4*z**6+z**4+z**3+6*z**2+4*z+3,
             "30/04/1777": 3*z**7+4*z**4+z**3+7*z**2+7*z+7,
             "20/11/1924": 2*z**7+z**5+z**4+z**3+9*z**2+2*z+4,
             "08/08/1900": z**2*(8*z**4+8*z**2+z+9),
             "01/01/2000": z**3*(z**3+z+2)}
for d, e in published.items():
    check(f"published polynomial of {d}", sp.expand(e - poly(d).as_expr()) == 0)

# 2. distinct roots and multiplicities (exact): mult(0) = final zeros, the rest squarefree
for d in dates:
    s = d.replace("/", "")
    m0 = len(s) - len(s.rstrip("0"))
    Q = sp.Poly(sp.quo(poly(d).as_expr(), z**m0), z)
    squarefree = sp.gcd(Q, Q.diff(z)).degree() == 0 and Q.eval(0) != 0
    k = Q.degree() + (1 if m0 > 0 else 0)
    k_num = len(distinct_roots(np.roots(coefficients(d)))[0])
    check(f"{d}: mult(0) = {m0}, other roots simple, k = {k}", squarefree and k_num == k)

# 3. isolated roots and the Girard heuristic
large = [d for d in dates if coefficients(d)[1] / coefficients(d)[0] >= 3.5]
small = [d for d in dates if d not in large]
errors = [abs(largest_real_root(d) + coefficients(d)[1] / coefficients(d)[0])
          / abs(largest_real_root(d)) for d in large]
check("14 dates with c1/c0 >= 3.5, all with an isolated root",
      len(large) == 14 and all(isolated(d) for d in large), f"{len(large)} dates")
ratios_small = [coefficients(d)[1] / coefficients(d)[0] for d in small]
check("no date with c1/c0 <= 3 has an isolated root",
      not any(isolated(d) for d in small) and max(ratios_small) <= 3,
      f"max c1/c0 among the others = {max(ratios_small)}")
check("relative errors of -c1/c0 between 1% and 10%", min(errors) >= 0.01 and max(errors) <= 0.10,
      f"min {min(errors):.4f}, max {max(errors):.4f}")
days = sorted({int(d[:2]) for d in large})
check("days of the 14 dates in 14-19 or 27-29", all(14 <= x <= 19 or 27 <= x <= 29 for x in days),
      str(days))
check("Lorenz: prediction -4.5, root -4.44", abs(largest_real_root("29/12/1972") + 4.44) < 0.005)
check("Euler: prediction -5, root -5.15", abs(largest_real_root("15/04/1707") + 5.15) < 0.005)
check("Neptune: prediction -1.5, root -2.38, not isolated",
      abs(largest_real_root("23/09/1846") + 2.38) < 0.005
      and len(framing_roots(np.roots(coefficients("23/09/1846")))) == 7)
check("15/08/1950: isolated root -5.29", abs(largest_real_root("15/08/1950") + 5.29) < 0.005)
check("03/07/1986 (Fig. 10b): degree 6 and no real roots (exact Sturm count)",
      poly("03/07/1986").degree() == 6 and poly("03/07/1986").count_roots() == 0)
check("15/08/1950 (Fig. 10d): simple root at z = 0",
      poly("15/08/1950").eval(0) == 0 and poly("15/08/1950").diff(z).eval(0) != 0)
rg = np.sort_complex(np.round(np.roots(coefficients("30/04/1777")), 3))
check("published roots of Gauss's date", set(np.round(rg, 3)) == {
      -1.302, 1.102 + 0.838j, 1.102 - 0.838j, 0.111 + 1.132j, 0.111 - 1.132j,
      -0.562 + 0.638j, -0.562 - 0.638j}, str(rg))

# 4. critical points: exact classification and numerical test of the regular free ones
res = {d: test_critical_points(coefficients(d)) for d in dates}
check("2k - 2 critical points for every date",
      all(res[d][0]["total"] == 2 * res[d][0]["k"] - 2 for d in dates))
with_pole = [d for d in dates if res[d][0]["critical_poles"]]
check("critical pole only for 11/02/2016", with_pole == ["11/02/2016"], str(with_pole))
P16 = poly("11/02/2016")
cp16 = res["11/02/2016"][0]["critical_poles"]
check("11/02/2016: p' = (z+1)^2(z^2-z+1)(7z^2-z+1), p(-1) = 5, double pole at -1",
      sp.expand(P16.diff(z).as_expr() - (z+1)**2*(z**2-z+1)*(7*z**2-z+1)) == 0
      and P16.eval(-1) == 5 and P16.diff(z).eval(-1) == 0 and P16.diff(z).diff(z).eval(-1) == 0
      and len(cp16) == 1 and cp16[0]["point"] == "-1" and cp16[0]["pole_order"] == 2)
converge = [d for d in dates if all(its is not None for _, its, _ in res[d][1])]
maxits = max(max(its for _, its, _ in res[d][1]) for d in converge)
check("37 dates whose regular free critical orbits all satisfy the criterion",
      len(converge) == 37, f"exception: {[d for d in dates if d not in converge]}")
check("at most 147 iterations (regular free critical points)", maxits == 147,
      f"{maxits} (in {[d for d in converge if max(its for _, its, _ in res[d][1]) == maxits]})")
lm = res["20/11/1924"][1]
check("Mandelbrot: 5 regular free points, the 4 complex converge, the real one does not",
      len(lm) == 5 and sum(1 for z0, its, _ in lm if abs(z0.imag) > 1e-9 and its) == 4
      and sum(1 for z0, its, _ in lm if abs(z0.imag) < 1e-9 and its is None) == 1
      and not res["20/11/1924"][0]["critical_poles"])
d2 = sp.factor(poly("01/01/2000").diff(z).diff(z).as_expr())
check("01/01/2000: p'' = 6z(5z^3+2z+2); 3 simple roots and 3 regular free points",
      sp.expand(d2 - 6*z*(5*z**3+2*z+2)) == 0 and len(res["01/01/2000"][0]["regular_free"]) == 3
      and res["01/01/2000"][0]["simple_roots"] == 3)

# 5. Lyapunov indicator and critical-orbit intervals
M = json.load(open("results/lyapunov.json"))
lam = {(r["precision"], r["x0"], r["n"]): r["lambda"] for r in M["runs"]}
check("lambda_n (double) ~ 0.19 / 0.18 / 0.19 for n = 1e4 / 1e5 / 1e6",
      round(lam[("float64", "c", 10**4)], 2) == 0.19 and round(lam[("float64", "c", 10**5)], 2) == 0.18
      and round(lam[("float64", "c", 10**6)], 2) == 0.19)
check("all runs between 0.182 and 0.192", min(lam.values()) >= 0.182 and max(lam.values()) <= 0.192,
      f"min {min(lam.values()):.4f}, max {max(lam.values()):.4f}")
check("same visited interval in all runs",
      len({(round(r["min"], 4), round(r["max"], 4)) for r in M["runs"]}) == 1)
C = json.load(open("results/orbit_intervals_certificate.json"))
check("critical-orbit intervals certified (Sturm + interval arithmetic)", C["certified"])
check("distance from I- U I+ to the real root > 0.55 and to the real poles > 0.21",
      C["min_distance_to_real_root"] > 0.55 and C["min_distance_to_real_poles"] > 0.21)

# 6. Table 1 and the constant 1/(2 sqrt 2)
xs = [Fraction(1)]
for _ in range(4):
    xs.append(xs[-1] - (xs[-1] ** 2 - 2) / (2 * xs[-1]))
check("fractions of Table 1", [str(v) for v in xs] == ["1", "3/2", "17/12", "577/408", "665857/470832"])
mp.mp.dps = 50
e3, e4 = abs(mp.mpf(577) / 408 - mp.sqrt(2)), abs(mp.mpf(665857) / 470832 - mp.sqrt(2))
check("0.354 x (2.1e-6)^2 ~ 1.6e-12", abs(e3 ** 2 / (2 * mp.sqrt(2)) - e4) / e4 < 1e-5,
      f"e3 = {mp.nstr(e3, 3)}, e4 = {mp.nstr(e4, 3)}")

# 7. Fig. 4: window of Gauss's date
cen, R = standard_window(coefficients("30/04/1777"))
check("Fig. 4: center 0 and half-widths 1.800 / 0.300 / 0.0375",
      abs(cen) < 1e-12 and round(R, 3) == 1.8 and round(R / 6, 3) == 0.3 and round(R / 48, 4) == 0.0375)

# 8. Propositions 1-3 and Corollary 1
r = Fraction(5, 2)
S = Fraction(9) / r**2 + Fraction(2) / r**3 + Fraction(9) / r**4 + Fraction(9) / r**5 + Fraction(9) / r**6
check("Prop. 1: S(5/2) = 1.927424", S == Fraction(1927424, 1000000), str(float(S)))
check("Prop. 1: c0 r + S < c1 in both cases", 1 * r + S < 5 and 2 * r + S < 7)
prop1_days = {15, 16, 17, 18, 19, 27, 28, 29}
n_total = n_bad = 0
for year in range(1000, 3000):
    for month in range(1, 13):
        for day in prop1_days:
            if day > calendar.monthrange(year, month)[1]:
                continue
            c = coefficients(f"{day:02d}/{month:02d}/{year}")
            rr = np.roots(c)
            outside = rr[np.abs(rr) > 2.5]
            ok = (len(outside) == 1 and abs(outside[0].imag) < 1e-9 and outside[0].real < 0
                  and abs(outside[0].real + c[1] / c[0]) < 2 / c[0])
            n_total += 1
            n_bad += not ok
check("Prop. 1 confirmed numerically for every valid date of the family", n_bad == 0,
      f"{n_total} dates, {n_bad} violations")
covered = [d for d in dates if int(d[:2]) in prop1_days]
check("14 dates with c1/c0 >= 3.5, all isolated; 13 covered by Prop. 1; exception 14/09/2015",
      len(large) == 14 and all(isolated(d) for d in large) and len(covered) == 13
      and set(large) - set(covered) == {"14/09/2015"})
quoted = {"15/04/1707": -5.15, "16/08/1821": -6.21, "29/04/1854": -4.60, "28/02/1878": -4.07,
          "29/05/1919": -4.62, "29/12/1972": -4.44, "16/05/1980": -6.13, "27/04/1952": -3.66,
          "15/08/1950": -5.29, "14/09/2015": -4.43}
all_outside = True
for d, v in quoted.items():
    cen, R = standard_window(coefficients(d))
    rr = np.roots(coefficients(d))
    x = rr[np.argmax(np.abs(rr))].real
    all_outside &= abs(x - cen.real) > R and round(x, 2) == v
check("isolated roots quoted in the captions and text, all outside the window", all_outside)
cen, R = standard_window(coefficients("23/09/1846"))
x = min(r_.real for r_ in np.roots(coefficients("23/09/1846")) if abs(r_.imag) < 1e-9)
check("Neptune: real root -2.38 inside the window", round(x, 2) == -2.38 and abs(x - cen.real) < R)
check("Prop. 2: 34 steps (m = 2) and 57 steps (m = 3) from 1 to 1e-10",
      math.floor(math.log(1e10) / math.log(2)) + 1 == 34
      and math.floor(math.log(1e10) / math.log(1.5)) + 1 == 57)
no_real = [d for d in dates if len(sp.real_roots(poly(d))) == 0]
check("six dates without real roots, all passing the critical-point test",
      len(no_real) == 6 and all(d in converge for d in no_real), ", ".join(no_real))
T = json.load(open("results/trapping_certificate.json"))
V = T["checks"]
check("Prop. 3: exact certificate for K- and K+", T["certified"]
      and T["K_minus"] == ["-3813/5000", "-78/125"] and T["K_plus"] == ["79/500", "503/2500"],
      f"{T['K_minus']} {T['K_plus']}")
check("Prop. 3: distances > 0.55 (real root) and > 0.21 (real poles)",
      T["min_distance_K_to_real_root"] > 0.55 and T["min_distance_K_to_real_poles"] > 0.21)
check("Prop. 3: N(alpha) ~ 0.20105, N(gamma) ~ -0.76248; poles -0.9730, -0.1130; root -1.3169",
      round(T["N(alpha)_approx"], 5) == 0.20105 and round(T["N(gamma)_approx"], 5) == -0.76248
      and [round(v, 4) for v in T["real_poles"]] == [-0.973, -0.113] and round(T["real_root"], 4) == -1.3169)
check("Prop. 3: distance bounds certified by Sturm counts (0.55 and 0.21)",
      sum("within" in k for k in V) == 4 and all(V[k] for k in V if "within" in k))
tc = T["two_cycle"]
mi = [float(v) for v in tc["multiplier_interval"].strip("[]").split(",")]
check("Cor. 1: unique 2-cycle in K-, x* ~ -0.7269194 <-> 0.1681738, multiplier in (-1.90446, -1.90445)",
      V["N o N has a single fixed point in K- (unique 2-cycle)"]
      and V["repelling 2-cycle: (N o N)'(x*) < -1"]
      and round(tc["x_star_approx"], 7) == -0.7269194 and round(tc["N(x_star)_approx"], 7) == 0.1681738
      and -1.90446 < mi[0] and mi[1] < -1.90445 and tc["degree_of_F"] == 49,
      f"{tc['x_star_approx']:.7f} <-> {tc['N(x_star)_approx']:.7f}, multiplier {mi[0]:.6f}")
cen, R = standard_window(coefficients("20/11/1924"))
check("grid rows nearest to the real axis at ~6e-4", 5e-4 < R / 3199 < 7e-4, f"{R / 3199:.2e}")
check("errors of 1e-16 grow to order 1 in ~200 iterations (lambda ~ 0.19)",
      150 < math.log(1e16) / 0.19 < 250, f"{math.log(1e16) / 0.19:.0f}")

# report
width = max(len(s) for s, _, _ in report)
with open("results/audit.txt", "w", encoding="utf-8") as f:
    for s, ok, detail in report:
        line = f"{'OK    ' if ok else 'FAILED'}  {s.ljust(width)}  {detail}".rstrip()
        print(line)
        f.write(line + "\n")
print(f"\n{sum(ok for _, ok, _ in report)}/{len(report)} statements confirmed")

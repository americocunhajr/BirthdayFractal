"""Computer-assisted verification of the critical-orbit intervals of 20/11/1924 (Table S5).

p(z) = 2z^7 + z^5 + z^4 + z^3 + 9z^2 + 2z + 4; c is the only real zero of p'' (near
-0.6888); x_j = N_p^j(c); I- = [x2, x4], I+ = [x1, x3].

Rigorous tools: exact rational arithmetic and Sturm sequences (SymPy) to isolate c and
count the zeros of p, p', p'' on intervals with rational endpoints; interval arithmetic
(mpmath.iv, 80 digits) to enclose x1, ..., x5. No inequality is decided in ordinary
floating point.

Checks (all must be True):
  A1  p'' is squarefree and c is a simple zero, isolated in [c_lo, c_hi];
  A2  on J- = [x2_lo, x4_hi] (which contains I-): p and p' have no zeros, p'' has exactly
      one zero, c;
  A3  on J+ = [x1_lo, x3_hi] (which contains I+): p, p', p'' have no zeros;
  A4  x2 < c < x4,  x1 < x3,  x2 < x4,  x5 < x3.
Consequences (argument in Online Resource 1): N_p(I-) = I+, N_p(I+) = I-, and the orbit
of c stays in I- U I+, which contains no roots and no poles of N_p.

Output: results/orbit_intervals_certificate.json
"""
import json
import os

import mpmath
import sympy as sp
from mpmath import iv

z = sp.symbols("z")
p = sp.Poly(2*z**7 + z**5 + z**4 + z**3 + 9*z**2 + 2*z + 4, z)
dp, d2p, d3p = p.diff(z), p.diff(z).diff(z), p.diff(z).diff(z).diff(z)

# --- A1: isolate c exactly
assert sp.gcd(d2p, d3p).degree() == 0, "p'' is not squarefree"
iso = [I for I in d2p.intervals() if I[0][0] <= sp.Rational(-6888, 10000) <= I[0][1]]
assert len(iso) == 1
(c_lo, c_hi), _ = iso[0]
c_lo, c_hi = d2p.refine_root(c_lo, c_hi, eps=sp.Rational(1, 10**40))
assert d2p.eval(c_lo) * d2p.eval(c_hi) < 0
n_real_d2p = d2p.count_roots()       # all real zeros of p''

# --- enclose x1..x5 with interval arithmetic
iv.dps = 80
coef = [int(a) for a in p.all_coeffs()]
dcoef = [int(a) for a in dp.all_coeffs()]


def horner(cs, X):
    acc = iv.mpf(0)
    for a in cs:
        acc = acc * X + a
    return acc


def N(X):
    d = horner(dcoef, X)
    assert not (d.a <= 0 <= d.b), "p' may vanish on the interval"
    return X - horner(coef, X) / d


def _rational(t):  # mpf tuple (sign, mantissa, exponent, bits) -> exact rational
    sign, man, exp, _ = t
    r = sp.Integer(int(man)) * sp.Integer(2) ** int(exp) if exp >= 0 else \
        sp.Rational(int(man), 2 ** int(-exp))
    return -r if sign else r


def lo_q(X):  # exact lower endpoint of an mpmath.iv interval
    return _rational(X._mpi_[0])


def hi_q(X):  # exact upper endpoint
    return _rational(X._mpi_[1])


A = iv.mpf(int(c_lo.p)) / iv.mpf(int(c_lo.q))   # encloses c_lo
B = iv.mpf(int(c_hi.p)) / iv.mpf(int(c_hi.q))   # encloses c_hi
C = iv.mpf([A.a, B.b])                            # encloses [c_lo, c_hi]
assert lo_q(C) <= c_lo and c_hi <= hi_q(C)
X = [C]
for j in range(5):
    X.append(N(X[-1]))
lo = [lo_q(Xi) for Xi in X]
hi = [hi_q(Xi) for Xi in X]
width = max(hi[j] - lo[j] for j in range(6))

# --- A2, A3: exact counts (Sturm) on rational intervals
Jm = (lo[2], hi[4])
Jp = (lo[1], hi[3])
A2 = (p.count_roots(*Jm) == 0 and dp.count_roots(*Jm) == 0 and d2p.count_roots(*Jm) == 1
      and Jm[0] < c_lo and c_hi < Jm[1])
A3 = (p.count_roots(*Jp) == 0 and dp.count_roots(*Jp) == 0 and d2p.count_roots(*Jp) == 0)

# --- A4: inequalities between rational endpoints
A4 = {"x2 < c": hi[2] < c_lo, "c < x4": c_hi < lo[4], "x1 < x3": hi[1] < lo[3],
      "x2 < x4": hi[2] < lo[4], "x5 < x3": hi[5] < lo[3]}

# distances from I- U I+ to the real root of p and to the real poles (exact lower bounds)
real_root = list(p.intervals())
assert len(real_root) == 1
(r_lo, r_hi), _ = real_root[0]
r_lo, r_hi = p.refine_root(r_lo, r_hi, eps=sp.Rational(1, 10**12))
poles = [dp.refine_root(a, b, eps=sp.Rational(1, 10**12)) for (a, b), _ in dp.intervals()]

out = {
    "c_interval": [str(sp.N(c_lo, 30)), str(sp.N(c_hi, 30))],
    "real_zeros_of_p2": int(n_real_d2p),
    "x": {f"x{j}": [str(sp.N(lo[j], 20)), str(sp.N(hi[j], 20))] for j in range(6)},
    "max_enclosure_width": str(sp.N(width, 3)),
    "A1_c_simple_isolated_zero": True,
    "A2_J_minus": bool(A2),
    "A3_J_plus": bool(A3),
    "A4": {k: bool(v) for k, v in A4.items()},
    "real_root_of_p": [float(r_lo), float(r_hi)],
    "real_poles": [[float(a), float(b)] for a, b in poles],
    "min_distance_to_real_root": float(lo[2] - r_hi),
    "min_distance_to_real_poles": float(min(min(abs(lo[2] - b), abs(a - hi[4]), abs(lo[1] - b),
                                                abs(a - hi[3])) for a, b in poles)),
    "versions": {"sympy": sp.__version__, "mpmath": mpmath.__version__},
}
out["certified"] = bool(A2 and A3 and all(A4.values()))

# --- exact record (rationals as "p/q"), for independent verification
R = lambda q_: str(sp.Rational(q_))
pairs = {"x2 < c": (hi[2], c_lo), "c < x4": (c_hi, lo[4]), "x1 < x3": (hi[1], lo[3]),
         "x2 < x4": (hi[2], lo[4]), "x5 < x3": (hi[5], lo[3])}
out["exact"] = {
    "polynomial": str(p.as_expr()),
    "gcd_p1_p2_degree": int(sp.gcd(dp, d2p).degree()),    # 0: no critical poles
    "gcd_p2_p3_degree": int(sp.gcd(d2p, d3p).degree()),   # 0: p'' squarefree
    "real_zeros_of_p2": int(n_real_d2p),
    "c_isolation": {"lo": R(c_lo), "hi": R(c_hi),
                    "sign_p2_lo": int(sp.sign(d2p.eval(c_lo))),
                    "sign_p2_hi": int(sp.sign(d2p.eval(c_hi)))},
    "x_enclosures": {f"x{j}": {"lo": R(lo[j]), "hi": R(hi[j])} for j in range(6)},
    "J_minus": {"lo": R(Jm[0]), "hi": R(Jm[1])},
    "J_plus": {"lo": R(Jp[0]), "hi": R(Jp[1])},
    "sturm_counts": {
        "J_minus": {"p": int(p.count_roots(*Jm)), "p1": int(dp.count_roots(*Jm)),
                    "p2": int(d2p.count_roots(*Jm))},
        "J_plus": {"p": int(p.count_roots(*Jp)), "p1": int(dp.count_roots(*Jp)),
                   "p2": int(d2p.count_roots(*Jp))}},
    "inequalities": {k: {"left_bound": R(a), "right_bound": R(b),
                         "positive_difference": bool(b - a > 0)} for k, (a, b) in pairs.items()},
    "real_root_of_p_isolation": {"lo": R(r_lo), "hi": R(r_hi)},
    "real_zeros_of_p1_isolation": [{"lo": R(a), "hi": R(b)} for a, b in poles],
    "real_root_distance_bound": R(lo[2] - r_hi),
}
assert out["exact"]["gcd_p1_p2_degree"] == 0 and out["exact"]["gcd_p2_p3_degree"] == 0
assert out["exact"]["c_isolation"]["sign_p2_lo"] * out["exact"]["c_isolation"]["sign_p2_hi"] < 0
assert all(v["positive_difference"] for v in out["exact"]["inequalities"].values())
os.makedirs("results", exist_ok=True)
json.dump(out, open("results/orbit_intervals_certificate.json", "w"), indent=1)
print("certified:", out["certified"], "| max enclosure width:", out["max_enclosure_width"])
print("Sturm counts:", out["exact"]["sturm_counts"])

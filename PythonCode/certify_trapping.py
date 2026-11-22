"""Exact certificate for Proposition 3 and Corollary 1 of the article (date 20/11/1924).

p(z) = 2z^7 + z^5 + z^4 + z^3 + 9z^2 + 2z + 4,   N(x) = x - p(x)/p'(x),
K- = [alpha, beta] = [-0.7626, -0.6240],   K+ = [gamma, delta] = [0.1580, 0.2012].

Checked in exact rational arithmetic (SymPy, Sturm sequences):
  (a) p and p' have no zeros in K- or in K+;
  (b) for t in {gamma, delta}, q_t(x) = (x - t) p'(x) - p(x) has no zeros in K-, and for
      t in {alpha, beta} it has no zeros in K+;
  (c) N(alpha) lies in (gamma, delta) and N(gamma) lies in (alpha, beta);
  (d) p'' has exactly one real zero c, and p''(alpha), p''(beta) have opposite signs;
  (e) p has no zeros within 0.55 of K- or K+, and p' has none within 0.21.
Hence N(K-) is in K+ and N(K+) is in K- (see the proof in the article), c lies in K-, and
no orbit starting in K- U K+ converges to a root.
Corollary 1: (f) the numerator F of N(N(x)) - x (degree 49) has exactly one real zero in
K-, so there is a unique 2-cycle meeting K- U K+; (g) interval arithmetic (mpmath.iv,
60 digits) on a rational isolating interval of x* gives (N o N)'(x*) < -1 (repelling).

Output: results/trapping_certificate.json
"""
import json
import os

import sympy as sp
from mpmath import iv

z = sp.symbols("z")
P = sp.Poly([2, 0, 1, 1, 1, 9, 2, 4], z)
P1 = P.diff(z)
P2 = P1.diff(z)

alpha, beta = sp.Rational("-0.7626"), sp.Rational("-0.6240")
gamma, delta = sp.Rational("0.1580"), sp.Rational("0.2012")
Km, Kp = (alpha, beta), (gamma, delta)


def q(t):
    return sp.Poly((z - t) * P1.as_expr() - P.as_expr(), z)


def N(x):
    return x - P.eval(x) / P1.eval(x)


def count(poly, K):
    return int(poly.count_roots(K[0], K[1]))


checks = {
    "p has no zeros in K-": count(P, Km) == 0,
    "p has no zeros in K+": count(P, Kp) == 0,
    "p' has no zeros in K-": count(P1, Km) == 0,
    "p' has no zeros in K+": count(P1, Kp) == 0,
    "q_gamma has no zeros in K-": count(q(gamma), Km) == 0,
    "q_delta has no zeros in K-": count(q(delta), Km) == 0,
    "q_alpha has no zeros in K+": count(q(alpha), Kp) == 0,
    "q_beta has no zeros in K+": count(q(beta), Kp) == 0,
    "N(alpha) in (gamma, delta)": gamma < N(alpha) < delta,
    "N(gamma) in (alpha, beta)": alpha < N(gamma) < beta,
    "p'' has a single real zero": int(P2.count_roots()) == 1,
    "p'' changes sign on K-": P2.eval(alpha) * P2.eval(beta) < 0,
}
d_root, d_pole = sp.Rational("0.55"), sp.Rational("0.21")
checks.update({
    "p has no zeros within 0.55 of K-": count(P, (alpha - d_root, beta + d_root)) == 0,
    "p has no zeros within 0.55 of K+": count(P, (gamma - d_root, delta + d_root)) == 0,
    "p' has no zeros within 0.21 of K-": count(P1, (alpha - d_pole, beta + d_pole)) == 0,
    "p' has no zeros within 0.21 of K+": count(P1, (gamma - d_pole, delta + d_pole)) == 0,
})

# ---- period-two orbit: N(N(x)) = x  <=>  F(x) = H_A(A, B) - x B H_B(A, B) = 0 on K-,
# where N = A/B, A = x p' - p, B = p', and H_f(A, B) = B^deg f * f(A/B)
A = sp.Poly(z * P1.as_expr() - P.as_expr(), z)
B = P1


def homogenize(f):
    n = f.degree()
    return sum((A ** i) * (B ** (n - i)) * c for i, c in zip(range(n, -1, -1), f.all_coeffs()))


F = homogenize(A) - sp.Poly(z, z) * B * homogenize(B)
checks["N o N has a single fixed point in K- (unique 2-cycle)"] = count(F, Km) == 1
iso = [I for I in F.intervals(eps=sp.Rational(1, 10**14)) if alpha <= I[0][0] and I[0][1] <= beta]
u, v = iso[0][0]

iv.dps = 60
cp = [iv.mpf(int(c)) for c in P.all_coeffs()]
cp1 = [iv.mpf(int(c)) for c in P1.all_coeffs()]
cp2 = [iv.mpf(int(c)) for c in P2.all_coeffs()]


def horner(cs, X):
    y = iv.mpf(0)
    for c in cs:
        y = y * X + c
    return y


def N_iv(X):
    return X - horner(cp, X) / horner(cp1, X)


def dN_iv(X):
    return horner(cp, X) * horner(cp2, X) / horner(cp1, X) ** 2


# enclosure of the rational interval [u, v]: 70-digit decimals, widened by 1e-50
# with outward-rounded interval arithmetic
X0 = iv.mpf([str(sp.N(u, 70)), str(sp.N(v, 70))])
X = iv.mpf([(X0.a - iv.mpf("1e-50")).a, (X0.b + iv.mpf("1e-50")).b])
assert X.a <= iv.mpf(str(sp.N(u, 70))).a and X.b >= iv.mpf(str(sp.N(v, 70))).b
multiplier = dN_iv(N_iv(X)) * dN_iv(X)
checks["repelling 2-cycle: (N o N)'(x*) < -1"] = bool(multiplier.b < -1)
checks = {k: bool(val) for k, val in checks.items()}
x_star = (u + v) / 2

real_root = sp.real_roots(P)
real_poles = sp.real_roots(P1)
dist_root = min(abs(sp.N(r, 30) - x) for r in real_root for x in (alpha, beta, gamma, delta))
dist_poles = min(min(abs(sp.N(r, 30) - a), abs(sp.N(r, 30) - b))
                 for r in real_poles for (a, b) in (Km, Kp))
out = {
    "polynomial": "2z^7+z^5+z^4+z^3+9z^2+2z+4",
    "K_minus": [str(alpha), str(beta)],
    "K_plus": [str(gamma), str(delta)],
    "N(alpha)": str(N(alpha)), "N(gamma)": str(N(gamma)),
    "N(alpha)_approx": float(N(alpha)), "N(gamma)_approx": float(N(gamma)),
    "checks": checks,
    "real_root": float(real_root[0]),
    "real_poles": [float(r) for r in real_poles],
    "min_distance_K_to_real_root": float(dist_root),
    "min_distance_K_to_real_poles": float(dist_poles),
    "two_cycle": {"x_star_isolating_interval": [str(u), str(v)],
                  "x_star_approx": float(x_star),
                  "N(x_star)_approx": float(N(x_star)),
                  "multiplier_interval": str(multiplier),
                  "degree_of_F": F.degree()},
    "certified": all(checks.values()),
}
os.makedirs("results", exist_ok=True)
json.dump(out, open("results/trapping_certificate.json", "w"), indent=1)
for k, val in checks.items():
    print(("OK    " if val else "FAIL  ") + k)
print("N(alpha) =", float(N(alpha)), " N(gamma) =", float(N(gamma)))
print("distance from K- U K+ to the real root:", round(float(dist_root), 4),
      " to the real poles:", round(float(dist_poles), 4))
print("2-cycle:", float(x_star), "<->", float(N(x_star)), " multiplier in", multiplier)
print("certified:", out["certified"])

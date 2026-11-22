"""
Exact certificate for Proposition 3 of the article (date 20/11/1924).

p(z) = 2z^7 + z^5 + z^4 + z^3 + 9z^2 + 2z + 4,   N(x) = x - p(x)/p'(x).
K- = [alpha, beta] = [-0.7626, -0.6240],   K+ = [gamma, delta] = [0.1580, 0.2012].

Claims, all checked in exact rational arithmetic (SymPy, Sturm sequences):
  (a) p and p' have no zeros in K- or in K+;
  (b) for t in {gamma, delta}, q_t(x) = (x - t) p'(x) - p(x) has no zeros in K-,
      and for t in {alpha, beta} it has no zeros in K+;
  (c) N(alpha) lies in (gamma, delta) and N(gamma) lies in (alpha, beta);
  (d) p'' has exactly one real zero c, and p''(alpha), p''(beta) have opposite signs;
  (e) p has no zeros in [alpha - 0.55, beta + 0.55] or in [gamma - 0.55, delta + 0.55],
      and p' has none in [alpha - 0.21, beta + 0.21] or in [gamma - 0.21, delta + 0.21],
      so every point of K- U K+ is farther than 0.55 from every real root and
      farther than 0.21 from every real pole of N.
Corollary (period-two orbit): N o N maps K- into itself, so it has a fixed point x*
in K- (intermediate value theorem), and {x*, N(x*)} is a 2-cycle because K- and K+
are disjoint. In addition, (f) the numerator F of N(N(x)) - x has exactly one real
zero in K- (Sturm), so this is the only 2-cycle meeting K- U K+, and (g) interval
arithmetic (mpmath.iv, 60 digits) on a rational isolating interval of x* gives
(N o N)'(x*) < -1: the cycle is repelling.
Since p' does not vanish on K-, N is continuous there, and N(x) = t is equivalent to
q_t(x) = 0. By (b) the interval N(K-) does not contain gamma or delta, and by (c) it
meets (gamma, delta); hence N(K-) is contained in K+. Likewise N(K+) is contained in K-.
By (d), c lies in K-. Every orbit that starts in K- or K+ therefore stays in K- U K+,
a compact set that contains no root of p, so it does not converge to a root.

Output: certificado_aprisionamento.json
"""
import json
import sympy as sp

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
    "p sem zeros em K-": count(P, Km) == 0,
    "p sem zeros em K+": count(P, Kp) == 0,
    "p' sem zeros em K-": count(P1, Km) == 0,
    "p' sem zeros em K+": count(P1, Kp) == 0,
    "q_gamma sem zeros em K-": count(q(gamma), Km) == 0,
    "q_delta sem zeros em K-": count(q(delta), Km) == 0,
    "q_alpha sem zeros em K+": count(q(alpha), Kp) == 0,
    "q_beta sem zeros em K+": count(q(beta), Kp) == 0,
    "N(alpha) em (gamma, delta)": gamma < N(alpha) < delta,
    "N(gamma) em (alpha, beta)": alpha < N(gamma) < beta,
    "p'' com um unico zero real": int(P2.count_roots()) == 1,
    "p'' muda de sinal em K-": P2.eval(alpha) * P2.eval(beta) < 0,
}
d_raiz, d_polo = sp.Rational("0.55"), sp.Rational("0.21")
checks.update({
    "p sem zeros a menos de 0,55 de K-": count(P, (alpha - d_raiz, beta + d_raiz)) == 0,
    "p sem zeros a menos de 0,55 de K+": count(P, (gamma - d_raiz, delta + d_raiz)) == 0,
    "p' sem zeros a menos de 0,21 de K-": count(P1, (alpha - d_polo, beta + d_polo)) == 0,
    "p' sem zeros a menos de 0,21 de K+": count(P1, (gamma - d_polo, delta + d_polo)) == 0,
})

# ---- period-two orbit: N(N(x)) = x  <=>  F(x) = H_A(A, B) - x B H_B(A, B) = 0 on K-,
# where N = A/B, A = x p' - p, B = p', and H_f(A, B) = B^deg f * f(A/B)
A = sp.Poly(z * P1.as_expr() - P.as_expr(), z)
B = P1


def homog(f):
    n = f.degree()
    return sum((A ** i) * (B ** (n - i)) * c for i, c in zip(range(n, -1, -1), f.all_coeffs()))


F = homog(A) - sp.Poly(z, z) * B * homog(B)
n_ciclo = count(F, Km)
checks["N o N com um unico ponto fixo em K- (2-ciclo unico)"] = n_ciclo == 1
iso = [iv for iv in F.intervals(eps=sp.Rational(1, 10**14)) if alpha <= iv[0][0] and iv[0][1] <= beta]
u, v = iso[0][0]

from mpmath import iv as I, mp
I.dps = 60
cp = [I.mpf(int(c)) for c in P.all_coeffs()]
cp1 = [I.mpf(int(c)) for c in P1.all_coeffs()]
cp2 = [I.mpf(int(c)) for c in P2.all_coeffs()]


def hv(cs, X):
    y = I.mpf(0)
    for c in cs:
        y = y * X + c
    return y


def Nint(X):
    return X - hv(cp, X) / hv(cp1, X)


def dNint(X):
    return hv(cp, X) * hv(cp2, X) / hv(cp1, X) ** 2


# enclosure of the rational interval [u, v]: 70-digit decimals, widened by 1e-50
# with outward-rounded interval arithmetic
X0 = I.mpf([str(sp.N(u, 70)), str(sp.N(v, 70))])
X = I.mpf([(X0.a - I.mpf("1e-50")).a, (X0.b + I.mpf("1e-50")).b])
assert X.a <= I.mpf(str(sp.N(u, 70))).a and X.b >= I.mpf(str(sp.N(v, 70))).b
mult = dNint(Nint(X)) * dNint(X)
checks["2-ciclo repulsor: (N o N)'(x*) < -1"] = bool(mult.b < -1)
checks = {k: bool(v) for k, v in checks.items()}
x_estrela = (u + v) / 2

raiz = sp.real_roots(P)
polos = sp.real_roots(P1)
dist_raiz = min(abs(sp.N(r, 30) - x) for r in raiz for x in (alpha, beta, gamma, delta))
dist_polos = min(min(abs(sp.N(r, 30) - a), abs(sp.N(r, 30) - b))
                 for r in polos for (a, b) in (Km, Kp))
out = {
    "polinomio": "2z^7+z^5+z^4+z^3+9z^2+2z+4",
    "K_menos": [str(alpha), str(beta)],
    "K_mais": [str(gamma), str(delta)],
    "N(alpha)": str(N(alpha)), "N(gamma)": str(N(gamma)),
    "N(alpha)_aprox": float(N(alpha)), "N(gamma)_aprox": float(N(gamma)),
    "verificacoes": checks,
    "raiz_real": float(raiz[0]),
    "polos_reais": [float(r) for r in polos],
    "dist_min_K_raiz_real": float(dist_raiz),
    "dist_min_K_polos_reais": float(dist_polos),
    "dois_ciclo": {"x_estrela_isolamento": [str(u), str(v)],
                   "x_estrela_aprox": float(x_estrela),
                   "N(x_estrela)_aprox": float(N(x_estrela)),
                   "multiplicador_intervalo": str(mult),
                   "grau_de_F": F.degree()},
    "certificado": all(checks.values()),
}
json.dump(out, open("certificado_aprisionamento.json", "w"), indent=1, ensure_ascii=False)
for k, v in checks.items():
    print(("OK   " if v else "FAIL ") + k)
print("N(alpha) =", float(N(alpha)), " N(gamma) =", float(N(gamma)))
print("distance from K-uK+ to the real root:", round(float(dist_raiz), 4),
      " to the real poles:", round(float(dist_polos), 4))
print("2-cycle:", float(x_estrela), "<->", float(N(x_estrela)), " multiplier in", mult)
print("certificate:", out["certificado"])

"""
Verificacao assistida por computador do confinamento da orbita critica real de
N_p para p(z) = 2z^7 + z^5 + z^4 + z^3 + 9z^2 + 2z + 4 (data 20/11/1924).

Notacao: c e o unico zero real de p'' perto de -0,6888; x_j = N_p^j(c);
I- = [x2, x4], I+ = [x1, x3].

Ferramentas rigorosas:
  - aritmetica racional exata e sequencias de Sturm (sympy) para isolar c e contar
    raizes de p, p' e p'' em intervalos com extremos racionais;
  - aritmetica intervalar (mpmath.iv, 80 digitos) para encerrar x1, ..., x5.
Nenhuma desigualdade e decidida em ponto flutuante comum.

Afirmacoes verificadas (todas devem dar True):
  A1  p'' e livre de quadrados e c e um zero simples, isolado em [c_lo, c_hi];
  A2  em J- = [x2_lo, x4_hi] (contem I-): p e p' nao tem zeros; p'' tem exatamente
      um zero, que e c;
  A3  em J+ = [x1_lo, x3_hi] (contem I+): p, p' e p'' nao tem zeros;
  A4  x2 < c < x4,  x1 < x3,  x2 < x4,  x5 < x3.
Consequencias (argumento no artigo): N_p(I-) = I+, N_p(I+) = I-, a orbita de c fica
em I- U I+, que nao contem raizes nem polos de N_p.
"""
import json
import sympy as sp
from mpmath import iv, mp

z = sp.symbols("z")
p = sp.Poly(2*z**7 + z**5 + z**4 + z**3 + 9*z**2 + 2*z + 4, z)
dp, d2p, d3p = p.diff(z), p.diff(z).diff(z), p.diff(z).diff(z).diff(z)

# --- A1: isolar c exatamente
assert sp.gcd(d2p, d3p).degree() == 0, "p'' nao e livre de quadrados"
iso = [iv_ for iv_ in d2p.intervals() if iv_[0][0] <= sp.Rational(-6888, 10000) <= iv_[0][1]]
assert len(iso) == 1
(c_lo, c_hi), _ = iso[0]
c_lo, c_hi = d2p.refine_root(c_lo, c_hi, eps=sp.Rational(1, 10**40))
assert d2p.eval(c_lo) * d2p.eval(c_hi) < 0
n_reais_d2p = d2p.count_roots()       # todos os zeros reais de p''

# --- encerrar x1..x5 com aritmetica intervalar
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
    assert not (d.a <= 0 <= d.b), "p' pode se anular no intervalo"
    return X - horner(coef, X) / d


def _rac(t):  # tupla mpf (sinal, mantissa, expoente, bits) -> racional exato
    sinal, man, exp, _ = t
    r = sp.Integer(int(man)) * sp.Integer(2) ** int(exp) if exp >= 0 else \
        sp.Rational(int(man), 2 ** int(-exp))
    return -r if sinal else r


def lo_q(X):  # extremo inferior exato de um intervalo mpmath.iv
    return _rac(X._mpi_[0])


def hi_q(X):  # extremo superior exato
    return _rac(X._mpi_[1])


A = iv.mpf(int(c_lo.p)) / iv.mpf(int(c_lo.q))   # encerra c_lo
B = iv.mpf(int(c_hi.p)) / iv.mpf(int(c_hi.q))   # encerra c_hi
C = iv.mpf([A.a, B.b])                            # encerra [c_lo, c_hi]
assert lo_q(C) <= c_lo and c_hi <= hi_q(C)
X = [C]
for j in range(5):
    X.append(N(X[-1]))
lo = [lo_q(Xi) for Xi in X]
hi = [hi_q(Xi) for Xi in X]
larg = max(hi[j] - lo[j] for j in range(6))

# --- A2, A3: contagens exatas (Sturm) em intervalos racionais
Jm = (lo[2], hi[4])
Jp = (lo[1], hi[3])
A2 = (p.count_roots(*Jm) == 0 and dp.count_roots(*Jm) == 0 and d2p.count_roots(*Jm) == 1
      and Jm[0] < c_lo and c_hi < Jm[1])
A3 = (p.count_roots(*Jp) == 0 and dp.count_roots(*Jp) == 0 and d2p.count_roots(*Jp) == 0)

# --- A4: desigualdades com extremos racionais
A4 = {"x2 < c": hi[2] < c_lo, "c < x4": c_hi < lo[4], "x1 < x3": hi[1] < lo[3],
      "x2 < x4": hi[2] < lo[4], "x5 < x3": hi[5] < lo[3]}

# distancias da uniao I- U I+ a raiz real de p e aos polos reais (cotas inferiores exatas)
raiz_real = [r for r in p.intervals()]
assert len(raiz_real) == 1
(r_lo, r_hi), _ = raiz_real[0]
r_lo, r_hi = p.refine_root(r_lo, r_hi, eps=sp.Rational(1, 10**12))
polos = [dp.refine_root(a, b, eps=sp.Rational(1, 10**12)) for (a, b), _ in dp.intervals()]

out = {
    "c_intervalo": [str(sp.N(c_lo, 30)), str(sp.N(c_hi, 30))],
    "zeros_reais_de_p2": int(n_reais_d2p),
    "x": {f"x{j}": [str(sp.N(lo[j], 20)), str(sp.N(hi[j], 20))] for j in range(6)},
    "largura_maxima_dos_encerramentos": str(sp.N(larg, 3)),
    "A1_c_zero_simples_isolado": True,
    "A2_J_menos": bool(A2),
    "A3_J_mais": bool(A3),
    "A4": {k: bool(v) for k, v in A4.items()},
    "raiz_real_de_p": [float(r_lo), float(r_hi)],
    "polos_reais": [[float(a), float(b)] for a, b in polos],
    "dist_min_a_raiz_real": float(lo[2] - r_hi),
    "dist_min_aos_polos_reais": float(min(min(abs(lo[2] - b), abs(a - hi[4]), abs(lo[1] - b),
                                              abs(a - hi[3])) for a, b in polos)),
    "versoes": {"sympy": sp.__version__, "mpmath": __import__("mpmath").__version__},
}
out["certificado"] = bool(A2 and A3 and all(A4.values()))

# --- registro exato (racionais como "p/q"), para verificacao independente
R = lambda q_: str(sp.Rational(q_))
pares = {"x2 < c": (hi[2], c_lo), "c < x4": (c_hi, lo[4]), "x1 < x3": (hi[1], lo[3]),
         "x2 < x4": (hi[2], lo[4]), "x5 < x3": (hi[5], lo[3])}
out["exato"] = {
    "polinomio": str(p.as_expr()),
    "gcd_p1_p2_grau": int(sp.gcd(dp, d2p).degree()),          # 0: nao ha polos criticos
    "gcd_p2_p3_grau": int(sp.gcd(d2p, d3p).degree()),         # 0: p'' livre de quadrados
    "zeros_reais_de_p2_em_R": int(n_reais_d2p),
    "c_isolamento": {"inf": R(c_lo), "sup": R(c_hi),
                     "sinal_p2_inf": int(sp.sign(d2p.eval(c_lo))),
                     "sinal_p2_sup": int(sp.sign(d2p.eval(c_hi)))},
    "encerramentos_x": {f"x{j}": {"inf": R(lo[j]), "sup": R(hi[j])} for j in range(6)},
    "J_menos": {"inf": R(Jm[0]), "sup": R(Jm[1])},
    "J_mais": {"inf": R(Jp[0]), "sup": R(Jp[1])},
    "contagens_de_Sturm": {
        "J_menos": {"p": int(p.count_roots(*Jm)), "p1": int(dp.count_roots(*Jm)),
                    "p2": int(d2p.count_roots(*Jm))},
        "J_mais": {"p": int(p.count_roots(*Jp)), "p1": int(dp.count_roots(*Jp)),
                   "p2": int(d2p.count_roots(*Jp))}},
    "desigualdades": {k: {"cota_esquerda": R(a), "cota_direita": R(b),
                          "diferenca_positiva": bool(b - a > 0)} for k, (a, b) in pares.items()},
    "isolamento_raiz_real_de_p": {"inf": R(r_lo), "sup": R(r_hi)},
    "isolamento_zeros_reais_de_p1": [{"inf": R(a), "sup": R(b)} for a, b in polos],
    "cota_dist_raiz_real": R(lo[2] - r_hi),
}
assert out["exato"]["gcd_p1_p2_grau"] == 0 and out["exato"]["gcd_p2_p3_grau"] == 0
assert out["exato"]["c_isolamento"]["sinal_p2_inf"] * out["exato"]["c_isolamento"]["sinal_p2_sup"] < 0
assert all(v["diferenca_positiva"] for v in out["exato"]["desigualdades"].values())
json.dump(out, open("certificado_mandelbrot.json", "w"), indent=1, ensure_ascii=False)
print("certificado:", out["certificado"], "| largura maxima:", out["largura_maxima_dos_encerramentos"])
print("Sturm:", out["exato"]["contagens_de_Sturm"])

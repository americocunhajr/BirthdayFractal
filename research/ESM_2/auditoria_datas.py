"""
Auditoria dos numeros citados no artigo sobre as datas. Cada afirmacao e recalculada
e comparada com o valor publicado; a saida lista OK/FALHA por afirmacao.
Usa: analise_suplementar.py (lista de datas e teste dos pontos criticos),
newton_aniversario.py, newton_arte.py, mandelbrot_orbita.json, certificado_mandelbrot.json.
"""
import json
from fractions import Fraction

import numpy as np
import sympy as sp

import analise_suplementar as S
from newton_aniversario import coeficientes, raizes_da_janela
from newton_arte import janela_padrao

z = sp.symbols("z")
rel = []


def check(nome, cond, detalhe=""):
    rel.append((nome, bool(cond), detalhe))


def poly(data):
    return sp.Poly([int(v) for v in coeficientes(data)], z)


# 1. lista de datas
datas = [d for d, _ in S.DATAS]
check("38 datas, todas distintas", len(datas) == 38 and len(set(datas)) == 38, f"{len(datas)} datas")

# 2. construcao dos polinomios: digitos de DDMMAAAA sem zeros iniciais
ok = all([int(ch) for ch in d.replace("/", "")] [len(d.replace("/", "")) - len(coeficientes(d)):]
         == [int(v) for v in coeficientes(d)] for d in datas)
check("polinomios = digitos de DDMMAAAA sem zeros iniciais", ok)
pub = {"25/12/1642": 2*z**7+5*z**6+z**5+2*z**4+z**3+6*z**2+4*z+2,
       "04/01/1643": 4*z**6+z**4+z**3+6*z**2+4*z+3,
       "30/04/1777": 3*z**7+4*z**4+z**3+7*z**2+7*z+7,
       "20/11/1924": 2*z**7+z**5+z**4+z**3+9*z**2+2*z+4,
       "08/08/1900": z**2*(8*z**4+8*z**2+z+9),
       "01/01/2000": z**3*(z**3+z+2)}
for d, e in pub.items():
    check(f"polinomio publicado de {d}", sp.expand(e - poly(d).as_expr()) == 0)

# 3. raizes distintas e multiplicidades (exato): mult. de 0 = zeros finais; resto livre de quadrados
for d in datas:
    P = poly(d)
    s = d.replace("/", "")
    m0 = len(s) - len(s.rstrip("0"))
    Q = sp.Poly(sp.quo(P.as_expr(), z**m0), z)
    sqf = sp.gcd(Q, Q.diff(z)).degree() == 0 and Q.eval(0) != 0
    k_exato = Q.degree() + (1 if m0 > 0 else 0)
    r = np.roots(coeficientes(d))
    k_num = len(S.distintas(r)[0])
    check(f"{d}: mult(0)={m0}, demais raizes simples, k={k_exato}", sqf and k_num == k_exato)

# 4. classificacao das raizes isoladas e heuristica de Girard
iso, nao_iso, erros = [], [], []
for d in datas:
    c = coeficientes(d)
    r = np.roots(c)
    rd = S.distintas(r)[0]
    isolada = len(raizes_da_janela(r)) < len(rd)
    razao = c[1] / c[0]
    if razao >= 3.5:
        iso.append((d, isolada))
        reais = rd[np.abs(rd.imag) < 1e-9].real
        rmax = reais[np.argmax(np.abs(reais))]
        erros.append(abs(rmax + razao) / abs(rmax))
    else:
        nao_iso.append((d, isolada, razao))
check("14 datas com c1/c0 >= 3,5, todas com raiz isolada",
      len(iso) == 14 and all(i for _, i in iso), f"{len(iso)} datas")
check("nenhuma data com c1/c0 <= 3 tem raiz isolada",
      all((not i) for _, i, rz in nao_iso) and all(rz <= 3 for _, _, rz in nao_iso),
      f"max c1/c0 entre as demais = {max(rz for _, _, rz in nao_iso)}")
check("erros relativos entre 1% e 10%", min(erros) >= 0.01 and max(erros) <= 0.10,
      f"min {min(erros):.4f}, max {max(erros):.4f}")
dias = sorted({int(d[:2]) for d, _ in iso})
check("dias das 14 datas entre 14-19 ou 27-29", all(14 <= x <= 19 or 27 <= x <= 29 for x in dias), str(dias))


def raiz_real(d):
    r = np.roots(coeficientes(d))
    re = r[np.abs(r.imag) < 1e-9].real
    return re[np.argmax(np.abs(re))]


check("Lorenz: previsao -4,5, raiz -4,44", abs(raiz_real("29/12/1972") + 4.44) < 0.005)
check("Euler: previsao -5, raiz -5,15", abs(raiz_real("15/04/1707") + 5.15) < 0.005)
check("Netuno: previsao -1,5, raiz -2,38, nao isolada", abs(raiz_real("23/09/1846") + 2.38) < 0.005
      and len(raizes_da_janela(np.roots(coeficientes("23/09/1846")))) == 7)
check("15/08/1950: raiz isolada -5,29", abs(raiz_real("15/08/1950") + 5.29) < 0.005)
check("03/07/1986 (Secao 10): grau 6 e nenhuma raiz real (Sturm exato)",
      poly("03/07/1986").degree() == 6 and poly("03/07/1986").count_roots() == 0)
check("15/08/1950 (Secao 10): raiz simples em z = 0",
      poly("15/08/1950").eval(0) == 0 and poly("15/08/1950").diff(z).eval(0) != 0)
rg = np.sort_complex(np.round(np.roots(coeficientes("30/04/1777")), 3))
check("raizes de Gauss publicadas", set(np.round(rg, 3)) == {-1.302, 1.102 + 0.838j, 1.102 - 0.838j,
      0.111 + 1.132j, 0.111 - 1.132j, -0.562 + 0.638j, -0.562 - 0.638j}, str(rg))

# 5. pontos criticos: classificacao exata e teste numerico dos livres regulares
res = {d: S.teste_criticos(coeficientes(d)) for d in datas}
check("contagem 2k-2 dos pontos criticos em todas as datas",
      all(res[d][0]["total"] == 2 * res[d][0]["k"] - 2 for d in datas))
com_polo = [d for d in datas if res[d][0]["polos_criticos"]]
check("polo critico so em 11/02/2016", com_polo == ["11/02/2016"], str(com_polo))
P16 = poly("11/02/2016")
pc16 = res["11/02/2016"][0]["polos_criticos"]
check("11/02/2016: p' = (z+1)^2(z^2-z+1)(7z^2-z+1), p(-1)=5, p'(-1)=p''(-1)=0, polo duplo em -1",
      sp.expand(P16.diff(z).as_expr() - (z+1)**2*(z**2-z+1)*(7*z**2-z+1)) == 0
      and P16.eval(-1) == 5 and P16.diff(z).eval(-1) == 0 and P16.diff(z).diff(z).eval(-1) == 0
      and len(pc16) == 1 and pc16[0]["ponto"] == "-1" and pc16[0]["ordem_do_polo"] == 2)
conv = [d for d in datas if all(it is not None for _, it, _ in res[d][1])]
falha = [d for d in datas if d not in conv]
itmax = max(max(it for _, it, _ in res[d][1]) for d in conv)
check("37 datas com todos os pontos criticos livres regulares convergindo", len(conv) == 37,
      f"excecao: {falha}")
check("maximo de iteracoes (livres regulares) = 147", itmax == 147,
      f"{itmax} (em {[d for d in conv if max(it for _, it, _ in res[d][1]) == itmax]})")
lm = res["20/11/1924"][1]
check("Mandelbrot: 5 livres regulares, 4 complexos convergem, o real nao; sem polos criticos",
      len(lm) == 5 and sum(1 for z0, it, _ in lm if abs(z0.imag) > 1e-9 and it) == 4
      and sum(1 for z0, it, _ in lm if abs(z0.imag) < 1e-9 and it is None) == 1
      and not res["20/11/1924"][0]["polos_criticos"])
# 01/01/2000: k=4, 2k-2 = 6 pontos criticos: 3 raizes simples + 3 zeros nao nulos de p''
P = poly("01/01/2000")
d2 = sp.factor(P.diff(z).diff(z).as_expr())
check("01/01/2000: p'' = 6z(5z^3+2z+2), 3 raizes simples e 3 livres regulares",
      sp.expand(d2 - 6*z*(5*z**3+2*z+2)) == 0 and len(res["01/01/2000"][0]["livres_regulares"]) == 3
      and res["01/01/2000"][0]["raizes_simples"] == 3)

# 6. indicador de Lyapunov e certificado
M = json.load(open("mandelbrot_orbita.json"))
lam = {(e["precisao"], e["x0"], e["n"]): e["lambda"] for e in M["execucoes"]}
vals = list(lam.values())
check("lambda_n (dupla) ~ 0,19 / 0,18 / 0,19 para n = 1e4/1e5/1e6",
      round(lam[("float64", "c", 10**4)], 2) == 0.19 and round(lam[("float64", "c", 10**5)], 2) == 0.18
      and round(lam[("float64", "c", 10**6)], 2) == 0.19)
check("todas as execucoes entre 0,182 e 0,192", min(vals) >= 0.182 and max(vals) <= 0.192,
      f"min {min(vals):.4f}, max {max(vals):.4f}")
check("intervalo visitado igual em todas as execucoes",
      len({(round(e["min"], 4), round(e["max"], 4)) for e in M["execucoes"]}) == 1)
C = json.load(open("certificado_mandelbrot.json"))
check("confinamento certificado (Sturm + intervalos)", C["certificado"])
check("distancia de I-uI+ a raiz real > 0,55 e aos polos reais > 0,21",
      C["dist_min_a_raiz_real"] > 0.55 and C["dist_min_aos_polos_reais"] > 0.21)

# 7. tabela de raiz de 2 e constante 1/(2 sqrt 2)
xs = [Fraction(1)]
for k in range(4):
    xs.append(xs[-1] - (xs[-1] ** 2 - 2) / (2 * xs[-1]))
check("frações da Tabela 1", [str(v) for v in xs] == ["1", "3/2", "17/12", "577/408", "665857/470832"])
import mpmath as mp
mp.mp.dps = 50
e3 = abs(mp.mpf(577) / 408 - mp.sqrt(2)); e4 = abs(mp.mpf(665857) / 470832 - mp.sqrt(2))
check("0,354 x (2,1e-6)^2 ~ 1,6e-12", abs(e3 ** 2 / (2 * mp.sqrt(2)) - e4) / e4 < 1e-5,
      f"e3={mp.nstr(e3, 3)}, e4={mp.nstr(e4, 3)}")

# 8. Figura 4: janela de Gauss e polo
cen, R = janela_padrao(np.array(coeficientes("30/04/1777"), float))
check("Figura 4: centro 0 e meia largura 1,800 / 0,300 / 0,0375",
      abs(cen) < 1e-12 and round(R, 3) == 1.8 and round(R / 6, 3) == 0.3 and round(R / 48, 4) == 0.0375)

# 9. revisao em ingles: Proposicoes 1-3 e numeros associados
import calendar, math
r = Fraction(5, 2)
Sr = Fraction(9) / r**2 + Fraction(2) / r**3 + Fraction(9) / r**4 + Fraction(9) / r**5 + Fraction(9) / r**6
check("Prop. 1: S(5/2) = 1,927424", Sr == Fraction(1927424, 1000000), str(float(Sr)))
check("Prop. 1: c0 r + S < c1 nos dois casos", 1 * r + Sr < 5 and 2 * r + Sr < 7)
dias = {15, 16, 17, 18, 19, 27, 28, 29}
n_tot = n_bad = 0
for ano in range(1000, 3000):
    for mes in range(1, 13):
        for dia in dias:
            if dia > calendar.monthrange(ano, mes)[1]:
                continue
            cc = np.array(coeficientes(f"{dia:02d}/{mes:02d}/{ano}"), float)
            rr = np.roots(cc); fora = rr[np.abs(rr) > 2.5]
            ok1 = (len(fora) == 1 and abs(fora[0].imag) < 1e-9 and fora[0].real < 0
                   and abs(fora[0].real + cc[1] / cc[0]) < 2 / cc[0])
            n_tot += 1; n_bad += (not ok1)
check("Prop. 1 conferida numericamente em todas as datas validas da familia", n_bad == 0,
      f"{n_tot} datas, {n_bad} violacoes")
grandes = [d for d in datas if np.array(coeficientes(d), float)[1] / np.array(coeficientes(d), float)[0] >= 3.5]
isol = {d: len(raizes_da_janela(np.roots(np.array(coeficientes(d), float))))
        < len(S.distintas(np.roots(np.array(coeficientes(d), float)))[0]) for d in datas}
cobertas = [d for d in datas if int(d[:2]) in dias]
check("14 datas com c1/c0 >= 3,5, todas isoladas; 13 cobertas pela Prop. 1; excecao 14/09/2015",
      len(grandes) == 14 and all(isol[d] for d in grandes) and len(cobertas) == 13
      and set(grandes) - set(cobertas) == {"14/09/2015"})
def raiz_isolada(d):
    rr = np.roots(np.array(coeficientes(d), float))
    return rr[np.argmax(np.abs(rr))].real
pub_r = {"15/04/1707": -5.15, "16/08/1821": -6.21, "29/04/1854": -4.60, "28/02/1878": -4.07,
         "29/05/1919": -4.62, "29/12/1972": -4.44, "16/05/1980": -6.13, "27/04/1952": -3.66,
         "15/08/1950": -5.29, "14/09/2015": -4.43}
fora_jan = True
for d, v in pub_r.items():
    cen, R = janela_padrao(np.array(coeficientes(d), float))
    x = raiz_isolada(d)
    fora_jan &= abs(x - cen.real) > R and round(x, 2) == v
check("raizes isoladas citadas nas legendas e no texto, todas fora da janela", fora_jan)
cen, R = janela_padrao(np.array(coeficientes("23/09/1846"), float))
x = min(r_.real for r_ in np.roots(np.array(coeficientes("23/09/1846"), float)) if abs(r_.imag) < 1e-9)
check("Netuno: raiz real -2,38 dentro da janela", round(x, 2) == -2.38 and abs(x - cen.real) < R)
check("Prop. 2: 34 passos (m=2) e 57 passos (m=3) de 1 a 1e-10",
      math.floor(math.log(1e10) / math.log(2)) + 1 == 34 and math.floor(math.log(1e10) / math.log(1.5)) + 1 == 57)
sem_reais = [d for d in datas if len(sp.real_roots(poly(d))) == 0]
check("seis datas sem raizes reais, todas aprovadas no teste dos pontos criticos",
      len(sem_reais) == 6 and all(d in conv for d in sem_reais), ", ".join(sem_reais))
A = json.load(open("certificado_aprisionamento.json"))
check("Prop. 3: certificado exato (K-, K+)", A["certificado"] and A["K_menos"] == ["-3813/5000", "-78/125"]
      and A["K_mais"] == ["79/500", "503/2500"], str(A["K_menos"]) + " " + str(A["K_mais"]))
check("Prop. 3: distancias > 0,55 (raiz real) e > 0,21 (polos reais)",
      A["dist_min_K_raiz_real"] > 0.55 and A["dist_min_K_polos_reais"] > 0.21)
check("Prop. 3: N(alpha) ~ 0,20105 e N(gamma) ~ -0,76248; polos -0,9730 e -0,1130; raiz -1,3169",
      round(A["N(alpha)_aprox"], 5) == 0.20105 and round(A["N(gamma)_aprox"], 5) == -0.76248
      and [round(v, 4) for v in A["polos_reais"]] == [-0.973, -0.113] and round(A["raiz_real"], 4) == -1.3169)
V = A["verificacoes"]
check("Prop. 3: distancias certificadas por contagens de Sturm (0,55 e 0,21)",
      all(V[k] for k in V if "a menos de" in k) and sum("a menos de" in k for k in V) == 4)
dc = A["dois_ciclo"]
mi = [float(v) for v in dc["multiplicador_intervalo"].strip("[]").split(",")]
check("Corolario 1: 2-ciclo unico em K-, x* ~ -0,7269194 <-> 0,1681738, multiplicador em (-1,90446; -1,90445)",
      V["N o N com um unico ponto fixo em K- (2-ciclo unico)"] and V["2-ciclo repulsor: (N o N)'(x*) < -1"]
      and round(dc["x_estrela_aprox"], 7) == -0.7269194 and round(dc["N(x_estrela)_aprox"], 7) == 0.1681738
      and -1.90446 < mi[0] and mi[1] < -1.90445 and dc["grau_de_F"] == 49,
      f"{dc['x_estrela_aprox']:.7f} <-> {dc['N(x_estrela)_aprox']:.7f}, mult {mi[0]:.6f}")
cen, R = janela_padrao(np.array(coeficientes("20/11/1924"), float))
check("linhas da grade mais proximas do eixo real a ~6e-4", 5e-4 < R / 3199 < 7e-4, f"{R / 3199:.2e}")
check("erros de 1e-16 crescem a ordem 1 em ~200 iteracoes (lambda ~ 0,19)",
      150 < math.log(1e16) / 0.19 < 250, f"{math.log(1e16) / 0.19:.0f}")

# relatorio
largura = max(len(n) for n, _, _ in rel)
with open("auditoria_datas.txt", "w", encoding="utf-8") as f:
    for n, okk, det in rel:
        linha = f"{'OK   ' if okk else 'FALHA'}  {n.ljust(largura)}  {det}"
        print(linha)
        f.write(linha + "\n")
print(f"\n{sum(o for _, o, _ in rel)}/{len(rel)} afirmacoes conferem")

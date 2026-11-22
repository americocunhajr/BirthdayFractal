"""
Material suplementar: para cada data usada no artigo,
  - polinomio, grau d e numero k de raizes distintas;
  - janela das figuras (centro e meia largura, regra da raiz isolada);
  - heuristica de Girard: -c1/c0 contra a raiz real de maior modulo;
  - classificacao EXATA (SymPy) dos pontos criticos de N_p na esfera de Riemann:
      raizes simples de p (pontos fixos superatratores);
      polos criticos: zeros de p' de multiplicidade mu >= 2 que nao sao raizes de p
        (polos de ordem mu, grau local mu; identificados por gcd(p', p''));
      pontos criticos livres regulares: zeros de p'' que nao sao zeros de p nem de p'
        (neles N_p' = 0); a soma, com multiplicidade, deve ser 2k - 2;
  - teste numerico dos pontos criticos livres regulares: iteracao de N_p em precisao
    dupla, a partir de aproximacoes com 30 digitos (sympy nroots), ate |passo| <
    1e-14 * escala e distancia a uma raiz < 1e-9, no maximo 10^4 iteracoes.
    Os polos criticos NAO sao iterados numericamente: sua orbita exata e polo -> oo -> oo.
Saidas: suplemento_datas.csv e tabelas LaTeX em suplemento_tabelas.tex.
"""
import csv
import numpy as np
import sympy as sp

from newton_aniversario import coeficientes
from newton_arte import janela_padrao
from newton_aniversario import raizes_da_janela

DATAS = [
    ("25/12/1642", "Nascimento de Isaac Newton (juliano)"),
    ("04/01/1643", "Nascimento de Isaac Newton (gregoriano)"),
    ("15/04/1707", "Nascimento de Leonhard Euler"),
    ("30/04/1777", "Nascimento de Carl Friedrich Gauss"),
    ("16/08/1821", "Nascimento de Arthur Cayley"),
    ("29/04/1854", "Nascimento de Henri Poincaré"),
    ("28/02/1878", "Nascimento de Pierre Fatou"),
    ("03/02/1893", "Nascimento de Gaston Julia"),
    ("20/11/1924", "Nascimento de Benoit Mandelbrot"),
    ("05/07/1687", "Publicação dos Principia (juliano)"),
    ("15/07/1687", "Publicação dos Principia (gregoriano)"),
    ("07/09/1822", "Independência do Brasil"),
    ("23/09/1846", "Galle observa Netuno"),
    ("13/05/1888", "Lei Áurea"),
    ("21/01/1889", "Anúncio do prêmio do rei Oscar II"),
    ("15/11/1889", "Proclamação da República"),
    ("30/11/1889", "Poincaré pede a interrupção da impressão"),
    ("08/08/1900", "Problemas de Hilbert"),
    ("23/10/1906", "Voo do 14-bis"),
    ("29/05/1919", "Eclipse de Sobral"),
    ("15/10/1952", "Fundação do IMPA"),
    ("29/12/1972", "Palestra de Lorenz"),
    ("16/05/1980", "Fundação do LNCC"),
    ("05/10/1988", "Constituição Federal"),
    ("23/06/1993", "Anúncio de Wiles"),
    ("22/04/1500", "Chegada de Cabral (juliano)"),
    ("02/05/1500", "Chegada de Cabral (gregoriano proléptico)"),
    ("04/10/1582", "Último dia juliano"),
    ("15/10/1582", "Primeiro dia gregoriano"),
    ("20/07/1969", "Apollo 11 (Brasília)"),
    ("21/07/1969", "Apollo 11 (UTC)"),
    ("14/09/2015", "Detecção de ondas gravitacionais"),
    ("11/02/2016", "Anúncio da detecção"),
    ("01/01/2000", "Virada do milênio"),
    ("30/08/1984", "Data pessoal (nascimento do autor)"),
    ("03/07/1986", "Data pessoal"),
    ("27/04/1952", "Data pessoal"),
    ("15/08/1950", "Data pessoal"),
]


def polinomio_tex(c):
    c = [int(v) for v in c]
    n = len(c) - 1
    termos = []
    for i, ci in enumerate(c):
        p = n - i
        if ci == 0:
            continue
        coef = "" if (ci == 1 and p > 0) else str(ci)
        var = "" if p == 0 else ("z" if p == 1 else f"z^{p}")
        termos.append(coef + var)
    return "+".join(termos)


def distintas(r, tol=1e-6):
    out, mult = [], []
    for z in r:
        for i, w in enumerate(out):
            if abs(z - w) < tol:
                mult[i] += 1
                break
        else:
            out.append(z); mult.append(1)
    return np.array(out), np.array(mult)


_z = sp.symbols("z")


def _remove_fatores(q, alvo):
    """Remove de q todos os fatores irredutiveis comuns com alvo (com multiplicidade)."""
    while True:
        h = sp.gcd(q, alvo)
        if h.degree() <= 0:
            return q
        q = sp.Poly(sp.quo(q, h), _z)


def pontos_criticos(c):
    """Classificacao exata dos pontos criticos de N_p (p com coeficientes inteiros c)."""
    P = sp.Poly([int(v) for v in c], _z)
    d1, d2 = P.diff(_z), P.diff(_z).diff(_z)
    _, fat = sp.factor_list(P.as_expr())
    simples = [f for f, e in fat if e == 1]
    k = sum(sp.Poly(f, _z).degree() for f, e in fat)
    s_ = sum(sp.Poly(f, _z).degree() for f in simples)
    Ps = sp.Poly(sp.Mul(*simples), _z) if simples else sp.Poly(1, _z)
    # multiplicidade extra de pontos criticos nas raizes simples (zeros de p'' em raizes simples)
    extra_simples = d2.degree() - _remove_fatores(d2, Ps).degree()
    # polos criticos: fatores de gcd(p', p'') que nao dividem p
    g = sp.gcd(d1, d2)
    polos = []
    if g.degree() > 0:
        for f, _e in sp.factor_list(g.as_expr())[1]:
            fp = sp.Poly(f, _z)
            if sp.rem(P, fp).is_zero:
                continue                        # zero de p: singularidade removivel
            mu = 0
            q = d1
            while sp.rem(q, fp).is_zero:
                q = sp.Poly(sp.quo(q, fp), _z); mu += 1
            for raiz in sp.Poly(f, _z).all_roots():
                polos.append({"ponto": str(raiz), "valor": complex(sp.N(raiz, 30)),
                              "ordem_do_polo": mu, "multiplicidade_critica": mu - 1})
    mult_polos = sum(pp["multiplicidade_critica"] for pp in polos)
    # livres regulares: zeros de p'' que nao sao zeros de p nem de p'
    q = _remove_fatores(d2, sp.Poly(P.as_expr() * d1.as_expr(), _z))
    livres = [complex(r) for r in q.nroots(n=30)] if q.degree() > 0 else []
    total = s_ + extra_simples + mult_polos + q.degree()
    assert total == 2 * k - 2, (c, total, k)
    return {"k": k, "raizes_simples": s_, "extra_simples": extra_simples,
            "polos_criticos": polos, "livres_regulares": livres, "total": total}


def teste_criticos(c, maxit=10_000):
    """Classificacao exata + iteracao numerica dos pontos criticos livres regulares."""
    pc = pontos_criticos(c)
    p, dp = np.poly1d(c), np.poly1d(c).deriv()
    rd, mult = distintas(np.roots(c))
    escala = max(1.0, np.max(np.abs(rd)))
    res = []
    for z0 in pc["livres_regulares"]:
        z = complex(z0)
        it = None
        for k in range(1, maxit + 1):
            d = dp(z)
            if d == 0:
                break
            passo = p(z) / d
            z -= passo
            if abs(passo) < 1e-14 * escala and np.min(np.abs(rd - z)) < 1e-9:
                it = k
                break
        res.append((complex(z0), it, complex(rd[np.argmin(np.abs(rd - z))]) if it else None))
    return pc, res


def fmt_c(z, nd=3):
    z = complex(z)
    if abs(z.imag) < 10 ** (-nd):
        return f"{z.real:.{nd}f}".replace(".", ",")
    s = "+" if z.imag >= 0 else "-"
    return f"{z.real:.{nd}f}{s}{abs(z.imag):.{nd}f}i".replace(".", ",")


def tnum(v, nd=3):
    """numero real em modo matematico, virgula decimal"""
    if v is None:
        return "---"
    s = f"{v:.{nd}f}"
    if s.lstrip("-").strip("0.") == "":
        s = s.lstrip("-")
    return "$" + s.replace(".", "{,}") + "$"


def tcomp(z, nd=3):
    z = complex(z)
    re_, im_ = round(z.real, nd), round(z.imag, nd)
    re_ = 0.0 if re_ == 0 else re_
    if im_ == 0:
        return tnum(re_, nd)
    sinal = "+" if im_ > 0 else "-"
    s = f"{re_:.{nd}f}{sinal}{abs(im_):.{nd}f}i"
    return "$" + s.replace(".", "{,}") + "$"


def calcular_linhas():
    linhas = []
    for data, oc in DATAS:
        c = np.array(coeficientes(data), float)
        r = np.roots(c)
        rd, mult = distintas(r)
        cen, R = janela_padrao(c)
        isolada = len(raizes_da_janela(r)) < len(rd)
        reais = rd[np.abs(rd.imag) < 1e-9].real
        rmax = reais[np.argmax(np.abs(reais))] if len(reais) else None
        razao = c[1] / c[0]
        erro = abs(rmax - (-razao)) / abs(rmax) if (rmax is not None and razao > 0) else None
        pc, res = teste_criticos(c)
        conv = all(it is not None for _, it, _ in res)
        linhas.append(dict(
            data=data, ocasiao=oc, grau=len(c) - 1, raizes_distintas=len(rd),
            mult_zero=len(data.replace("/", "")) - len(data.replace("/", "").rstrip("0")),
            mult_max=int(mult.max()), polinomio=polinomio_tex(c),
            centro=round(cen.real, 4), meia_largura=round(R, 4), raiz_fora_da_janela=isolada,
            c1_c0=round(razao, 3), previsao=round(-razao, 3),
            raiz_real_max=None if rmax is None else round(float(rmax), 4),
            erro_rel=None if erro is None else round(erro, 3),
            pontos_criticos_total=pc["total"],
            n_livres_regulares=len(pc["livres_regulares"]),
            n_polos_criticos=len(pc["polos_criticos"]),
            polos_criticos="; ".join(f"{pp['ponto']} (polo de ordem {pp['ordem_do_polo']})"
                                     for pp in pc["polos_criticos"]),
            livres_regulares_convergem=conv,
            iter_max=max((it or 0) for _, it, _ in res) if res else 0,
            _res=res, _polos=pc["polos_criticos"],
            detalhe="; ".join(f"{fmt_c(z0)} -> {fmt_c(rz) if rz is not None else 'nao converge'}"
                              f" ({it if it else '>10^4'} it.)" for z0, it, rz in res)))
    return linhas


def escrever(linhas):
    with open("suplemento_tabelas.tex", "w", encoding="utf-8") as f:
        f.write("% gerado por analise_suplementar.py\n")
        f.write("\\begin{longtable}{@{}l p{5.0cm} l c c c@{}}\n"
                "\\caption{Datas, polinômios, grau $d$, número $k$ de raízes distintas e "
                "multiplicidade $m_0$ da raiz $z=0$ (0 se $z=0$ não é raiz). As demais raízes "
                "são simples (verificação exata).}"
                "\\label{tab:s1}\\\\\n\\toprule\n"
                "Data & Ocasião & $p(z)$ & $d$ & $k$ & $m_0$\\\\\n\\midrule\\endfirsthead\n"
                "\\toprule Data & Ocasião & $p(z)$ & $d$ & $k$ & $m_0$\\\\\\midrule\\endhead\n")
        for L in linhas:
            f.write(f"{L['data']} & {L['ocasiao']} & ${L['polinomio']}$ & {L['grau']} & "
                    f"{L['raizes_distintas']} & {L['mult_zero']}\\\\\n")
        f.write("\\bottomrule\n\\end{longtable}\n\n")

        f.write("\\begin{longtable}{@{}l r r c r r r@{}}\n"
                "\\caption{Janela (centro real e meia largura $R$), raiz isolada fora da janela, "
                "razão $c_1/c_0$, raiz real de maior módulo e erro relativo da previsão "
                "$-c_1/c_0$.}\\label{tab:s2}\\\\\n\\toprule\n"
                "Data & centro & $R$ & isolada & $c_1/c_0$ & raiz real & erro rel.\\\\\n"
                "\\midrule\\endfirsthead\n"
                "\\toprule Data & centro & $R$ & isolada & $c_1/c_0$ & raiz real & erro rel.\\\\"
                "\\midrule\\endhead\n")
        for L in linhas:
            f.write(f"{L['data']} & {tnum(L['centro'])} & {tnum(L['meia_largura'])} & "
                    f"{'sim' if L['raiz_fora_da_janela'] else 'não'} & {tnum(L['c1_c0'], 2)} & "
                    f"{tnum(L['raiz_real_max'])} & {tnum(L['erro_rel'])}\\\\\n")
        f.write("\\bottomrule\n\\end{longtable}\n\n")

        f.write("\\begin{longtable}{@{}l p{12.3cm}@{}}\n"
                "\\caption{Pontos críticos livres de $N_p$. Para os regulares $z_0$ (zeros de "
                "$p''$ que não são zeros de $p$ nem de $p'$), a raiz alcançada e o número de "
                "iterações, com o critério: passo menor que $10^{-14}\\max(1,\\max|\\zeta|)$ e "
                "distância a uma raiz menor que $10^{-9}$, no máximo $10^4$ iterações. Os polos "
                "críticos (zeros múltiplos de $p'$ que não são raízes de $p$) foram identificados "
                "em aritmética exata e não foram iterados numericamente.}\\label{tab:s3}\\\\\n"
                "\\toprule\n"
                "Data & $z_0\\to$ raiz alcançada (iterações)\\\\\n\\midrule\\endfirsthead\n"
                "\\toprule Data & $z_0\\to$ raiz alcançada (iterações)\\\\\\midrule\\endhead\n")
        for L in linhas:
            partes = []
            for z0, it, rz in L["_res"]:
                alvo = tcomp(rz) if rz is not None else "não converge"
                n = str(it) if it else "$>10^4$"
                partes.append(f"{tcomp(z0)} $\\to$ {alvo} ({n})")
            for pp in L["_polos"]:
                partes.append(f"polo crítico ${sp.latex(sp.sympify(pp['ponto']))}$, de ordem "
                              f"{pp['ordem_do_polo']}: órbita exata "
                              f"${sp.latex(sp.sympify(pp['ponto']))}\\mapsto\\infty\\mapsto\\infty$")
            f.write(f"{L['data']} & {{\\footnotesize " + "; ".join(partes) + "}\\\\\n")
        f.write("\\bottomrule\n\\end{longtable}\n")

    campos = [k for k in linhas[0] if not k.startswith("_")]
    with open("suplemento_datas.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader(); w.writerows(linhas)


if __name__ == "__main__":
    linhas = calcular_linhas()
    escrever(linhas)
    for L in linhas:
        print(L["data"], "k=", L["raizes_distintas"], "crit=", L["pontos_criticos_total"],
              "livres_reg=", L["n_livres_regulares"], "polos_crit=", L["n_polos_criticos"],
              "converge=", L["livres_regulares_convergem"], "itmax=", L["iter_max"])

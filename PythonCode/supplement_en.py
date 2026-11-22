"""
English Tables S1-S6 (Online Resource 1 of the article) and an English per-date CSV.

Tables S1-S3 and the CSV come from the same computation as analise_suplementar.py
(exact classification of critical points with SymPy, numerical iteration of the
regular free critical points). Tables S4-S6 are read from mandelbrot_orbita.json,
certificado_mandelbrot.json and auditoria_figuras.json, written by
analise_mandelbrot.py, certifica_mandelbrot.py and auditoria_figuras.py.

Output: supplementary_dates.csv; the LaTeX code of Tables S1-S6 is written into
        esm/ESM_1.tex (source of Online Resource 1), between the lines
        "% BEGIN TABLES S1-S6" and "% END TABLES S1-S6", or, if that file is not
        present, into tables_s1_s6.tex.
Usage:  python supplement_en.py [ESM_1.tex]   (default: esm/ESM_1.tex)
"""
import csv
import json
import math
import os
import sys

import sympy as sp

from analise_suplementar import calcular_linhas

MANUSCRIPT = "esm/ESM_1.tex"

OCCASION = {
    "25/12/1642": "Birth of Isaac Newton (Julian)",
    "04/01/1643": "Birth of Isaac Newton (Gregorian)",
    "15/04/1707": "Birth of Leonhard Euler",
    "30/04/1777": "Birth of Carl Friedrich Gauss",
    "16/08/1821": "Birth of Arthur Cayley",
    "29/04/1854": "Birth of Henri Poincar\u00e9",
    "28/02/1878": "Birth of Pierre Fatou",
    "03/02/1893": "Birth of Gaston Julia",
    "20/11/1924": "Birth of Benoit Mandelbrot",
    "05/07/1687": "Publication of the Principia (Julian)",
    "15/07/1687": "Publication of the Principia (Gregorian)",
    "07/09/1822": "Independence of Brazil",
    "23/09/1846": "Galle observes Neptune",
    "13/05/1888": "Golden Law (abolition of slavery in Brazil)",
    "21/01/1889": "Announcement of King Oscar II's prize",
    "15/11/1889": "Proclamation of the Republic of Brazil",
    "30/11/1889": "Poincar\u00e9 asks to stop the printing",
    "08/08/1900": "Hilbert's problems",
    "23/10/1906": "Flight of the 14-bis",
    "29/05/1919": "Sobral eclipse",
    "15/10/1952": "Founding of IMPA",
    "29/12/1972": "Lorenz's talk",
    "16/05/1980": "Founding of LNCC",
    "05/10/1988": "Brazilian Federal Constitution",
    "23/06/1993": "Wiles's announcement",
    "22/04/1500": "Arrival of Cabral (Julian)",
    "02/05/1500": "Arrival of Cabral (proleptic Gregorian)",
    "04/10/1582": "Last Julian day",
    "15/10/1582": "First Gregorian day",
    "20/07/1969": "Apollo 11 (Bras\u00edlia time)",
    "21/07/1969": "Apollo 11 (UTC)",
    "14/09/2015": "Detection of gravitational waves",
    "11/02/2016": "Announcement of the detection",
    "01/01/2000": "New Year's Day 2000",
    "30/08/1984": "Personal date (author's birth)",
    "03/07/1986": "Personal date",
    "27/04/1952": "Personal date",
    "15/08/1950": "Personal date",
}

# figure panels of the article, for Table S6
FIGURE = {
    "newton_25-12-1642": "25/12/1642 (Fig. 1a)", "newton_04-01-1643": "04/01/1643 (Fig. 1b)",
    "newton_30-04-1777": "30/04/1777 (Figs. 4a, 5)",
    "newton_15-04-1707": "15/04/1707 (Fig. 7a)", "newton_16-08-1821": "16/08/1821 (Fig. 7b)",
    "newton_29-04-1854": "29/04/1854 (Fig. 7c)", "newton_28-02-1878": "28/02/1878 (Fig. 7d)",
    "newton_03-02-1893": "03/02/1893 (Fig. 7e)", "newton_20-11-1924": "20/11/1924 (Fig. 7f)",
    "newton_05-07-1687": "05/07/1687 (Fig. 8a)", "newton_23-09-1846": "23/09/1846 (Fig. 8b)",
    "newton_30-11-1889": "30/11/1889 (Fig. 8c)", "newton_29-05-1919": "29/05/1919 (Fig. 8d)",
    "newton_29-12-1972": "29/12/1972 (Fig. 8e)", "newton_16-05-1980": "16/05/1980 (Fig. 8f)",
    "newton_08-08-1900": "08/08/1900 (Fig. 9)",
    "newton_30-08-1984": "30/08/1984 (Fig. 10a)", "newton_03-07-1986": "03/07/1986 (Fig. 10b)",
    "newton_27-04-1952": "27/04/1952 (Fig. 10c)", "newton_15-08-1950": "15/08/1950 (Fig. 10d)",
    "cayley_z2-1": "$z^2-1$ (Fig. 3a)", "cayley_z3-1": "$z^3-1$ (Fig. 3b)",
    "falha_z3-2z+2": "$z^3-2z+2$ (Fig. 11)",
    "zoom_2": "Gauss, magnification (Fig. 4b)", "zoom_3": "Gauss, magnification (Fig. 4c)",
}


def tnum(v, nd=3):
    """real number in math mode, decimal point"""
    if v is None:
        return "---"
    s = f"{v:.{nd}f}"
    if s.lstrip("-").strip("0.") == "":
        s = s.lstrip("-")
    return "$" + s + "$"


def ccomp(z, nd=3):
    """complex number as plain text, decimal point"""
    z = complex(z)
    re_, im_ = round(z.real, nd), round(z.imag, nd)
    re_ = 0.0 if re_ == 0 else re_
    if im_ == 0:
        return f"{re_:.{nd}f}"
    return f"{re_:.{nd}f}{'+' if im_ > 0 else '-'}{abs(im_):.{nd}f}i"


def tcomp(z, nd=3):
    return "$" + ccomp(z, nd) + "$"


def thousands(x):
    return f"{x:,}"


def tables_s1_s3(lines, w_occasion="3.4cm", w_detail="11.2cm"):
    out = []
    out.append("\\begin{longtable}{@{}l >{\\raggedright\\arraybackslash}p{" + w_occasion
               + "} l c c c@{}}\n"
               "\\caption{Dates, polynomials, degree $d$, number $k$ of distinct roots, and "
               "multiplicity $m_0$ of the root $z=0$ (0 if $z=0$ is not a root). All other roots "
               "are simple (exact verification)}\\label{tab:s1}\\\\\n\\toprule\n"
               "Date & Occasion & $p(z)$ & $d$ & $k$ & $m_0$\\\\\n\\midrule\\endfirsthead\n"
               "\\toprule Date & Occasion & $p(z)$ & $d$ & $k$ & $m_0$\\\\\\midrule\\endhead\n")
    for L in lines:
        out.append(f"{L['data']} & {OCCASION[L['data']]} & ${L['polinomio']}$ & {L['grau']} & "
                   f"{L['raizes_distintas']} & {L['mult_zero']}\\\\\n")
    out.append("\\bottomrule\n\\end{longtable}\n\n")

    out.append("\\begin{longtable}{@{}l r r c r r r@{}}\n"
               "\\caption{Window (real center and half-width $R$), isolated root outside the "
               "window, ratio $c_1/c_0$, real root of largest modulus, and relative error of the "
               "prediction $-c_1/c_0$}\\label{tab:s2}\\\\\n\\toprule\n"
               "Date & center & $R$ & isolated & $c_1/c_0$ & real root & rel.\\ error\\\\\n"
               "\\midrule\\endfirsthead\n"
               "\\toprule Date & center & $R$ & isolated & $c_1/c_0$ & real root & rel.\\ error\\\\"
               "\\midrule\\endhead\n")
    for L in lines:
        out.append(f"{L['data']} & {tnum(L['centro'])} & {tnum(L['meia_largura'])} & "
                   f"{'yes' if L['raiz_fora_da_janela'] else 'no'} & {tnum(L['c1_c0'], 2)} & "
                   f"{tnum(L['raiz_real_max'])} & {tnum(L['erro_rel'])}\\\\\n")
    out.append("\\bottomrule\n\\end{longtable}\n\n")

    out.append("\\begin{longtable}{@{}l >{\\raggedright\\arraybackslash}p{" + w_detail + "}@{}}\n"
               "\\caption{Free critical points of $N_p$. For each regular one, $z_0$ (a zero of "
               "$p''$ that is not a zero of $p$ or of $p'$), the root reached and the number of "
               "iterations, with the criterion: step smaller than $10^{-14}\\max(1,\\max|\\zeta|)$ "
               "and distance to a root smaller than $10^{-9}$, at most $10^4$ iterations. Critical "
               "poles (multiple zeros of $p'$ that are not roots of $p$) were identified in exact "
               "arithmetic and were not iterated numerically}\\label{tab:s3}\\\\\n\\toprule\n"
               "Date & $z_0\\to$ root reached (iterations)\\\\\n\\midrule\\endfirsthead\n"
               "\\toprule Date & $z_0\\to$ root reached (iterations)\\\\\\midrule\\endhead\n")
    for L in lines:
        parts = []
        for z0, it, rz in L["_res"]:
            if rz is not None:
                target = tcomp(rz)
            elif L["data"] == "20/11/1924" and abs(complex(z0).imag) < 1e-12:
                target = ("no convergence in $10^4$ iterations; nonconvergence proved "
                          "(Proposition~3 of the article)")
            else:
                target = "no convergence in $10^4$ iterations"
            if it:
                parts.append(f"{tcomp(z0)} $\\to$ {target} ({it})")
            else:
                parts.append(f"{tcomp(z0)}: {target}")
        for pp in L["_polos"]:
            pt = sp.latex(sp.sympify(pp["ponto"]))
            parts.append(f"critical pole ${pt}$, of order {pp['ordem_do_polo']}: exact orbit "
                         f"${pt}\\mapsto\\infty\\mapsto\\infty$")
        out.append(f"{L['data']} & " + "; ".join(parts) + "\\\\\n")
    out.append("\\bottomrule\n\\end{longtable}\n")
    return "".join(out)


def table_s4(M):
    def n4(v):
        return f"${v:.4f}$"

    def power(x):
        k = int(round(math.log10(x)))
        return f"$10^{{{k}}}$" if abs(x - 10 ** k) < 1e-9 else f"${x / 10 ** k:.0f}\\times10^{{{k}}}$"

    x0 = {"c": "$c$", "c+1e-08": "$c+10^{-8}$", "c-1e-08": "$c-10^{-8}$", "c+1e-04": "$c+10^{-4}$"}
    prec = {"float64": "double", "long double": "extended (\\texttt{long double})",
            "mpmath 40 digitos": "40 digits (\\texttt{mpmath})"}
    rows = [f"{prec[r['precisao']]} & {x0[r['x0']]} & {power(r['n'])} & {n4(r['lambda'])} & "
            f"[{n4(r['min'])}, {n4(r['max'])}]\\\\" for r in M["execucoes"]]
    return (r"""\begin{table}[htbp]
\centering
\caption{Orbit of the real free critical point $c\approx""" + f"{M['ponto_critico']:.6f}"
            + r"""$ of $N_p$ for 20/11/1924: indicator $\lambda_n=\frac1n\sum_{j=1}^n\log|N_p'(x_j)|$,
with $x_j=N_p^j(c)$, and interval visited by the orbit (minimum and maximum of the iterates)}
\label{tab:s4}
\begin{tabular}{@{}l l r r l@{}}
\toprule
Precision & starting point & $n$ & $\lambda_n$ & visited interval\\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
""")


def table_s5(C):
    from fractions import Fraction
    enc = C["exato"]["encerramentos_x"]

    def mid_rad(j):
        a, b = Fraction(enc[f"x{j}"]["inf"]), Fraction(enc[f"x{j}"]["sup"])
        m, rad = (a + b) / 2, (b - a) / 2
        e = int(math.floor(math.log10(float(rad))))
        return f"${float(m):.15f}$", f"${float(rad) / 10 ** e:.1f}\\times10^{{{e}}}$"

    rows = "\n".join(f"$x_{j}$ & {mid_rad(j)[0]} & {mid_rad(j)[1]}\\\\" for j in range(6))
    A4 = C["A4"]
    yes = lambda b: "yes" if b else "no"
    return (r"""\begin{table}[htbp]
\centering
\caption{Critical-orbit intervals for 20/11/1924: $I_-=[x_2,x_4]$ and $I_+=[x_1,x_3]$, the
smallest closed intervals that contain the even and the odd iterates of $c$, respectively.
Validated enclosures of $x_j=N_p^j(c)$ ($x_0=c$), given by midpoint and radius, and result of
each exact check. The exact rational endpoints, the isolating intervals, the Sturm counts, and
the bounds for each inequality are in the file \texttt{certificado\_mandelbrot.json} of
Online Resource~2}
\label{tab:s5}
\begin{tabular}{@{}c@{}}
\begin{tabular}{@{}l r r@{}}
\toprule
 & midpoint & radius\\
\midrule
""" + rows + r"""
\bottomrule
\end{tabular}\\[8pt]
\begin{tabular}{@{}l l@{}}
\toprule
Check & result\\
\midrule
$\gcd(p',p'')=1$: $N_p$ has no critical poles & """ + yes(C["exato"]["gcd_p1_p2_grau"] == 0) + r"""\\
$p''$ squarefree, a single real zero, $c$ isolated & """
            + yes(C["A1_c_zero_simples_isolado"] and C["zeros_reais_de_p2"] == 1) + r"""\\
on $J_-\supset I_-$: $p$ and $p'$ without zeros, $p''$ with one zero ($c$) & """
            + yes(C["A2_J_menos"]) + r"""\\
on $J_+\supset I_+$: $p$, $p'$, and $p''$ without zeros & """ + yes(C["A3_J_mais"]) + r"""\\
$x_2<c<x_4$, $x_1<x_3$, $x_2<x_4$ & """
            + yes(A4["x2 < c"] and A4["c < x4"] and A4["x1 < x3"] and A4["x2 < x4"]) + r"""\\
$x_5<x_3$ & """ + yes(A4["x5 < x3"]) + r"""\\
\bottomrule
\end{tabular}
\end{tabular}
\end{table}
""")


def table_s6(F):
    rows = []
    near = None
    for k, v in F.items():
        c = v["contagem"]
        if "nao_classificados_perto_do_ciclo_1e-6" in v:
            near = (v["nao_classificados_perto_do_ciclo_1e-6"], c["2"])
        rows.append(f"{FIGURE[k]} & {thousands(c['0'])} & {thousands(c['1'])} & "
                    f"{thousands(c['2'])} & {thousands(c['3'])}\\\\")
    quant = (f"all {thousands(near[1])}" if near[0] == near[1]
             else f"{thousands(near[0])} of the {thousands(near[1])}")
    note = (f"For Fig.~11, {quant} points of class (c) have their last iterate within "
            f"$10^{{-6}}$ of $0$ or $1$")
    return (r"""\begin{longtable}{@{}l r r r r@{}}
\caption{Classification of the $3200\times3200$ points of each fractal picture of the article:
(a) classified; (b) non-finite value in floating point; (c) unclassified after 300 iterations;
(d) stopped by the step criterion farther than $10^{-4}R$ from every root. """ + note + r"""}
\label{tab:s6}\\
\toprule
Picture & (a) & (b) & (c) & (d)\\
\midrule\endfirsthead
\toprule Picture & (a) & (b) & (c) & (d)\\\midrule\endhead
""" + "\n".join(rows) + "\n\\bottomrule\n\\end{longtable}\n")


def write_csv(lines, path="supplementary_dates.csv"):
    fields = ["date", "occasion", "degree", "distinct_roots", "multiplicity_of_zero",
              "max_multiplicity", "polynomial", "window_center", "window_half_width",
              "isolated_root_outside_window", "c1_over_c0", "girard_prediction",
              "real_root_largest_modulus", "relative_error", "critical_points_total",
              "regular_free_critical_points", "critical_poles_count", "critical_poles",
              "regular_free_all_converge", "max_iterations", "regular_free_detail"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(fields)
        for L in lines:
            poles = "; ".join(f"{pp['ponto']} (pole of order {pp['ordem_do_polo']})"
                              for pp in L["_polos"])
            detail = "; ".join(
                f"{ccomp(z0)} -> {ccomp(rz) if rz is not None else 'does not converge'}"
                f" ({it if it else '>10^4'} it.)" for z0, it, rz in L["_res"])
            w.writerow([L["data"], OCCASION[L["data"]], L["grau"], L["raizes_distintas"],
                        L["mult_zero"], L["mult_max"], L["polinomio"], L["centro"],
                        L["meia_largura"], L["raiz_fora_da_janela"], L["c1_c0"], L["previsao"] + 0.0,
                        L["raiz_real_max"], L["erro_rel"], L["pontos_criticos_total"],
                        L["n_livres_regulares"], L["n_polos_criticos"], poles,
                        L["livres_regulares_convergem"], L["iter_max"], detail])


def appendix_block(lines, M, C, F):
    return ("% BEGIN TABLES S1-S6 (written by supplement_en.py; do not edit by hand)\n"
            "{\\small\n"
            + tables_s1_s3(lines, "5.6cm", "12.4cm") + "\n" + table_s4(M) + "\n" + table_s5(C)
            + "\n" + table_s6(F) + "}\n% END TABLES S1-S6\n")


def inline(path, block):
    s = open(path, encoding="utf-8").read()
    a, b = "% BEGIN TABLES S1-S6", "% END TABLES S1-S6\n"
    i, j = s.index(a), s.index(b) + len(b)
    open(path, "w", encoding="utf-8").write(s[:i] + block + s[j:])


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else MANUSCRIPT
    lines = calcular_linhas()
    assert set(OCCASION) == {L["data"] for L in lines}
    F = json.load(open("auditoria_figuras.json"))
    assert set(F) == set(FIGURE)
    block = appendix_block(lines, json.load(open("mandelbrot_orbita.json")),
                           json.load(open("certificado_mandelbrot.json")), F)
    if os.path.exists(path):
        inline(path, block)
    else:
        path = "tables_s1_s6.tex"
        open(path, "w", encoding="utf-8").write(block)
    write_csv(lines)
    print("Tables S1-S6 written into", path, "and supplementary_dates.csv;",
          sum(L["livres_regulares_convergem"] for L in lines), "of", len(lines),
          "dates with all regular free critical points converging; max iterations",
          max(L["iter_max"] for L in lines if L["livres_regulares_convergem"]))

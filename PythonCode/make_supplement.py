"""Tables S1-S6 of Online Resource 1 and the per-date table results/dates.csv.

Tables S1-S3 and the CSV are computed here, for the 38 dates of the article: polynomial,
degree, multiplicities, window, Girard heuristic, and the exact classification and
numerical test of the critical points (birthday.critical). Tables S4-S6 are read from
results/lyapunov.json, results/orbit_intervals_certificate.json and
results/picture_classification.json.

The LaTeX code of the tables is written into supplement/ESM_1.tex, between the lines
"% BEGIN TABLES S1-S6" and "% END TABLES S1-S6". Usage: python make_supplement.py
"""
import csv
import json
import math
from fractions import Fraction

import numpy as np
import sympy as sp

from birthday.critical import distinct_roots, test_critical_points
from birthday.dates import DATES, coefficients, polynomial_latex
from birthday.newton import framing_roots, standard_window

SUPPLEMENT = "supplement/ESM_1.tex"

# pictures of Table S6 (keys of results/picture_classification.json)
PICTURE_LABEL = {
    "Fig1a": "25/12/1642 (Fig. 1a)", "Fig1b": "04/01/1643 (Fig. 1b)",
    "Fig4a": "30/04/1777 (Figs. 4a, 5)",
    "Fig7a": "15/04/1707 (Fig. 7a)", "Fig7b": "16/08/1821 (Fig. 7b)",
    "Fig7c": "29/04/1854 (Fig. 7c)", "Fig7d": "28/02/1878 (Fig. 7d)",
    "Fig7e": "03/02/1893 (Fig. 7e)", "Fig7f": "20/11/1924 (Fig. 7f)",
    "Fig8a": "05/07/1687 (Fig. 8a)", "Fig8b": "23/09/1846 (Fig. 8b)",
    "Fig8c": "30/11/1889 (Fig. 8c)", "Fig8d": "29/05/1919 (Fig. 8d)",
    "Fig8e": "29/12/1972 (Fig. 8e)", "Fig8f": "16/05/1980 (Fig. 8f)",
    "Fig9a": "08/08/1900 (Fig. 9)",
    "Fig10a": "30/08/1984 (Fig. 10a)", "Fig10b": "03/07/1986 (Fig. 10b)",
    "Fig10c": "27/04/1952 (Fig. 10c)", "Fig10d": "15/08/1950 (Fig. 10d)",
    "Fig3a": "$z^2-1$ (Fig. 3a)", "Fig3b": "$z^3-1$ (Fig. 3b)",
    "Fig11": "$z^3-2z+2$ (Fig. 11)",
    "Fig4b": "Gauss, magnification (Fig. 4b)", "Fig4c": "Gauss, magnification (Fig. 4c)",
}


# ---------------------------------------------------------------- per-date data
def compute_rows():
    rows = []
    for date, occasion in DATES:
        c = coefficients(date)
        r = np.roots(c)
        rd, mult = distinct_roots(r)
        center, R = standard_window(c)
        isolated = len(framing_roots(r)) < len(rd)
        real = rd[np.abs(rd.imag) < 1e-9].real
        rmax = real[np.argmax(np.abs(real))] if len(real) else None
        ratio = c[1] / c[0]
        error = abs(rmax - (-ratio)) / abs(rmax) if (rmax is not None and ratio > 0) else None
        cp, res = test_critical_points(c)
        digits = date.replace("/", "")
        rows.append(dict(
            date=date, occasion=occasion, degree=len(c) - 1, distinct_roots=len(rd),
            multiplicity_of_zero=len(digits) - len(digits.rstrip("0")),
            max_multiplicity=int(mult.max()), polynomial=polynomial_latex(c),
            window_center=round(center.real, 4), window_half_width=round(R, 4),
            isolated_root_outside_window=isolated, c1_over_c0=round(ratio, 3),
            girard_prediction=round(-ratio, 3) + 0.0,
            real_root_largest_modulus=None if rmax is None else round(float(rmax), 4),
            relative_error=None if error is None else round(error, 3),
            critical_points_total=cp["total"],
            regular_free_critical_points=len(cp["regular_free"]),
            critical_poles_count=len(cp["critical_poles"]),
            regular_free_all_converge=all(its is not None for _, its, _ in res),
            max_iterations=max((its or 0) for _, its, _ in res) if res else 0,
            _results=res, _poles=cp["critical_poles"]))
    return rows


# ---------------------------------------------------------------- formatting
def tnum(v, nd=3):
    """real number in math mode"""
    if v is None:
        return "---"
    s = f"{v:.{nd}f}"
    if s.lstrip("-").strip("0.") == "":
        s = s.lstrip("-")
    return "$" + s + "$"


def ctext(z, nd=3):
    """complex number as plain text"""
    z = complex(z)
    re_, im_ = round(z.real, nd), round(z.imag, nd)
    re_ = 0.0 if re_ == 0 else re_
    if im_ == 0:
        return f"{re_:.{nd}f}"
    return f"{re_:.{nd}f}{'+' if im_ > 0 else '-'}{abs(im_):.{nd}f}i"


def tcomp(z, nd=3):
    return "$" + ctext(z, nd) + "$"


# ---------------------------------------------------------------- tables
def tables_s1_s3(rows, w_occasion="5.6cm", w_detail="12.4cm"):
    out = ["\\begin{longtable}{@{}l >{\\raggedright\\arraybackslash}p{" + w_occasion
           + "} l c c c@{}}\n"
           "\\caption{Dates, polynomials, degree $d$, number $k$ of distinct roots, and "
           "multiplicity $m_0$ of the root $z=0$ (0 if $z=0$ is not a root). All other roots "
           "are simple (exact verification)}\\label{tab:s1}\\\\\n\\toprule\n"
           "Date & Occasion & $p(z)$ & $d$ & $k$ & $m_0$\\\\\n\\midrule\\endfirsthead\n"
           "\\toprule Date & Occasion & $p(z)$ & $d$ & $k$ & $m_0$\\\\\\midrule\\endhead\n"]
    for L in rows:
        out.append(f"{L['date']} & {L['occasion']} & ${L['polynomial']}$ & {L['degree']} & "
                   f"{L['distinct_roots']} & {L['multiplicity_of_zero']}\\\\\n")
    out.append("\\bottomrule\n\\end{longtable}\n\n")

    out.append("\\begin{longtable}{@{}l r r c r r r@{}}\n"
               "\\caption{Window (real center and half-width $R$), isolated root outside the "
               "window, ratio $c_1/c_0$, real root of largest modulus, and relative error of the "
               "prediction $-c_1/c_0$}\\label{tab:s2}\\\\\n\\toprule\n"
               "Date & center & $R$ & isolated & $c_1/c_0$ & real root & rel.\\ error\\\\\n"
               "\\midrule\\endfirsthead\n"
               "\\toprule Date & center & $R$ & isolated & $c_1/c_0$ & real root & rel.\\ error\\\\"
               "\\midrule\\endhead\n")
    for L in rows:
        out.append(f"{L['date']} & {tnum(L['window_center'])} & {tnum(L['window_half_width'])} & "
                   f"{'yes' if L['isolated_root_outside_window'] else 'no'} & "
                   f"{tnum(L['c1_over_c0'], 2)} & {tnum(L['real_root_largest_modulus'])} & "
                   f"{tnum(L['relative_error'])}\\\\\n")
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
    for L in rows:
        parts = []
        for z0, its, rz in L["_results"]:
            if rz is not None:
                target = tcomp(rz)
            elif L["date"] == "20/11/1924" and abs(complex(z0).imag) < 1e-12:
                target = ("no convergence in $10^4$ iterations; nonconvergence proved "
                          "(Proposition~3 of the article)")
            else:
                target = "no convergence in $10^4$ iterations"
            parts.append(f"{tcomp(z0)} $\\to$ {target} ({its})" if its else f"{tcomp(z0)}: {target}")
        for pp in L["_poles"]:
            pt = sp.latex(sp.sympify(pp["point"]))
            parts.append(f"critical pole ${pt}$, of order {pp['pole_order']}: exact orbit "
                         f"${pt}\\mapsto\\infty\\mapsto\\infty$")
        out.append(f"{L['date']} & " + "; ".join(parts) + "\\\\\n")
    out.append("\\bottomrule\n\\end{longtable}\n")
    return "".join(out)


def table_s4(M):
    def n4(v):
        return f"${v:.4f}$"

    def power(x):
        k = int(round(math.log10(x)))
        return f"$10^{{{k}}}$" if abs(x - 10 ** k) < 1e-9 else f"${x / 10 ** k:.0f}\\times10^{{{k}}}$"

    x0 = {"c": "$c$", "c+1e-08": "$c+10^{-8}$", "c-1e-08": "$c-10^{-8}$", "c+1e-04": "$c+10^{-4}$"}
    precision = {"float64": "double", "long double": "extended (\\texttt{long double})",
                 "mpmath 40 digits": "40 digits (\\texttt{mpmath})"}
    rows = [f"{precision[r['precision']]} & {x0[r['x0']]} & {power(r['n'])} & {n4(r['lambda'])} & "
            f"[{n4(r['min'])}, {n4(r['max'])}]\\\\" for r in M["runs"]]
    return (r"""\begin{table}[htbp]
\centering
\caption{Orbit of the real free critical point $c\approx""" + f"{M['critical_point']:.6f}"
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
    enc = C["exact"]["x_enclosures"]

    def mid_rad(j):
        a, b = Fraction(enc[f"x{j}"]["lo"]), Fraction(enc[f"x{j}"]["hi"])
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
the bounds for each inequality are in the file
\texttt{results/orbit\_intervals\_certificate.json} of the code repository}
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
$\gcd(p',p'')=1$: $N_p$ has no critical poles & """ + yes(C["exact"]["gcd_p1_p2_degree"] == 0) + r"""\\
$p''$ squarefree, a single real zero, $c$ isolated & """
            + yes(C["A1_c_simple_isolated_zero"] and C["real_zeros_of_p2"] == 1) + r"""\\
on $J_-\supset I_-$: $p$ and $p'$ without zeros, $p''$ with one zero ($c$) & """
            + yes(C["A2_J_minus"]) + r"""\\
on $J_+\supset I_+$: $p$, $p'$, and $p''$ without zeros & """ + yes(C["A3_J_plus"]) + r"""\\
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
    for name, v in F.items():
        n = v["counts"]
        if "unclassified_near_cycle_1e-6" in v:
            near = (v["unclassified_near_cycle_1e-6"], n["unclassified"])
        rows.append(f"{PICTURE_LABEL[name]} & {n['classified']:,} & {n['non_finite']:,} & "
                    f"{n['unclassified']:,} & {n['stopped_far']:,}\\\\")
    quantity = f"all {near[1]:,}" if near[0] == near[1] else f"{near[0]:,} of the {near[1]:,}"
    note = (f"For Fig.~11, {quantity} points of class (c) have their last iterate within "
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


CSV_FIELDS = ["date", "occasion", "degree", "distinct_roots", "multiplicity_of_zero",
              "max_multiplicity", "polynomial", "window_center", "window_half_width",
              "isolated_root_outside_window", "c1_over_c0", "girard_prediction",
              "real_root_largest_modulus", "relative_error", "critical_points_total",
              "regular_free_critical_points", "critical_poles_count", "critical_poles",
              "regular_free_all_converge", "max_iterations", "regular_free_detail"]


def write_csv(rows, path="results/dates.csv"):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(CSV_FIELDS)
        for L in rows:
            row = dict(L)
            row["critical_poles"] = "; ".join(f"{pp['point']} (pole of order {pp['pole_order']})"
                                              for pp in L["_poles"])
            row["regular_free_detail"] = "; ".join(
                f"{ctext(z0)} -> {ctext(rz) if rz is not None else 'does not converge'}"
                f" ({its if its else '>10^4'} it.)" for z0, its, rz in L["_results"])
            w.writerow([row[k] for k in CSV_FIELDS])


def write_tables(block, path=SUPPLEMENT):
    s = open(path, encoding="utf-8").read()
    a, b = "% BEGIN TABLES S1-S6", "% END TABLES S1-S6\n"
    i, j = s.index(a), s.index(b) + len(b)
    open(path, "w", encoding="utf-8").write(s[:i] + block + s[j:])


if __name__ == "__main__":
    rows = compute_rows()
    F = json.load(open("results/picture_classification.json"))
    assert set(F) == set(PICTURE_LABEL)
    block = ("% BEGIN TABLES S1-S6 (written by make_supplement.py; do not edit by hand)\n"
             "{\\small\n" + tables_s1_s3(rows) + "\n"
             + table_s4(json.load(open("results/lyapunov.json"))) + "\n"
             + table_s5(json.load(open("results/orbit_intervals_certificate.json"))) + "\n"
             + table_s6(F) + "}\n% END TABLES S1-S6\n")
    write_tables(block)
    write_csv(rows)
    converge = [L for L in rows if L["regular_free_all_converge"]]
    print(f"Tables S1-S6 written into {SUPPLEMENT}, and results/dates.csv; {len(converge)} of "
          f"{len(rows)} dates with all regular free critical points converging; max iterations "
          f"{max(L['max_iterations'] for L in converge)}")

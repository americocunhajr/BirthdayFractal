"""The 38 dates of the article and the rule that turns a date into a polynomial."""
import numpy as np

# (date, occasion), in the order of Tables S1-S3 of Online Resource 1
DATES = [
    ("25/12/1642", "Birth of Isaac Newton (Julian)"),
    ("04/01/1643", "Birth of Isaac Newton (Gregorian)"),
    ("15/04/1707", "Birth of Leonhard Euler"),
    ("30/04/1777", "Birth of Carl Friedrich Gauss"),
    ("16/08/1821", "Birth of Arthur Cayley"),
    ("29/04/1854", "Birth of Henri Poincaré"),
    ("28/02/1878", "Birth of Pierre Fatou"),
    ("03/02/1893", "Birth of Gaston Julia"),
    ("20/11/1924", "Birth of Benoit Mandelbrot"),
    ("05/07/1687", "Publication of the Principia (Julian)"),
    ("15/07/1687", "Publication of the Principia (Gregorian)"),
    ("07/09/1822", "Independence of Brazil"),
    ("23/09/1846", "Galle observes Neptune"),
    ("13/05/1888", "Golden Law (abolition of slavery in Brazil)"),
    ("21/01/1889", "Announcement of King Oscar II's prize"),
    ("15/11/1889", "Proclamation of the Republic of Brazil"),
    ("30/11/1889", "Poincaré asks to stop the printing"),
    ("08/08/1900", "Hilbert's problems"),
    ("23/10/1906", "Flight of the 14-bis"),
    ("29/05/1919", "Sobral eclipse"),
    ("15/10/1952", "Founding of IMPA"),
    ("29/12/1972", "Lorenz's talk"),
    ("16/05/1980", "Founding of LNCC"),
    ("05/10/1988", "Brazilian Federal Constitution"),
    ("23/06/1993", "Wiles's announcement"),
    ("22/04/1500", "Arrival of Cabral (Julian)"),
    ("02/05/1500", "Arrival of Cabral (proleptic Gregorian)"),
    ("04/10/1582", "Last Julian day"),
    ("15/10/1582", "First Gregorian day"),
    ("20/07/1969", "Apollo 11 (Brasília time)"),
    ("21/07/1969", "Apollo 11 (UTC)"),
    ("14/09/2015", "Detection of gravitational waves"),
    ("11/02/2016", "Announcement of the detection"),
    ("01/01/2000", "New Year's Day 2000"),
    ("30/08/1984", "Birth of the author"),
    ("03/07/1986", "Birth of the author's wife"),
    ("27/04/1952", "Birth of the author's father"),
    ("15/08/1950", "Birth of the author's mother"),
]


def coefficients(date, order="DMY"):
    """Coefficients of the polynomial of a date DD/MM/YYYY, highest degree first:
    the digits of DDMMYYYY (order="DMY") or YYYYMMDD (order="YMD"), leading zeros removed."""
    d, m, y = date.split("/")
    digits = {"DMY": d + m + y, "YMD": y + m + d}[order]
    c = np.trim_zeros(np.array([int(ch) for ch in digits], float), "f")
    if c.size < 3:
        raise ValueError("polynomial of degree < 2")
    return c


def polynomial_latex(c):
    """LaTeX of the polynomial with integer coefficients c, e.g. 2z^7+z^5+4."""
    c = [int(v) for v in c]
    n = len(c) - 1
    terms = []
    for i, ci in enumerate(c):
        power = n - i
        if ci == 0:
            continue
        coef = "" if (ci == 1 and power > 0) else str(ci)
        var = "" if power == 0 else ("z" if power == 1 else f"z^{power}")
        terms.append(coef + var)
    return "+".join(terms)

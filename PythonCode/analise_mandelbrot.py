"""
Orbita do ponto critico livre real de N_p para 20/11/1924 (Mandelbrot).

Indicador de Lyapunov de tempo finito ao longo da orbita, excluindo o ponto critico
(onde N_p' = 0):
    lambda_n = (1/n) * sum_{k=1}^{n} log |N_p'(x_k)|,   x_k = N_p^k(c).
Sensibilidade: numero de iteracoes, precisao aritmetica (float64, long double,
mpmath com 40 digitos) e ponto inicial perturbado.
"""
import json
import numpy as np
import mpmath as mp

from newton_aniversario import coeficientes

c = np.array(coeficientes("20/11/1924"), float)          # 2 0 1 1 1 9 2 4
dc, d2c = np.polyder(c), np.polyder(c, 2)
raizes = np.roots(c)
raiz_real = float(raizes[np.argmin(np.abs(raizes.imag))].real)
crit = [z for z in np.roots(d2c) if abs(z.imag) < 1e-12]
c0 = float(crit[0].real)
polos = sorted(float(z.real) for z in np.roots(dc) if abs(z.imag) < 1e-12)


def orbita(x0, n, dtype=np.float64):
    cc = c.astype(dtype); dd = dc.astype(dtype); ee = d2c.astype(dtype)
    x = dtype(x0)
    soma = dtype(0)
    xmin = xmax = None
    dmin_raiz = np.inf
    for k in range(n):
        x = x - np.polyval(cc, x) / np.polyval(dd, x)      # x_{k+1}
        px, dpx, d2px = np.polyval(cc, x), np.polyval(dd, x), np.polyval(ee, x)
        soma += np.log(abs(px * d2px / dpx ** 2))
        xmin = x if xmin is None else min(xmin, x)
        xmax = x if xmax is None else max(xmax, x)
        dmin_raiz = min(dmin_raiz, abs(float(x) - raiz_real))
    return float(soma / n), float(xmin), float(xmax), dmin_raiz


def orbita_mp(x0, n, dps=40):
    mp.mp.dps = dps
    C = [mp.mpf(int(v)) for v in c]
    D = [mp.mpf(int(v)) * (len(c) - 1 - i) for i, v in enumerate(c[:-1])]
    E = [D[i] * (len(D) - 1 - i) for i in range(len(D) - 1)]
    x = mp.mpf(x0)
    soma = mp.mpf(0)
    xmin, xmax = mp.inf, -mp.inf
    for k in range(n):
        x = x - mp.polyval(C, x) / mp.polyval(D, x)
        soma += mp.log(abs(mp.polyval(C, x) * mp.polyval(E, x) / mp.polyval(D, x) ** 2))
        xmin, xmax = min(xmin, x), max(xmax, x)
    return float(soma / n), float(xmin), float(xmax)


if __name__ == "__main__":
    out = {"polinomio": c.tolist(), "ponto_critico": c0, "raiz_real": raiz_real,
           "polos_reais": polos, "execucoes": []}
    print("ponto critico c =", c0, " raiz real =", raiz_real, " polos reais =", polos)
    for n in (10**4, 10**5, 10**6):
        lam, a, b, dr = orbita(c0, n)
        out["execucoes"].append({"precisao": "float64", "x0": "c", "n": n, "lambda": lam,
                                 "min": a, "max": b, "dist_min_raiz": dr})
        print(f"float64  n={n:>8}: lambda={lam:.4f}  intervalo=[{a:.4f}, {b:.4f}]  dist min a raiz={dr:.3f}")
    lam, a, b, dr = orbita(c0, 10**5, np.longdouble)
    out["execucoes"].append({"precisao": "long double", "x0": "c", "n": 10**5, "lambda": lam,
                             "min": a, "max": b, "dist_min_raiz": dr})
    print(f"longdbl  n=  100000: lambda={lam:.4f}  intervalo=[{a:.4f}, {b:.4f}]")
    for dx in (1e-8, -1e-8, 1e-4):
        lam, a, b, dr = orbita(c0 + dx, 10**5)
        out["execucoes"].append({"precisao": "float64", "x0": f"c{dx:+.0e}", "n": 10**5,
                                 "lambda": lam, "min": a, "max": b, "dist_min_raiz": dr})
        print(f"float64  x0=c{dx:+.0e} n=100000: lambda={lam:.4f}  intervalo=[{a:.4f}, {b:.4f}]")
    lam, a, b = orbita_mp(c0, 2 * 10**4)
    out["execucoes"].append({"precisao": "mpmath 40 digitos", "x0": "c", "n": 2 * 10**4,
                             "lambda": lam, "min": a, "max": b})
    print(f"mpmath40 n=   20000: lambda={lam:.4f}  intervalo=[{a:.4f}, {b:.4f}]")
    json.dump(out, open("mandelbrot_orbita.json", "w"), indent=1)

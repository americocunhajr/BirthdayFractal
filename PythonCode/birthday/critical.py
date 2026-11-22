"""Critical points of the Newton map N_p = z - p/p' on the Riemann sphere.

Exact classification (SymPy), for p with integer coefficients:
  simple roots of p          superattracting fixed points;
  critical poles             zeros of p' of multiplicity mu >= 2 that are not roots of p
                             (poles of order mu, local degree mu), found with gcd(p', p'');
  regular free critical      zeros of p'' that are not zeros of p or of p' (N_p' = 0 there).
  points
With multiplicity, the total must be 2k - 2, where k is the number of distinct roots.

Numerical test: each regular free critical point is iterated in double precision, from
a 30-digit approximation, until |step| < 1e-14 * max(1, max|root|) and the distance to a
root is < 1e-9, with at most 10^4 iterations. Critical poles are not iterated: their
exact orbit is pole -> infinity -> infinity.
"""
import numpy as np
import sympy as sp

_z = sp.symbols("z")


def distinct_roots(r, tol=1e-6):
    """Distinct values among the numerical roots r, and their multiplicities."""
    out, mult = [], []
    for z in r:
        for i, w in enumerate(out):
            if abs(z - w) < tol:
                mult[i] += 1
                break
        else:
            out.append(z)
            mult.append(1)
    return np.array(out), np.array(mult)


def _remove_factors(q, target):
    """Remove from q every irreducible factor it shares with target, with multiplicity."""
    while True:
        h = sp.gcd(q, target)
        if h.degree() <= 0:
            return q
        q = sp.Poly(sp.quo(q, h), _z)


def critical_points(c):
    """Exact classification of the critical points of N_p (p with integer coefficients c)."""
    P = sp.Poly([int(v) for v in c], _z)
    d1, d2 = P.diff(_z), P.diff(_z).diff(_z)
    _, factors = sp.factor_list(P.as_expr())
    simple = [f for f, e in factors if e == 1]
    k = sum(sp.Poly(f, _z).degree() for f, e in factors)
    n_simple = sum(sp.Poly(f, _z).degree() for f in simple)
    Ps = sp.Poly(sp.Mul(*simple), _z) if simple else sp.Poly(1, _z)
    # extra critical multiplicity at simple roots (zeros of p'' at simple roots)
    extra = d2.degree() - _remove_factors(d2, Ps).degree()
    # critical poles: factors of gcd(p', p'') that do not divide p
    g = sp.gcd(d1, d2)
    poles = []
    if g.degree() > 0:
        for f, _e in sp.factor_list(g.as_expr())[1]:
            fp = sp.Poly(f, _z)
            if sp.rem(P, fp).is_zero:
                continue                      # zero of p: removable singularity
            mu = 0
            q = d1
            while sp.rem(q, fp).is_zero:
                q = sp.Poly(sp.quo(q, fp), _z)
                mu += 1
            for root in sp.Poly(f, _z).all_roots():
                poles.append({"point": str(root), "value": complex(sp.N(root, 30)),
                              "pole_order": mu, "critical_multiplicity": mu - 1})
    n_poles = sum(pp["critical_multiplicity"] for pp in poles)
    # regular free critical points: zeros of p'' that are not zeros of p or of p'
    q = _remove_factors(d2, sp.Poly(P.as_expr() * d1.as_expr(), _z))
    free = [complex(r) for r in q.nroots(n=30)] if q.degree() > 0 else []
    total = n_simple + extra + n_poles + q.degree()
    assert total == 2 * k - 2, (c, total, k)
    return {"k": k, "simple_roots": n_simple, "extra_at_simple_roots": extra,
            "critical_poles": poles, "regular_free": free, "total": total}


def test_critical_points(c, maxit=10_000):
    """Exact classification, then numerical iteration of the regular free critical points.
    Returns (classification, [(z0, iterations or None, root reached or None), ...])."""
    cp = critical_points(c)
    p, dp = np.poly1d(c), np.poly1d(c).deriv()
    rd, _ = distinct_roots(np.roots(c))
    scale = max(1.0, np.max(np.abs(rd)))
    results = []
    for z0 in cp["regular_free"]:
        z = complex(z0)
        its = None
        for k in range(1, maxit + 1):
            d = dp(z)
            if d == 0:
                break
            step = p(z) / d
            z -= step
            if abs(step) < 1e-14 * scale and np.min(np.abs(rd - z)) < 1e-9:
                its = k
                break
        results.append((complex(z0), its, complex(rd[np.argmin(np.abs(rd - z))]) if its else None))
    return cp, results

"""Newton's method on a grid of starting points, and the colored pictures of the article.

Each starting point z0 is iterated by z <- z - a p(z)/p'(z) (a = 1: Newton's method)
until the step is smaller than tol * R, where R is the half-width of the window, or
for at most maxit iterations. A point is classified by the nearest root if it ends
within 1e-4 R of it. Brightness comes from a smoothed iteration count, normalized in
each basin; colors are defined in OKLCH and mixed in OKLab.
"""
import numpy as np
from PIL import Image, ImageDraw


# ---------------------------------------------------------------- window
def framing_roots(roots, kappa=3.0):
    """Distinct roots used to frame the picture: excludes isolated roots, those farther
    than kappa times the median distance from the robust center (medians of the real
    and imaginary parts of the distinct roots)."""
    distinct = []
    for z in roots:
        if all(abs(z - w) > 1e-6 for w in distinct):
            distinct.append(z)
    distinct = np.array(distinct)
    center = np.median(distinct.real) + 1j * np.median(distinct.imag)
    d = np.abs(distinct - center)
    return distinct[d <= kappa * np.median(d)]


def standard_window(c, margin=1.3, focus=True):
    """Center and half-width of the standard window (isolated roots left out if focus)."""
    roots = np.roots(c)
    r = framing_roots(roots) if focus else roots
    center = r.mean()
    return center, margin * np.max(np.abs(r - center))


# ---------------------------------------------------------------- iteration
def iterate(c, n=800, ss=2, maxit=300, tol=1e-10, margin=1.3, focus=True,
            center=None, R=None, Rref=None, relax=1.0, with_status=False):
    """Newton's method from an (n*ss) x (n*ss) grid.

    center, R: window (default: standard_window). Rref: scale used in the tolerance and
    in the classification (default: R; magnifications pass the R of the full view).
    relax: step factor a in z <- z - a p/p'. Returns (basin, nu, roots, center, R), where
    basin is the index of the root reached (-1: unclassified) and nu the smoothed
    iteration count. with_status=True also returns a dictionary with the class of each
    point (0 classified; 1 non-finite step or iterate; 2 unclassified after maxit
    iterations; 3 stopped farther than 1e-4 Rref from every root) and the last iterate.
    """
    roots = np.roots(c)
    dc = np.polyder(c)
    if center is None or R is None:
        center, R = standard_window(c, margin, focus)
    Rref = R if Rref is None else Rref
    N = n * ss
    x = np.linspace(center.real - R, center.real + R, N)
    y = np.linspace(center.imag + R, center.imag - R, N)  # row 0 = top
    Z = (x[None, :] + 1j * y[:, None]).ravel()

    kstop = np.full(Z.size, np.nan)   # iteration at which the point stopped
    logstep = np.full(Z.size, np.nan)  # log of the last step
    active = np.arange(Z.size)
    tolR = tol * Rref
    logtol = np.log(tolR)
    with np.errstate(all="ignore"):
        for k in range(maxit):
            za = Z[active]
            step = np.polyval(c, za) / np.polyval(dc, za)
            if relax != 1:
                step = relax * step
            Z[active] = za - step
            astep = np.abs(step)
            done = ~np.isfinite(step) | (astep < tolR)
            kstop[active[done]] = k
            logstep[active[done]] = np.log(np.maximum(astep[done], 1e-300))
            active = active[~done]
            if active.size == 0:
                break

    # classification by the nearest root
    best = np.full(Z.size, np.inf)
    basin = np.full(Z.size, -1)
    with np.errstate(invalid="ignore"):
        for j, r in enumerate(roots):
            d = np.abs(Z - r)
            closer = d < best
            best[closer] = d[closer]
            basin[closer] = j
    basin[~np.isfinite(Z) | ~(best < 1e-4 * Rref) | np.isnan(kstop)] = -1

    # smoothed count: quadratic convergence (simple root) or linear with ratio
    # (m-1)/m (root of multiplicity m)
    mult = np.array([np.sum(np.abs(roots - r) < 1e-6) for r in roots])
    m = np.where(basin >= 0, mult[np.maximum(basin, 0)], 1)
    with np.errstate(all="ignore"):
        frac_quadratic = np.log2(np.clip(logstep / logtol, 1.0, 2.0 ** 1.5))
        rho = np.where(m > 1, (m - 1) / m, 0.5)
        frac_linear = np.clip((logtol - logstep) / -np.log(rho), 0.0, 1.0)
    nu = kstop - np.where(m > 1, frac_linear, frac_quadratic)
    if with_status:
        status = np.zeros(Z.size, dtype=np.int8)
        stopped = ~np.isnan(kstop)
        with np.errstate(invalid="ignore"):
            status[stopped & ~np.isfinite(Z)] = 1
            status[~stopped] = 2
            status[stopped & np.isfinite(Z) & ~(best < 1e-4 * Rref)] = 3
        info = {"status": status.reshape(N, N), "Zlast": Z.reshape(N, N)}
        return basin.reshape(N, N), nu.reshape(N, N), roots, center, R, info
    return basin.reshape(N, N), nu.reshape(N, N), roots, center, R


# ---------------------------------------------------------------- colors
def oklch_to_oklab(L, C, h):
    h = np.deg2rad(h)
    return np.stack([L, C * np.cos(h), C * np.sin(h)], axis=-1)


def oklab_to_srgb(lab):
    L, a, b = lab[..., 0], lab[..., 1], lab[..., 2]
    l_ = L + 0.3963377774 * a + 0.2158037573 * b
    m_ = L - 0.1055613458 * a - 0.0638541728 * b
    s_ = L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    rgb = np.stack([
        4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s], axis=-1)
    rgb = np.clip(rgb, 0, 1)
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)


def hex_to_oklab(hx):
    v = np.array([int(hx[i:i + 2], 16) / 255 for i in (1, 3, 5)])
    lin = np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)
    l = 0.4122214708 * lin[0] + 0.5363325363 * lin[1] + 0.0514459929 * lin[2]
    m = 0.2119034982 * lin[0] + 0.6806995451 * lin[1] + 0.1073969566 * lin[2]
    s = 0.0883024619 * lin[0] + 0.2817188376 * lin[1] + 0.6299787005 * lin[2]
    l_, m_, s_ = np.cbrt([l, m, s])
    return np.array([
        0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_,
        1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_,
        0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_])


# palette of the article: basin colors with fixed OKLCH lightness and chroma, hues
# spread around H0; light tone near the roots, ink near the boundaries
LIGHTNESS, CHROMA, H0 = 0.88, 0.075, 220
INK = "#2b2d42"
LIGHTEN = 0.40
INK_RGB = (43, 45, 66)
RING_RGB = (21, 23, 38)


def hues(roots, center, h0=H0):
    """Hue of each root: distinct roots sorted by angle around the center and spread
    evenly around h0, so that conjugate pairs get hues symmetric about h0."""
    distinct = []
    for z in roots:
        if all(abs(z - w) > 1e-6 for w in distinct):
            distinct.append(z)
    distinct = np.array(distinct)
    order = np.argsort(np.angle(distinct - center))
    d = len(distinct)
    h = np.empty(d)
    for pos, i in enumerate(order):
        h[i] = h0 + 360 * (pos + 0.5) / d - 180
    return np.array([h[np.argmin(np.abs(distinct - z))] for z in roots])


def smoothstep(t):
    return t * t * (3 - 2 * t)


def colorize(basin, nu, roots, center, scale=6.0, normalize=True, floors=None):
    """sRGB image of a computed grid. center: reference point for the order of the hues.
    normalize: subtract the iteration floor (2nd percentile) of each basin, so that slow
    basins are not dark; floors: {basin: floor} from another image (magnifications use
    the floors of the full view). Returns (rgb, floors)."""
    hs = hues(roots, center)
    base = oklch_to_oklab(np.full(len(roots), LIGHTNESS), np.full(len(roots), CHROMA), hs)
    ink = hex_to_oklab(INK)
    white = np.array([1.0, 0.0, 0.0])

    t = np.zeros(basin.shape)
    used = {}
    valid = basin >= 0
    global_floor = np.percentile(nu[valid], 2) if valid.any() else 0.0
    for j in np.unique(basin[valid]):
        m = basin == j
        if not normalize:
            floor = global_floor
        elif floors is not None and j in floors:
            floor = floors[j]
        else:
            floor = np.percentile(nu[m], 2)
        used[int(j)] = floor
        t[m] = 1 - np.exp(-np.maximum(nu[m] - floor, 0) / scale)
    t[~valid] = 1.0

    b = base[np.maximum(basin, 0)]
    light = b + LIGHTEN * (white - b)
    u = smoothstep(np.clip(t / 0.35, 0, 1))[..., None]
    v = smoothstep(np.clip((t - 0.35) / 0.65, 0, 1))[..., None] ** 1.3
    lab = light + u * (b - light)      # near the root: light tone; then the base color
    lab = lab + v * (ink - lab)        # near the boundary: ink
    lab[~valid] = ink
    return oklab_to_srgb(lab), used


# ---------------------------------------------------------------- picture
def render(c, n=800, ss=2, show_roots=True, focus=True, center=None, R=None, Rref=None,
           color_center=None, floors=None, normalize=True, relax=1.0):
    """Picture of Newton's method for the polynomial with coefficients c, n x n pixels
    (computed on an (n*ss)^2 grid and reduced). Returns (image, fraction of unclassified
    points, floors, (center, R))."""
    c = np.asarray(c, float)
    basin, nu, roots, center, R = iterate(c, n=n, ss=ss, focus=focus, center=center, R=R,
                                          Rref=Rref, relax=relax)
    cc = center if color_center is None else color_center
    rgb, floors = colorize(basin, nu, roots, cc, normalize=normalize, floors=floors)
    N = n * ss
    img = Image.fromarray((rgb * 255 + 0.5).astype(np.uint8))
    if show_roots:
        draw = ImageDraw.Draw(img)
        rad = 0.0055 * N
        for z in roots:
            px = (z.real - (center.real - R)) / (2 * R) * (N - 1)
            py = ((center.imag + R) - z.imag) / (2 * R) * (N - 1)
            if -rad <= px <= N + rad and -rad <= py <= N + rad:
                draw.ellipse([px - rad, py - rad, px + rad, py + rad],
                             fill=(255, 255, 255), outline=RING_RGB, width=max(1, int(rad / 3)))
    return img.resize((n, n), Image.LANCZOS), np.mean(basin < 0), floors, (center, R)

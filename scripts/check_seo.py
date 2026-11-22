#!/usr/bin/env python3
"""Validate search-facing static assets without external dependencies."""
from __future__ import annotations

import json
from pathlib import Path
import re
import struct
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent.parent
SITE = "https://americocunhajr.github.io/BirthdayFractal/"
SITEMAP = ROOT / "sitemap.xml"
HTML_FILES = [ROOT / "index.html", *sorted((ROOT / "pages").glob("*.html"))]
assert len(HTML_FILES) == 6, f"Expected homepage and 5 articles; found {len(HTML_FILES)}"

urls = [
    el.text
    for el in ET.parse(SITEMAP).findall(
        ".//{http://www.sitemaps.org/schemas/sitemap/0.9}loc"
    )
]
expected = {
    SITE + ("" if page.name == "index.html" else page.relative_to(ROOT).as_posix())
    for page in HTML_FILES
}
assert set(urls) == expected, f"Sitemap mismatch: missing {expected - set(urls)}, extra {set(urls) - expected}"
assert len(urls) == len(set(urls)), "Duplicate URLs in sitemap"

for page in HTML_FILES:
    html = page.read_text(encoding="utf-8")
    url = SITE + ("" if page.name == "index.html" else page.relative_to(ROOT).as_posix())
    canonical = re.findall(r'<link rel="canonical" href="([^"]+)"', html)
    assert canonical == [url], f"Missing/incorrect canonical for {page}"
    assert '<meta name="description" content="' in html, f"Missing description: {page}"
    assert 'name="robots" content="index, follow' in html, f"Unexpected robots rules: {page}"
    assert f'<meta property="og:url" content="{url}">' in html, f"Incorrect Open Graph URL: {page}"
    assert f'<meta property="og:image" content="{SITE}og-preview.png">' in html, f"Incorrect preview: {page}"
    scripts = re.findall(r'<script type="application/ld\+json">([\s\S]*?)</script>', html)
    assert len(scripts) == 1, f"Missing/duplicate JSON-LD: {page}"
    schema = json.loads(scripts[0])
    assert schema["@context"] == "https://schema.org"
    assert schema["@graph"], f"Empty graph: {page}"
    assert "<h1>" in html or re.search(r"<h1[\s>]", html), f"Missing H1: {page}"
    if page.name != "index.html":
        assert 'href="../?date=' in html, f"No fractal-generator link: {page}"
        assert 'href="../styles.css"' in html, f"Broken relative styles: {page}"

home = (ROOT / "index.html").read_text(encoding="utf-8")
for p in HTML_FILES[1:]:
    assert f'href="{p.relative_to(ROOT).as_posix()}"' in home, f"Not linked from homepage: {p}"
assert "applicationCategory" in home

image = (ROOT / "og-preview.png").read_bytes()
assert image[:8] == b"\x89PNG\r\n\x1a\n", "Invalid PNG signature"
assert struct.unpack(">II", image[16:24]) == (1200, 630), "Preview dimensions are not 1200x630"

print("SEO validation passed: 6 canonical/indexable HTML URLs, valid JSON-LD, sitemap, internal links, and 1200x630 PNG.")

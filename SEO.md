# Search indexing for BirthdayFractal

Canonical URL: https://americocunhajr.github.io/BirthdayFractal/

The static [XML sitemap](sitemap.xml) includes the homepage and five stand-alone mathematical date articles. All articles have unique canonical URLs, descriptions, Open Graph and social-image metadata, and JSON-LD. The homepage has WebSite, WebApplication, and Person metadata.

The [social preview](og-preview.png) is an illustrative cubic Newton fractal, not a rendering of a specific birthday polynomial.

## Steps requiring the site owner's action

1. In https://search.google.com/search-console, add the **URL-prefix property** `https://americocunhajr.github.io/BirthdayFractal/`.
2. Follow Google's ownership verification method. If given an HTML verification file, add it unchanged to the repository root; if given an HTML tag, place it in `index.html`'s head. Do not fabricate a verification token.
3. Submit the sitemap URL `https://americocunhajr.github.io/BirthdayFractal/sitemap.xml` under **Sitemaps**.
4. Use **URL Inspection** on the homepage and a date page, select **Test live URL**, then **Request indexing** when offered.
5. Monitor the Pages and Search results reports for actual indexing/impressions.

For a GitHub Pages project at `/BirthdayFractal/`, search crawler robots rules are taken from the **host-root** `https://americocunhajr.github.io/robots.txt`, *not* a project-level `/BirthdayFractal/robots.txt`. We do not add a misleading project-level robots.txt or promise special AI search rankings. The content is already static HTML, so crawlers can read these notes without running JavaScript.

To add a date page, update the sitemap and the homepage's linked date story collection. The academic software listing is in https://americocunha.org/software/ and the source and documentation are in the public GitHub repository.

## Publication and validation

In **Settings → Pages**, select **Deploy from a branch** → `main` → `/ (root)`. GitHub's built-in Pages workflow publishes the site. The independent **Validate BirthdayFractal website** workflow checks JavaScript syntax, canonical URLs, JSON-LD, the sitemap, internal links, and the social image. Keeping these jobs separate avoids duplicate deployment artifacts.

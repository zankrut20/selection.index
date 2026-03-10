/* Force default dark mode if no user preference is stored. */
window.addEventListener('DOMContentLoaded', function () {
  if (localStorage.getItem('theme') === null) {
    document.documentElement.setAttribute('data-bs-theme', 'dark');
    var darkButton = document.querySelector('[data-bs-theme-value="dark"]');
    if (darkButton) {
      // Simulate click to update button state and lightswitch.js internal state
      darkButton.click();
      // Remove it from localStorage so it doesn't force this preference later
      localStorage.removeItem('theme');
    }
  }

  /* ── SEO: canonical link ─────────────────────────────────────────── */
  if (!document.querySelector('link[rel="canonical"]')) {
    var canonical = document.createElement('link');
    canonical.rel = 'canonical';
    canonical.href = window.location.href.split('?')[0].split('#')[0];
    document.head.appendChild(canonical);
  }

  /* ── SEO: meta keywords ──────────────────────────────────────────── */
  if (!document.querySelector('meta[name="keywords"]')) {
    var kw = document.createElement('meta');
    kw.name = 'keywords';
    kw.content = [
      'selection index', 'plant breeding', 'quantitative genetics',
      'R package', 'CRAN', 'genomic selection', 'marker-assisted selection',
      'genetic advance', 'phenotypic selection', 'LPSI', 'LGSI', 'LMSI',
      'eigen selection index', 'multistage selection', 'stochastic simulation',
      'heritability', 'variance-covariance matrix', 'economic weights',
      'breeding pipeline', 'selection.index'
    ].join(', ');
    document.head.appendChild(kw);
  }

  /* ── SEO: Open Graph type (pkgdown sets title/desc, we fix type) ─── */
  var ogType = document.querySelector('meta[property="og:type"]');
  if (!ogType) {
    var meta = document.createElement('meta');
    meta.setAttribute('property', 'og:type');
    meta.content = 'website';
    document.head.appendChild(meta);
  }

  /* ── SEO: Twitter Card tags ──────────────────────────────────────── */
  var twitterCard = document.querySelector('meta[name="twitter:card"]');
  if (!twitterCard) {
    [
      { name: 'twitter:card',        content: 'summary' },
      { name: 'twitter:title',       content: document.title },
      { name: 'twitter:description', content: (document.querySelector('meta[name="description"]') || {}).content || '' },
      { name: 'twitter:site',        content: '@zankrut20' }
    ].forEach(function (attrs) {
      var m = document.createElement('meta');
      m.name    = attrs.name;
      m.content = attrs.content;
      document.head.appendChild(m);
    });
  }

  /* ── SEO: JSON-LD SoftwareApplication structured data ────────────── */
  if (!document.querySelector('script[type="application/ld+json"]')) {
    var schema = {
      '@context': 'https://schema.org',
      '@type': 'SoftwareApplication',
      'name': 'selection.index',
      'alternateName': 'selection.index R Package',
      'description': 'A production-ready R package for plant breeders and quantitative geneticists to compute classical phenotypic, genomic, marker-assisted, restricted, constrained, and eigen selection indices. Maximizes genetic advance with multi-trait LPSI, LGSI, LMSI, and stochastic simulation tools.',
      'applicationCategory': 'ScientificApplication',
      'applicationSubCategory': 'Plant Breeding & Quantitative Genetics',
      'operatingSystem': 'Windows, macOS, Linux',
      'programmingLanguage': 'R',
      'softwareVersion': '2.0.1',
      'dateModified': '2026-03-06',
      'license': 'https://www.gnu.org/licenses/gpl-3.0',
      'url': 'https://zankrut20.github.io/selection.index/',
      'downloadUrl': 'https://cran.r-project.org/package=selection.index',
      'installUrl': 'https://cran.r-project.org/package=selection.index',
      'author': {
        '@type': 'Person',
        'name': 'Zankrut Goyani',
        'email': 'zankrut20@gmail.com',
        'url': 'https://github.com/zankrut20'
      },
      'maintainer': {
        '@type': 'Person',
        'name': 'Zankrut Goyani',
        'email': 'zankrut20@gmail.com'
      },
      'publisher': {
        '@type': 'Person',
        'name': 'Zankrut Goyani'
      },
      'sameAs': [
        'https://github.com/zankrut20/selection.index',
        'https://cran.r-project.org/package=selection.index',
        'https://cloud.r-project.org/package=selection.index'
      ],
      'offers': {
        '@type': 'Offer',
        'price': '0',
        'priceCurrency': 'USD'
      },
      'keywords': [
        'selection index', 'plant breeding', 'quantitative genetics',
        'genomic selection', 'marker-assisted selection', 'LPSI', 'LGSI',
        'eigen selection index', 'multistage selection', 'stochastic simulation',
        'genetic advance', 'heritability', 'R package'
      ]
    };
    var ld = document.createElement('script');
    ld.type = 'application/ld+json';
    ld.textContent = JSON.stringify(schema, null, 2);
    document.head.appendChild(ld);
  }
});


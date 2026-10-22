"""Math macros for the OOFEM manuals, defined once for every builder.

The manuals inherit a set of math macros from their original LaTeX sources
(``\\mbf``, ``\\del``, ``\\eps`` and the rest).  Those have to reach two
different consumers: MathJax, which renders the HTML, and the LaTeX preamble,
which renders the PDF.

Declaring them twice is what broke the PDF build: ``matlibmanual`` and
``elementlibmanual`` configured ``mathjax3_config`` only, so the HTML had the
definitions and LaTeX had none, and ``make all-pdf`` stopped at the first
``Undefined control sequence``.  This extension derives both from a single
dictionary, so the two can no longer drift apart.

Usage -- in a manual's ``conf.py``::

    extensions = [..., 'oofemmath']

    oofem_math_macros = {
        'eps':  r'\\mbf{\\varepsilon}',      # no argument
        'del':  [r'\\displaystyle\\frac{#1}{#2}', 2],   # two arguments
    }

The macros common to all the manuals (:data:`COMMON_MACROS`) are added
automatically; a manual only declares what is specific to it, and may override
a common one by redefining it.

The dictionary uses the MathJax spelling -- a string for a macro without
arguments, or ``[body, nargs]`` with them -- because that is the form the HTML
already renders.  The LaTeX definitions are generated from it, so the PDF comes
out matching the HTML rather than the other way round.
"""

from sphinx.util import logging

logger = logging.getLogger(__name__)

#: Macros every OOFEM manual may use.  These come from the shared part of the
#: original LaTeX preambles, where ``\mbf`` was ``\mbox{\boldmath$#1$}``;
#: ``\boldsymbol`` is the equivalent that both MathJax and amsmath understand.
COMMON_MACROS = {
    'mbf':  [r'\boldsymbol{#1}', 1],
    'del':  [r'\displaystyle\frac{#1}{#2}', 2],
    'der':  [r'\frac{{\rm d}{#1}}{{\rm d}{#2}}', 2],
    'pard': [r'\frac{\partial{#1}}{\partial{#2}}', 2],
}


def _split(spec):
    """Return ``(body, nargs)`` for a MathJax macro specification."""
    if isinstance(spec, (list, tuple)):
        body, nargs = spec[0], int(spec[1]) if len(spec) > 1 else 0
    else:
        body, nargs = spec, 0
    return body, nargs


def collect_macros(config):
    """Merge the common macros with the ones the manual declares."""
    macros = dict(COMMON_MACROS)
    macros.update(getattr(config, 'oofem_math_macros', None) or {})
    return macros


def latex_preamble(macros):
    """Render the macros as LaTeX definitions.

    ``\\providecommand`` followed by ``\\renewcommand`` rather than a plain
    ``\\newcommand``: a few of these names (``\\sym``, ``\\ud``, ``\\e``) could
    well be taken by some package in the future, and that combination ends up
    with our definition either way instead of failing the build.
    """
    lines = ['% -- math macros shared by the OOFEM manuals (doc/_ext/oofemmath.py)']
    for name in sorted(macros):
        body, nargs = _split(macros[name])
        args = '[%d]' % nargs if nargs else ''
        lines.append(r'\providecommand{\%s}{}' % name)
        lines.append(r'\renewcommand{\%s}%s{%s}' % (name, args, body))
    return '\n'.join(lines)


def _on_config_inited(app, config):
    macros = collect_macros(config)
    if not macros:
        return

    # MathJax, for the HTML builders.
    mathjax = dict(getattr(config, 'mathjax3_config', None) or {})
    tex = dict(mathjax.get('tex') or {})
    declared = dict(tex.get('macros') or {})
    # A manual that still configures mathjax3_config by hand keeps winning, so
    # that this extension can be adopted without touching existing settings.
    merged = dict(macros)
    merged.update(declared)
    tex['macros'] = merged
    mathjax['tex'] = tex
    config.mathjax3_config = mathjax

    # LaTeX preamble, for the PDF builder.
    elements = dict(getattr(config, 'latex_elements', None) or {})
    preamble = elements.get('preamble', '')
    block = latex_preamble(merged)
    if block not in preamble:
        elements['preamble'] = (preamble + '\n' + block) if preamble else block
    config.latex_elements = elements

    logger.info('oofemmath: %d math macros defined for html and latex',
                len(merged))


def setup(app):
    app.add_config_value('oofem_math_macros', {}, 'env')
    app.connect('config-inited', _on_config_inited)

    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

"""Settings of the combined documentation site, read from the manuals.

The top-level ``doc/conf.py`` builds one site holding every manual.  Each manual
also builds on its own, from its own ``conf.py``, and the combined site needs
what those declare -- the extensions providing ``:param:``, ``record`` or
``tikz``, the math macros, the files to skip.  Repeating them in the top-level
``conf.py`` is how the two drift apart, so they are read from the manuals::

    from oofemcombined import collect
    extensions, oofem_math_macros, exclude_patterns = collect(MANUALS)

This is a helper for ``conf.py``, not a Sphinx extension.

The manuals' ``conf.py`` files are parsed, not executed: executing one would
rerun its ``sys.path`` and ``global_conf`` setup, which is written relative to
the manual's own directory.  Only settings assigned a literal value are read;
anything else is taken as unset.
"""

import ast
import os

#: Extensions a manual loads that have no business in the combined site.
#: oofemlatex only fixes up the LaTeX output; the combined site is HTML only.
SKIP_EXTENSIONS = {'oofemlatex'}


def _literals(path):
    """Return the literal top-level assignments of a ``conf.py`` file."""
    with open(path, encoding='utf-8') as f:
        tree = ast.parse(f.read(), path)
    values = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name):
                try:
                    values[target.id] = ast.literal_eval(node.value)
                except ValueError:
                    values.pop(target.id, None)
    return values


def collect(manuals, root='.'):
    """Return ``(extensions, oofem_math_macros, exclude_patterns)`` for the
    combined site, read from ``<root>/<manual>/conf.py`` of each manual.

    The extensions are the union of the manuals' lists, in first-seen order.
    The math macros are merged; two manuals defining the same macro
    differently is an error, as the combined site can only have one of them.
    The exclude patterns are each manual's own, prefixed with its directory.
    """
    extensions, macros, excludes = [], {}, []
    origin = {}
    for manual in manuals:
        conf = _literals(os.path.join(root, manual, 'conf.py'))

        for ext in conf.get('extensions', []):
            if ext not in extensions and ext not in SKIP_EXTENSIONS:
                extensions.append(ext)

        for name, spec in conf.get('oofem_math_macros', {}).items():
            if name in macros and macros[name] != spec:
                raise ValueError(
                    'math macro \\%s is defined differently in %s and %s'
                    % (name, origin[name], manual))
            macros[name] = spec
            origin.setdefault(name, manual)

        for pattern in conf.get('exclude_patterns', []):
            excludes.append(manual + '/' + pattern)

    return extensions, macros, excludes


def connect_substitutions(app, manuals, author, copyright, root='.'):
    """Give the pages of each manual its own ``|author|`` and ``|copyright|``.

    A standalone build appends them through ``rst_epilog`` (see
    ``global_conf.setup``), but the combined site has a single ``rst_epilog``,
    which would show one author list in every manual.  Instead the manual's
    substitutions are appended to each of its ``.rst`` pages as they are read.
    ``author`` and ``copyright`` are used for a manual that does not set its
    own.
    """
    from global_conf import substitutions

    text = {}
    for manual in manuals:
        conf = _literals(os.path.join(root, manual, 'conf.py'))
        text[manual] = substitutions(conf.get('author', author),
                                     conf.get('copyright', copyright))

    def _source_read(app, docname, source):
        manual = docname.split('/', 1)[0]
        if manual in text and str(app.env.doc2path(docname)).endswith('.rst'):
            source[0] += text[manual]

    app.connect('source-read', _source_read)

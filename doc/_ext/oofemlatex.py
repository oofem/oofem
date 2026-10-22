"""Shared LaTeX/PDF presentation fixes for the OOFEM manuals.

Things that have to be right in every manual's PDF and that Sphinx does not
get right on its own.  Load it from a manual's ``conf.py``::

    extensions = [..., 'oofemlatex']

Author list on the title page
-----------------------------

``sphinxmanual.cls`` typesets the author inside a one-column ``tabular``::

    \\begin{tabular}[t]{c}
      \\@author
    \\end{tabular}

A ``tabular`` row never line-breaks, so the comma-separated author string runs
off the page: the material and element library manuals, with ten names, were
overfull by 763pt -- more than twice the text width -- and even the five-name
manuals overflowed by 36pt.  The class expects the names to be separated by
``\\\\`` so that each gets its own row, which is what this extension does for
the LaTeX builder only; the HTML keeps the comma-separated form.
"""

from sphinx.util import logging

logger = logging.getLogger(__name__)


def latex_author(author):
    """Turn a comma-separated author string into ``tabular`` rows."""
    parts = [a.strip() for a in author.split(',') if a.strip()]
    return r' \\ '.join(parts)


def _on_config_inited(app, config):
    # latex_documents entries are (startdocname, targetname, title, author,
    # theme, toctree_only).  Sphinx has already filled in its default from
    # config.author by now, and has LaTeX-escaped it on the way -- "et al."
    # becomes "et al.\@{}" -- so the field is rewritten in place rather than
    # compared against config.author, which would no longer match.
    documents, rows = [], 0
    for entry in (config.latex_documents or []):
        entry = list(entry)
        if len(entry) > 3 and isinstance(entry[3], str):
            wrapped = latex_author(entry[3])
            if wrapped != entry[3]:
                entry[3] = wrapped
                rows = max(rows, wrapped.count(r'\\') + 1)
        documents.append(tuple(entry))
    if rows:
        config.latex_documents = documents
        logger.info('oofemlatex: author split over %d title-page rows', rows)


def setup(app):
    # 'config-inited' fires before the LaTeX builder reads latex_documents, but
    # after Sphinx has filled in its default from project/author.
    app.connect('config-inited', _on_config_inited, priority=800)

    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

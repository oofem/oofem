"""A ``tikz`` directive for the OOFEM manuals.

The element library manual draws its element figures with TikZ, and the
conversion keeps that source rather than pre-rendering it, so the pictures stay
editable.  When `sphinxcontrib-tikz`_ is installed we hand the directive over to
it and the pictures are rendered properly.  When it is not, the documentation
must still build, so we fall back to showing the TikZ source as a labelled,
captioned literal block.

.. _sphinxcontrib-tikz: https://pypi.org/project/sphinxcontrib-tikz/
"""

from docutils import nodes
from docutils.parsers.rst import Directive, directives


class TikzFallback(Directive):
    """Render a TikZ picture as its source, with the caption beneath it."""

    has_content = True
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = {
        # Accepted so that sources stay compatible with sphinxcontrib-tikz.
        "libs": directives.unchanged,
        "libraries": directives.unchanged,
        "stringsubst": directives.flag,
        "xscale": directives.unchanged,
        "align": directives.unchanged,
        "include": directives.unchanged,
    }

    def run(self):
        source = "\n".join(self.content)
        container = nodes.container(classes=["oofem-tikz"])

        literal = nodes.literal_block(source, source)
        literal["language"] = "latex"
        container += literal

        if self.arguments:
            caption = self.arguments[0]
            node = nodes.caption(caption, "")
            self.state.nested_parse(
                self.content.__class__([caption], source=caption), 0, node)
            container += node

        self.add_name(container)
        return [container]


def setup(app):
    try:
        app.setup_extension("sphinxcontrib.tikz")
    except Exception:
        # sphinxcontrib-tikz is unavailable; show the source instead.
        app.add_directive("tikz", TikzFallback)

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }

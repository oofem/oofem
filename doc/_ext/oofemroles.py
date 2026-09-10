"""Sphinx roles mirroring the OOFEM manual LaTeX macros.

The LaTeX sources define a handful of markup commands for describing input
records (see the preamble of matlibmanual.tex and ../include/include.tex)::

    \\descitem{name}            bold keyword introducing a record
    \\elemparam{name}{type}     mandatory parameter, rendered as "name (type)"
    \\optelemparam{name}{type}  optional parameter, rendered as "[name (type)]"
    \\elemstring{name}          verbatim string
    \\optelemstring{name}       optional verbatim string
    \\param{name}               parameter name in running text
    \\optparam{name}            optional parameter name in running text

The LaTeX->RST conversion emitted these as interpreted text roles, so we define
matching roles here rather than flattening them in the sources.  Keeping them as
roles preserves the mandatory/optional distinction and gives a single place to
restyle every parameter reference.

The parameter type (the second LaTeX argument) is carried inside the role text
in braces, e.g. :elemparam:`E{rn}`.
"""

import re

from docutils import nodes

#: ``name{type}`` as produced by the normalisation pass over the RST sources.
_TYPED = re.compile(r"^(?P<name>.*?)\s*\{(?P<type>[^{}]*)\}$", re.DOTALL)


def _split_type(text):
    """Split ``name{type}`` into ``(name, type)``; ``type`` is None if absent."""
    m = _TYPED.match(text)
    if m:
        return m.group("name"), m.group("type")
    return text, None


def _param_nodes(text, optional):
    """Build the nodes for a parameter reference, optionally bracketed."""
    name, ptype = _split_type(text)
    result = [nodes.literal(name, name, classes=["oofem-param"])]
    if ptype:
        result.append(nodes.inline(ptype, "(%s)" % ptype,
                                   classes=["oofem-paramtype"]))
    if optional:
        result = ([nodes.Text("[")] + result + [nodes.Text("]")])
    return result


def _make_param_role(optional):
    def role(name, rawtext, text, lineno, inliner, options=None, content=None):
        return _param_nodes(nodes.unescape(text), optional), []
    return role


def _make_wrapper_role(node_class, prefix="", suffix="", classes=None):
    def role(name, rawtext, text, lineno, inliner, options=None, content=None):
        value = nodes.unescape(text)
        node = node_class(value, prefix + value + suffix,
                          classes=list(classes or []))
        return [node], []
    return role


def setup(app):
    # Parameters of an input record: name plus optional type in braces.
    app.add_role("elemparam", _make_param_role(optional=False))
    app.add_role("optelemparam", _make_param_role(optional=True))
    app.add_role("param", _make_param_role(optional=False))
    app.add_role("optparam", _make_param_role(optional=True))

    # Record keyword and verbatim strings.
    app.add_role("descitem", _make_wrapper_role(
        nodes.strong, classes=["oofem-descitem"]))
    app.add_role("elemstring", _make_wrapper_role(nodes.literal))
    app.add_role("optelemstring", _make_wrapper_role(
        nodes.literal, prefix="[", suffix="]"))

    # Leftovers from LaTeX font-switching commands.
    for alias in ("emph", "em", "it", "sl"):
        app.add_role(alias, _make_wrapper_role(nodes.emphasis))
    app.add_role("bf", _make_wrapper_role(nodes.strong))
    app.add_role("tt", _make_wrapper_role(nodes.literal))
    app.add_role("mmt", _make_wrapper_role(nodes.literal))

    return {
        "version": "1.0",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }

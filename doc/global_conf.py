# Configuration file for the Sphinx documentation builder.
#
# This file only contains global oofem documentation settings. Individual
# manuals may override these settings by importing this file.
#
copyright = '2025, Bořek Patzák, www.oofem.org'
#copyright = '%Y, Bořek Patzák'
author = 'Bořek Patzák, Martin Horák, Vít Šmilauer, Milan Jirásek, et al.'

# The full version, including alpha/beta/rc tags
release = '3.0'
version = '3.0'


def substitutions(author, copyright):
    """The |author| and |copyright| substitutions, as reST source."""
    return (f"\n.. |author| replace:: {author}\n"
            f".. |copyright| replace:: {copyright}\n")


# |author| and |copyright| are appended to every .rst file once the
# configuration is final, rather than fixed here: a manual that sets its own
# author after "from global_conf import *" (the library manuals do) would
# otherwise show the generic list above.  The "from global_conf import *" in a
# manual's conf.py brings this setup() along, and Sphinx calls it.
def setup(app):
    def _epilog(app, config):
        config.rst_epilog = ((config.rst_epilog or '')
                             + substitutions(config.author, config.copyright))
    app.connect('config-inited', _epilog)

# -- Shared HTML presentation -------------------------------------------------
#
# Every manual uses the alabaster theme, so they look and navigate alike.  Do
# not set html_theme or html_theme_options in a manual; anything a manual needs
# on top of these should extend the dict rather than replace it.
html_theme = 'alabaster'
html_theme_options = {
    'sidebar_width': '250px',
    'show_relbar_bottom': True,
}
# alabaster's default sidebar, listed explicitly so that every manual shows the
# same one: title, global table of contents, prev/next links and search.
html_sidebars = {
    '**': ['about.html', 'navigation.html', 'relations.html', 'searchbox.html'],
}

# doc/_static/oofem.css is the single place where the page width and the figure
# sizing are defined, for every manual and every theme.  A manual that needs
# static directories of its own must keep '../_static' among them, so that the
# shared stylesheet is still copied.
html_static_path = ['../_static']
html_css_files = ['oofem.css']

# Number figures and tables, so :numref: works in all manuals.
numfig = True

# -- Shared figure conversion -------------------------------------------------
#
# Manuals loading 'sphinx.ext.imgconverter' rasterise their PDF figures for the
# HTML builders.  The extension defaults to 96 dpi, which leaves a vector
# drawing visibly soft once it is scaled up to the text column; 200 dpi is
# still a small file and stays sharp.  Ignored by manuals that do not load the
# extension.
image_converter_args = ['-density', '200', '-define', 'pdf:use-cropbox=true']
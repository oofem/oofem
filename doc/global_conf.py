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

# Inject into all .rst files as a substitution
rst_epilog = f"""
.. |author| replace:: {author}
.. |copyright| replace:: {copyright}
"""

# -- Shared HTML presentation -------------------------------------------------
#
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
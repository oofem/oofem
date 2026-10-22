# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# -- Project information -----------------------------------------------------

project = 'OOFEM Element Library Manual'
copyright = '2026, Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
author = 'Bořek Patzák, Martin Horák, Mikael Öhman, Milan Jirásek, Vít Šmilauer, Peter Grassl, Petr Havlásek, Edita Dvořáková, Jim Brouzoulis, Carl Sandström'
release = '3.0'

# -- General configuration ---------------------------------------------------

# doc/_ext holds the extensions shared by the OOFEM manuals: oofemroles defines
# the roles mirroring the \param, \field and \descitem macros, and oofemtikz
# provides the "tikz" directive used by the element figures.
sys.path.insert(0, os.path.abspath(os.path.join('..', '_ext')))

extensions = ['oofemroles', 'oofemtikz', 'oofemmath', 'oofemlatex']

# The source directory also holds the LaTeX sources and the latex2html output.
exclude_patterns = [
    '_build',
    'html',
    'auto',
    'Thumbs.db',
    '.DS_Store',
]

# ":numref:" is used throughout to reference tables and figures by number.
numfig = True

# The equations of this manual use \mbf, \del, \der and \pard, inherited from
# elementlibmanual.tex.  All four are common to the OOFEM manuals and come from
# oofemmath, which defines them for MathJax (HTML) and in the LaTeX preamble
# (PDF) alike, so no macros need to be declared here.  Anything specific to this
# manual would go in an "oofem_math_macros" dict.

# -- Options for HTML output -------------------------------------------------

html_theme = 'alabaster'
# doc/_static is shared by the OOFEM manuals.
html_static_path = ['../_static']
html_css_files = ['oofem.css']

# The element summary tables are wide, so the text column is widened to 80% of
# the browser window.  That, and the figure sizing, is set once for every
# manual in doc/_static/oofem.css -- do not add page_width or body_max_width
# here, as the stylesheet loads after the theme and would override them.
html_theme_options = {
    'sidebar_width': '250px',
}

master_doc = 'index'

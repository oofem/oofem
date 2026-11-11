# Configuration file for the Sphinx documentation builder.
#
# This builds the combined site holding every OOFEM manual.  Each manual still
# builds on its own from its own conf.py; the settings this site needs from
# them are read from those files rather than repeated here (see
# _ext/oofemcombined.py).
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('.'))
# doc/_ext holds the extensions shared by the OOFEM manuals.
sys.path.insert(0, os.path.abspath('_ext'))

# Shared settings: the alabaster theme, the stylesheet, numfig.
from global_conf import *
import global_conf
from oofemcombined import collect, connect_substitutions

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OOFEM'
copyright = '2025, Bořek Patzák &al.'
author = 'Bořek Patzák &al.'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# The manuals of the combined site.  Keep in step with the toctree in index.md.
MANUALS = ['usermanual', 'oofemInput', 'theorymanual2', 'developerGuide',
           'python', 'extractorInput', 'elementlibmanual', 'matlibmanual']

_extensions, oofem_math_macros, _excludes = collect(MANUALS)

extensions = ['myst_parser'] + [e for e in _extensions if e != 'myst_parser']
myst_enable_extensions=['attrs_inline','linkify','colon_fence']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store'] + _excludes


# Replaces the setup() imported from global_conf: its single rst_epilog would
# give every manual the same |author|.  Each manual's pages get that manual's
# own author instead, or global_conf's for a manual that sets none.
def setup(app):
    connect_substitutions(app, MANUALS, global_conf.author,
                          global_conf.copyright)

#breathe_projects={'oofem':'refman/xml/'}
#breathe_default_project='oofem'

# -- Options for HTML output -------------------------------------------------

# This directory holds doc/_static itself, hence not global_conf's
# '../_static'.
html_static_path = ['_static']
html_theme_options = dict(
    html_theme_options,
    # Shown under the table of contents in the sidebar.
    extra_nav_links={
        'C++ reference': 'http://www.oofem.org/resources/doc/oofemrefman/index.html',
        'GitHub': 'https://github.com/oofem/oofem',
        'Homepage': 'https://oofem.org',
    },
)

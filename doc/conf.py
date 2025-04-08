# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OOFEM'
copyright = '2025, Bořek Patzák &al.'
author = 'Bořek Patzák &al.'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['myst_parser']
myst_enable_extensions=['attrs_inline','linkify','colon_fence']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

#breathe_projects={'oofem':'refman/xml/'}
#breathe_default_project='oofem'

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_theme_options=dict(
    external_links=[
        {'name':'Element Library','url':'http://www.oofem.org/resources/doc/elementlibmanual/html/elementlibmanual.html'},
        {'name':'Material Library','url':'http://www.oofem.org/resources/doc/matlibmanual/html/matlibmanual.html'},
        {'name':'Extractor','url':'http://www.oofem.org/resources/doc/extractorInput/html/extractorInput.html'},
        {'name':'C++ reference','url':'http://www.oofem.org/resources/doc/oofemrefman/index.html'},
    ],
    icon_links=[
        {'name':'GitHub','url':'/https://github.com/oofem/oofem','icon':'fa-brands fa-square-github','type':'fontawesome'},
        {'name':'Homepage','url':'https://oofem.org','icon':'fa-solid fa-house','type':'fontawesome'},
    ],
    # global navigation (left sidebar)
    show_nav_level=4,
    navigation_depth=4,
    # page TOC (right sidebar)
    show_toc_level=2,
    secondary_sidebar_items=['page-toc'] # no sourcelink
)

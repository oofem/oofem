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
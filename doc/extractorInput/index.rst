OOFEM Extractor Manual
===============================================

Authors: |author|

Version: |version|, Date: |today|

``extractor.py`` is a small Python tool, distributed with
`OOFEM <https://oofem.org>`_ in the ``tools`` directory, that reads selected
values from an OOFEM output file.  It can print them step by step, which is
handy for plotting a load--displacement curve or a time history.  It can also
check them against reference values, which is what it was originally written
for: regression testing.

This manual describes how to run the tool and the format of its input file.
The same record syntax is understood, in a slightly stricter dialect, by the
``errorcheck`` export module built into OOFEM itself, which runs the regression
tests; :doc:`errorcheck` lists the differences.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   inputformat
   errorcheck


.. The index and search pages exist in the HTML only; in the PDF this would be
   an empty chapter.

.. only:: html

   .. rubric:: Indices and tables

   * :ref:`genindex`
   * :ref:`search`

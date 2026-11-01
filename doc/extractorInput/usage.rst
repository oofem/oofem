Running the extractor
=====================

The extractor locates user-defined values in an OOFEM output file.  The values
are described in an *extractor input file*, whose format is given in
:doc:`inputformat`.  The tool works in one of two modes.

Extraction mode
    The default.  For every solution step found in the output file, the tool
    prints one row to the standard output.  The columns are the extracted
    values, in the order in which their records appear in the input file.

Checking mode
    Each value is located in the output file and compared with the reference
    value given in its record.  If the difference exceeds the tolerance, an
    error is reported.

Command line
------------

.. code-block:: none

   extractor.py [-c] [-p] -f inputFileNamePath

``-f inputFileNamePath``
    Path to the extractor input file.  The option is compulsory; without it the
    tool prints usage information and exits.

``-c``
    Turns on the checking mode.  Without it, the extraction mode is used.

``-p``
    Lets the check pass even when some values cannot be found in the output
    file; each missing value is then reported as a warning instead of an
    error.  The same can be set in the input file with ``error_when_missing``,
    see :ref:`value-records-section`.

The tool is a plain Python 3 script, so it can also be started as
``python3 tools/extractor.py ...`` from the OOFEM source tree.

Output
------

In extraction mode every value is printed in ``%15e`` format; a value that was
not found in a given step is printed as ``---``.  A row is printed for each
solution step, that is, for each ``Output for time`` header in the output file.

In checking mode, one line is printed per input file::

   Checker.py: patch100.in                                            [  OK  ]

When a check fails, the line ends with ``[FAILED]`` and each failing record is
listed above it, by its index in the input file, with the found value, the
expected value and their difference.

Exit status
-----------

``0``
    Extraction finished, or all the checks passed.

``1``
    At least one check failed, or a record is missing a compulsory parameter.

``2``
    The input file is malformed: the ``#%BEGIN_CHECK%`` or ``#%END_CHECK%``
    record was not found, or ``#INCLUDE`` directives are nested too deeply.

When files are included with ``#INCLUDE``, the statuses of the included files
are added to that of the including file, so any non-zero status means failure.

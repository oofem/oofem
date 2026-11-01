The errorcheck export module
============================

OOFEM has a built-in counterpart of the checking mode: the ``errorcheck``
export module, which is what the regression tests in ``tests/regression`` run.
It reads its rules from a value records section like the one described in
:ref:`value-records-section`, taken from the input file itself or from the file
given by its ``filename`` parameter.  The module and its parameters
(``filename``, ``extract``, ``writeist``) are described in the
`OOFEM input manual
<https://www.oofem.org/resources/doc/oofemInput/html/pages/app-export-modules.html#error-checking-module>`_.

Unlike ``extractor.py``, the module does not parse the output file.  It takes
the values straight from the analysis while it runs, so they are not rounded
to the printed precision, and a check can use a tolerance as tight as
``1.e-12``.  Its record syntax differs from that of ``extractor.py`` in the
following points.

Differences from extractor.py
-----------------------------

Field order
    The rules in a text input file are read with fixed patterns, so the
    parameters **must be given in the order shown below**, with no others in
    between.  Optional parameters may be left out.

Header
    The tolerance is read from a fixed position, so the header must be exactly
    ``#%BEGIN_CHECK% tolerance <value>``, or just ``#%BEGIN_CHECK%``, in which
    case the default tolerance is ``1.e-6``.  ``error_when_missing`` is not
    supported.

Solution step
    ``tStep`` is the solution step *number*, an integer, not its time.  The
    optional ``tStepVer`` selects the version of the step, for an analysis that
    solves a step more than once; it defaults to ``0``.

Per-record tolerance
    Every record may end with ``tolerance #``, which overrides the tolerance
    from the header for that record.

Records
    ``#DOFMAN`` is not accepted, only ``#NODE``; ``unknown`` must be ``d``,
    ``v`` or ``a`` (total value, velocity, acceleration), and ``type`` is not
    accepted.  ``#ERRORESTIMATE`` and ``#INCLUDE`` are not supported, and any
    other record starting with a single ``#`` is an error.  Lines starting with
    ``##`` are comments.  ``#EIGVAL`` and ``#ELEMENTNODE`` are additional
    records.

Keywords
    The ``keyword`` of ``#ELEMENT`` is the integer number of an
    ``InternalStateType`` (``1`` for stresses, ``4`` for strains, and so on;
    see ``src/core/internalstatetype.h``).  The ``keyword`` of
    ``#BEAM_ELEMENT`` is ``0`` for the local end displacements and ``1`` for
    the local end forces.  ``extractor.py`` reads the latter as stresses, so
    ``#BEAM_ELEMENT`` rules written for the module cannot be checked with
    ``extractor.py``.

Record syntax
-------------

.. code-block:: none

   #%BEGIN_CHECK% [tolerance #]
   #NODE         tStep # [tStepVer #] number # dof # unknown # value #
   #ELEMENTNODE  tStep # number # dofman # dof # unknown # value #
   #ELEMENT      tStep # [tStepVer #] number # [irule #] gp # keyword # component # value #
   #BEAM_ELEMENT tStep # [tStepVer #] number # keyword # component # value #
   #REACTION     tStep # [tStepVer #] number # dof # value #
   #LOADLEVEL    tStep # value #
   #EIGVAL       tStep # EigNum # value #
   #TIME
   #%END_CHECK%

Every record except ``#TIME`` may also end with ``tolerance #``, as described
above.  The records not found in ``extractor.py`` are:

``#ELEMENTNODE``
    An unknown of an internal dof manager of an element, which is not a node of
    the mesh.  ``number`` is the element number, ``dofman`` the index of the
    internal dof manager within the element; the other parameters are as for
    ``#NODE``.

``#EIGVAL``
    The eigenvalue number ``EigNum`` of an eigenvalue analysis.

``#TIME``
    Always passes.  It is useful in extraction mode, where it writes the time
    of the step as a column.

``#REACTION`` is supported for structural problems only.

Running the checks
------------------

A failed rule is reported as a warning in the step it applies to, with the
computed and the expected value.  After the last step, OOFEM stops with an
error if any rule has failed, so the run ends with a non-zero exit status.

With the ``extract`` parameter, the module checks nothing.  In every step it
writes the values of all the rules as one row to the file
``<output base name>.m<n>``, where ``<n>`` is the number of the export module,
ignoring ``tStep``.  With ``writeist``, it prints a rule set for the current
results to the standard output: ``#NODE`` rules for all the DOFs, and
``#ELEMENT`` rules for the listed internal state types.  This is the usual way
to create the rules for a new regression test.

Example
-------

An excerpt from ``tests/regression/sm/patch100.in``, the same patch test as in
the examples of :doc:`inputformat`:

.. code-block:: none

   #%BEGIN_CHECK% tolerance 1.e-4
   ## check reactions
   #REACTION tStep 1 number 1 dof 1 value 2.5
   #REACTION tStep 1 number 1 dof 2 value 1.40625
   ## check all nodes
   #NODE tStep 1 number 4 dof 1 unknown d value -1.5625
   ## check element no. 1 strain and stress vectors
   #ELEMENT tStep 1 number 1 gp 1 keyword 4 component 1  value -0.520833333
   #ELEMENT tStep 1 number 1 gp 1 keyword 1 component 1  value -8.333333333
   #%END_CHECK%

The same rules can also be checked with ``extractor.py -c``: since the solution
step at time 1 is step number 1, ``tStep 1`` means the same for both tools.

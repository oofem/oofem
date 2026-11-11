Input file format
=================

An extractor input file consists of two compulsory sections: the **source file
section**, followed by the **value records section**, which determines the
values to be extracted or checked.

Source file section
-------------------

The first line of the file that does not start with ``#`` holds the path to
the OOFEM output file the values are read from.  The path is opened as given,
so a relative path is resolved against the current working directory, not
against the directory of the input file.  The line should contain nothing but
the path.

Any lines between the source file section and the value records section are
ignored.  The first line of an OOFEM input file is the name of its output file,
and every extractor record starts with ``#``, which OOFEM treats as a comment.
An OOFEM input file and an extractor input file can therefore be merged into
one: the regression tests in ``tests/regression`` do exactly that.

.. _value-records-section:

Value records section
---------------------

The value records section describes the values to be extracted or checked.  Its
general format is

.. code-block:: none

   #%BEGIN_CHECK% [tolerance #] [error_when_missing #]
   #(DOFMAN|NODE)  {tStep #} number # dof # (unknown|type) # {value #}
   #ELEMENT        {tStep #} number # [irule #] gp # keyword # component # {value #}
   #BEAM_ELEMENT   {tStep #} number # keyword # component # {value #}
   #REACTION       {tStep #} number # dof # {value #}
   #LOADLEVEL      {tStep #} {value #}
   #ERRORESTIMATE  {tStep #} {value #}
   <#TIME>
   #INCLUDE fileName
   #%END_CHECK%

The notation is as follows:

* ``#`` standing alone after a parameter name is the value to be substituted,
  an integer, a real number or a string depending on the parameter;
* parameters in ``[ ]`` are optional;
* parameters in ``{ }`` are required in checking mode only, and ignored in
  extraction mode;
* records in ``< >`` are available in extraction mode only;
* ``(a|b)`` means ``a`` or ``b``.

The ``#%BEGIN_CHECK%`` and ``#%END_CHECK%`` records are compulsory.  Between
them, the records may come in any order, and so may the parameters within a
record.  In extraction mode the values are printed in columns in the order of
their records.  A string value, such as a ``keyword``, may be enclosed in double
quotes.  Lines that do not match any record are ignored, so a line starting
with ``##`` can serve as a comment.

The ``#%BEGIN_CHECK%`` record takes two options.

``tolerance``
    The absolute tolerance used in checking mode: a check fails when the
    found value differs from the expected one by more than this.  The default
    is ``1.e-4``.  Write the number with a decimal point (``1.e-3``, not
    ``1e-3``); otherwise it is not recognised and the default is used.

``error_when_missing``
    With ``0``, a value that cannot be found in the output file produces only a
    warning, and the check can still pass.  With ``1``, the default, it is an
    error.  This is the same as the ``-p`` command line option.

Common parameters
~~~~~~~~~~~~~~~~~

``tStep``
    In checking mode, the solution step to which the expected value belongs,
    given by its *time*: the value printed in the ``Output for time`` header of
    the step in the output file.  The two are compared for equality, so the
    time must be given to the precision with which it is printed.

``value``
    In checking mode, the expected value.

DOFMAN, NODE
~~~~~~~~~~~~

.. code-block:: none

   #(DOFMAN|NODE) {tStep #} number # dof # (unknown|type) # {value #}

An unknown associated with a degree of freedom (DOF) of a node or another dof
manager.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``number``
     - Number (label) of the dof manager.
   * - ``dof``
     - DOF number.  The DOF numbering is determined by the domain type, or by
       the DOF ordering defined for the dof manager (see the OOFEM input
       manual).
   * - ``unknown``
     - A single character selecting the physical meaning of the value, as
       printed in front of the value in the DOF lines of the output file.
       ``d`` is the primary unknown: the displacement in structural analysis,
       the temperature or other field in transport problems.  The dynamic
       analyses also print ``v`` for velocity and ``a`` for acceleration.  The
       tool further accepts ``t``, ``f`` and ``p``.  ``type`` is accepted as a
       synonym of ``unknown``.

ELEMENT
~~~~~~~

.. code-block:: none

   #ELEMENT {tStep #} number # [irule #] gp # keyword # component # {value #}

A value associated with an integration point of an element.  Usually many
values are printed for each integration point, for example the strain and
stress vectors in structural analysis, and which ones depends on the problem
being solved and the constitutive model used.  Each is identified by its
keyword and by its position in the array printed after that keyword.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``number``
     - Number (label) of the element.
   * - ``irule``
     - Number of the integration rule of the element.  Defaults to ``1``, the
       first integration rule.
   * - ``gp``
     - Number of the integration point within the rule.  In the output file,
       the integration point is the one headed ``GP irule.gp :``.
   * - ``keyword``
     - Name of the quantity, the word preceding its values in the
       integration point record, for example ``stresses`` or ``strains``.  Any
       word followed by numbers on a line of the record can be used.  The
       integer codes below are accepted as well; they are the
       ``InternalStateType`` numbers that the ``errorcheck`` module uses (see
       :doc:`errorcheck`).
   * - ``component``
     - Position of the value in the array following the keyword, starting
       from 1.

.. list-table:: Integer keyword codes
   :widths: 18 82
   :header-rows: 1

   * - Code
     - Keyword
   * - ``1``
     - ``stresses``
   * - ``4``
     - ``strains``
   * - ``13``
     - ``damage``
   * - ``39``
     - ``DoH``
   * - ``56``
     - ``flow``
   * - ``91``
     - ``F``
   * - ``92``
     - ``P``
   * - ``98``
     - ``jump``
   * - ``99``
     - ``traction``
   * - ``135``
     - ``equivalentTime``

BEAM_ELEMENT
~~~~~~~~~~~~

.. code-block:: none

   #BEAM_ELEMENT {tStep #} number # keyword # component # {value #}

A value associated with a beam element formulated without integration
points, with exact analytical integration, for which the results are printed
as end forces and end displacements.  At present these are the ``Beam2d`` and
``Beam3d`` elements.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``number``
     - Number (label) of the element.
   * - ``keyword``
     - ``displacements`` for the local end displacements, ``forces`` for the
       local end forces.
   * - ``component``
     - Position of the value in the printed vector, starting from 1.

REACTION
~~~~~~~~

.. code-block:: none

   #REACTION {tStep #} number # dof # {value #}

The reaction force associated with a DOF subjected to a Dirichlet boundary
condition, as printed in the reaction table of the output file.

.. list-table::
   :widths: 18 82
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``number``
     - Number (label) of the dof manager.
   * - ``dof``
     - DOF number.

LOADLEVEL
~~~~~~~~~

.. code-block:: none

   #LOADLEVEL {tStep #} {value #}

The load level reached in the solution step, from the ``Reached load level :``
line printed by the nonlinear static and nonlinear explicit dynamic analyses.

ERRORESTIMATE
~~~~~~~~~~~~~

.. code-block:: none

   #ERRORESTIMATE {tStep #} {value #}

The relative error estimate of an adaptive analysis, read from a line starting
with ``Relative error estimate:``.

.. note::

   The error estimators of the current OOFEM version print the estimate in a
   different form (the ZZ estimator as an indented percentage, the Huerta
   estimator to the log only), so this record does not find a value in their
   output.

TIME
~~~~

.. code-block:: none

   #TIME

Available in extraction mode only.  Prints the time of the solution step as a
column, which is useful as the abscissa of a time history.

INCLUDE
~~~~~~~

.. code-block:: none

   #INCLUDE fileName

Processes another extractor input file, with its own source file and value
records sections, as part of this one.  ``#include`` is accepted too.  The file
name may contain only letters, digits, underscores and dots, so the file must
be in the current working directory.  Included files may not include further
files.  In checking mode, the including file fails if any included file fails,
so a single master file can run a whole set of checks.

Examples
--------

The following input file checks two values in the integration point of
element 1, and one displacement, in the step at time 1:

.. code-block:: none

   patch100.out
   #%BEGIN_CHECK% tolerance 1.e-4
   #ELEMENT tStep 1 number 1 gp 1 component 1 keyword "stresses" value -8.3333e+00
   #ELEMENT tStep 1 number 1 gp 1 component 1 keyword "strains"  value -5.2083e-01
   #DOFMAN  tStep 1 number 4 dof 1 type d value -1.5625
   #%END_CHECK%

Run with ``extractor.py -c -f`` on this file, the tool prints a single line
ending in ``[  OK  ]`` when all three values match.

Without ``-c`` and with the checking-only parameters left out, the same kind of
input extracts a history, one row per solution step: here the time, the
displacement of node 4 and the reaction in DOF 1 of node 1.

.. code-block:: none

   patch100.out
   #%BEGIN_CHECK%
   #TIME
   #NODE number 4 dof 1 unknown d
   #REACTION number 1 dof 1
   #%END_CHECK%

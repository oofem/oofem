.. _eigensolverssection:

Eigenvalue solver parameters
============================

Wherever a record syntax shows :optelemparam:`eigensolverparams{...}`, that
stands for the following group of attributes, written directly in the same
record.

.. record:: Syntax

   :elemparam:`stype{in}` :optelemparam:`smtype{in}`
   :elemparam:`solverParams{s}`

**Parameters**

:param:`stype{in}`
    Eigenvalue solver type; see :numref:`eigenvaluesolverparamtable`.

:optparam:`smtype{in}`
    Sparse matrix storage scheme.  It must be compatible with the solver type;
    see :ref:`sparselinsolver`.

:param:`solverParams{s}`
    The parameters specific to the selected solver.

.. _eigenvaluesolverparamtable:

.. table:: Eigenvalue solver parameters.

   +--------------------+-----------------+-------------------------------------+
   | Solver type        | ``stype`` id    | Solver parameters                   |
   +====================+=================+=====================================+
   | Subspace Iteration | 0 (default)     | None.                               |
   +--------------------+-----------------+-------------------------------------+
   | Inverse Iteration  | 1               | None.                               |
   +--------------------+-----------------+-------------------------------------+
   | SLEPc solver       | 2               | Requires ``smtype 7``.              |
   |                    |                 | See also the SLEPc manual.          |
   +--------------------+-----------------+-------------------------------------+

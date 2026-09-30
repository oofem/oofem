.. _sparselinsolver:

Sparse linear solver parameters
===============================

Wherever a record syntax shows :optelemparam:`sparselinsolverparams{...}`, that
stands for the following group of attributes, written directly in the same
record.

.. record:: Syntax

   :optelemparam:`lstype{in}` :optelemparam:`smtype{in}`
   :elemparam:`solverParams{s}`

**Parameters**

:optparam:`lstype{in}`
    Solver for the linear system of equations.  Default ``0``.  See
    :numref:`linsolvstoragecompattable` for the identifiers and
    :numref:`sparsesolverparams` for the parameters each one takes.

:optparam:`smtype{in}`
    Sparse matrix storage scheme.  Default ``0``, the symmetric skyline
    (``SMT_Skyline``).  The scheme must be compatible with the solver type; the
    permitted combinations are in :numref:`linsolvstoragecompattable`.

:param:`solverParams{s}`
    The parameters specific to the selected solver, listed in
    :numref:`sparsesolverparams`.

The available storage formats are the symmetric skyline (``SMT_Skyline``),
unsymmetric skyline (``SMT_SkylineU``), compressed column
(``SMT_CompCol``), dynamically growing compressed column
(``SMT_DynCompCol``), symmetric compressed column (``SMT_SymCompCol``),
dynamically growing compressed row (``SMT_DynCompRow``), the Spooles library
format (``SMT_SpoolesMtrx``), the PETSc library matrix representation
(``SMT_PetscMtrx``, a sparse serial or parallel matrix in AIJ format), and the
DSS-compatible representations (``SMT_DSS_*``).

Solver identifiers
------------------

.. _lstypetable:

.. list-table:: Linear solver types.
   :header-rows: 1
   :widths: 14 40 46

   * - ``lstype``
     - Solver
     - Compatible storage schemes (``smtype``)
   * - 0
     - Direct (``ST_Direct``), the default
     - 0, 1
   * - 1
     - Iterative Method Library (``ST_IML``)
     - 0, 1, 2, 3, 4, 5; recommended 4 if symmetric, otherwise 2
   * - 2
     - Spooles (``ST_Spooles``)
     - 6
   * - 3
     - PETSc (``ST_Petsc``)
     - 7
   * - 4
     - DirectSparseSolver (``ST_DSS``)
     - 8, 9, 10; recommended 8 if symmetric, otherwise 10
   * - 6
     - Intel MKL Pardiso (``ST_MKLPardiso``)
     - 2
   * - 7
     - SuperLU_MT (``ST_SuperLU_MT``)
     - 2
   * - 8
     - Pardiso Project (``ST_PardisoProjectOrg``)
     - 2

.. _linsolvstoragecompattable:

.. list-table:: Sparse matrix storage schemes.
   :header-rows: 1
   :widths: 14 44 42

   * - ``smtype``
     - Storage format
     - Compatible solvers (``lstype``)
   * - 0
     - ``SMT_Skyline`` — symmetric skyline, the default
     - 0, 1
   * - 1
     - ``SMT_SkylineU`` — unsymmetric skyline
     - 0, 1
   * - 2
     - ``SMT_CompCol`` — compressed column
     - 1, 6, 7, 8
   * - 3
     - ``SMT_DynCompCol`` — dynamically growing compressed column
     - 1
   * - 4
     - ``SMT_SymCompCol`` — symmetric compressed column
     - 1
   * - 5
     - ``SMT_DynCompRow`` — dynamically growing compressed row
     - 1
   * - 6
     - ``SMT_SpoolesMtrx`` — Spooles library format
     - 2
   * - 7
     - ``SMT_PetscMtrx`` — PETSc sparse serial/parallel matrix, AIJ format
     - 3
   * - 8
     - ``SMT_DSS_sym_LDL`` — DSS, symmetric LDL
     - 4
   * - 9
     - ``SMT_DSS_sym_LL`` — DSS, symmetric LL
     - 4
   * - 10
     - ``SMT_DSS_unsym_LU`` — DSS, unsymmetric LU
     - 4

.. _sparsesolverparams:

.. table:: Solver parameters.

   +------------------------+----+-----------------------------------------------------------+
   | Solver type            | id | Solver parameters and notes                               |
   +========================+====+===========================================================+
   | ST_Direct              | 0  | None.                                                     |
   +------------------------+----+-----------------------------------------------------------+
   | ST_IML                 | 1  | ``lstol`` (rn), ``lsiter`` (in), ``lsprecond`` (in),      |
   |                        |    | optional ``stype`` (in) and ``precondattributes`` (s).    |
   |                        |    | Included in OOFEM; compile with ``USE_IML``.              |
   +------------------------+----+-----------------------------------------------------------+
   | ST_Spooles             | 2  | Optional ``msglvl`` (in) and ``msgfile`` (s).             |
   |                        |    | See the Spooles documentation.                            |
   +------------------------+----+-----------------------------------------------------------+
   | ST_Petsc               | 3  | See the PETSc manual, and the run-time options below.     |
   +------------------------+----+-----------------------------------------------------------+
   | ST_DSS                 | 4  | Sparse direct solver included in OOFEM;                   |
   |                        |    | compile with ``USE_DSS``.                                 |
   +------------------------+----+-----------------------------------------------------------+
   | ST_MKLPardiso          | 6  | Requires Intel MKL Pardiso.                               |
   +------------------------+----+-----------------------------------------------------------+
   | ST_SuperLU_MT          | 7  | SuperLU for shared memory machines.                       |
   +------------------------+----+-----------------------------------------------------------+
   | ST_PardisoProjectOrg   | 8  | Requires the Pardiso solver from pardiso-project.org.     |
   +------------------------+----+-----------------------------------------------------------+

Iterative solvers from the IML library
--------------------------------------

:optparam:`stype{in}`
    Which iterative solver from the IML library to use: ``0`` (default) the
    conjugate gradient solver, ``1`` GMRES.

:param:`lstol{rn}`
    Maximum value of the residual after the final iteration.

:param:`lsiter{in}`
    Maximum number of iterations.

:param:`lsprecond{in}`
    Type of preconditioner; see :numref:`precondtable`.

:optparam:`precondattributes{s}`
    Optional preconditioner parameters; see :numref:`precondtable`.

.. _precondtable:

.. table:: Preconditioning summary.

   +--------------+----+--------------------+--------------------------------------------+
   | Precond type | id | Compatible storage | Description and parameters                 |
   +==============+====+====================+============================================+
   | IML_VoidPrec | 0  | all                | No preconditioning.                        |
   +--------------+----+--------------------+--------------------------------------------+
   | IML_DiagPrec | 1  | all                | Diagonal preconditioning.                  |
   +--------------+----+--------------------+--------------------------------------------+
   | IML_ILUPrec  | 2  | SMT_CompCol,       | Incomplete LU decomposition with no        |
   |              |    | SMT_DynCompCol     | fill-up.                                   |
   +--------------+----+--------------------+--------------------------------------------+
   | IML_ILUPrec  | 3  | SMT_DynCompRow     | Incomplete LU (ILUT) with fill-up.  The    |
   |              |    |                    | ``precondattributes`` are                  |
   |              |    |                    | ``droptol`` (rn), the dropping tolerance,  |
   |              |    |                    | and ``partfill`` (in), the level of        |
   |              |    |                    | fill-up.                                   |
   +--------------+----+--------------------+--------------------------------------------+
   | IML_ICPrec   | 4  | SMT_SymCompCol,    | Incomplete Cholesky with no fill-up.       |
   |              |    | SMT_CompCol        |                                            |
   +--------------+----+--------------------+--------------------------------------------+

PETSc run-time options
----------------------

With ``ST_PETSC`` several run-time options can be set on the command line, for
example:

.. code-block:: none

   -ksp_type [cg, gmres, bicg, bcgs]
   -pc_type  [jacobi, bjacobi, none, ilu, ...]
   -ksp_monitor
   -ksp_rtol <value>
   -ksp_view
   -ksp_converged_reason

These options override the defaults, because PETSc's ``KSPSetFromOptions()``
routine is called after any other customization routine.

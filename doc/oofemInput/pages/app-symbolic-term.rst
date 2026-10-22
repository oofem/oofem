.. _SymbolicTermSec:

SymbolicTerm expression language
================================

``SymbolicTerm`` lets a term of the weak form be written directly in the input
file, in a small matrix language, instead of existing as a C++ class.  The
expressions are parsed and compiled into bytecode at initialization and then
executed from a pre-allocated pool during the assembly loops, so a term written
this way costs little more than a compiled one.

This page is the reference for what may appear inside those expressions.  For
the ``Variable``, ``Term`` and ``Integral`` records that frame them, see
:ref:`TermsSec`; for the formulation the module implements, see the
**Theory Manual**.

Left- and right-hand side expressions
-------------------------------------

A term has to supply both a left-hand side matrix and a right-hand side vector,
and ``SymbolicTerm`` keeps these as two separate expressions.

:elemparam:`lexpression{s}`
    The linearization: the tangent or stiffness matrix contribution of the
    term.  It must evaluate to a two-dimensional matrix.

:elemparam:`rexpression{s}`
    The internal force or residual vector contribution.  It evaluates to a
    column matrix, which the term converts to a vector.

A pure source term has no unknown to depend on, so its
:param:`lexpression` is zero.

Record syntax
-------------

``SymbolicTerm`` inherits from ``GenericCellTerm`` and adds the two expression
strings.

.. record::

   :descitem:`SymbolicTerm` :elemparam:`num{in}` :elemparam:`variable{s}`
   :elemparam:`testvariable{s}` :elemparam:`mmode{et}`
   :elemparam:`lexpression{s}` :elemparam:`rexpression{s}`
   :optelemparam:`ctype{in}` :optelemparam:`atype{in}`
   :optelemparam:`uvmt{in}`

**Parameters**

:param:`variable{s}`
    Name of the primary unknown field of the term, for example ``"u"`` or
    ``"p"``.  It must match a ``Variable`` record of the problem; see
    :ref:`VariablesSec`.

:param:`testvariable{s}`
    Name of the test (weighting) function, for example ``"w"`` or ``"dp"``.

:param:`mmode{et}`
    Material mode used to evaluate the term, a ``MaterialMode`` value.  The
    name may be written instead of the number, with the leading underscore
    dropped, so ``mmode "PlaneStrain"`` and ``mmode 7`` are the same thing.
    Prefer the name: the integer values are an implementation detail of
    ``src/core/materialmode.h``.

    .. important::

       An enum name has to be **quoted**.  Names are matched case-sensitively,
       and the parser lower-cases everything outside double quotes (see
       :ref:`anatomy-of-an-input-file`), so a bare ``mmode PlaneStrain``
       arrives as ``planestrain`` and is not recognised.

:param:`lexpression{s}`, :param:`rexpression{s}`
    The expressions described above.  Both are strings, so they must be
    enclosed in double quotes — which also stops the parser from lower-casing
    them, and matters because the functor names are mixed case.

:optparam:`ctype{in}`, :optparam:`atype{in}`, :optparam:`uvmt{in}`
    Attributes of the parent ``GenericCellTerm``, accepted here as well.

Naming a test function as an unknown
------------------------------------

A ``Variable`` record that is a test function should say so with
:param:`dualto`, naming the unknown it weights:

.. code-block:: none

   Variable name "pw"  interpolation "feilin" type 0 quantity "Pressure" \
       size 1 dofs 1 11
   Variable name "dpw" interpolation "feilin" type 0 quantity "Pressure" \
       size 1 dofs 1 11 dualto "pw"

Without it the two records are indistinguishable — they may even carry
different interpolations, which is what a non-symmetric Petrov-Galerkin
formulation looks like — and the assembly cannot tell which field to read the
nodal unknowns from when it pushes state to the material.

:param:`dualto` is required only where a deck names a test function as some
term's :param:`variable`, which happens for pure source terms, since those
have no unknown to depend on.

.. note::

   :param:`quantity` is what the assembly matches against the material's
   declared state layout, so it must say what the field *is*;
   :param:`dofs` says only where that field's unknowns are numbered.  The two
   are checked for agreement at initialization, and a mismatch warns rather
   than errors, since a deck may legitimately place a field on an
   unconventional DOF.  See :ref:`VariablesSec`.

Expression syntax
-----------------

Several statements may be chained with semicolons, the last one giving the
value of the expression.

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Feature
     - Syntax
     - Example
   * - Assignment
     - ``variable = expression``
     - ``K = MDer(gp, ts, MatResponseMode::TangentStiffness);``
   * - Arithmetic
     - ``+``, ``-``, ``*``, unary ``-``
     - ``A + B``, ``A * 2.0``, ``-A``
   * - Matrix transpose
     - ``.T``
     - ``Grad_s(w,gp).T``
   * - Scalar literal
     - a plain number
     - ``1.0``, ``5``, ``-2.4e-7``
   * - Matrix literal
     - ``[[r1c1, r1c2], [r2c1, r2c2]]``
     - ``I = [[1, 0], [0, 1]];``
   * - Comparison
     - ``==``, ``!=``, ``>``, ``<``, ``>=``, ``<=``
     - ``a > 5``
   * - Logical operators
     - ``&&``, ``||``
     - ``(a > 5) && (b < 2)``
   * - Ternary if
     - ``if(condition, then, else)``
     - ``val = if(a > 0, 1, -1);``
   * - Slicing
     - ``matrix[row, col]``
     - ``val = K[0, 1];``, ``row = K[0, :];``

Context variables
-----------------

Within a term evaluated on an element, the following names are always in scope.

``gp``
    The current Gauss integration point (``GaussPoint*``).

``ts``
    The current time step (``TimeStep*``).

``cell``
    The element being evaluated (``MPElement*``).

The field variables of the problem are also in scope under their own names —
``u``, ``p``, ``w``, ``dp`` and so on — as ``Variable*`` objects.

Built-in functions
------------------

.. list-table::
   :header-rows: 1
   :widths: 27 52 21

   * - Signature
     - Description
     - Returns
   * - ``Grad_s(var, gp)``
     - Symmetric gradient (strain) operator of a vector field; its shape
       depends on the material mode.
     - matrix
   * - ``Grad(var, gp)``
     - Gradient operator of a scalar field, one row per spatial direction.
     - matrix
   * - ``Div(var, gp)``
     - Divergence operator of a vector field.
     - row vector
   * - ``N(var, gp)``
     - Shape function interpolation matrix.
     - matrix
   * - ``MDer(gp, ts, mode)``
     - Characteristic **matrix** of the material for the given
       ``MatResponseMode`` — tangent, conductivity, permeability and so on.
     - matrix
   * - ``MVec(gp, ts, mode)``
     - Characteristic **vector** of the material for the given mode — stress,
       flux and so on.
     - column vector
   * - ``MProp(gp, ts, mode)``
     - Characteristic **scalar** of the material for the given mode —
       capacity, Biot constant, density and so on.
     - scalar
   * - ``Sig(var, gp, ts)``
     - Stress vector; shorthand for
       ``MVec(gp, ts, MatResponseMode::Stress)``.
     - column vector
   * - ``Sig_dev(var, gp, ts)``
     - Deviatoric stress vector; shorthand for
       ``MVec(gp, ts, MatResponseMode::DeviatoricStress)``.
     - column vector
   * - ``ru(var, cell, ts)``
     - "Reads unknowns": nodal values of a field on the current element.
     - column vector
   * - ``rv(var, cell, ts)``
     - Nodal velocities (rates) of a field on the current element.
     - column vector
   * - ``eval(var, gp, ts)``
     - Value of a field interpolated at the integration point.
     - scalar
   * - ``vcat(a, b, ...)``
     - Vertical concatenation of matrices or vectors.
     - matrix
   * - ``LumpMatrix(m)``
     - HRZ lumping of a consistent mass matrix.
     - matrix
   * - ``print(x)``
     - Debugging aid: prints its argument and returns it.
     - its argument

Response modes
~~~~~~~~~~~~~~

A ``MatResponseMode`` may be named — ``MatResponseMode::TangentStiffness``,
``MatResponseMode::Permeability``, and so on for every enumerator — or given as
its integer value.

.. important::

   Prefer the names.  The integer values are an implementation detail, and a
   deck that encodes them would change meaning silently if the enumeration were
   ever renumbered.

The material queries are reads, not evaluations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``MDer``, ``MVec``, ``MProp``, ``Sig`` and ``Sig_dev`` return values the
material has **already** computed.  Before any term is evaluated, the solver
pushes the generalized state of every integration point to its material once
per equilibrium iteration, and the material does all of its constitutive work
there and caches the results.

Two consequences when writing a term:

* A term neither assembles nor passes the state.  Earlier decks did this by
  hand, building the state with ``vcat`` and handing it to ``MVec``; that is no
  longer needed, and ``MVec`` no longer takes a state argument.  What the
  material expects, and in what order, is declared by the material itself.
* The result does not depend on evaluation order.  The tangent and the residual
  are assembled in separate sweeps over the integration points; previously only
  the residual sweep supplied a state, so a tangent query consumed whatever the
  previous sweep had left behind.

The ``var`` argument of ``Sig`` and ``Sig_dev`` is retained for compatibility.
It is checked against the field that actually supplies the strain on the cell,
but it no longer selects what is read.

For a field to be available to be pushed at all, some term on the cell must
name it as its **unknown** :param:`variable`; test functions are excluded
through :param:`dualto`.

.. seealso::

   The state push and the declared state layout — what a material advertises
   through ``giveStateVariableIDs``, the ``(field, operator)`` pairs and the
   available state operators — are described in the **Theory Manual**, in the
   symbolic MPM chapter.  The design rationale for the split and for the
   push/pull contract is in ``doc/unified_material_interface.md``.

Examples
--------

Deviatoric stress term, solid mechanics
    The standard deviatoric stiffness and internal force vector:

    .. code-block:: none

       SymbolicTerm 1 variable "u" testvariable "w" mmode 7 \
           lexpression "Grad_s(w,gp).T * MDer(gp,ts,MatResponseMode::DeviatoricStiffness) * Grad_s(u,gp)" \
           rexpression "Grad_s(w,gp).T * Sig_dev(u,gp,ts)"

Incompressibility constraint, mixed u-p formulation
    The volumetric constraint relating the displacement or velocity divergence
    to the pressure:

    .. code-block:: none

       SymbolicTerm 2 variable "p" testvariable "w" mmode 7 \
           lexpression "Div(w,gp).T * N(p,gp)" \
           rexpression "Div(w,gp).T * N(p,gp) * ru(p, cell, ts)"

Pressure stabilization term
    A stabilization term for the pressure field, with a small scaling factor:

    .. code-block:: none

       SymbolicTerm 4 variable "p" testvariable "dp" mmode 7 ctype 27 uvmt 1 \
           lexpression "N(dp,gp).T * 2.4e-7 * N(p,gp)" \
           rexpression "N(dp,gp).T * 2.4e-7 * N(p,gp) * ru(p, cell, ts)"

A complete weak form written this way — the five terms of the Cook membrane
u-p formulation — is in
``tests/regression/mpm/mpms_cook2_u2p1.in``, and the same formulation is
derived in the **Theory Manual**.

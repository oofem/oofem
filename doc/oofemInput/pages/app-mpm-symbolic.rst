Symbolic MPM records
====================

The symbolic ``mpm`` module lets a problem be defined by writing its weak form
in terms of *variables*, *terms* and *integrals*.  The counts of these records
are declared by :param:`nvariables`, :param:`nterms` and :param:`nintegrals` in
the analysis record, and the records themselves follow the export module
records; see :ref:`AnalysisRecord` and :ref:`StationaryMPMSProblem`.

The three record kinds build on one another:

* a **Variable** is an unknown field, or the test (weighting) function of one;
* a **Term** is an integrand, written in terms of a variable and a test
  variable;
* an **Integral** integrates one term over a set of elements.

.. _VariablesSec:

Variable records
----------------

Describes a variable, that is a field appearing in the weak form.

.. record::

   :descitem:`Variable` :elemparam:`name{s}`
   :elemparam:`interpolation{s}` :elemparam:`type{et}`
   :elemparam:`quantity{et}` :elemparam:`size{in}`
   :elemparam:`dofs{ia}` :optelemparam:`dualto{s}`

**Parameters**

:param:`name{s}`
    Name of the variable.  Terms refer to variables by this name.

:param:`interpolation{s}`
    Interpolation of the variable.  The supported values are:

    ``feiconst``
        Constant interpolation.
    ``feilin``
        Linear interpolation.
    ``feiquad``
        Quadratic interpolation.
    ``fei1dlin_x``, ``fei1dlin_y``, ``fei1dlin_z``
        Linear interpolation on a 1D cell aligned with the given global axis.

:param:`type{et}`
    Rank of the variable: ``scalar`` (``0``) or ``vector`` (``1``).  Either the
    name or the number may be written.

:param:`quantity{et}`
    What the variable physically is, as a ``FieldType`` value.  Either the name
    or the number may be written, and the ``FT_`` prefix may be dropped — so
    ``quantity "Temperature"`` and ``quantity 5`` are the same thing.  The
    values are defined in ``src/core/field.h``:

    .. list-table::
       :header-rows: 1
       :widths: 12 34 54

       * - Value
         - Name
         - Meaning
       * - 0
         - ``FT_Unknown``
         - Unspecified
       * - 1
         - ``FT_Velocity``
         - Velocity field
       * - 2
         - ``FT_Displacements``
         - Displacement field
       * - 3
         - ``FT_VelocityPressure``
         - Combined velocity and pressure field
       * - 4
         - ``FT_Pressure``
         - Pressure field
       * - 5
         - ``FT_Temperature``
         - Temperature field
       * - 6
         - ``FT_HumidityConcentration``
         - Humidity or concentration field
       * - 7
         - ``FT_TransportProblemUnknowns``
         - Generic transport problem unknowns
       * - 8
         - ``FT_TemperatureAmbient``
         - Ambient temperature
       * - 9
         - ``FT_EigenStrain``
         - Eigenstrain field
       * - 10
         - ``FT_VOF``
         - Volume of fluid
       * - 11
         - ``FT_Pressure2``
         - Pressure of the second phase
       * - 12
         - ``FT_Concentration1``
         - Mass concentration of the first constituent
       * - 13
         - ``FT_Concentration2``
         - Mass concentration of the second constituent

:param:`size{in}`
    Size, that is the dimension, of the variable.

:param:`dofs{ia}`
    Physical meaning of the variable's DOFs, as ``DofIDItem`` values; see the
    table under :ref:`NodeElementSideRecords`.  The size of the array should
    equal :param:`size`.

:optparam:`dualto{s}`
    Name of the unknown field that this variable is the test (weighting)
    function of.  This is the *only* thing that distinguishes an unknown from
    its weighting function: the two are declared as two Variable records that
    differ by name, and they may even carry different interpolations, which is
    what a non-symmetric Petrov-Galerkin formulation looks like.

    Declaring it matters because assembly has to read the nodal unknowns
    through the unknown field's interpolation rather than through its weighting
    function's.

**Example**

.. code-block:: none

   Variable name "u" interpolation "feilin" type 1 quantity "Temperature" \
       size 1 dofs 1 10
   Variable name "w" interpolation "feilin" type 1 quantity "Temperature" \
       size 1 dofs 1 10 dualto "u"

Here ``u`` is the temperature unknown and ``w`` its test function.  See
``tests/regression/mpm/mpms06.in`` for a complete input file.

.. _TermsSec:

Term records
------------

A term is what an integral in the weak form integrates.

.. record:: Record syntax common to all terms

   :descitem:`TermType` :elemparam:`num{in}` :elemparam:`variable{s}`
   :elemparam:`testvariable{s}` :elemparam:`mmode{et}`

**Parameters**

:param:`variable{s}`
    Name of the unknown variable (field) of the term.

:param:`testvariable{s}`
    Name of the test variable (field) of the term.

:param:`mmode{et}`
    Material mode used to evaluate the term, as a ``MaterialMode`` value; see
    ``src/core/materialmode.h``.

Individual term types may introduce further parameters of their own.

Notation
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Symbol
     - Meaning
   * - :math:`q`, :math:`p`
     - Scalar valued functions
   * - :math:`\mathbf{w}`, :math:`\mathbf{u}`
     - Vector valued functions

Supported term types
~~~~~~~~~~~~~~~~~~~~

``BTSigmaTerm``
    :math:`\int_\Omega \nabla^s \mathbf{w}\ \mathbf{\sigma}(\nabla^s
    \mathbf{u})`, where :math:`\mathbf{\sigma}` is the (nonlinear) operator
    evaluated by the constitutive model.

    Supported material modes: ``_3dMat``, ``_3dUP``, ``_2dUP``,
    ``_PlaneStress``.

    Optional parameter: :optelemparam:`lhsmatmode{in}`.

``BTamNTerm``
    :math:`\int_\Omega \nabla^s \mathbf{w}\ a\mathbf{m}\ p`, where ``a`` is a
    material parameter defined by the constitutive model whose meaning is fixed
    by the :param:`atype` parameter, and
    :math:`\mathbf{m}^T=[1,1,1,0,0,0]^T`.

    Note that :math:`\mathbf{\sigma}_{D}=\mathbf{\sigma}-\mathbf{m}p`.

``NTamTBTerm``
    The transposed version of ``BTamNTerm``:
    :math:`\int_\Omega w\ \alpha\mathbf{m}\ \nabla^s \mathbf{u}`.

``NTcN``
    :math:`\int_\Omega q c p`, where ``c`` is a constant defined by the
    material model whose meaning is fixed by the :param:`ctype` parameter.

``NTfTerm``
    :math:`\int_\Omega \mathbf{w}\cdot\bar{\mathbf{t}}`, where
    :math:`\bar{\mathbf{t}}` is a given flux vector supplied by the
    :elemparam:`flux{ra}` parameter.

``SymbolicTerm``
    A term whose integrand is written directly in the input file, in a symbolic
    mathematical language that is compiled to bytecode during initialization.

    It adds two mandatory string parameters: :elemparam:`lexpression{s}`, the
    left-hand side expression evaluated to the tangent/stiffness matrix
    contribution, and :elemparam:`rexpression{s}`, the right-hand side
    expression evaluated to the internal force or residual vector
    contribution.

    The expression language — its functions ``Grad_s``, ``Grad``, ``Div``,
    ``N``, ``Sig``, ``Sig_dev``, ``MDer``, ``MVec``, ``MProp`` and the rest, and
    the full list of attributes — is documented separately in
    ``doc/symbolic_term.md``.

.. _IntegralsSec:

Integral records
----------------

An integral integrates one term over a set of elements.

.. record::

   :descitem:`Integral` :elemparam:`num{in}` :elemparam:`domain{in}`
   :elemparam:`set{in}` :elemparam:`term{in}`

**Parameters**

:param:`domain{in}`
    Number of the domain in which the integral is defined.

:param:`set{in}`
    Number of the set of elements over which the integration is performed; see
    :ref:`SetRecords`.

    .. note::

       Boundary integrals require the boundary elements to be present.

:param:`term{in}`
    Number of the term to integrate.

The integrals contributing to the left- and right-hand sides of the problem are
selected by the :param:`lhsterms` and :param:`rhsterms` arrays of the analysis
record; see :ref:`StationaryMPMSProblem`.

.. note::

   Despite their names, :param:`lhsterms` and :param:`rhsterms` list
   *integral* numbers, not term numbers.

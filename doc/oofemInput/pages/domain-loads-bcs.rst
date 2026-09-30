.. _LoadBoundaryInitialConditions:

Load and boundary condition records
===================================

These records describe the boundary conditions, including loads and contact
conditions.  They are all counted by :param:`nbc` in the
:ref:`ComponentsSizeRecord`.

.. record:: Record syntax common to all boundary conditions

   :descitem:`EntType` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :optelemparam:`set{in}`
   :optelemparam:`valType{in}` :optelemparam:`dofs{ia}`
   :optelemparam:`isImposedTimeFunction{in}`

The order of the records is arbitrary and the number is given by :param:`num`.
The numbering must start at one and end at *n*, where *n* is the number of
records.

**Common parameters**

:param:`loadTimeFunction{in}`
    Number of the time function whose value multiplies every component of the
    condition, so describing its variation in time.  See
    :ref:`TimeFunctionsRecords`.

:optparam:`set{in}`
    Number of the set (see :ref:`SetRecords`) of entities the condition is
    applied to.  This is the usual way of attaching a condition; the
    alternative is to reference the condition from a dof manager's ``bc`` or
    ``load`` array, or from an element's ``bodyLoads`` or ``boundaryLoads``
    array.

:optparam:`dofs{ia}`
    Which DOFs the condition acts upon, using the ``DofIDItem`` values listed
    under :ref:`NodeElementSideRecords`.  Not relevant to every condition
    type.

:optparam:`valType{in}`
    Physical meaning of the condition's value, which some formulations
    require: ``1`` temperature, ``2`` force/traction, ``3`` pressure,
    ``4`` humidity, ``5`` velocity, ``6`` displacement.

:optparam:`isImposedTimeFunction{in}`
    Number of a time function controlling whether the condition exists at all.
    A nonzero value means the condition is active, zero means it is inactive —
    it does not exist — at that time.  By default the condition applies at all
    times.

Dirichlet boundary conditions
-----------------------------

.. _BoundaryConditionRecord:

BoundaryCondition
~~~~~~~~~~~~~~~~~

Prescribes a value on the selected DOFs.  The physical meaning of the value
follows from the DOF it is applied to.

.. record::

   :descitem:`BoundaryCondition` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :optelemparam:`set{in}`
   :optelemparam:`dofs{ia}`
   (:elemparam:`values{ra}` | :elemparam:`prescribedvalue{rn}` |
   :elemparam:`d{rn}`)

**Parameters**

:param:`values{ra}`
    Vector of prescribed values, one per DOF listed in :param:`dofs`.  This is
    the preferred form.

:param:`prescribedvalue{rn}`
    A single value prescribed on *all* the selected DOFs.  Used only when
    :param:`values` is absent.

:param:`d{rn}`
    An alternative spelling of :param:`prescribedvalue`, retained for
    compatibility.  It is read only when neither :param:`values` nor
    :param:`prescribedvalue` is present.

**Example**

.. code-block:: none

   BoundaryCondition 1 loadTimeFunction 1 dofs 2 1 2 values 2 0.0 0.0 set 2

.. _PrescribedGradientRecord:

PrescribedGradient
~~~~~~~~~~~~~~~~~~

Prescribes :math:`v_i = d_{ij}(x_j-\bar{x}_j)`, or
:math:`s = d_{1j}(x_j - \bar{x}_j)` for a scalar field, where :math:`v_i` are
the primary unknowns, :math:`x_j` the coordinate of the node, :math:`\bar x` is
:param:`ccoord` and :math:`d` is :param:`gradient`.

This is the typical boundary condition of multiscale analysis, where
:math:`d = \partial_x s` would be a macroscopic gradient at the integration
point — that is, it is a boundary condition for prolongation.  It is also
convenient for testing an arbitrary specimen in shear.

.. record::

   :descitem:`PrescribedGradient` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`gradient{rm}`
   :optelemparam:`set{in}` :optelemparam:`ccoord{ra}`

**Parameters**

:param:`gradient{rm}`
    The prescribed gradient :math:`d`.

:optparam:`ccoord{ra}`
    Centre coordinate :math:`\bar x`.  Defaults to zero.

Mixed gradient / pressure conditions
------------------------------------

These conditions ensure that the deviatoric gradient and the pressure are at
least weakly fulfilled on the prescribed domain.  They are used for the
computational homogenization of incompressible flow or elasticity problems.

MixedGradientPressureDirichlet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prescribes
:math:`v_i = d_{\mathrm{dev},ij}(x_j-\bar{x}_j) + d_\mathrm{vol}(x_i-\bar{x}_i)`
together with a pressure :math:`p`, where :math:`v_i` are the primary unknowns,
:math:`x_j` the coordinate of the node, :math:`\bar x` is :param:`ccoord` and
:math:`d_\mathrm{dev}` is :param:`devGradient`.

An additional unknown appears, :math:`d_\mathrm{vol}`, which the boundary
condition handles itself with no control from the input file; it is related to
the applied pressure.  This is useful in multiscale computations of RVEs that
behave incompressibly, typically fluid problems, and is convenient for testing
an arbitrary specimen in shear with a free volumetric part — in which case the
pressure is set to zero.

.. record::

   :descitem:`MixedGradientPressureDirichlet` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`devGradient{ra}`
   :elemparam:`pressure{rn}` :optelemparam:`set{in}`
   :optelemparam:`ccoord{ra}`

**Parameters**

:param:`devGradient{ra}`
    The prescribed deviatoric gradient.

:param:`pressure{rn}`
    The prescribed pressure.

:optparam:`ccoord{ra}`
    Centre coordinate :math:`\bar x`.  Defaults to zero.

MixedGradientPressureNeumann
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prescribes a constant, unknown deviatoric stress tensor along the specified
boundaries.  The additional unknown
:math:`\boldsymbol{\sigma}_\mathrm{dev}` is handled by the boundary condition
itself, with no control from the input file, and the input
:param:`devGradient` is fulfilled weakly, homogenized over the element sides.

As with the Dirichlet variant the volumetric gradient is free.  This is useful
in multiscale computations of RVEs that behave incompressibly, typically fluid
problems, in which case the element sides should cover the entire RVE boundary.
Symmetry is not assumed, so rigid body rotations are removed, but translations
must be prescribed separately.

.. record::

   :descitem:`MixedGradientPressureNeumann` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`devGradient{ra}`
   :elemparam:`pressure{rn}` :optelemparam:`set{in}`

MixedGradientPressureWeakPeriodic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prescribes a periodic constant, unknown stress tensor along the specified
boundaries.  With :param:`order` set to 1 the results are the same as for the
Neumann boundary condition.

.. note::

   The record keyword is ``mixedgradientpressureweakperiodic``.

.. record::

   :descitem:`MixedGradientPressureWeakPeriodic` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`devGradient{ra}`
   :elemparam:`pressure{rn}` :elemparam:`order{in}`
   :optelemparam:`set{in}`

**Parameters**

:param:`order{in}`
    Order of the periodic stress approximation.

Nodal loads
-----------

.. _NodalLoadRecord:

NodalLoad
~~~~~~~~~

A concentrated nodal load.

.. record::

   :descitem:`NodalLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`components{ra}`
   :optelemparam:`set{in}` :optelemparam:`dofs{ia}`
   :optelemparam:`cstype{in}`

**Parameters**

:param:`components{ra}`
    Components of the nodal load vector.  Its size corresponds to the total
    number of nodal DOFs, and the *i*-th value to the *i*-th DOF of the
    associated dof manager.

:optparam:`cstype{in}`
    Coordinate system of the load: ``0`` global, ``1`` the entity-specific
    local system (the default).

PrescribedTractionPressureBC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A pressure boundary condition of Dirichlet type arising from prescribed
tractions.  In the CBS algorithm formulation a prescribed traction boundary
condition leads indirectly to a pressure boundary condition in the
corresponding nodes, and this condition implements that pressure condition.
Its value follows from the applied tractions, which should be specified on
element edges or surfaces using suitable boundary loads.

.. record::

   :descitem:`PrescribedTractionPressureBC` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :optelemparam:`set{in}`

Body loads
----------

DeadWeight
~~~~~~~~~~

Dead weight loading applied to the element volume, for structural elements.
For transport problems it represents the internal source, that is the rate of
heat generated per unit volume.

The magnitude of the load for the *i*-th DOF is the product of the material
density, the corresponding volume, and the *i*-th member of
:param:`components`.

.. record::

   :descitem:`DeadWeight` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`components{ra}`
   :optelemparam:`set{in}`

StructTemperatureLoad
~~~~~~~~~~~~~~~~~~~~~

Temperature loading imposed on elements.  The members of
:param:`components` represent the change of temperature, or the change of the
temperature gradient, corresponding to the individual element strain
components.  See the Element Library Manual for details.

.. record::

   :descitem:`StructTemperatureLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`components{ra}`
   :optelemparam:`set{in}`

StructEigenstrainLoad
~~~~~~~~~~~~~~~~~~~~~

Prescribes an eigenstrain — a stress-free strain — on a structural element.

:param:`components` is defined in the global coordinate system, and the number
of components corresponds to the material mode: plane stress has three
components, 3D has six.

.. record::

   :descitem:`StructEigenstrainLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`components{ra}`
   :optelemparam:`set{in}`

**Periodic boundary conditions through eigenstrains**

Periodic boundary conditions can be imposed using eigenstrains together with
master-slave nodes.  Consider the decomposition of strain into an average and a
fluctuating part,

.. math::

   \boldsymbol{\varepsilon}(\boldsymbol{x}) =
   \langle \boldsymbol{\varepsilon} \rangle +
   \boldsymbol{\varepsilon}^*(\boldsymbol{x}),

where :math:`\langle \boldsymbol{\varepsilon} \rangle` can be imposed as an
eigenstrain over the domain, and the solution then gives the fluctuating part
:math:`\boldsymbol{\varepsilon}^*(\boldsymbol{x})`.  The master-slave nodes
have to interconnect the opposing boundary nodes of the unit cell.

.. _BoundaryLoads:

Boundary loads
--------------

.. _loadtype-values:

Load types
~~~~~~~~~~

The edge and surface loads share the :param:`loadType` parameter, which
distinguishes the kind of boundary condition.  The supported values are defined
in ``src/core/bctype.h``:

``loadType`` = 2
    Prescribed flux input, that is a Neumann boundary condition.

``loadType`` = 3
    Uniform distributed load, or the convection (Newton) boundary condition.
    In the latter case :param:`components` holds the environmental values — the
    temperature of the environment — corresponding to the element unknowns, and
    the :param:`properties` dictionary should hold the transfer (convection)
    coefficient, assumed constant, under the key ``a``.

``loadType`` = 7
    Radiative boundary condition (Stefan-Boltzmann).  It requires the
    emissivity :math:`\varepsilon\in\langle 0,1\rangle`, and
    :param:`components` holds the environmental values — the temperature of the
    environment.  Default units are degrees Celsius; set
    :param:`temperOffset` to 0 to compute in Kelvin.

If the boundary condition corresponds to a distributed force load,
:param:`components` holds the components of the distributed load corresponding
to the element unknowns.  If it corresponds to a prescribed flux input,
:param:`components` holds the components of the prescribed input flux.

Shared parameters
~~~~~~~~~~~~~~~~~

:param:`components{ra}`
    The load components.  The load is specified for all the DOFs of the object
    it is associated with.

:optparam:`dofexcludemask{ia}`
    For some boundary condition types a zero load value does not mean that the
    load is not applied — Newton's type of condition, for instance — so a mask
    is needed to exclude specific DOFs.  It should have the same size as
    :param:`components` and is filled with zeroes by default; where a value is
    nonzero, the corresponding component is set to zero and the load is not
    applied to that DOF.

:optparam:`csType{in}`
    Coordinate system of the load: ``0`` global (the default), ``1`` the
    entity-specific local system.

:optparam:`properties{dc}`
    Dictionary of additional properties required by the load type, for example
    the convection coefficient under the key ``a``.

:optparam:`propertytf{dc}`
    Dictionary giving, for selected properties, the number of a time function
    describing the property's variation in time.  A time function must be
    registered under the same key as in :param:`properties`.  The property
    value is then the product of the value from :param:`properties` and the
    time function evaluated at the given time.  Where no time function is
    provided for a property, a unit constant function is assumed.

:optparam:`temperOffset{rn}`
    Temperature offset used by the radiative boundary condition; set it to 0 to
    compute in Kelvin rather than degrees Celsius.

:optparam:`reference{}`
    A flag marking the load as part of the non-scaling reference load vector;
    see :param:`refloadmode` under :ref:`NonLinearStatic`.

ConstantEdgeLoad
~~~~~~~~~~~~~~~~

A constant edge load or boundary condition.

.. record::

   :descitem:`ConstantEdgeLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`loadType{in}`
   :elemparam:`components{ra}` :optelemparam:`set{in}`
   :optelemparam:`dofexcludemask{ia}` :optelemparam:`csType{in}`
   :optelemparam:`properties{dc}` :optelemparam:`propertytf{dc}`
   :optelemparam:`temperOffset{rn}`

ConstantSurfaceLoad
~~~~~~~~~~~~~~~~~~~

A constant surface load or boundary condition.

.. record::

   :descitem:`ConstantSurfaceLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`loadType{in}`
   :elemparam:`components{ra}` :optelemparam:`set{in}`
   :optelemparam:`dofexcludemask{ia}` :optelemparam:`csType{in}`
   :optelemparam:`properties{dc}` :optelemparam:`propertytf{dc}`
   :optelemparam:`temperOffset{rn}` :optelemparam:`loadoffset{rn}`

**Additional parameter**

:optparam:`loadoffset{rn}`
    Offset of the load application point from the surface, along the surface
    normal.

LinearEdgeLoad
~~~~~~~~~~~~~~

A linearly varying edge load.  :param:`components` holds the load components
for the unknowns at the beginning of the edge, followed by the values valid at
the end of the edge.  :param:`csType` and :param:`loadType` mean the same as
for ``ConstantEdgeLoad``.

.. record::

   :descitem:`LinearEdgeLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`loadType{in}`
   :elemparam:`components{ra}` :optelemparam:`set{in}`
   :optelemparam:`dofexcludemask{ia}` :optelemparam:`csType{in}`
   :optelemparam:`formulation{in}` :optelemparam:`sc{ra}`
   :optelemparam:`ec{ra}`

**Additional parameters**

:optparam:`formulation{in}`
    Whether the load is given in terms of the entity's local parametric
    coordinates or of global coordinates.

:optparam:`sc{ra}`, :optparam:`ec{ra}`
    Start and end coordinates used by the global formulation.

UserDefinedBoundaryLoad
~~~~~~~~~~~~~~~~~~~~~~~

A user-defined boundary flux, whose intensity is given by a configured function.
``PythonExpression`` (see :ref:`TimeFunctionsRecords`) makes this completely
configurable.

.. note::

   The record keyword is ``usrdefboundaryload``.

.. record::

   :descitem:`usrdefboundaryload` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`dofs{ia}`
   :elemparam:`set{in}` :elemparam:`intensityfunction{in}`
   :elemparam:`components{ra}` :optelemparam:`geomtype{in}`
   :optelemparam:`approxorder{in}`

**Parameters**

:param:`intensityfunction{in}`
    Number of the function defining the flux intensity.  The function is
    defined in the global coordinate system as a function of position — the
    ``x`` array variable — and time — the ``t`` variable.

:optparam:`geomtype{in}`
    Type of geometry the flux is applied to: ``4`` for a surface flux (the
    default), ``3`` for an edge flux.

:optparam:`approxorder{in}`
    Approximate order of the intensity function, used to set up the
    integration rule on the boundary.

:param:`dofs{ia}`
    The DOFs the flux is applied to.

:param:`set{in}`
    The set of element boundary entities the flux is applied to.

.. note::

   :param:`loadTimeFunction` and :param:`components` have no effect here, but
   they still have to be provided.  Write ``components 0`` for an empty array.

InteractionLoad
~~~~~~~~~~~~~~~

A fluid-pressure-induced load in the solid part of the fluid-structure problem
(see :ref:`fluidstructureproblem`).  :param:`ndofs`, :param:`csType` and
:param:`loadType` mean the same as for ``LinearEdgeLoad``.

:param:`components` holds the load components for the unknowns at the beginning
of the edge — :param:`ndofs` values — followed by the values valid at the end of
the edge, another :param:`ndofs` values.  The load should be defined in the
global coordinate system (:param:`csType` = 0), since it acts in the direction
normal to the edge.

.. record::

   :descitem:`InteractionLoad` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`ndofs{in}`
   :elemparam:`loadType{in}` :elemparam:`components{ra}`
   :elemparam:`coupledparticles{ia}` :optelemparam:`set{in}`
   :optelemparam:`csType{in}`

**Parameters**

:param:`coupledparticles{ia}`
    The PFEM particles from the fluid part of the problem that supply the fluid
    pressure.

Constraints and couplings
-------------------------

LinearConstraintBC
~~~~~~~~~~~~~~~~~~

Implements a linear constraint of the form :math:`\sum_i w_i r_i = c`, where
:math:`r_i` are the unknowns of the DOFs selected by :param:`dofmans` and
:param:`dofs`.

The condition is introduced as an additional stationary condition through a
Lagrange multiplier, which is itself an additional degree of freedom introduced
by this boundary condition.

.. record::

   :descitem:`LinearConstraintBC` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`weights{ra}`
   :elemparam:`dofmans{ia}` :elemparam:`dofs{ia}`
   :elemparam:`rhs{rn}` :elemparam:`lhstype{ia}`
   :elemparam:`rhstype{ia}` :optelemparam:`weightsLtf{ia}`
   :optelemparam:`rhsLtf{in}`

**Parameters**

:param:`dofmans{ia}`, :param:`dofs{ia}`
    Dof manager numbers and the corresponding DOF indices identifying the
    participating DOFs.

:param:`weights{ra}`
    Weights of the participating DOFs.  By default all weights are 1.

:optparam:`weightsLtf{ia}`
    Number of a time function associated with each individual weight; each
    weight is multiplied by the value that function returns.

:param:`rhs{rn}`
    The constant :math:`c`.

:optparam:`rhsLtf{in}`
    Number of the time function multiplying :param:`rhs`.  Defaults to 1.

:param:`lhstype{ia}`, :param:`rhstype{ia}`
    The characteristic components this condition contributes to, as
    ``CharType`` enumeration values.  The left-hand side contribution is
    assembled into the terms identified by :param:`lhstype`, the right-hand
    side contribution into the term identified by :param:`rhstype`.  Multiple
    values are allowed, which makes it possible to select all variants of the
    stiffness matrix, for example.

.. important::

   The :param:`dofmans`, :param:`dofs`, :param:`weights` and
   :param:`weightsLtf` arrays must all have the same size.

.. _interactionbc:

InteractionBoundaryCondition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A special boundary condition prescribed on interaction PFEM particles (see
:ref:`interactionparticle`) in the PFEM part of the fluid-structure problem.
Particles of this sort are treated as if they had prescribed velocities, but
the values change dynamically as the solid part deforms; the velocities are
obtained from the coupled structural nodes.

.. record::

   :descitem:`InteractionBoundaryCondition` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}`

LatticeNeumannCoupling
~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Part of the lattice module; compile with ``USE_LM``.

An active boundary condition coupling a transport (pore-pressure) lattice to a
mechanical lattice in a staggered analysis.  For each mechanical node listed in
:param:`smnodes` it reads the fluid pressure :math:`P_f` from the corresponding
transport node in :param:`tmnodes`, in the coupled transport slave problem of a
``StaggeredProblem``, and applies a nodal force
:math:`f = P_f\, l\, \mathbf{n}`, where :math:`l` is the distance between the
mechanical node and its transport counterpart and :math:`\mathbf{n}` is the
unit :param:`direction`.

This implements Approach 1 of Grassl, Fahy, Gallipoli and Wheeler (2015).  The
transport slave problem is identified by the ``coupling`` field of the driving
``StaggeredProblem``.

.. record::

   :descitem:`LatticeNeumannCoupling` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`smnodes{ia}`
   :elemparam:`tmnodes{ia}` :elemparam:`direction{ra}`

**Parameters**

:param:`smnodes{ia}`, :param:`tmnodes{ia}`
    Mechanical and transport node numbers.  The two arrays must have the same
    length.

:param:`direction{ra}`
    Unit direction vector :math:`\mathbf{n}` of the applied force.

LatticeDirichletCoupling
~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Part of the lattice module; compile with ``USE_LM``.

Prescribes the pore pressure :math:`P_f` at a transport lattice node as the
distance-weighted average of the compression-only normal stress of the
mechanical lattice elements listed in :param:`couplingelements`, read from the
coupled mechanical slave problem of a ``StaggeredProblem``.  Tensile normal
stress is clamped to zero.

The coupling is explicit, lagged by one step: the pressure prescribed at step
:math:`n` uses the mechanical stress of step :math:`n-1`, and the first step
prescribes zero.  This implements Approach 2 of Grassl, Fahy, Gallipoli and
Wheeler (2015).  The mechanical slave problem is identified by the ``coupling``
field of the driving ``StaggeredProblem``.

The node the condition acts on is given through the associated set and
:param:`dofs` — DOF 11, :math:`P_f` — as for a standard ``BoundaryCondition``.

.. record::

   :descitem:`LatticeDirichletCoupling` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`couplingelements{ia}`
   :optelemparam:`set{in}` :optelemparam:`dofs{ia}`

.. _ContactBoundaryConditions:

Contact boundary conditions
---------------------------

Contact boundary conditions pair two contact surfaces (see
:ref:`ContactSurfaceRecords`), which must both be defined before the condition
that references them.

.. _structural-penalty-contact:

Structural penalty contact
~~~~~~~~~~~~~~~~~~~~~~~~~~

A penalty-based contact boundary condition for modelling contact between
deformable bodies.  The formulation enforces normal and tangential constraints
between the surfaces defined by ``StructuralFEContactSurface`` records and their
underlying ``StructuralContactElement_*`` elements.

It follows Konyukhov and Schweizerhof, *Computational Contact Mechanics:
Geometrically Exact Theory for Arbitrary Shaped Bodies* (Springer, LNACM 67,
2013): closest-point projection (Ch. 3), penalty/friction evolution and return
mapping (Sec. 6.1), and consistent tangents (Sec. 7.1).

The contact is enforced by the penalty method and contributes to the residual
and to the tangent system of equations at each iteration.

.. record::

   :descitem:`structuralpenaltycontactbc` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`dofs{ia}`
   :elemparam:`pn{rn}` :elemparam:`pt{rn}` :elemparam:`friction{rn}`
   :elemparam:`mastersurface{in}` :elemparam:`slavesurface{in}`
   :elemparam:`nsd{in}` :optelemparam:`frictiontransition{rn}`
   :optelemparam:`frictionhardening{rn}` :optelemparam:`algo{in}`
   :optelemparam:`searchpadding{rn}` :optelemparam:`searchtol{rn}`
   :optelemparam:`facethysteresis{rn}`
   :optelemparam:`generalizedfeatures{in}`
   :optelemparam:`directionalprojection{in}`
   :optelemparam:`autopenalty{in}` :optelemparam:`tangentmode{in}`
   :optelemparam:`fdcheck{in}` :optelemparam:`fdperturbation{rn}`
   :optelemparam:`fdoutputprefix{s}` :optelemparam:`fdtolerance{rn}`

**Required parameters**

:param:`pn{rn}`
    Normal penalty stiffness, controlling the resistance against penetration
    between the contacting surfaces.  Still required as a placeholder when
    :param:`autopenalty` is ``1``, in which case its value is ignored.

:param:`pt{rn}`
    Tangential penalty stiffness, controlling the tangential response.  The
    same placeholder rule as :param:`pn` applies when :param:`autopenalty` is
    ``1``.

:param:`friction{rn}`
    Coulomb coefficient of friction.  Set it to ``0.0`` for frictionless
    contact.

    .. warning::

       The frictional branch of this boundary condition is an experimental,
       unverified development feature.  Use ``friction 0`` for production
       analyses until it has been validated.

:param:`mastersurface{in}`, :param:`slavesurface{in}`
    Identifiers of the master and slave contact surfaces, defined by the
    corresponding ``StructuralFEContactSurface`` records; see
    :ref:`ContactSurfaceRecords`.

:param:`dofs{ia}`
    The affected degrees of freedom, typically ``2 1 2`` in 2D problems or
    ``3 1 2 3`` in 3D.

:param:`nsd{in}`
    Number of spatial dimensions: ``2`` for plane strain or plane stress, ``3``
    for 3D problems.

**Friction parameters**

:optparam:`frictiontransition{rn}`
    Smoothness of the differentiable stick/slip projection used by the
    experimental friction model.  Default ``0.0``, range ``[0, 1)``; zero
    recovers the sharp Coulomb return map.  Requires :param:`tangentmode`
    ``0``, ``2`` or ``3``.

:optparam:`frictionhardening{rn}`
    Dimensionless post-yield tangential-slip hardening ratio used by the
    experimental friction model.  Default ``0.0``, range ``[0, 1)``; zero is
    perfect Coulomb friction.  Positive values require a positive :param:`pt`
    and :param:`tangentmode` ``0``, ``2`` or ``3``.

**Contact search parameters**

:optparam:`algo{in}`
    Contact search algorithm.  Default ``0``, a plain surface-to-surface
    search; ``1`` selects a sweep-and-prune broad-phase search, 3D only
    (``nsd 3``), recommended for larger numbers of contact elements.

:optparam:`searchpadding{rn}`
    Absolute broad-phase bounding-box padding used by the contact search.  A
    negative value, or omitting the field, selects the automatic,
    geometry-based default.

:optparam:`searchtol{rn}`
    Parametric-domain margin used by the "is this point still inside this
    facet" test.  Default ``1.e-10``.  Comparable to other codes'
    sliding-elastic-interface search tolerance, whose typical default is
    ``0.01``.

:optparam:`facethysteresis{rn}`
    Relative distance-squared margin a competing master facet must exceed
    before it may replace the facet currently owned by a contact pair.  Default
    ``0.0``, which disables the hysteresis.  A small positive value, for
    example ``1.e-6``, suppresses Newton chattering caused by a slave point
    sitting near a shared edge between two adjacent facets, where the true
    closest facet would otherwise flip every iteration.

:optparam:`generalizedfeatures{in}`
    Default ``0``.  Set to ``1`` to extend the closest-point search from plain
    facet (surface) projection to also consider the edge and vertex features of
    the master surface.  Mutually exclusive with
    :param:`directionalprojection`.

:optparam:`directionalprojection{in}`
    Default ``0``.  Set to ``1`` to project each slave point onto the master
    surface along the slave surface's own normal direction, instead of a
    general closest-point search.  Currently requires ``nsd 3``.  Mutually
    exclusive with :param:`generalizedfeatures`.

**Penalty and tangent parameters**

:optparam:`autopenalty{in}`
    Default ``0``.  Set to ``1`` to ignore the :param:`pn` and :param:`pt`
    values and compute the normal and tangential penalty stiffnesses
    automatically for each slave contact element, from the contacting
    materials' initial Young's modulus and the element geometry, through an
    :math:`E_n A/V` factor.

:optparam:`tangentmode{in}`
    Contact tangent formulation.  Default ``0``, automatic selection; ``1``
    rate-form analytical, for diagnostic use; ``2`` branch-frozen
    finite-difference; ``3`` exact finite-step analytical tangent, including
    facet-history columns, with a finite-difference fallback for the projection
    features it does not yet support.

**Finite-difference verification**

:optparam:`fdcheck{in}`
    Default ``0``.  Set to ``1`` to enable an internal finite-difference
    verification of the analytical contact tangent on every call.  A
    development and debugging aid only; it adds significant runtime cost.

:optparam:`fdperturbation{rn}`
    Relative perturbation size used by :param:`fdcheck`.  Default ``1.e-7``.

:optparam:`fdoutputprefix{s}`
    Filename prefix for the diagnostic output written by :param:`fdcheck`.
    Default ``"contact_fd"``.

:optparam:`fdtolerance{rn}`
    Relative tolerance used by :param:`fdcheck` when comparing the analytical
    and finite-difference tangents.  Default ``0.0``.

**Example**

.. code-block:: none

   structuralpenaltycontactbc 3 loadTimeFunction 1 dofs 2 1 2 \
       pn 1.e8 pt 1.e8 friction 0.0 mastersurface 1 slavesurface 2 nsd 2

This defines a frictionless penalty contact condition between master surface 1
and slave surface 2, with normal and tangential stiffness equal to 1.e8.  The
condition acts on degrees of freedom 1 and 2 — displacement in the *x* and *y*
directions — in a 2D plane strain or plane stress domain.

For bidirectional, two-pass contact, two such conditions can be defined with
the master and slave surfaces swapped:

.. code-block:: none

   structuralpenaltycontactbc 6 loadTimeFunction 3 dofs 3 1 2 3 \
       pn 0.01 pt 0.01 friction 0 autopenalty 1 \
       mastersurface 2 slavesurface 1 nsd 3 algo 1 directionalprojection 1
   structuralpenaltycontactbc 7 loadTimeFunction 3 dofs 3 1 2 3 \
       pn 0.01 pt 0.01 friction 0 autopenalty 1 \
       mastersurface 1 slavesurface 2 nsd 3 algo 1 directionalprojection 1

This frictionless two-pass example uses automatic penalty stiffness and the
sweep-and-prune search with directional projection, each surface acting as
master for the other.

.. _thermal-s2s-contact:

Thermal surface-to-surface contact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A surface-to-surface thermal contact boundary condition modelling heat transfer
across the interface between two bodies.  It couples the temperature degrees of
freedom on the surfaces defined by the corresponding thermal contact surface
records, and contributes to the residual vector and the tangent matrix of the
thermal problem.

.. record::

   :descitem:`S2SThermalContact` :elemparam:`num{in}`
   :elemparam:`loadTimeFunction{in}` :elemparam:`dofs{ia}`
   :elemparam:`mastersurface{in}` :elemparam:`slavesurface{in}`
   :elemparam:`gapConductivityFunction{in}`

**Parameters**

:param:`loadTimeFunction{in}`
    Identifier of the time/load scaling function applied to the thermal contact
    contribution.

:param:`dofs{ia}`
    The affected degrees of freedom.  In a standard thermal analysis this is
    typically the single temperature DOF, number ``10``.

:param:`mastersurface{in}`, :param:`slavesurface{in}`
    Identifiers of the master and slave thermal contact surfaces.

:param:`gapConductivityFunction{in}`
    Identifier of a function giving the effective interface conductivity, or
    conductance, as a function of the gap.

**Example**

.. code-block:: none

   S2SThermalContact 2 loadTimeFunction 1 dofs 1 10 \
       mastersurface 1 slavesurface 2 gapConductivityFunction 2

This defines a thermal surface-to-surface contact condition between master
surface 1 and slave surface 2.  The condition acts on the temperature degree of
freedom, numbered ``10``, and uses ``gapConductivityFunction 2`` to prescribe a
gap-dependent thermal conductivity across the interface.

The gap conductivity function itself may be defined as, for example:

.. code-block:: none

   PiecewiseLinFunction 2 npoints 2 \
       t 2 0 0.1 \
       f(t) 2 10000 0

Here the variable ``t`` is the distance — the gap — between the master and slave
surfaces, and ``f(t)`` is the gap-dependent thermal conductivity.  A high
conductivity is prescribed for zero gap, that is perfect contact, decreasing
linearly to zero as the gap grows to ``0.1``.

.. _other-boundary-conditions:

Other registered boundary conditions
------------------------------------

OOFEM registers considerably more boundary condition types than are documented
above.  The following are available in input files but not yet described here;
their attributes can be read off the ``_IFT_`` definitions in the corresponding
header.

**Multiscale and homogenization**

``prescribedgradientbcneumann``, ``prescribedgradientperiodic``,
``prescribedgradientbcweakdirichlet``, ``prescribedgradientbcweakperiodic``,
``prescribedgradientmultiple``, ``weakperiodicbc``, ``prescribedmean``,
``solutionbasedshapefunction``, ``transversereinforcementconstraint``,
``prescribeddispslipbcdirichletrc``, ``prescribeddispslipbcneumannrc``,
``prescribeddispslipmultiple``

**Transport gradient conditions**

``tmgraddirichlet``, ``tmgradneumann``, ``tmgradperiodic``

**Loads**

``pointload``, ``momentload``, ``constantpressureload``,
``freeconstantsurfaceload``, ``gravitypressure``, ``surfacetension``,
``reinforcement``, ``depositedheatsource``

**Temperature fields**

``externaltempfieldload``, ``usrdeftempfield``, ``tf1``

**Other**

``rotatingboundary``, ``userdefdirichletbc``,
``PrescribedGenStrainShell7``

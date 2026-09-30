.. _NodeElementSideRecords:

Dof manager records
===================

These records describe the individual dof managers: nodes, or element sides if
they manage any DOFs.

.. record::

   :descitem:`DofManagerType` :elemparam:`num{in}`
   :optelemparam:`load{ia}` :optelemparam:`DofIDMask{ia}`
   :optelemparam:`bc{ia}` :optelemparam:`ic{ia}`
   :optelemparam:`doftype{ia}` :optelemparam:`masterMask{ia}`
   :optelemparam:`boundary{}`
   ``<``\ :optelemparam:`globnum{in}`\ ``>``
   ``<``\ :optelemparam:`shared{}` | :optelemparam:`remote{}` |
   :optelemparam:`null{}`\ ``>``
   ``<``\ :optelemparam:`partitions{ia}`\ ``>``

The order of the records is arbitrary and the dof manager number is given by
:param:`num`.  Numbering may be non-contiguous: think of the number as a label
assigned to an individual dof manager, by which it is referenced.

Degrees of freedom
------------------

By default the nodal DOFs are determined by asking all the connected elements.
Additional DOFs can be requested with the :param:`DofIDMask` array, which also
fixes their physical interpretation.

:optparam:`DofIDMask{ia}`
    Each item describes the physical meaning of the corresponding DOF in the
    dof manager, using the ``DofIDItem`` values in the table below.

.. important::

   Two DOFs with the same physical meaning may not appear in the same dof
   manager.

.. list-table:: ``DofIDItem`` values accepted by ``DofIDMask``
   :header-rows: 1
   :widths: 12 20 68

   * - Value
     - Name
     - Meaning
   * - 1, 2, 3
     - ``D_u``, ``D_v``, ``D_w``
     - Displacement along the *x*, *y*, *z* axis
   * - 4, 5, 6
     - ``R_u``, ``R_v``, ``R_w``
     - Rotation about the *x*, *y*, *z* axis (right-hand rule)
   * - 7, 8, 9
     - ``V_u``, ``V_v``, ``V_w``
     - Velocity along the *x*, *y*, *z* axis
   * - 10
     - ``T_f``
     - Temperature field
   * - 11
     - ``P_f``
     - Pressure field
   * - 12, 13
     - ``G_0``, ``G_1``
     - DOFs of gradient-type constitutive formulations
   * - 14
     - ``C_1``
     - Mass concentration of the first constituent
   * - 15, 16, 17
     - ``W_u``, ``W_v``, ``W_w``
     - Components of the change in the director field
   * - 18
     - ``Gamma``
     - Inhomogeneous thickness strain along the director field
   * - 19–22
     - ``D_u_edge_const`` …
     - Constant and linear parts of the boundary displacement, used by the
       Trefftz element
   * - 23, 24
     - ``Warp_PsiTheta``, ``Warp_Theta``
     - Relative twist times the deplanation function, and the relative twist,
       used by the ``TrWarp`` element
   * - 25, 26, 27
     - ``LMP_u``, ``LMP_v``, ``LMP_w``
     - Lagrange multipliers in the *x*, *y*, *z* direction
   * - 28, 29, 30
     - ``Trac_u``, ``Trac_v``, ``Trac_w``
     - Independent traction field in the *x*, *y*, *z* direction
   * - 31–39
     - ``E_xx`` … ``E_yx``
     - Macroscopic strain components
   * - 40, 41, 42
     - ``G_yz``, ``G_xz``, ``G_xy``
     - Macroscopic shear strain components
   * - 43–50
     - ``K_xx`` … ``K_yx``
     - Macroscopic curvature components
   * - 51, 52, 53
     - ``S_u``, ``S_v``, ``S_w``
     - Macroscopic reinforcement slip field
   * - 54
     - ``VF``
     - Volume fraction
   * - 55
     - ``P_f2``
     - Pressure field of the second phase
   * - 56
     - ``C_2``
     - Mass concentration of the second constituent

The authoritative list is the ``DofIDItem`` enumeration in
``src/core/dofiditem.h``.

.. note::

   XFEM enrichment DOFs are **not** in this range.  They are allocated
   dynamically by the enrichment items, starting above ``MaxDofID`` (500), so
   they are never written in a ``DofIDMask``.

Boundary and initial conditions
-------------------------------

Primary (Dirichlet) boundary conditions are applied through :param:`bc`,
natural boundary conditions through :param:`load`.

:optparam:`bc{ia}`
    Primary boundary conditions.  The size must equal the number of DOFs in the
    dof manager, and the *i*-th value relates to the *i*-th DOF; the ordering
    and physical meaning of the DOFs follow from the domain record, and may be
    set per dof manager with :param:`DofIDMask`.  The values are boundary
    condition record numbers, or zero where no primary condition applies to
    that DOF.  A compatible type is required: primary conditions need
    ``BoundaryCondition`` records.  When omitted, no primary condition is
    applied.

:optparam:`load{ia}`
    Record numbers of the natural boundary conditions that are applied; the
    required record type is ``NodalLoad``.  Where more than one natural
    condition applies, the actual value is the sum of all contributions.  Note
    that the *values* of a natural condition for the individual DOFs are given
    in its own record, not here.  When omitted, no natural condition is
    applied.

:optparam:`ic{ia}`
    Initial conditions, in the same layout as :param:`bc`.  The size should
    equal the number of DOFs in the dof manager, and the values are initial
    condition record numbers, or zero where no initial condition applies — in
    which case a zero initial value is assumed.

:optparam:`boundary{}`
    A flag marking the dof manager as lying on the boundary of the domain.

Master-slave DOFs
-----------------

:param:`doftype` and :param:`masterMask` connect some of a dof manager's DOFs —
the *slave* DOFs — to the corresponding DOF, by physical meaning, of another dof
manager — the *master* DOF.  The master-slave principle makes it easy to model,
for example, a structural hinge, where several elements are connected by
introducing several nodes with the same coordinates that share the same
displacement DOFs while each keeps its own rotational DOFs.

:optparam:`doftype{ia}`
    Type of each DOF to create: ``0`` a master DOF, ``1`` a simple slave DOF
    linked to a single master DOF, ``2`` a general slave DOF that may depend on
    different DOFs belonging to different dof managers.  When omitted, all DOFs
    are created as master DOFs.

:optparam:`masterMask{ia}`
    Required whenever :param:`doftype` is given.  Its meaning depends on the
    type of the dof manager and is described in the sections below.

.. _parallel-dofmanagers:

Parallel mode
-------------

.. note::

   This section applies only to the parallel build.  See
   :ref:`parallel-model` for the overall picture.

Exactly one of :param:`shared`, :param:`remote` and :param:`null` may be used
for a given dof manager.  If none is used, the dof manager is a local one for
the partition.

:optparam:`shared{}`
    The dof manager is shared by neighbouring partitions, and the contributions
    from all contributing domains are summed.  Typical of the node-cut
    algorithm; see :numref:`fig-nodecut-lm` and :numref:`fig-nodecut-nlm`.

:optparam:`remote{}`
    The dof manager in the active domain is only a mirror of a remote dof
    manager, so the remote values have to be copied into the local ones.
    Typical of the element cut; see :numref:`fig-elementcut-lm`.

:optparam:`null{}`
    A so-called null dof manager, which should be shared only by remote
    elements.  Remote elements exist only so that a nonlocal constitutive model
    can average efficiently, transferring just the local material values to be
    averaged.  Null nodes therefore serve only to compute the real integration
    point coordinates of remote elements; there is no reason to maintain their
    unknowns, so they are assigned no equation number and contribute nothing to
    the governing equations of the local partition.  See
    :numref:`fig-nodecut-nlm`.

:optparam:`partitions{ia}`
    For a shared dof manager, the list of remote partitions sharing it; for a
    remote one, the partition holding its remote counterpart.  The local
    partition must not be included in the list.

:optparam:`globnum{in}`
    Global number of the dof manager, unique across all partitions.  The
    parallel version requires a unique global numbering of dof managers to link
    the partitions together.

Slave DOFs are allowed, but their masters have to be in the same partition.
The masters may themselves be remote copies.

.. _parallel-partitioning-figures:

Partitioning strategies
~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../figs/nodecut0cb.pdf
   :alt: Node-cut partitioning.
   :name: fig-nodecut

   Node-cut partitioning.

.. figure:: ../figs/nodecut1cb.pdf
   :alt: Node-cut partitioning, local constitutive mode.
   :name: fig-nodecut-lm

   Node-cut partitioning, local constitutive mode.

.. figure:: ../figs/nodecutnonloc1.pdf
   :alt: Node-cut partitioning, nonlocal constitutive mode.
   :name: fig-nodecut-nlm

   Node-cut partitioning, nonlocal constitutive mode.

.. figure:: ../figs/elementcut0.pdf
   :alt: Element-cut partitioning.
   :name: fig-elementcut

   Element-cut partitioning.

.. figure:: ../figs/elementcut1.pdf
   :alt: Element-cut partitioning, local constitutive mode.
   :name: fig-elementcut-lm

   Element-cut partitioning, local constitutive mode.

Supported dof manager types
---------------------------

Node
~~~~

An abstraction for a finite element node.

.. record::

   :descitem:`Node` :elemparam:`num{in}` :elemparam:`coords{ra}`
   :optelemparam:`lcs{ra}`

**Parameters**

:param:`coords{ra}`
    Coordinates of the node in the global coordinate system: the *x*, *y* and,
    depending on the problem, *z* coordinate.

:optparam:`lcs{ra}`
    A user-defined local coordinate system in the node.  The array holds six
    numbers: the first three are a direction vector of the local *x*-axis, the
    next three a direction vector of the local *y*-axis, and the local *z*-axis
    follows from the vector product.  A right-handed system is assumed.  When a
    local system is given, the boundary conditions and the applied loading are
    specified in it, and the reactions and displacements are also reported in
    it.  By default the coordinate system in a node is the global one.

**Master-slave support**

A node can create only master and simple slave DOFs, so the permitted values of
:param:`doftype` are ``0`` and ``1``.  For a node, :param:`masterMask` is an
array whose size equals the number of DOFs, the *i*-th value giving the master
dof manager to which the *i*-th DOF is directly linked; DOFs with the same
physical meaning are linked together.  A local coordinate system in a node with
linked DOFs is supported, but it must be exactly the same as on the master.

Rigid arm node
~~~~~~~~~~~~~~

A node connected to another node — the master — by a rigid arm.  Its DOFs may be
linked to the master through the rigid arm transformation, or be independent.

The rigid arm node avoids the very stiff elements otherwise used to model a
rigid-arm connection: it maps its DOFs onto the master DOFs by simple
transformations, assuming small rotations, so contributions to the rigid arm
node can be localised directly into the master's equations.  A rigid arm node
can be loaded independently of its master, but cannot have boundary or initial
conditions of its own — these follow entirely from the master's conditions.

.. record::

   :descitem:`RigidArmNode` :elemparam:`num{in}` :elemparam:`coords{ra}`
   :elemparam:`master{in}` :optelemparam:`masterMask{ia}`
   :optelemparam:`lcs{ra}`

**Parameters**

:param:`coords{ra}`
    Coordinates of the node in the global coordinate system.

:param:`master{in}`
    Number of the master node onto which the rigid arm node's DOFs are mapped.

:optparam:`masterMask{ia}`
    How each mapped DOF depends on the master DOFs.  Its size should equal the
    number of DOFs, and for every linked DOF — one whose :param:`doftype` value
    is ``2`` — the corresponding value should be ``1``.

:optparam:`lcs{ra}`
    Local coordinate system, as for ``Node``.  The rigid arm node and its
    master may have different local coordinate systems; if none is given, the
    global system applies.

Currently only certain DOFs can be mapped: see :param:`doftype`.  Linked DOFs
should have a :param:`doftype` value of ``2``, non-linked (primary) DOFs ``0``.

Hanging node
~~~~~~~~~~~~

A node connected to a master element by generalized interpolation.  A hanging
node has no degrees of freedom of its own, apart from unlinked ones: all values
are interpolated from the master element and its DOFs.  This is useful, for
example, for an arbitrary FE mesh of a concrete specimen, or to allow local
refinement of an FE mesh.  Hanging nodes may form a chain.

Contributions of a hanging node are localised directly into the master's
equations.  A hanging node may have boundary or initial conditions of its own,
but only for primary, unlinked DOFs; for linked DOFs these conditions follow
entirely from the master's.  The local coordinate system must be the same for
all master nodes.  A hanging node can be loaded independently of its master.

.. record::

   :descitem:`HangingNode` :elemparam:`num{in}` :elemparam:`coords{ra}`
   :elemparam:`doftype{ia}` :optelemparam:`masterElement{in}`
   :optelemparam:`masterRegion{in}`

**Parameters**

:param:`coords{ra}`
    Coordinates of the node.

:param:`doftype{ia}`
    Permitted values are ``0`` for a primary DOF and ``2`` for a linked DOF.

:optparam:`masterElement{in}`
    Number of the element to which the hanging node is attached.  The node may
    sit at any coordinate inside the master element, which must support the
    necessary interpolation classes; the same interpolation is assumed for
    unknowns and for geometry.  When omitted, or set to ``-1``, the node
    locates the element closest to its coordinates.

:optparam:`masterRegion{in}`
    When omitted or zero, all regions are searched; otherwise only the elements
    in the cross section with this number.  Unused when
    :param:`masterElement` is given directly.

Continuum frame node
~~~~~~~~~~~~~~~~~~~~

A frame or beam node embedded in a continuum (solid) mesh.  Its input is
identical to ``HangingNode`` and its translational DOFs are interpolated from
the master element in exactly the same way.

Its linked rotational DOFs — those with :param:`doftype` ``2``, that is
``R_u``, ``R_v``, ``R_w`` — are instead constrained to the infinitesimal
rotation of the master element,
:math:`\boldsymbol{\omega} = \tfrac12\,\nabla\times\boldsymbol{u}`, evaluated
from the shape-function gradients of the translations.  This lets a frame node
embedded in a solid mesh, which carries no rotational DOFs, inherit the local
continuum rotation.  All three rotations, including the torsional one about the
beam axis, are determined automatically: no torsional restraint has to be
applied, and the result does not depend on the node ordering of the frame
element.

Rotational DOFs left as primary (:param:`doftype` ``0``) or fixed are not
touched.  To bond only the translations to the matrix and let the frame element
itself carry the rotations, use a plain ``HangingNode`` instead.

.. note::

   The continuum rotational constraint is currently implemented for the linear
   tetrahedron only; a non-tetrahedral master element raises an error.

.. record::

   :descitem:`ContinuumFrameNode` :elemparam:`num{in}`
   :elemparam:`coords{ra}` :elemparam:`doftype{ia}`
   :optelemparam:`masterElement{in}` :optelemparam:`masterRegion{in}`

Slave node
~~~~~~~~~~

Works exactly like a hanging node, except that the weights are not computed
from an element but given explicitly, along with the connected dof managers.

.. record::

   :descitem:`SlaveNode` :elemparam:`num{in}` :elemparam:`coords{ra}`
   :elemparam:`doftype{ia}` :elemparam:`masterDofMan{ia}`
   :elemparam:`weights{ra}`

**Parameters**

:param:`masterDofMan{ia}`
    Numbers of the connected master dof managers.

:param:`weights{ra}`
    Interpolation weights, in the same order as :param:`masterDofMan`.

General slave node
~~~~~~~~~~~~~~~~~~

A generalization of the slave node, in which each slave DOF may depend on a
different number of master DOFs.

.. record::

   :descitem:`GeneralSlaveNode` :elemparam:`num{in}`
   :elemparam:`coords{ra}` :elemparam:`doftype{ia}`
   :elemparam:`masterSizes{ia}` :elemparam:`masterList{ia}`
   :elemparam:`masterWeights{ra}`

**Parameters**

:param:`doftype{ia}`
    As for ``SlaveNode`` and ``HangingNode``: ``0`` a primary DOF, ``2`` a
    linked DOF.

:param:`masterSizes{ia}`
    Number of master DOFs for each slave DOF.

:param:`masterList{ia}`
    Master nodes and their DOFs, for all the slave DOFs.

:param:`masterWeights{ra}`
    Weights associated with the entries of :param:`masterList`.

Element side
~~~~~~~~~~~~

An abstraction for an element side that holds some unknowns.

.. record::

   :descitem:`ElementSide` :elemparam:`num{in}`

.. _pfemparticles:

PFEM particle
~~~~~~~~~~~~~

The particle used in PFEM analysis; see :ref:`pfemIncomp`.

.. record::

   :descitem:`PFEMParticle` :elemparam:`num{in}` :elemparam:`coords{ra}`

.. _interactionparticle:

Interaction PFEM particle
~~~~~~~~~~~~~~~~~~~~~~~~~

A special particle used in the PFEM part of the fluid-structure problem (see
:ref:`fluidstructureproblem`).  The particle is attached to
:param:`coupledNode` from the structural counterpart.

.. record::

   :descitem:`InteractionPFEMParticle` :elemparam:`num{in}`
   :elemparam:`coords{ra}` :elemparam:`bc{ia}`
   :elemparam:`coupledNode{in}`

**Parameters**

:param:`coupledNode{in}`
    Number of the structural node to which this particle is attached.

:param:`bc{ia}`
    Must prescribe an ``InteractionBoundaryCondition`` (see
    :ref:`interactionbc`) so that the velocities can be read from the solid
    nodes.

Other registered dof manager types
----------------------------------

The following dof manager types are registered but not documented here:

``particle``
    A particle used by the discrete element / particle models;
    ``src/sm/particle.C``.

``qcnode``
    A node of the quasicontinuum formulation, used with the
    ``qclinearstatic`` engineering model; ``src/core/qcnode.C``.

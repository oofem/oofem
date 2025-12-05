	
.. _DomainRecord:

Domain record(s)
================

This set of records describes the whole domain and its type. Depending
on the type of problem, there may be one or several domain records. If
not indicated, one domain record is default for all problem types.

| The domain type is used to resolve the default number of DOFs in node
  and their physical meaning. Format is following
| ``domain`` ``domainType``
| The ``domainType`` can be one from the following

-  The ``2dPlaneStress`` and ``2d-Truss`` modes declare two default dofs
   per node (u-displacement, v-displacement),

-  The ``3d`` mode declares three default dofs per node (u-displacement,
   v-displacement, w-displacement),

-  The ``2dMindlinPlate`` mode declares three default dofs per node
   (w-displacent, u-rotation, v-rotation). Strain vector contains
   :math:`\kappa_{xx}`, :math:`\kappa_{yy}`, :math:`\kappa_{xy}`,
   :math:`\gamma_{xz}`, :math:`\gamma_{yz}`. Stress vector contains
   :math:`m_{xx}`, :math:`m_{yy}`, :math:`m_{xy}`, :math:`q_{xz}`,
   :math:`q_{yz}`.

-  The ``3dShell`` mode declares six default dofs per node (displacement
   and rotation along each axis).

-  The ``2dBeam`` mode declares three default dofs per node
   (u-displacement, w-displacement, v-rotation).

-  The ``2dIncompFlow`` mode declares three default dofs per node
   (u-velocity, v-velocity, and pressure). The default number of dofs
   per node as well as their physical meaning can be overloaded in
   particular dof manager record (see section
   :ref:`NodeElementSideRecords`).

   The further records describe particular domain components -
   OutputManagers, DofManagers, Elements, CrossSection models, Material
   Models, Boundary and Initial Conditions and Load time functions.

.. _OutputManagerRecord:

Output manager record
---------------------

| The output manager controls output. It can filter output to specific
  solution steps, and within these selected steps allows also to filter
  output only to specific dof managers and elements. The format of
  output manager record is
| [``tstep_all``] [``tstep_step #(in)``] [``tsteps_out #(rl)``]
  [``dofman_all``] [``dofman_output #(rl)``]
  [``dofman_except #(rl)``] [``element_all``]
  [``element_output #(rl)``] [``element_except #(rl)``]

| To select all solution steps, in which output will be performed, use
  ``tstep_all``. To select each ``tstep_step``-nth step, use
  ``tstep_step`` parameter. In order to select only specific solution
  steps, the ``tsteps_out``\ list can be specified, supplying solution
  step number list in which output will be done. The combination of
  ``tstep_step`` and ``tsteps_out`` parameters is allowed.

| Output manager allows also to filter output to only specific dof
  managers and elements. If these specific members are selected, the
  output happens only in selected solution steps. The ``dofman_all`` and
  ``element_all`` parameters select all dof managers or elements
  respectively. Parameter arrays ``dofman_output`` and
  ``element_output`` allow to select only specific members. Numbers of
  selected members are then contained in ``dofman_output`` or
  ``element_output`` lists respectively. The previously selected members
  can be explicitly de-selected by specifying their component numbers in
  ``dofman_except`` or ``element_except`` lists. A few examples:
| ``dofman_output {1 3}``  prints nodes 1,3
| ``dofman_output {(1 3)}``  prints nodes 1,2,3
| ``element_output {1 3}``  prints elements 1,3
| ``element_output {(1 3)}``  prints elements 1,2,3
| ``element_output {(1 3) 5 6}``  prints elements 1,2,3,5,6

.. _ComponentsSizeRecord:

Components size record
----------------------

| This record describes the number of components in related domain. The
  particular records will follow immediately in input file. The general
  format is:

| ``ndofman #(in)`` ``nelem #(in)``
  ``ncrosssect #(in)`` ``nmat #(in)`` ``nbc #(in)``
  ``nic #(in)`` ``nltf #(in)`` [``nset #(in)``] [``ncontactsurf #(in)``] [``nbarrier #(in)``] 
  
| where ``ndofman`` represents number of dof managers (e.g. nodes) and their
  associated records, ``nelem`` represents number of elements and their
  associated records, ``ncrosssect`` is number of cross sections and
  their records, ``nmatdnMat`` is number of material models and their
  records, ``nbc`` represents number of boundary conditions (including
  loads and contact conditions) and their records, ``nic`` parameter determines the number of
  initial conditions, and ``nltf`` represents number of time functions
  and their associated records. The optional ``nset`` parameter specifies the
  number of set records (see :ref:`SetRecords`), and the optional ``ncontactsurf``
  parameter specifies the number of contact surface records (see
  :ref:`ContactSurfaceRecords`).The optional parameter ``nbarrier``
  represents the number of nonlocal barriers and their records. If not
  specified, no optional entities (sets, contact surfaces, or barriers)
  are assumed.
  
.. _NodeElementSideRecords:

Dof manager records
-------------------

These records describe individual DofManager records (i.e. nodes or
element sides (if they manage some DOFs)). The general format is
following:

``DofManagerType`` ``#(in)`` [``load #(ra)``] [``DofIDMask #(ia)``]
[``bc #(ia)``] [``ic #(ia)``] [``doftype #(ia)``
``masterMask #(ia)``]  <[``shared``]>
:math:`|` <[``remote``]> :math:`|` <[``null``]>
<[``partitions #(ia)``]>

The order
of particular records is optional, the dof manager number is determined
by (``#(in)`` parameter. The numbering of individual dof managers
is arbitrary, it could be even non-continuous. In this context, one
could think of dof manager number as a label that is assigned to
individual dof manager and by which the dof manager is referenced.

By default, the nodal DOFs are determined by asking all the connected
elements. Specifying additional dofs can be done using the using the
``DofIDMask`` array which determines their physical interpretation. Each
item of ``DofIDMask`` array describes the physical meaning of
corresponding DOF in dof manager. Currently the following values are
supported: {u-displacement=1, v-displacement=2, w-displacement=3,
u-rotation=4, v-rotation=5, w-rotation=6, u-velocity=7, v-velocity=8,
w-velocity=9, temperature=10, pressure=11, special dofs for
gradient-type constitutive models=12 and 13, mass concentration=14,
special dofs for extended finite elements (XFEM)=15–30}. **It is not
allowed to have two DOFs with the same physical meaning in the same
DofManager.**

The applied primary (Dirichlet) boundary conditions are specified using
"bc" record, while natural boundary conditions using "load" parameter.

-  The size of "bc" array (primary bc) should be equal to number of DOFs
   in dof manager and i-th value relates to i-th DOF - the ordering and
   physical meaning of DOFs is determined by domain record and can be
   optionally specified for each dof manager individually (see next
   paragraph). The values of this array are corresponding boundary
   condition record numbers or zero, if no primary bc is applied to
   corresponding DOF. The compatible boundary condition type are
   required: primary conditions require "BoundaryCondition" records.

-  The load "array" contains record numbers of natural boundary
   conditions that are applied. The required record type for natural
   condition is "NodalLoad". The actual value is the summation of all
   contributions, if more than one natural bc is applied. See section on
   boundary conditions for the syntax. Please note, that the values of
   natural bc for individual DOFs are specified in its record, not in
   dofmanager record.

By default, if "bc" and/or "load" parameters are omitted, no primary
and/or natural bc are applied. Analogously, initial conditions are
represented using ``ic`` array. The size of ``ic`` array should be equal
to number of DOFs in dof manager. The values of this array are
corresponding initial condition record numbers or zero, if no initial
condition is applied to corresponding DOF (in this case zero value is
assumed as value of initial condition).

Parameters ``dofType`` and ``masterMask`` allows to connect some dof
manager’s dofs (so-called “slave” dofs) to corresponding dof (according
to their physical meaning) of another dof manager (so-called “master”
dof). The master slave principle allows for example simple modeling of
structure hinges, where multiple elements are connected by introducing
multiple nodes (with same coordinates) sharing the same displacement
dofs and each one possessing their own rotational dofs. Parameter
``dofType`` determines the type of (slave) dof to create. Currently
supported values are 0 for master DOF, 1 for simpleSlave DOF (linked to
another single master DOF), and 2 for general slave dof, that can depend
on different DOFs belonging to different dof managers. If ``dofType`` is
not specified, then by default all DOFs are created as master DOFs. If
provided, masterMask is also required. The meaning of ``masterMask``
parameter is depending on type of particular dofManager, and will be
described in corresponding sections.

Supported DofManagerType keywords are

-  Node record

   | ``Node`` ``coords #(ra)`` [``lcs #(ra)``]
   | Represent an abstraction for finite element node. The node
     coordinates in space (given by global coordinate system) are
     described using ``coords`` attribute. This array contains x, y and
     possibly z (depends on problem under consideration) coordinate of
     node. By default, the coordinate system in node is global
     coordinate system. User defined local coordinate system in node is
     described using ``lcs`` array. This array contains six numbers,
     where the first three numbers represent a directional vector of the
     local x-axis, and the next three numbers represent a directional
     vector of the local y-axis. The local z-axis is determined using a
     vector product. A right-hand coordinate system is assumed. If user
     defined local coordinate system in node is specified, then the
     boundary conditions and applied loading are specified in this local
     coordinate system. The reactions and displacements are also in
     ``lcs`` system at the output.

   The node can create only master DOFs and SimpleSlave DOFs, so the
   allowable values of ``dofType`` array are in range 0,1. For the Node
   dof manager, the ``masterMask`` is the array of size equal to number
   of DOFs, and the i-th value determines the master dof manager, to
   which i-th dof is directly linked (the dof with same physical meaning
   are linked together). The local coordinate system in node with same
   linked dofs is supported, but it should be exactly the same as on
   master.

-  Rigid arm record

   | ``RigidArmNode`` ``coords #(ra)`` ``master #(in)``
     [``masterMask #(ia)``] [``lcs #(ra)``]
   | Represent node connected to other node (called master) using rigid
     arm. Rigid arm node DOFs can be linked to master (via rigid arm
     transformation) or can be independent. The rigid arm node allows to
     avoid very stiff elements used for modelling the rigid-arm
     connection. The rigid arm node maps its dofs to master dofs using
     simple transformations (small rotations are assumed). Therefore,
     the contribution to rigid arm node can be localized directly to
     master related equations. The rigid arm node can not have its own
     boundary or initial conditions, they are determined completely from
     master dof conditions. Currently it is possible to map only certain
     dofs - see ``dofType``. Linked DOFs should have dofType value equal
     to 2, non-linked (primary) DOFs 0.

   Rigid arm node can be loaded independently of master. The node
   coordinates in space (given by global coordinate system) are
   described using ``coords`` field. This array contains x, y and
   possibly z (depends on problem under consideration) coordinate of
   node. The ``master`` parameter is the master node number, to which
   rigid arm node dofs are mapped. The rigid arm node and master can
   have arbitrary local coordinate systems (if not specified, global
   one is assumed).

   The optional parameter ``masterMask`` allows to specify how
   particular mapped DOF depends on master DOFs. The size of
   ``masterMask`` array should be equal to number of DOFs. For all
   linked DOFs (with corresponding dofType value equal to 2) the
   corresponding value of ``masterMask`` array should be 1.

   The local coordinate system in rigid arm node is supported, the
   coordinate system in master and slave can be different. If no lcs is
   set, global one is assumed.the global cs applies.

-  Hanging node

   ``HangingNode`` ``coords #(ra)`` ``dofType #(in)``
   [``masterElement #(in)``] [``masterRegion #(in)``]

   Hanging node is connected to an a master element using generalized
   interpolation. Hanging node posses no degrees of freedom (except
   unlined dofs) - all values are interpolated from corresponding master
   elements and its DOFs. arbitrary FE mesh of concrete specimen or to
   facilitate the local refinement of FE mesh. The hanging nodes can be
   in a chain.

   The contributions of hanging node are localized directly to master
   related equations. The hanging node can have its own boundary or
   initial conditions, but only for primary unlinked DOFs. For linked
   DOFs, these conditions are determined completely from master DOF
   conditions. The local coordinate system should be same for all master
   nodes. The hanging node can be loaded independently of its master.

   Values of array ``dofType`` can have following values: 0-primary DOF,
   2-linked DOF.

   The value of ``masterElement`` specifies the element number to which
   the hanging node is attached. The node can be attached to any
   arbitrary coordinate within the master element. The element must
   support the necessary interpolation classes. The same interpolation
   for unknowns and geometry is assumed.

   The no (or -1) value for ``masterElement`` is supplied, then the node
   will locate the element closest to its coordinate. If no (or zero)
   value for ``masterRegion`` is supplied, then all regions will be
   searched, otherwise only the elements in cross section with number
   ``masterRegion``. If ``masterElement`` is directly supplied
   ``masterRegion`` is unused.

-  Slave node

   ``SlaveNode``  ``coords #(ra)``  ``dofType #(in)``
   ``masterDofMan #(ia)``  ``weights #(ra)``

   Works identical to hanging node, but the weights (``weights``) are
   not computed from any element, but given explicitly, as well as the
   connected dof managers (``masterDMan``).

-  General Slave node

   ``GeneralSlaveNode``  ``coords #(ra)``  ``dofType #(ia)``
   ``masterSizes #(ia)``  ``masterList #(ia)``  ``masterWeights #(ra)``
   
   Generalization of SlaveNode; Meaning of ``dofType`` is the same as for SlaveNode or HangingNode, i.e., 0 means primary DOF, while
   2 means linked DOF.  Array ``masteSizes`` specifies number of master dofs for each slave dof. Array ``masterList`` specifies master nodes and their dofs for all the slave dofs. Finally, ``masterList`` provides the associated weights. 

   
-  Element side

   | ``ElementSide``
   | Represents an abstraction for element side, which holds some
     unknowns.

.. _pfemparticles:  

-  PFEMParticle

   | ``PFEMParticle`` ``coords #(ra)``
   | Represent the particle used in PFEM analysis.

.. _interactionparticle:
   
-  InteractionPFEMParticle 

   | ``InteractionPFEMParticle`` ``coords #(ra)`` ``bc #(ia)``
     ``coupledNode #(in)``
   | Represent a special particle used in the PFEM-part of the
     FluidStructureProblem. The particle is attached to ``coupledNode``
     from the structural counter part.
     InteractionBoundaryCondition (see  :ref:`interactionbc`)
     must be prescribed under ``bc`` to access the velocities from solid
     nodes.

.. _ElementsRecords:

Element records
---------------

These records specify a description of particular elements. The general
format is following:

``ElementType`` ``#(in)`` ``mat #(in)`` ``crossSect #(in)``
``nodes #(ia)`` [``bodyLoads #(ia)``] [``boundaryLoads #(ia)``]
[``activityltf #(in)``] [``lcs #(ra)``]
<[``partitions #(ia)``]> <[``remote``]>

The order of element records is optional, the element number is
determined by ``#(in)`` parameter. The numbering of individual
elements is arbitrary, it could be even non-continuous. In this context,
one could think of element number as a label that is assigned to
individual elements and by which the element is referenced.

Element material is described by parameter ``mat``, which contains
corresponding material record number. Element cross section is
determined by cross section with ``crossSect`` record number. Element
dof managers (nodes, sides, etc.) defining element geometry are
specified using ``nodes`` array.

Body load acting on element is specified using ``bodyLoads`` array.
Components of this array are corresponding load record numbers. The
loads should have the proper type (body load type), otherwise error will
be generated.

Boundary load acting on element boundary is specified using
``boundaryLoads`` array. The format of this array is

.. math:: 2\cdot size \; lnum(1)~id(1)~\dots~lnum(size)~id(size),

where :math:`size` is total number of loadings applied to element,
:math:`lnum(i)` is the applied load number, and :math:`id(i)` is the
corresponding entity number, to which the load is applied (for example a
side or a surface number). The entity numbering is element dependent and
is described in element specific sections. The applied loads must be of
proper type (boundary load type), otherwise error is generated.

The support for element insertion and removal during the analysis is
provided. One can specify optional time function (identified by its id
using ``activityltf`` parameter). The nonzero value of this time
function indicates, whether the element is active (nonzero value, the
default) or inactive (zero value) at particulat solution step. Tested
for structural and transport elements. This feature allows considering
temperature evolution of layered casting of concrete, where certain
layers needs to be inactive before they are cast. See a corresponding
example in oofem tests how to enforce hydrating material model, boundary
conditions and element activity acting concurrently.

Orientation of local coordinates can be specified using ``lcs`` array.
This array contains six numbers, where the first three numbers represent
a directional vector of local x-axis, and the next three numbers
represent a directional vector of local y-axis. The local z-axis is
determined using the vector product. The ``lcs`` array on the element is
particularly useful for modeling of orthotropic materials which follow
the element orientation. On a beam or truss element, the ``lcs`` array
has no effect and the 1D element orientation is aligned with the global
:math:`xx` component.

Available material models, their outline and corresponding parameters
are described in separate **Element Library Manual.**

.. _SetRecords:

Set records
-----------

Sets specify regions of the geometry as a combination of volumes,
surfaces, edges, and nodes. The main usage of sets are to connect
regions of elements to a given cross section or apply a boundary
condition, though sets can be used for many other things as well.

``Set`` ``#(in)`` [``elements #(ia)``] [``elementranges #(rl)``]
[``allElements``] [``nodes #(ia)``] [``noderanges #(rl)``] [``allNodes``]
[``elementboundaries #(ia)``] [``elementedges #(ia)``] 
<``ver 1.6`` [``dofmanprops #(s)``] [``elemprops #(s)``]>

Volumes (elements) and nodes can be specified using either a list,
``elements``, ``nodes``, or with a range list ``elementranges``,
``noderanges``. Edges ``elementedges``, and surfaces
``elementboundaries``, are specified in a interleaved list, every other
number specifying the element, and edge/surface number (the total length
of the list being twice the number of surfaces/edges). The internal
numbering of edges/surfaces is available in the **Element Library
Manual**.

Note that edge loads (singular loads given in “newton per length” (or
equivalent), should be applied to ``elementedges``, surface loads
“newton per area” on ``elementboundaries``, and bulk loads “newton per
volume” on ``elements``.

Example 1: A deadweight (gravity) load would be applied to the
``elements`` in a set, while a distributed line load would be applied to
the midline “edge” of the beam element, thus should be applied to a
``elementedges`` set. In the latter case, the midline of the beam is
defined as the first (and only) “edge” of the beam.

Example 2: Axisymmetric structural element analysis: A deadweight load
would be applied to ``elements`` in a set. A external pressure would be
defined as a surface load an be applied to the ``elementboundaries`` in
a set. The element integrates the load (analytically) around the axis,
so the load would still count as a surface load.

Arbitrary dof manager and element properties can be specified using
``dofmanprops`` and ``elemprops`` parameters (supported by ver. 1.6 and higher). The optional ``dofmanprops``
parameter is a string, which contains any dof manager parameters that are subsequently set for all 
dof managers in the set. The optional ``elemprops`` parameter is a string, which contains any element
parameters that are subsequently set for all elements in the set. The syntax for setting the individual parameters is the syntax of input file records.
Note that parameters defined using this mechanism have the lowest priority, i.e. 
they can be overridden by the parameters defined in the dof manager or element records (see tests/sm/setprops01.in for an example).




.. _ContactSurfaceRecords:

Contact surface records
-----------------------

Contact surfaces define geometric entities that participate in contact
interactions. Each surface is formed by a set of contact elements and can
be assigned the role of *master* or *slave* within a contact boundary
condition. Contact surfaces provide the link between the geometric
description of the contacting interfaces and the boundary conditions that
enforce their interaction.

**Syntax**

``FEContactSurfaceType`` ``#(in)`` ``ce_set #(in)``

**Parameters**

- ``#(in)`` — unique surface identifier.
- ``ce_set`` — identifier of the element set that contains the contact
  elements forming this surface (see :ref:`SetRecords`).

**Description**

Each contact surface groups one or more contact elements of type
``ContactElementType_*`` into a single logical entity. These surfaces
are referenced by contact boundary conditions (see
:ref:`ContactBoundaryConditions`) as ``mastersurface`` or ``slavesurface``.
Both surfaces must exist before defining any contact boundary condition.

**Example**

::

   StructuralFEContactSurface 1 ce_set 3
   StructuralFEContactSurface 2 ce_set 4

In this example, structural finite element surface 1 is created from contact elements included in set
3, and surface 2 from structural contact elements in set 4. These surfaces can then be
paired using a ``structuralpenaltycontactbc`` record to define a contact
interaction.

**Notes**

- Each contact surface can serve as a master or slave surface depending on
  the definition in the contact boundary condition.
- The number of contact surface records is given by ``ncontactsurf`` in the
  Components size record.
- Contact surfaces do not introduce new DOFs; they provide only geometric
  organization for contact evaluation.









.. _CrossSectionRecords:

Cross section records
---------------------

These records specify a cross section model descriptions. The general
format is following:

``CrossSectType`` ``#(in)``

The order of particular cross section records is optional, cross section
model number is determined by ``#(in)`` parameter. The numbering
should start from one and should end at n, where n is the number of
records.

The crossSectType keyword can be one from following possibilities

-  | Integral cross section with constant properties
   | ``SimpleCS`` [``thick #(rn)``] [``width #(rn)``] [``area #(rn)``]
     [``iy #(rn)``] [``iz #(rn)``] [``ik #(rn)``]
     [``shearareay #(rn)``] [``shearareaz #(rn)``]
     ``beamshearcoeff #(rn)``
   | Represents integral type of cross section model. In current
     implementation, such cross section is described using cross section
     thick (``thickVal``) and width (``widthVal``). For some problems
     (for example 3d), the corresponding volume and cross section
     dimensions are determined using element geometry, and then you can
     omit some (or all) parameters (refer to documentation of individual
     elements for required cross section properties). Parameter ``area``
     allows to set cross section area, parameters ``iz``, ``iz``, and
     ``ik`` represent inertia moment along y and z axis and Saint-Venant torsional
     constant. Parameter ``beamshearcoeff`` allows to set shear
     correction factor, or equivalent shear areas (``shearareay`` and
     ``shearareaz`` parameters) can be provided. These cross section
     properties are assumed to be defined in local coordinate system of
     element.

-  | Integral cross section with variable properties
   | ``VariableCS`` [``thick #(expr)``] [``width #(expr)``] [``area #(expr)``]
     [``iy #(expr)``] [``iz #(expr)``] [``ik #(expr)``]
     [``shearareay #(expr)``] [``shearareaz #(expr)``]
   | Represents integral type of cross section model, where individual
     cross section parameters can be expressed as an arbitrary function
     of global coordinates x,y,z. Similar to SimpleCS, for some problems
     (for example 3d), the corresponding volume and cross section
     dimensions are determined using element geometry, then you can omit
     many (or some) parameters (refer to documentation of individual
     elements for required cross section properties). Parameter ``area``
     allows to set cross section area, parameters ``iz``, ``iz``, and
     ``ik`` represent inertia moment along y and z axis and Saint-Venant torsional
     constant. Parameters (``shearareay`` and ``shearareaz``
     determine shear area, which is required by beam and plate elements.
     All cross section properties are assumed to be defined in local
     coordinate system of element.

-  | Layered cross section
   | ``LayeredCS`` ``nLayers #(in)`` ``LayerMaterials #(ia)``
     ``Thicks #(ra)`` ``Widths #(ra)``  [``midSurf #(rn)``] [``nintegrationpoints #(in)``] [``layerintegrationpoints #(ia)``] [``beamshearcoeffxz #(rn)``]
   | Represents the layered cross section model, based on geometrical
     hypothesis, that cross sections remain planar after deformation.
     Number of layers is determined by ``nLayers`` parameter. Materials
     for each layer are specified by ``LayerMaterials`` array. For each
     layer is necessary to input geometrical characteristic, thick -
     using ``Thicks`` array, and width - using ``Widths`` array.
     Position of mid surface is determined by its distance from bottom
     of cross section using ``midSurf`` parameter (normal and momentum
     forces are then computed with regard to it’s position, by default it is located at average thickness position). 
     The number of integration points per layer can be specified using ``nintegrationpoints`` parameter, default is one integration point. It is also possible to set up different number of integration points per individual layer using ``layerintegrationpoints`` array, where its size should be equal to number of layers configured. The ``layerintegrationspoints`` parameter overrides the ``nitengrationpoints`` setting. The Gauss integration rule is used for setting up integration points in each layer.
   | The optional parameter ``beamshearcoeffxz`` allows to set shear correction factor for 2D beam sections, 
     controlling shear effective area used to evaluate shear force (default value is 1.0).
     Elements using this cross section model must implement layered cross section
     extension. For information see element library manual.

-  | Fibered cross section
   | ``FiberedCS`` ``nfibers #(in)`` ``fibermaterials #(ia)``
     ``thicks #(ra)`` ``widths #(ra)`` ``thick #(rn)``
     ``width #(rn)`` ``fiberycentrecoords #(ra)``
     ``fiberzcentrecoords #(ra)``
   | Cross section represented as a set of rectangular fibers. It is
     based on geometrical hypothesis, that cross sections remain planar
     after deformation (3d generalization of layered approach for
     beams). Paramater ``nfibers`` determines the number of fibers that
     together form the overall cross section. The model requires to
     specify a material model corresponding to particular fiber using
     ``fibermaterials`` array. This array should contain for each fibre
     corresponding material model number (the material model specified
     on element level has no meaning in this particular case). **The
     geometry of cross section is determined from fiber dimensions and
     fiber positions, all input in local coordinate system of the beam
     (yz plane).** The thick and width of each fiber are determined
     using ``thicks`` and ``widths`` arrays. The overall thick and width
     are specified using parameters ``thick`` and ``width``. Positions
     of particular fibers are specified by providing coordinates of
     center of each fiber using ``fiberycentrecoords`` array for
     y-coordinates and ``fiberzcentrecoords`` array for z-coordinates.

-  | Warping cross section
   | ``WarpingCS`` ``WarpingNode #(in)``
   | Represents the cross section for Free warping analysis, see section
     :ref:`FreeWarping`. The ``WarpingNode`` parametr defines the
     number of external node with prescribed boundary condition which
     corresponds to the relative twist of warping cross section.

.. _MaterialTypeRecords:

Material type records
---------------------

These records specify a material model description. The general format
is following:

``MaterialType`` ``#(in)`` ``d #(rn)``

The order of particular material records is optional, the material
number is determined by ``#(in)`` parameter. The numbering should
start from one and should end at n, where n is the number of records.
Material density is compulsory parameter and it’s value is given by
``d`` parameter.

Available material models, their outline and corresponding parameters
are described in separate **Material Library Manual**.

.. _NonlocalBarrierRecords:

Nonlocal barrier records
------------------------

Nonlocal material models of integral type are based on replacement of
certain suitable local quantity in local constitutive law by their
nonlocal counterparts, that are obtained as weighted average over some
characteristic volume. The weighted average is computed as a sum of a
remote value multiplied by weight function value. The weight function
typically depend on a distance between remote and receiver points and
decreases with increasing distance. In some cases, it is necessary to
disregard mutual interaction between some points (for example if they
are on the opposite sides of a thin notch, which prevents the nonlocal
interactions to take place). The barriers are the way how to introduce
these constrains. The barrier represent a curve (in 2D) or surface (in
3D). When the line connecting receiver and remote point intersects a
barrier, the barriers is activated and the corresponding interaction is
not taken into account.

Currently, the supported barrier types are following:

-  Polyline barrier

   | ``polylinebarrier`` ``#(in)`` ``vertexnodes #(ia)`` [``xcoordindx #(in)``]
     [``ycoordindx #(in)``]
   | This represents a polyline barrier for 2D problems. Barrier is a
     polyline, defined as a sequence of nodes representing vertices. The
     vertices are specified using parameter ``vertexnodes`` array, which
     contains the node numbers. The optional parameters ``xcoordindx``
     and ``ycoordindx`` allow to select the plane (xy, yz, or xz), where
     the barrier is defined. The ``xcoordindx`` is the first coordinate
     index, ``ycoordindx`` is the second. The default values are 1 for
     ``xcoordindx`` and 2 for ``ycoordindx``, representing barrier in xy
     plane.

-  Symmetry barrier

   | ``symmetrybarrier`` ``#(in)`` ``origin #(ra)`` ``normals #(ra)``
     ``activemask #(ia)``
   | Implementation of symmetry barier, that allows to specify up to
     three planes (orthogonal ones) of symmetry. This barrier allows to
     model the symmetry of the averaged field on the boundary without
     the need of modeling the other part of structure across the plane
     of symmetry. It is based on modifying the integration weights of
     source points to take into account the symmetry. The potential
     symmetry planes are determined by specifying orthogonal
     right-handed coordinate system, where axes represent the normals of
     corresponding symmetry planes. Parameter ``origin`` determines the
     origin of the coordinate system, the ``normals`` array contains
     three components of x-axis direction vector, followed by three
     components of y-axis direction vector (expressed in global
     coordinate system). The z-axis is determined from the orthogonality
     conditions. Parameter ``activemask`` allows to specify active
     symmetry planes; i-th nonzero value activates the symmetry barrier
     for plane with normal determined by corresponding coordinate axis
     (x=1, y=2, z=3).

.. _LoadBoundaryInitialConditions:

Load and boundary conditions
----------------------------

These records specify description of boundary conditions. The general
format is following:

``EntType`` ``#(in)`` ``loadTimeFunction #(in)`` [``valType #(in)``]
[``dofs #(ia)``] [``isImposedTimeFunction #(in)``]

The order of particular records is optional, boundary condition number
is determined by ``#(in)`` parameter. The numbering should start
from one and should end at n, where n is the number of records. Time
function value (given by ``loadTimeFunction`` parameter) is a
multiplier, using which each component (value of loading or value of
boundary condition) describes its time variation. The optional parameter
``valType`` allows to determine the physical meaning of bc value, which
is sometimes required. Supported values are (1 - temperature, 2 -
force/traction, 3 - pressure, 4 - humudity, 5 - velocity, 6 -
displacement). Another optional parameter ``dofs`` is used to determine
which dofs the boundary condition should act upon. It is not relevant
for all BCs..

The nonzero value of ``isImposedTimeFunction`` time function indicates
that given boundary condition is active, zero value indicates not active
boundary condition in given time (the bc does not exist). By default,
the boundary condition applies at any time.

Currently, EntType keyword can be one from

-  Dirichlet boundary condition

   ``BoundaryCondition`` ``prescribedvalue #(rn)`` [``d #(rn)``]

   Represents boundary condition. Prescribed value is specified using
   ``prescribedvalue`` parameter. The physical meaning of value is fully
   determined by corresponding DOF. Optionally, the prescribed value can
   be specified using ``d`` parameter. It is introduced for
   compatibility reasons. If ``prescribedvalue`` is specified, then
   ``d`` is ignored.

-  Prescribed gradient boundary condition (Dirichlet type)

   ``PrescribedGradient`` ``gradient #(rm)`` [``cCoords #(ra)``]

   Prescribes :math:`v_i = d_{ij}(x_j-\bar{x}_j)` or
   :math:`s = d_{1j}(x_j - \bar{x}_j)` where :math:`v_i` are primary
   unknowns, :math:`x_j` is the coordinate of the node, :math:`\bar x`
   is ``cCoords`` and :math:`d` is ``gradient``. The parameter
   ``cCoords`` defaults to zero. This is typical boundary condition in
   multiscale analysis where :math:`d = \partial_x s` would a
   macroscopic gradient at the integration point, i.e. this is a
   boundary condition for prolongation. It is also convenient to use
   when one wants to test a arbitrary specimen for shear.

-  Mixed prescribed gradient / pressure boundary condition (Active type)

   ``MixedGradientPressure`` ``devGradient #(ra)`` ``pressure #(rn)``
   [``cCoord #(ra)``]

   All boundary conditions of ensures that the deviatoric gradient and
   pressure is at least weakly fullfilled on the prescribed domain. They
   are used for computational homogenization of incompressible flow or
   elasticity problems.

-  Mixed prescribed gradient / pressure boundary condition (Weakly
   periodic type)

   ``MixedGradientPressureWeaklyPeriodic`` ``order #(rn)``

   Prescribes a periodic constant (unknown) stress tensor along the
   specified boundaries. For ``order`` set to 1, one obtains the same
   results as the Neumann boundary condition.

-  Mixed prescribed gradient / pressure boundary condition (Neumann
   type)

   ``MixedGradientPressureNeumann``

   Prescribes a constant (unknown) deviatoric stress tensor along the
   specified boundaries. Additional unknowns appears,
   :math:`\boldsymbol{\sigma}_\mathrm{dev}`, which is handled by the
   boundary condition itself (no control from the input file). The input
   devGradient is weakly fulfilled (homogenized over the elementsides).
   As with the the Dirichlet type, the volumetric gradient is free. This
   is useful in multiscale computations of RVE’s that experience
   incompressible behavior, typically fluid problems. In that case, the
   element sides should cover the entire RVE boundary. It is also
   convenient to use when one wants to test a arbitrary specimen for
   shear, with a free volumetric part (in which case the pressure is set
   to zero). Symmetry is not assumed, so rigid body rotations are
   removed, but translations need to be prescribed separately.

-  Mixed prescribed gradient / pressure boundary condition (Dirichlet
   type)

   ``MixedGradientPressureDirichlet``

   Prescribes
   :math:`v_i = d_{\mathrm{dev},ij}(x_j-\bar{x}_j) + d_\mathrm{vol}(x_i-\bar{x}_i)`,
   and a pressure :math:`p`. where :math:`v_i` are primary unknowns,
   :math:`x_j` is the coordinate of the node, :math:`\bar x` is
   ``cCoords`` and :math:`d_\mathrm{dev}` is ``devGradient``. The
   parameter ``cCoords`` defaults to zero. An additional unknown
   appears, :math:`d_\mathrm{vol}`, which is handled by the boundary
   condition itself (no control from the input file). This unknown is in
   a way related to the applied pressure. This is useful in multiscale
   computations of RVE’s that experience incompressible behavior,
   typically fluid problems. It is also convenient to use when one wants
   to test a arbitrary specimen for shear, with a free volumetric part
   (in which case the pressure is set to zero).

-  Nodal fluxes (loads)
   ``NodalLoad`` ``components #(ra)`` [``cstype #(in)``]
   Concentrated nodal load. The components of nodal load vector are
   given by ``components`` parameter. The size of this vector
   corresponds to a total number of nodal DOFs, and i-th value
   corresponds to i-th DOF in associated dof manager. The load can be
   defined in global coordinate system (``cstype`` = 0) or in entity -
   specific local coordinate system (``cstype`` = 1, default).

-  ``PrescribedTractionPressureBC``

   Represents pressure boundary condition (of Dirichlet type) due to
   prescribed tractions. In CBS algorithm formulation the prescribed
   traction boundary condition leads indirectly to pressure boundary
   condition in corresponding nodes. This boundary condition implements
   this pressure bc. The value of bc is determined from applied
   tractions, that should be specified on element edges/surfaces using
   suitable boundary loads.

- Linear constraint boundary condition

   ``LinearConstraintBC`` ``weights #(ra)`` [``weightsLtf #(ia)``]
   ``dofmans #(in)`` ``dofs #(in)`` ``rhs #(rn)``
   [``rhsLtf #(in)``] ``lhstype #(ia)`` ``rhsType #(ia)``

   This boundary condition implements a linear constraint in the form
   :math:`\sum_i w_ir_i = c`, where :math:`r_i` are unknowns related to
   DOFs determined by ``dofmans`` and ``dofs``, the weights are
   determined by ``weights`` and ``weightsLtf``. The constant is
   determined by ``rhs`` and ``rhsLtf`` parameters. This boundary
   condition is introduced as additional stationary condition using
   Lagrange multiplier, which is an additional degree of freedom
   introduced by this boundary condition.

   The individual DOFs are determined using dof manager numbers
   (``dofmans`` array) and corresponding DOF indices (``dofs``). The
   weights corresponding to participating DOFs are specified using
   ``weights`` array. The weights are multiplied by value returned by
   load time function, associated to individual weight using optional
   ``weightsLtf`` array. By default, all weights are set to 1. The
   constant :math:`c` is determined by ``rhs`` parameter and it is
   multiplied by the value of load time function, specified using
   ``rhsLtf`` parameter, or by 1 by default. The characteristic
   component, to which this boundary condition contributes must be
   identified using ``lhstype`` and ``rhsType`` parameters, values of
   which are corresponding to CharType enum. The left hand side
   contribution is assembled into terms identified by ``lhstype``. The
   rhs contribution is assembled into the term identified by ``rhsType``
   parameter. Note, that multiple values are allowed, this allows to
   select all variants of stifness matrix, for example. Note, that the
   size of ``dofmans``, ``dofs``, ``weights``, ``weightsLtf`` arrays
   should be equal.

.. _interactionbc:

-  InteractionBoundaryCondition 

   ``InteractionBoundaryCondition``

   Is a special boundary condition prescribed on
   InteractionPFEMParticles (see interactionparticle_
   in the PFEM part of the FluidStructureProblem. This sort of particles
   is regarded as it would have prescribed velocities, but the values
   change dynamically, as the solid part deforms. The velocities are
   obtained from coupled structural nodes.

- Body loads

   - Volume flux (load)

     | ``DeadWeight`` ``components #(ra)``
     | Represents dead weight loading applied on element volume (for
       structural elements). For transport problems, it represents the
       internal source, i.e. the rate of (heat) generated per unit volume.
       The magnitude of load for specific i-th DOF is computed as product
       of material density, corresponding volume and i-th member of
       ``components`` array.

   - Structural temperature load

     ``StructTemperatureLoad`` ``components #(ra)``

     Represents temperature loading imposed to some elements. The members
     of ``components`` array represent the change of temperature (or
     change of temperature gradient) corresponding to specific element
     strain components. See element library manual for details.

   - Structural eigenstrain load
     ``StructEigenstrainLoad`` ``components #(ra)``

     Prescribes eigenstrain (or stress-free
     strain) to a structural element. The array of ``components`` is
     defined in the global coordinate system. The number of components
     corresponds to a material mode, e.g. plane stress has three
     components and 3D six. Periodic boundary conditions can be imposed
     using eigenstrains and master-slave nodes. Consider decomposition of
     strain into average and fluctuating par
     .. math:: \boldsymbol{\varepsilon}(\boldsymbol{x}) = \langle \boldsymbol{\varepsilon} \rangle + \boldsymbol{\varepsilon}^*(\boldsymbol{x})

     where :math:`\langle \boldsymbol{\varepsilon} \rangle` can be imposed
     as eigenstrain over the domain and the solution gives the fluctuating
     part :math:`\boldsymbol{\varepsilon}^*(\boldsymbol{x})`. Master-slave
     nodes have to interconnect opposing boundary nodes of a unit cell.

- Boundary loads
   -  Constant edge fluxes (load)

     ``ConstantEdgeLoad`` ``loadType #(in)`` ``components #(ra)``
     [``dofexcludemask #(ia)``] [``csType #(in)``]
     [``properties #(dc)``] [``propertytf #(dc)``]

   - Constant surface fluxes (load)

     ``ConstantSurfaceLoad`` ``loadType #(in)`` ``components #(ra)``
     [``dofexcludemask #(ia)``] [``csType #(in)``]
     [``properties #(dc)``] [``propertytf #(dc)``]

     Represent constant edge/surface loads or boundary conditions.
     Parameter ``loadType`` distinguishes the type of boundary condition.
     Supported values are specified in bctype.h:

     *  ``loadType`` = 2 prescribed flux input (Neumann boundary
        condition),

     *  ``loadType`` = 3 uniform distributed load or the convection
        (Newton) BC. Parameter ``components`` contains the environmental
        values (temperature of the environment) corresponding to element
        unknowns, and ``properties`` dictionary should contain value of
        transfer (convection) coefficient (assumed to be a constant) under
        the key ’a’,

      * ``loadType`` = 7 specifies radiative boundary condition
        (Stefan-Boltzmann). It requires to specify emmisivity
        :math:`\varepsilon\in\langle 0,1\rangle`, the ``components`` array
        contains the environmental values (temperature of the
        environment). Default units are Celsius. Optional parameter
        ``temperOffset`` = 0 can be used to calculate in Kelvin.


     If the boundary condition corresponds to distributed force load, the
     ``components`` array contains components of distributed load
     corresponding to element unknowns. The load is specified for all DOFs
     of object to which is associated. For some types of boundary
     conditions the zero value of load does not mean that the load is not
     applied (Newton’s type of bc, for example). Then some mask, which
     allows to exclude specific dofs is necessary. The ``dofexcludemask``
     parameter is introduced to alow this. It should have the same size as
     ``components`` array, and by default is filled with zeroes. If some
     value of dofExcludeMask is set to nonzero, then the corresponding
     componentArray is set to zero and load is not applied for this DOF.
     If the boundary condition corresponds to prescribed flux input, then
     the ``components`` array contains the components of prescribed input
     flux corresponding to element unknowns.

     The properties can vary in time. Each property can have associated
     time function which determines its time variation. The time functions
     are set up using optional ``propertytf`` dictionary, containing for
     selected properties the corresponding time function number. The time
     function must be registered under the same key as in ``properties``
     dictionary. The property value is then computed by product of
     property value (determined by ``properties``) and corresponding time
     function evaluated at given time. If no time function provided for
     particula property, a unit constant function is assumed.

     The load can be defined in global coordinate system (``csType`` = 0,
     default) or in entity - specific local coordinate system (``csType``
     = 1).

   -  Linear edge flux (load)

     ``LinearEdgeLoad`` ``loadType #(in)`` ``components #(ra)``
     [``dofexcludemask #(ia)``] [``csType #(in)``]

     Represents linear edge load. The meanings of parameters ``csType``
     and ``loadType`` are the same as for **ConstantEdgeLoad**. In
     ``components`` array are stored load components for corresponding
     unknowns at the beginning of edge, followed by values valid for end
     of edge. The load can be defined in global coordinate system
     (``csType`` = 0, default) or in entity - specific local coordinate
     system (``csType`` = 1).

   -  User Defined Boundary Flux (load) 
      ``UserDefinedBoundaryLoad`` ``loadTimeFunction #(in)`` ``dofs #(ia)`` ``set #(in)`` ``intensityfunction #(in)`` ``geomtype #(in)`` ``approxorder #(in)`` ``components 0``

      Represents user defined boundary flux, where the flux is specified using configured function (see :ref:`TimeFunctionsRecords`). Note that PythonExpression provides complete configurability. 
      The flux intensity is defined by the ``intensityfunction``, which is a reference to specific function. 
      The function is defined in the global coordinate system as a function of position (x array variable) and time (t variable). 
      The ``geomtype`` parameter determines the type of geometry, where the flux is applied (4 for surface flux (defaut), 3 for edge flux). 
      The ``approxorder`` parameter determines the (approximate) order of intensity function. This is used to set up integration rule on the boundary.
      The ``dofs`` array determines the DOFs, to which the flux is applied. 
      The ``set`` parameter determines the set of element boundary entities, to which the flux is applied. 
      The ``loadTimeFunction`` and ``components`` parameters have no effect, but needs to be provided.

   -  InteractionLoad

     ``InteractionLoad`` ``ndofs #(in)`` ``loadType #(in)``
     ``Components #(ra)`` [``csType #(in)``]
     ``coupledparticles #(ia)``

     Represents a fluid pressure induced load in the solid part of the
     FluidStructureProblem. The meanings of parameters ``ndofs``,
     ``csType``, and ``loadType`` are the same as for **LinearEdgeLoad**.
     In ``Components`` array are stored load components for corresponding
     unknowns at the beginning of edge (``ndofs`` values), followed by
     values valid for end of edge (``ndofs`` values). The load should be
     defined in global coordinate system (``csType`` = 0) as it acts in
     normal direction of the edge. Array ``coupledparticles`` assign
     PFEMParticles from the fluid part of the problem providing fluid
     pressure.
     
          
     -  Structural penalty contact boundary condition

   ``structuralpenaltycontactbc`` ``loadTimeFunction #(in)`` ``dofs #(ia)``
   ``pn #(rn)`` ``pt #(rn)`` ``friction #(rn)`` ``mastersurface #(in)``
   ``slavesurface #(in)`` ``nsd #(in)``

   Represents a penalty-based contact boundary condition used to model
   contact interactions between deformable bodies. The contact formulation
   enforces normal and tangential constraints between surfaces defined by
   ``StructuralFEContactSurface`` records and their underlying
   ``StructuralContactElement_*`` elements.

   The parameters have the following meaning:

   - ``pn`` — normal penalty stiffness, controlling resistance against
     penetration between contacting surfaces.
   - ``pt`` — tangential penalty stiffness, controlling tangential
     response.
   - ``friction`` — coefficient of friction (currently experimental and
     under development). Set to ``0.0`` for frictionless contact.
   - ``mastersurface`` and ``slavesurface`` — identifiers of the master
     and slave contact surfaces defined by corresponding
     ``StructuralFEContactSurface`` records (see :ref:`ContactSurfaceRecords`).
   - ``dofs`` — list of affected degrees of freedom (typically ``1 2`` in 2D
     problems or ``1 2 3`` in 3D problems).
   - ``nsd`` — number of spatial dimensions (2 for plane strain or plane stress, or 3 for
     3D problems).

   The contact is enforced via the penalty method and contributes to the
   residual and tangent system of equations at each iteration. This
   boundary condition supports frictionless contact and a preliminary
   version of frictional contact, which is currently experimental and
   under development.

   **Example:**

   ::

      structuralpenaltycontactbc 3 loadTimeFunction 1 dofs 2 1 2 \
          pn 1.e8 pt 1.e8 friction 0.0 mastersurface 1 slavesurface 2 nsd 2

   This example defines a frictionless penalty contact condition between
   master surface 1 and slave surface 2 with normal and tangential
   stiffness equal to 1.e8. The condition acts on degrees of freedom 1 and
   2 (displacement is X and Y direction) in a 2D plane strain/stress domain.

   For bidirectional contact, two such boundary conditions can be defined
   with swapped master and slave surfaces.

.. _InitialConditions:

Initial conditions
------------------

These records specify description of initial conditions. The general
format is following:

``InitialCondition`` ``#(in)`` ``conditions #(dc)`` (deprecated)
   
``InitialCondition`` ``#(in)`` ``f #(expr)`` ``dfdt #(expr)`` ``d2fdt2 #(expr)`` ``set #(in)`` ``dofs #(ia)`` 

The order of particular
records is optional, load, boundary or initial condition number is
determined by (``num``\ #)(in) parameter. The numbering should start
from one and should end at n, where n is the number of records. Initial
parameters are listed in ``conditions`` dictionary using keys followed
by their initial values. Now ’v’ key represents velocity and ’a’ key
represents acceleration.

The second alternative allows to specify initial conditions providing expressions for unknown value, velocity and acceleration using expressions (``f``, ``dfdt``, ``d2fdt2`` parameters). 
The expressions can depend on position (’x’, ’y’, and ’z’ variables corresponding to x,y, and z coordinates of given node (dof manager)). 
The individual DOFs are determined using ``dofs`` array. 
The ``set`` parameter determines the set of nodes, to which the initial condition is applied.

.. _TimeFunctionsRecords:

(Time) Functions records
----------------------------

These records specify description of time functions, which generally
describe time variation of components during solution. The general
format is following:

``TimeFunctType`` ``#(in)`` [``initialValue #(rn)``]

The order of these records is optional, time function number is
determined by ``#(in)`` parameter. The ``initialValue`` parameter
allows to control the way, how increment of receiver is evaluated for
the first solution step. This first solution step increment is evaluated
as the difference of value of receiver at this first step and given
initial value (which is by default set to zero). The increments of
receiver in subsequent steps are computed as a difference between
receiver evaluated at given solution step and in previous step.

The numbering should start from one and should end at n, where n is the
number of records.

Currently, TimeFunctType keyword can be one from

-  Constant function

   ``ConstantFunction`` ``f(t) #(rn)``

   Represents the constant time function, with value ``f(t)``.

-  Peak function

   ``PeakFunction`` ``t #(rn)`` ``f(t) #(rn)``

   Represents peak time function. If time is equal to ``t`` value, then
   the value of time function is given by ``f(t)`` value, otherwise zero
   value is returned.

-  Piecewise function

   ``PiecewiseLinFunction`` [``nPoints #(in)`` ``t #(ra)`` ``f(t) #(ra)``] [
   ``datafile #("string")``]

   Represents the piecewise time function. The particular time values in
   ``t`` array should be sorted according to time scale. Corresponding
   time function values are in ``f(t)`` array. Value for time, which is
   not present in ``t`` is computed using liner interpolation scheme.
   Number of time-value pairs is in ``nPoints`` parameter.

   The second alternative allows reading input data from an external
   ASCII file. A hash commented line (#) is skipped during reading. File
   name should be eclosed with " ".

-  Heaviside-like time function

   ``HeavisideLTF`` ``origin #(rn)`` ``value #(rn)``

   Up to time, given by parameter ``origin``, the value of time function
   is zero. If time greater than ``origin`` parameter, the value is
   equal to parameter ``value`` value.

-  User defined

   ``UsrDefLTF`` ``f(t) #(expr)`` [``dfdt(t) #(expr)``]
   [``d2fdt2(t) #(expr)``]

   Represents user defined time function. The expressions can depend on
   “t” parameter, for which actual time will be substituted and
   expression evaluated. The function is defined using ``f(t)``
   parameter, and optionally, its first and second time derivatives
   using ``dfdt(t)`` and ``d2fdt2(t)`` parameters. The first and second
   derivatives may be required, this depend on type of analysis.

-  User defined Python expression <requires to compile with  USE_PYTHON_EXTENSION = ON >

   ``PythonExpression`` (``f #(s)`` | ``fpath #(s)``) [``dfdt #(s)`` | ``dfdtpath #(s)``] [``d2fdt2 #(s)`` | ``d2fdt2path #(s)``]


   Represents user defined function defined by Python script. The expressions can depend on
   ``t`` parameter, for which actual time will be substituted and ``x`` array containing the position.
   The python expression or script should define ``ret`` variable, which will be returned as the value of the function.

   The parameter ``f`` is an string containing Python expression, alternativaly, the parameter ``fpath`` is a string containing the path to the Python script file defining the python code for the function. 
   Optional parameters ``dfdt``, ``dfdtpath``, ``d2fdt2``, and ``d2fdt2path`` allow to prowide expressions or path to python scripts for evaluating derivative and second derivative of a function, 
   that may be required, depending on the context of use.


.. _XFEMManagerRecords:

Xfem manager record and associated records
------------------------------------------

| This record specifies the number of enrichment items and simulation
  options common for all enrichment items. Functions used for enrichment
  (e.g. Heaviside, abs or branch functions) are not specified here, they
  are specified for each enrichment item separately. The same holds for
  the geometrical representation of each enrichment item (e.g. a polygon
  line or a circle). Currently, OOFEM supports XFEM simulations of
  cracks and material interfaces in 2D. The input format for the XFEM
  manager is:
| ``XfemManager`` ``numberofenrichmentitems #(in)``
  ``numberofgppertri #(in)`` ``debugvtk #(in)``
  ``vtkexport #(in)`` ``exportfields #(in)``

  where
  ``numberofenrichmentitems`` represents number of enrichment items,
  ``numberofgppertri`` denotes the number of Gauss points in each
  subtriangle of a cut element (default 12) and ``debugvtk`` controls if
  additional debug vtk files should be written (1 activates the option,
  0 is default).

| The specification of an enrichment item may consist of several lines,
  see e.g. the test *sm/xFemCrackValBranch.in*. First, the enrichment
  item type is specified together with some optional parameters
  according to
| ``EntType`` ``#(in)`` ``enrichmentfront #(in)`` ``propagationlaw #(in)``

  where ``enrichmentfront`` specifies an
  enrichment front (we may for example employ branch functions at a
  crack tip and Heaviside enrichment along the rest of the crack, hence
  the “front” of the enrichment is treated separately) and
  ``propagationlaw`` specifies a rule for crack propagation (this
  feature is still highly experimental though). Specification of an
  ``enrichmentfront`` and a ``propagationlaw`` is optional.

| The next line specifies the enrichment function to be used:
| ``EntType`` ``#(in)``

| This is followed by a line specifying the geometric description (e.g.
  a polygon line or a circle) according to
| ``EntType`` ``#(in)`` ``extra attributes``
  where the number and type of extra
  attributes to specify will vary depending on the geometry chosen, e.g.
  center and radius for a circle or a number of points for a polygon
  line.

| If an enrichment front was specified previously, the type and
  properties of the enrichment front are specified on the next line
  according to
| ``EntType`` ``#(in)`` ``extra attributes``

| If a propagation law was specified previously, it’s type and
  properties are also specified on a separate line according to
| ``EntType`` ``#(in)`` ``extra attributes``


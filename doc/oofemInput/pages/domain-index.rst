.. _DomainRecord:

Domain record
=============

The domain record describes the whole domain and its type.  Depending on the
problem type there may be one or several domain records; unless stated
otherwise, one domain record is the default for all problem types.

The domain type resolves the default number of DOFs in a node and their
physical meaning.

.. record::

   :descitem:`domain` :elemparam:`domainType{et}`

Domain types
------------

The keyword is matched on a prefix, so any trailing text after the name below
is ignored.  Since the whole input line is lower-cased before parsing, case
does not matter — but punctuation does: write ``2dtruss``, not ``2d-Truss``.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - ``domainType``
     - Default DOFs
     - Physical meaning
   * - ``2dplanestress``
     - 2
     - *u*-displacement, *v*-displacement
   * - ``2dplanestressrot``
     - 3
     - *u*-, *v*-displacement and in-plane rotation
   * - ``planestrain``
     - 2
     - *u*-displacement, *v*-displacement
   * - ``3d``
     - 3
     - *u*-, *v*-, *w*-displacement
   * - ``3daxisymm``
     - 2
     - radial and axial displacement
   * - ``2dmindlinplate``
     - 3
     - *w*-displacement, *u*-rotation, *v*-rotation
   * - ``3dshell``
     - 6
     - displacement and rotation about each axis
   * - ``3ddirshell``
     - 7
     - seven-parameter shell based on director fields
   * - ``1dtruss``
     - 1
     - *u*-displacement
   * - ``2dtruss``
     - 2
     - *u*-displacement, *v*-displacement
   * - ``2dbeam``
     - 3
     - *u*-displacement, *w*-displacement, *v*-rotation
   * - ``warping``
     - 1
     - warping DOF, for :ref:`FreeWarping`
   * - ``heattransfer``
     - 1
     - temperature
   * - ``mass1transfer``
     - 1
     - mass concentration of one matter
   * - ``hema1``
     - 2
     - coupled heat and mass transfer of one matter
   * - ``2dincompflow``
     - 3
     - *u*-velocity, *v*-velocity, pressure
   * - ``3dincompflow``
     - 4
     - *u*-, *v*-, *w*-velocity and pressure
   * - ``2dlattice``
     - 3
     - 2D lattice structural mode
   * - ``3dlattice``
     - 6
     - 3D lattice structural mode
   * - ``2dmasslatticetransport``
     - 1
     - 2D lattice mass transport
   * - ``3dmasslatticetransport``
     - 1
     - 3D lattice mass transport
   * - ``2dheatlattice``, ``3dheatlattice``
     - 1
     - lattice heat transfer

.. note::

   In the current implementation both ``2dheatlattice`` and ``3dheatlattice``
   resolve to the *mass* transport lattice mode in
   ``Domain::resolveDomainDofsDefaults``.  Since both modes declare a single
   DOF per node this rarely shows, but do not rely on the heat-lattice
   spellings selecting a distinct mode.

For the ``2dmindlinplate`` mode the strain vector contains
:math:`\kappa_{xx}`, :math:`\kappa_{yy}`, :math:`\kappa_{xy}`,
:math:`\gamma_{xz}`, :math:`\gamma_{yz}`, and the stress vector contains
:math:`m_{xx}`, :math:`m_{yy}`, :math:`m_{xy}`, :math:`q_{xz}`,
:math:`q_{yz}`.

The default number of DOFs per node, and their physical meaning, can be
overridden in an individual dof manager record; see
:ref:`NodeElementSideRecords`.

.. _OutputManagerRecord:

Output manager record
---------------------

The output manager controls the text output.  It can filter output to specific
solution steps and, within those steps, to specific dof managers and elements.

.. record::

   :descitem:`OutputManager` :optelemparam:`tstep_all{}`
   :optelemparam:`tstep_step{in}` :optelemparam:`tsteps_out{rl}`
   :optelemparam:`dofman_all{}` :optelemparam:`dofman_output{rl}`
   :optelemparam:`dofman_except{rl}` :optelemparam:`element_all{}`
   :optelemparam:`element_output{rl}` :optelemparam:`element_except{rl}`

**Selecting solution steps**

:optparam:`tstep_all{}`
    Produce output in all solution steps.

:optparam:`tstep_step{in}`
    Produce output in every :param:`tstep_step`-th step.

:optparam:`tsteps_out{rl}`
    Produce output only in the listed solution steps.  May be combined with
    :param:`tstep_step`.

**Selecting components**

Output can also be filtered to specific dof managers and elements.  When such
members are selected, the output happens only in the selected solution steps.

:optparam:`dofman_all{}`, :optparam:`element_all{}`
    Select all dof managers, respectively all elements.

:optparam:`dofman_output{rl}`, :optparam:`element_output{rl}`
    Select only the listed members.

:optparam:`dofman_except{rl}`, :optparam:`element_except{rl}`
    De-select members previously selected, by listing their component numbers.

**Examples**

.. code-block:: none

   dofman_output {1 3}          prints nodes 1, 3
   dofman_output {(1 3)}        prints nodes 1, 2, 3
   element_output {1 3}         prints elements 1, 3
   element_output {(1 3)}       prints elements 1, 2, 3
   element_output {(1 3) 5 6}   prints elements 1, 2, 3, 5, 6

See :ref:`attribute-types` for the range list (**rl**) syntax.

.. _ComponentsSizeRecord:

Components size record
----------------------

This record declares how many components of each kind the domain contains.  The
corresponding records follow immediately in the input file, and the counts must
match the records actually present.

.. record::

   :elemparam:`ndofman{in}` :elemparam:`nelem{in}`
   :elemparam:`ncrosssect{in}` :elemparam:`nmat{in}`
   :elemparam:`nbc{in}` :elemparam:`nic{in}` :elemparam:`nltf{in}`
   :optelemparam:`nset{in}` :optelemparam:`ncontactsurf{in}`
   :optelemparam:`nbarrier{in}` :optelemparam:`nxfemman{in}`
   :optelemparam:`nfracman{in}` :optelemparam:`nsd{in}`
   :optelemparam:`axisymm{}` :optelemparam:`topology{s}`

**Parameters**

:param:`ndofman{in}`
    Number of dof managers, for example nodes, and their records.

:param:`nelem{in}`
    Number of elements and their records.

:param:`ncrosssect{in}`
    Number of cross sections and their records.

:param:`nmat{in}`
    Number of material models and their records.

:param:`nbc{in}`
    Number of boundary conditions — including loads and contact conditions —
    and their records.

:param:`nic{in}`
    Number of initial conditions and their records.

:param:`nltf{in}`
    Number of time functions and their records.

:optparam:`nset{in}`
    Number of set records; see :ref:`SetRecords`.

:optparam:`ncontactsurf{in}`
    Number of contact surface records; see :ref:`ContactSurfaceRecords`.

:optparam:`nbarrier{in}`
    Number of nonlocal barriers and their records; see
    :ref:`NonlocalBarrierRecords`.

:optparam:`nxfemman{in}`
    Number of XFEM managers, ``0`` or ``1``; see :ref:`XFEMManagerRecords`.

:optparam:`nfracman{in}`
    Number of fracture managers, ``0`` or ``1``.

:optparam:`nsd{in}`
    Number of spatial dimensions of the domain.  When omitted it is deduced
    from the domain type.

:optparam:`axisymm{}`
    A flag; when present, the problem is axisymmetric.

:optparam:`topology{s}`
    Name of a topology description to create.

If an optional count is omitted, none of the corresponding entities is assumed
to be present.

.. warning::

   A count larger than the number of records actually present leaves unread
   lines in the input file, and OOFEM warns about them.  The most common cause
   is a missing ``nset`` in this record.

Domain component records
------------------------

The remaining records describe the individual domain components.

.. toctree::
    :maxdepth: 2

    domain-dofmanagers
    domain-elements
    domain-sets
    domain-contact-surfaces
    domain-cross-sections
    domain-materials
    domain-loads-bcs
    domain-initial-conditions
    domain-time-functions
    domain-xfem

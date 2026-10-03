.. _anatomy-of-an-input-file:

Anatomy of an input file
========================

Before reading the record reference it helps to see a whole input file at once.
The file below is ``tests/regression/sm/patch100.in``: a linear static analysis
of five plane-stress quadrilaterals in pure compression.  It is complete — this
is everything OOFEM needs.

.. code-block:: none
   :linenos:

   patch100.out
   Patch test of PlaneStress2d elements -> pure compression
   LinearStatic nsteps 1
   domain 2dPlaneStress
   OutputManager tstep_all dofman_all element_all
   ndofman 8 nelem 5 ncrosssect 1 nmat 1 nbc 3 nic 0 nltf 1 nset 3
   node 1 coords 3  0.0   0.0   0.0
   node 2 coords 3  0.0   4.0   0.0
   node 3 coords 3  2.0   2.0   0.0
   node 4 coords 3  3.0   1.0   0.0
   node 5 coords 3  8.0   0.8   0.0
   node 6 coords 3  7.0   3.0   0.0
   node 7 coords 3  9.0   0.0   0.0
   node 8 coords 3  9.0   4.0   0.0
   PlaneStress2d 1 nodes 4 1 4 3 2  NIP 1
   PlaneStress2d 2 nodes 4 1 7 5 4  NIP 1
   PlaneStress2d 3 nodes 4 4 5 6 3  NIP 1
   PlaneStress2d 4 nodes 4 3 6 8 2  NIP 1
   PlaneStress2d 5 nodes 4 5 7 8 6  NIP 1
   Set 1 elementranges {(1 5)}
   Set 2 nodes 2 1 2
   Set 3 nodes 2 7 8
   SimpleCS 1 thick 1.0 width 1.0 material 1 set 1
   IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0
   BoundaryCondition 1 loadTimeFunction 1 dofs 2 1 2 values 1 0.0 set 2
   BoundaryCondition 2 loadTimeFunction 1 dofs 1 2 values 1 0.0 set 3
   NodalLoad 3 loadTimeFunction 1 dofs 2 1 2 components 2 2.5 0.0 set 3
   ConstantFunction 1 f(t) 1.0

Reading it line by line
-----------------------

The first six lines are the *header*.  They are positional: each one must appear
exactly once, in exactly this order.

Line 1 — output file
    The path of the output file, with no keyword.  An existing file of the same
    name is overwritten.  See :ref:`OutputFileRecord`.

Line 2 — job description
    Free text describing the job, with no keyword.  It is copied into the output
    file.  See :ref:`JobDescriptionRecord`.

Line 3 — analysis
    ``LinearStatic nsteps 1`` selects the analysis type and its parameters.  The
    record keyword is the name of the engineering model; here one solution step
    is requested, i.e. one load case.  See :ref:`AnalysisRecord`.

Line 4 — domain
    ``domain 2dPlaneStress`` fixes the default number of degrees of freedom per
    node and their physical meaning — here two displacements, *u* and *v*.  See
    :ref:`DomainRecord`.

Line 5 — output manager
    ``OutputManager tstep_all dofman_all element_all`` asks for all nodes and
    all elements to be printed in every solution step.  See
    :ref:`OutputManagerRecord`.

Line 6 — components size
    ``ndofman 8 nelem 5 ...`` declares how many records of each kind follow.
    These counts must match the records actually present, or the analysis stops.
    See :ref:`ComponentsSizeRecord`.

Everything after the header is a *component record*.  Component records are
grouped by kind, and within each group the number following the keyword is the
component's label:

Lines 7–14 — dof managers
    ``node 8 coords 3 ...`` — eight nodes, each with three coordinates.  The
    ``3`` after ``coords`` is the array size, not a coordinate; see
    :ref:`attribute-types`.  See :ref:`NodeElementSideRecords`.

Lines 15–19 — elements
    ``PlaneStress2d 1 nodes 4 1 4 3 2 NIP 1`` — element 1 is a four-node
    quadrilateral on nodes 1, 4, 3, 2 with one integration point.  Element types
    and their parameters are documented in the **Element Library Manual**.  See
    :ref:`ElementsRecords`.

Lines 20–22 — sets
    A set names a group of nodes, elements, edges or surfaces so that cross
    sections and boundary conditions can be attached to it.  ``Set 1`` is all
    five elements, ``Set 2`` the two nodes on the left edge, ``Set 3`` the two
    on the right.  See :ref:`SetRecords`.

Line 23 — cross section
    ``SimpleCS 1 thick 1.0 width 1.0 material 1 set 1`` gives the section
    properties and binds material 1 to the elements of set 1.  See
    :ref:`CrossSectionRecords`.

Line 24 — material
    ``IsoLE 1 d 0. E 15.0 n 0.25 talpha 1.0`` — isotropic linear elastic
    material with zero density, :math:`E = 15`, :math:`\nu = 0.25`.  Material
    models are
    documented in the **Material Library Manual**.  See
    :ref:`MaterialTypeRecords`.

Lines 25–27 — boundary conditions and loads
    ``BoundaryCondition`` prescribes a value (here zero displacement) on the
    DOFs listed in ``dofs`` for the nodes of the given set; ``NodalLoad``
    applies a concentrated force.  Each refers to a time function through
    ``loadTimeFunction``.  See :ref:`LoadBoundaryInitialConditions`.

Line 28 — time function
    ``ConstantFunction 1 f(t) 1.0`` — the multiplier used by the boundary
    conditions above, constant and equal to one.  See
    :ref:`TimeFunctionsRecords`.

Things worth knowing early
--------------------------

Input is case-insensitive, except inside quotes
    Every input line is converted to lower case before parsing, so ``node``,
    ``Node`` and ``NODE`` are the same keyword.  This manual capitalises record
    keywords for readability only.  Text between double quotes is left
    untouched, which matters for file paths and for Python expressions (see
    :ref:`TimeFunctionsRecords`).

Arrays carry their own length
    Every integer or real array starts with its size.  ``dofs 2 1 2`` is "two
    DOFs, numbers 1 and 2" — a frequent source of confusion.  See
    :ref:`attribute-types`.

Component numbers are labels
    Nodes and elements may be numbered arbitrarily, even non-contiguously; the
    number is a label by which the component is referenced.  Cross sections,
    materials, boundary conditions, initial conditions and time functions,
    however, must be numbered from 1 to *n* without gaps.

Attribute order within a record is free
    Only the order of the *records* is fixed.

One record per line, unless continued
    End a line with ``\`` to continue a record onto the next line, and start a
    line with ``#`` to write a comment.

Where to go next
----------------

* :ref:`syntax-and-general-rules` — the full record order and the notation used
  in the reference sections.
* :ref:`AnalysisRecord` — the analysis types OOFEM can run.
* :doc:`examples` — further complete input files, including parallel ones.
* The **Element Library Manual** and **Material Library Manual** document the
  element and material records themselves, which this manual only frames.

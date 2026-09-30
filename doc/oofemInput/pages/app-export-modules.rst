.. _ExportModulesSec:

Export modules
==============

Export modules write the computed data for post-processing by external
software.  The number of export module records is declared by :param:`nmodules`
in the analysis record; see :ref:`AnalysisRecord`.

.. record:: Record syntax common to all export modules

   :descitem:`EntType` :optelemparam:`tstep_all{}`
   :optelemparam:`tstep_step{in}` :optelemparam:`tsteps_out{rl}`
   :optelemparam:`subtsteps_out{in}` :optelemparam:`domain_all{}`
   :optelemparam:`domain_mask{ia}` :optelemparam:`regionsets{ia}`
   :optelemparam:`timescale{rn}` :optelemparam:`pythonexport{}`

Selecting solution steps
------------------------

:optparam:`tstep_all{}`
    Produce output in all solution steps.

:optparam:`tstep_step{in}`
    Produce output in every :param:`tstep_step`-th step.

:optparam:`tsteps_out{rl}`
    Produce output only in the listed solution steps.

:optparam:`subtsteps_out{in}`
    Set to ``1`` to export intermediate results as well, for example during
    substepping or the individual equilibrium iterations.  Requires support
    from the solver.

Selecting domains
-----------------

:optparam:`domain_all{}`
    Export all domains of the problem.

:optparam:`domain_mask{ia}`
    Export only the domains whose numbers are listed.

.. _exportregions:

Regions
-------

The export is done on a region basis: on each region the nodal recovery is
performed independently and the results are exported as a separate piece.  This
makes it possible to account for discontinuities, and to export variables
defined only by a particular material model.

Region volumes are defined by sets containing the individual elements.  By
default one region is created, containing every element in the problem domain.

:optparam:`regionsets{ia}`
    Numbers (ids) of the domain sets to use as regions instead.  Note that
    regions are determined solely from elements.

**Example**

.. code-block:: none

   vtkxml tstep_all cellvars 1 46 vars 1 1 primvars 1 1 stype 2 \
       regionsets 2 1 2

Other common parameters
-----------------------

:optparam:`timescale{rn}`
    Scales the time in the output.  In a transport problem the basic time unit
    is the second, so ``timescale 2.777777e-4`` — that is 1/3600 — converts all
    the time data in the VTK XML output from seconds to hours.

:optparam:`pythonexport{}`
    A flag enabling the export of the module's data to the Python interface.

VTK export modules
------------------

.. _vtkxml:

vtkxml, vtkhdf5 and vtk
~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   The ``vtk`` module is **deprecated**.  Use ``vtkxml`` or ``vtkhdf5``
   instead: ``vtkxml`` exports results recovered region by region and has many
   more features.

``vtkxml``
    Exports results into ``vtu`` files — VTK unstructured grids — one file per
    solution step.

``vtkhdf5``
    Exports results into a VTK HDF5 file, a single HDF file for all time steps.

    .. note::

       Available from version 1.6.  Requires the HDF5 library and a build
       configured with ``USE_HDF5=ON``; HDF5 version 1.14 or later is
       recommended.

Both output formats can be visualized with ParaView.

.. record::

   :descitem:`vtkxml` :optelemparam:`vars{ia}`
   :optelemparam:`primvars{ia}` :optelemparam:`cellvars{ia}`
   :optelemparam:`ipvars{ia}` :optelemparam:`externalforces{ia}`
   :optelemparam:`stype{in}` :optelemparam:`setmembership{}`

   :descitem:`vtkhdf5` :optelemparam:`vars{ia}`
   :optelemparam:`primvars{ia}` :optelemparam:`cellvars{ia}`
   :optelemparam:`ipvars{ia}` :optelemparam:`externalforces{ia}`
   :optelemparam:`stype{in}`

   :descitem:`vtk` :optelemparam:`vars{ia}`
   :optelemparam:`primvars{ia}` :optelemparam:`cellvars{ia}`
   :optelemparam:`stype{in}` :optelemparam:`regionstoskip{ia}`

**Parameters**

:optparam:`vars{ia}`
    Identifiers of the internal variables to export.  These variables are
    smoothed and transferred to the nodes.  The values are those of the
    ``InternalStateType`` enumeration in ``src/core/internalstatetype.h``.

:optparam:`primvars{ia}`
    Identifiers of the primary variables to export.  The values correspond to
    the ``UnknownType`` enumeration in ``src/core/unknowntype.h``.  Note that
    the enumeration values start from zero unless specified directly, and that
    not every value is supported by every material model or analysis type.

:optparam:`cellvars{ia}`
    Identifiers of constant variables defined on an element (cell), such as a
    material number.  The identifiers are in
    ``src/core/internalstatetype.h``.

:optparam:`ipvars{ia}`
    Identifiers of internal variables to export **without** smoothing,
    directly as a point dataset in which each point corresponds to an
    individual integration point.  A separate ``vtu`` file is created for this
    raw point data.

:optparam:`externalforces{ia}`
    Identifiers of external force components to export.

:optparam:`stype{in}`
    Smoothing procedure used to compute nodal values from the values at the
    integration points: ``0`` simple nodal averaging, generally supported only
    by triangular and tetrahedral elements; ``1`` Zienkiewicz-Zhu recovery, the
    default; ``2`` Superconvergent Patch Recovery (SPR), based on least-squares
    fitting.

:optparam:`setmembership{}`
    Triggers the export of set membership information for dof managers and
    cells, as a point and cell variable.  Since the VTK format does not support
    sparse sets, the membership is encoded into datasets of byte values
    (``UINT8``, named ``VertexSetMembership`` and ``CellSetMembership``) in
    which the *i*-th bit is 1 if the component is a member of set *i* and zero
    otherwise.

:optparam:`regionstoskip{ia}`
    ``vtk`` only: regions to leave out of the export.

vtkpfem
~~~~~~~

VTK PFEM (particle FEM) export.  Exports the particle positions to VTK as a
point dataset.

.. record::

   :descitem:`vtkpfem` :optelemparam:`vars{ia}`
   :optelemparam:`primvars{ia}` :optelemparam:`cellvars{ia}`
   :optelemparam:`ipvars{ia}` :optelemparam:`stype{in}`

vtkmemory
~~~~~~~~~

VTK memory export.  This module produces no output of its own; it prepares the
data structures needed to support VTK export or VTK visualization, and is used
by the Python interface to access the VTK datasets.

.. record::

   :descitem:`vtkmemory` :optelemparam:`vars{ia}`
   :optelemparam:`primvars{ia}` :optelemparam:`cellvars{ia}`
   :optelemparam:`ipvars{ia}`

vtkxmlxfem
~~~~~~~~~~

VTK XFEM export.  Exports XFEM-related data; what is exported is determined by
the XFEM manager's :param:`exportfields` parameter, see
:ref:`XFEMManagerRecords`.

.. record::

   :descitem:`vtkxmlxfem`

Homogenization export
---------------------

Homogenizes integration point quantities in the global coordinate system —
stress, strain, damage, heat flow and so on — by summing and averaging them over
the volume.  The region sets the averaging runs over can be selected.

Averaging works for all domains, with an extension for trusses: a truss is
treated as a volume element with oriented stress and strain components along
the truss axis, and the transformation to global components happens before the
averaging.

.. record::

   :descitem:`hom` :elemparam:`ists{ia}` :optelemparam:`scale{rn}`
   :optelemparam:`regionSets{ia}` :optelemparam:`reactions{}`
   :optelemparam:`strain_energy{}`

**Parameters**

:param:`ists{ia}`
    Internal state types to export, as defined in
    ``src/core/internalstatetype.h``.

:optparam:`scale{rn}`
    Multiplies all the averaged integration point quantities.  Default ``1``.

:optparam:`regionSets{ia}`
    Region sets to average over.  The whole domain is averaged by default.

:optparam:`reactions{}`
    A flag; when present, the reactions are exported as well.

:optparam:`strain_energy{}`
    Calculates the strain energy over the selected elements as

    .. math::

       W^*=\int_V \int \sigma \,\mathrm{d}
       (\varepsilon-\varepsilon_{eig}) \,\mathrm{d} V,

    where :math:`\sigma` is the stress tensor, :math:`\varepsilon` the strain
    tensor and :math:`\varepsilon_{eig}` the eigenstrain tensor, which
    originates from a temperature load or a prescribed eigenstrain.  The strain
    energy increment and the total strain energy are reported in each step.
    The integration uses the mid-point rule for the stress and is exact for
    linear elastic materials.

Gauss point export
------------------

Useful when a variable such as damage has to be plotted against a spatial
coordinate with a tool like gnuplot.  It generates files with the data organised
in columns, one row per Gauss point, so that for instance the damage
distribution along a one-dimensional bar can be plotted directly.

.. record::

   :descitem:`gpexportmodule` :optelemparam:`vars{ia}`
   :optelemparam:`ncoords{in}`

**Parameters**

:optparam:`vars{ia}`
    Identifiers of the internal variables to export, as defined in
    ``src/core/internalstatetype.h``.

:optparam:`ncoords{in}`
    Number of spatial coordinates to export at each Gauss point.  Depending on
    the spatial dimension of the domain a point may have one, two or three
    coordinates.

    ``-1``
        Export only the coordinates that are actually used.
    ``0``
        Export no coordinates.
    a positive integer
        Export exactly that many coordinates.  If it exceeds the actual number
        of coordinates, the actual ones are padded with zeros.

    For instance, in a 2D problem the actual number of coordinates is 2: with
    ``ncoords 3`` the two actual coordinates followed by a 0 are exported,
    while with ``ncoords 1`` only the first coordinate is exported.

Output file format
~~~~~~~~~~~~~~~~~~

The module creates a file with the extension ``gp`` after each step for which
output is performed.  The file starts with a header whose lines begin with
``#``, followed by the data section.  Each data line corresponds to one Gauss
point and contains:

#. element number,
#. material number,
#. Gauss point number,
#. contributing volume around the Gauss point,
#. Gauss point global coordinates, written as a real array of length
   :param:`ncoords`,
#. internal variables according to :param:`vars`, each written as a real array
   of the corresponding length.

**Example**

.. code-block:: none

   GPExportModule 1 tstep_step 100 domain_all ncoords 2 vars 5 4 13 31 64 65

The ``*.gp`` file is written after every 100 steps and contains, for each Gauss
point in the entire domain, its 2 coordinates and the internal variables of
type 4, 13, 31, 64 and 65 — the strain tensor, the damage tensor, the maximum
equivalent strain level, the stress work density and the dissipated work
density.  The material model must of course be able to deliver these.

The size of the strain tensor depends on the spatial dimension, and the size of
the damage tensor on the spatial dimension and on the type of model: a simple
isotropic damage model has just one component, while an anisotropic damage
model may have more.  The other variables in this example are scalars, but they
are still written as arrays of length 1, so the value is always preceded by a
``1`` giving the array length.

Since certain internal variables have the meaning of densities — per unit
volume or area, again depending on the spatial dimension — it is useful to have
the contributing volume of the Gauss point.  The product of that volume and the
density is an additive contribution to the total value of the corresponding
variable, which can be exploited, for example, to evaluate the total dissipated
energy over the entire domain.

Solution status monitor
-----------------------

.. note::

   Available from version 1.6.

Creates a configurable report on the progress of the solution.

.. record::

   :descitem:`solutionstatus` :optelemparam:`fmt{s}`

**Parameters**

:optparam:`fmt{s}`
    Controls what is reported for each solution step.  The string contains data
    codes separated by colons.  The default is ``"m:s:a:nite:t:dt:st:cr"``.

.. list-table:: Data codes accepted by ``fmt``
   :header-rows: 1
   :widths: 16 84

   * - Data code
     - Description
   * - ``m``
     - Meta step number
   * - ``s``
     - Solution step number
   * - ``a``
     - Attempt number
   * - ``nite``
     - Number of iterations
   * - ``t``
     - Solution step target time
   * - ``dt``
     - Solution step time increment
   * - ``st``
     - Time spent solving the solution step, in seconds
   * - ``cr``
     - Convergence reason status: *Converged*, *Diverged_I* (convergence not
       reached within the maximum iteration limit), *Diverged_T* (diverged
       solution), or *Failed*
   * - ``-``
     - Placeholder; prints a literal ``-``

**Example**

.. code-block:: none

   solutionstatus tstep_all fmt m:s:a:nite:t:dt:st:cr

Error checking module
---------------------

Implements the error checking rules used by the regression tests: it compares
selected computed results against reference values.

The individual rules are defined in a separate section with a dedicated syntax,
which can be part of the input file itself or live in an external file.  The
syntax follows that of the Extractor package, specifically its
`Value records section
<https://www.oofem.org/resources/doc/extractorInput/html/node2.html>`_.

.. record::

   :descitem:`errorcheck` :optelemparam:`filename{s}`
   :optelemparam:`extract{}` :optelemparam:`writeist{ia}`

**Parameters**

:optparam:`filename{s}`
    Path to the file in which the rules are defined.  By default the rules are
    read from the input file itself.

:optparam:`extract{}`
    A flag switching the module from testing to extraction: the selected
    quantities are written out rather than compared against expected values.
    In this mode output is produced in every time step, ignoring the ``tstep``
    selection.

:optparam:`writeist{ia}`
    Internal state types to write rules for, when generating a rule set.

By default the module performs the regression test and raises an error when one
or more checks fail.

.. _other-export-modules:

Other registered export modules
-------------------------------

The following export modules are registered but not documented here:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Record keyword
     - Purpose
   * - ``output``
     - The standard text output module that writes the ``.out`` file.
   * - ``dm``
     - Dof manager export, tabulating selected dof manager values.
   * - ``poi``
     - Export at user-defined points of interest.
   * - ``gnuplot``
     - Output prepared for direct plotting with gnuplot.
   * - ``matlab``
     - Output prepared for reading into MATLAB.
   * - ``crackvectorexport``
     - Exports crack vectors.
   * - ``vtkxmllattice``
     - VTK XML export specialised for lattice models.
   * - ``vtkxmlperiodic``
     - VTK XML export for periodic problems.
   * - ``qcvtkxml``
     - VTK XML export for the quasicontinuum formulation.

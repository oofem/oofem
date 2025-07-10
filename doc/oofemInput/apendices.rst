

Appendix
========

.. _sparselinsolver:

Sparse linear solver parameters
-------------------------------

| The sparselinsolverparams field has the following general syntax:
| [``lstype #(in)``] [``smtype #(in)``] ``solverParams #(string)``
  
where parameter ``lstype`` allows to select the solver for the linear system of equations. Currently
supported values are 0 (default) for a direct solver (ST_Direct), 1
for an Iterative Method Library (IML) solver (ST_IML), 2 for a Spooles
direct solver, 3 for Petsc library family of solvers, and 4 for a
DirectSparseSolver (ST_DSS). Parameter ``smtype`` allows to select the
sparse matrix storage scheme. The scheme should be compatible with the
solver type. Currently supported values (marked as “id”) are
summarized in table :ref:`linsolvstoragecompattable`. The 0
value is default and selects the symmetric skyline (SMT_Skyline).
Possible storage formats include unsymmetric skyline (SMT_SkylineU),
compressed column (SMT_CompCol), dynamically growing compressed column
(SMT_DynCompCol), symmetric compressed column (SMT_SymCompCol),
spooles library storage format (SMT_SpoolesMtrx), PETSc library matrix
representation (SMT_PetscMtrx, a sparse serial/parallel matrix in AIJ
format), and DSS compatible matrix representations (SMT_DSS). The
allowed ``lstype`` and ``smtype`` combinations are summarized in the
table :ref:`linsolvstoragecompattable`,
together with solver parameters related to specific solver.



.. raw:: latex

   \begin{landscape}

.. _linsolvstoragecompattable:
.. table:: Solver and storage scheme compatibility.

   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |Storage format  |``smtype id``| ``lstype (id)``                                                              |
   +================+=============+==========+=======+===========+=========+=======+==============+==============+
   |                | id          |Direct (0)|IML (1)|Spooles (2)|Petsc (3)|DSS (4)|MKLPardiso (6)|SuperLU_MT (7)|
   |                |             |          |       |           |         |       |Pardiso.org(8)|              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_Skyline     | 0           | +        |  +    |           |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_SkylineU    | 1           | +        |  +    |           |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+  
   |SMT_CompCol     | 2           |          |  +    |           |         |       |  +           |   +          |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_DynCompCol  | 3           |          |  +    |           |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_SymCompCol  | 4           |          |  +    |           |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_DynCompRow  | 5           |          |  +    |           |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_SpoolesMtrx | 6           |          |       |   +       |         |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_PetscMtrx   | 7           |          |       |           |  +      |       |              |              |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_DSS_sym_LDL | 8           |          |       |           |         |       |              |   +          |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_DSS_sym_LL  | 9           |          |       |           |         |       |              |   +          |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+
   |SMT_DSS_unsym_LU| 10          |          |       |           |         |       |              |   +          |
   +----------------+-------------+----------+-------+-----------+---------+-------+--------------+--------------+

.. raw:: latex

   \end{landscape}
   
The solver parameters in ``solverParams`` depend on the solver type and
are summarized in table
(sparsesolverparams_).

.. _sparsesolverparams:

.. table:: Solver parameters.

   ==================== == =======================================================================
   Solver type          id Solver parameters/notes
   ==================== == =======================================================================
   ST_Direct            0
   ST_IML               1  [``stype`` #(in)] ``lstol`` #(rn) ``lsiter`` #(in)\ ``lsprecond`` #(in)
   \                       [``precondattributes`` #(string)]
   \                       Included in OOFEM, requires to compile with USE_IML
   ST_Spooles           2  [``msglvl`` #(in)] [``msgfile`` #(s)]
   \                       http://www.netlib.org/linalg/spooles/spooles.2.2.html
   ST_Petsc             3  See Petsc manual, for details
   ST_DSS               4  Sparse direct solver, included in OOFEM
   \                       Requires to compile with USE_DSS
   ST_MKLPardiso        6  Requires Intel MKL Pardiso
   ST_SuperLU_MT        7  SuperLU for shared memory machines
   \                       http://crd-legacy.lbl.gov/ xiaoye/SuperLU/
   ST_PardisoProjectOrg 8  Requires Pardiso solver(http://www.pardiso-project.org/)
   ==================== == =======================================================================

In case of ``ST_PETSC``, the user can set several run-time options, e.g.,
   ``-ksp\_type`` ``[cg, gmres, bicg, bcgs]``
   ``-pc\_type`` ``[jacobi, bjacobi,none,ilu,...]``
   ``-ksp\_monitor`` ``-ksp\_rtol #`` ``-ksp\_view`` ``-ksp\_converged_reason``.
   These options will override those that are default (PETSC KSPSetFromOptions() routine is called after any other customization
   routines).}

The ``stype`` allows to select particular iterative solver from IML
library, currently supported values are 0 (default) for
Conjugate-Gradient solver, 1 for GMRES solver. Parameter ``lstol``
represents the maximum value of residual after the final iteration and
the ``lsiter`` is maximum number of iteration for iterative solver. The
``precondattributes`` parameters contains the optional preconditioner
parameters. The ``lsprecond`` parameter determines the type of
preconditioner to be used. The possible values of ``lsprecond`` together
with supported storage schemes and their descriptions are summarized in
table :ref:`precondtable`.

.. _precondtable:

.. table:: Preconditioning summary.

   ============ == ================== =========================================
   Precond type id Compatible storage Description and parameters
   ============ == ================== =========================================
   IML_VoidPrec 0  all                No preconditioning
   IML_DiagPrec 1  all                Diagonal preconditioning
   IML_ILUPrec  2  SMT_CompCol        Incomplete LU Decomposition
   \               SMT_DynCompCol     with no fill up
   IML_ILUPrec  3  SMT_DynCompRow     Incomplete LU (ILUT) with
   \                                  fillup.
   \                                  The ``precondattributes`` are:
   \                                  [``droptol`` #(rn)] [``partfill`` #(in)].
   \                                  ``droptol`` dropping tolerance
   \                                  ``partfill`` level of fill-up
   IML_ICPrec   4  SMT_SymCompCol     Incomplete Cholesky
   \               SMT_CompCol        with no fill up
   ============ == ================== =========================================

.. _eigensolverssection:

Eigen value solvers
-------------------

| The eigensolverparams field has the following general syntax:
| ``stype #(in)`` [``smtype #(in)``] ``solverParams #(string)``
  where parameter ``stype`` allows to
  select solver type. Parameter ``smtype`` allows to select sparse
  matrix storage scheme. The scheme should be compatible with solver
  type. Currently supported values of ``stype`` are summarized in
  table :ref:`eigenvaluesolverparamtable`.

.. _eigenvaluesolverparamtable:

.. table:: Eigen Solver parameters.

   ================== =========== =====================
   Solver type        stype id    solver parameters
   Subspace Iteration 0 (default)
   Inverse Iteration  1
   SLEPc solver       2           requires “smtype 7”
                                  see also SLEPc manual
   ================== =========== =====================

.. _dynamicloadbalancing:

| There are in general two basic factors causing load imbalance between
  individual subdomains: (i) one comming from application nature, such
  as switching from linear to nonlinear response in certain regions or
  local adaptive refinment, and (ii) external factors, caused by
  resourse realocation, typical for nondedicated cluster environments,
  where indivudual processors are shared by different applications and
  users, leading to time variation in allocated processing power. The
  load balance recovery is achieved by repartitioning of the problem
  domain and transferring the work (represented typically by finite
  elements) from one subdomain to another. This section describes the
  structure and syntax of parameters related to dynamic load balancing.
  The corresponding part of analysis record has the following general
  syntax:
| [``lbflag #(in)``] [``forcelb1 #(in)``] [``wtp #(ia)``]
  [``lbstep #(in)``] [``relwct #(rn)``] [``abswct #(rn)``]
  [``minwct #(rn)``]
  
where the parameters have following meaning:

-  ``lbflag``, when set to nonzero value activates the dynamic load
   balancing. Default value is zero.

-  ``forcelb1`` forces the load rebalancing after the first solution
   step, when set to nonzero value.

-  ``wtp`` allows to activate optional load balancing plugins. At
   present, the only supported value is 1, that activates nonlocal
   plugin, necessary for nonlocal averaging to work properly when
   dynamic load balancing is active.

-  ``lbstep`` rebalancing, if needed, is performed only every lbstep
   solution step. Default value is 1 (recover balance after every step,
   if necessary).

-  ``relwcr`` sets relative wall-clock imbalance treshold. When achieved
   relative imbalance between wall clock solution time of individual
   processors is greater than provided treshold, the rebalancing
   procedure will be activated.

-  ``abswct`` sets absolute wall-clock imbalance treshold. When achieved
   absolute imbalance between wall clock solution time of individual
   processors is greater than provided treshold, the rebalancing
   procedure will be activated.

-  ``minwct`` minimum absolute imbalance to perform relative imbalance
   check using ``relwcr`` parameter, otherwise only absolute check is
   done. Default value is 0.

At present, the load balancing support requires ParMETIS module to be
configured and compiled.

.. _errorestimators:

Error estimators and indicators
-------------------------------

The currently supported values of ``eetype`` are in table
:ref:`eetypestable`.

-  EET_SEI - Represents scalar error indicator. It indicates element
   error based on the value of some suitable scalar value (for example
   damage level, plastic strain level) obtained from the element
   integration points and corresponding material model.

-  EET_ZZEE - The implementation of Zienkiewicz Zhu Error Estimator. It
   requires the special element algorithms, which may not be available
   for all element types.

   Please note, that in the actual version, the error on the element
   level is evaluated using default integration rule. For example, in
   case of ZZ error estimator, the error (L2 or energy norm) is
   evaluated from the difference of computed and “recovered” stresses,
   which are approximated using the same interpolation functions as
   displacements). Therefore, in many cases, the default integration
   rule order is not sufficient and higher integration must be used on
   elements (consult element library manual and related NIP parameter).

-  EET_CZZSI - The implementation of combined criteria: Zienkiewicz Zhu
   Error Estimator for elastic regime and scalar error indicator in
   non-linear regime.

.. _eetypestable:

.. table:: Supported error estimators and indicators.

   ========================= ==========
   Error estimator/indicator ``eetype``
   ========================= ==========
   EET_SEI                   0
   EET_ZZEE                  1
   EET_CZZSI                 2
   ========================= ==========

The sets of parameters (``errorestimatorparams`` field) used to
configure each error estimator are different

-  | ``EET_SEI``
   | [``regionskipmap #(ia)``]  ``vartype #(in)``
     ``minlim #(rn)`` ``maxlim #(rn)`` ``mindens #(rn)``
     ``maxdens #(rn)`` ``defdens #(rn)``
    [``remeshingdensityratio #(rn)``]

   -  ``regionskipmap`` parameter allows to skip some regions. The error
      is not evaluated in these regions and default mesh density is
      used. The size of this array should be equal to number of regions
      and nonzero entry indicates region to skip.

   -  ``vartype`` parameter determines the type of internal variable to
      be used as error indicator. Currently supported value is 1,
      representing damage based indicator.

   -  If the indicator value is in range given by parameters
      (``minlim``, ``maxlim``) then the proposed mesh density is
      linearly interpolated within range given by parameters
      (``mindens``, ``maxdens``). If indicator value is less than value
      of ``minlim`` parameter then value of ``defdens`` parameter is
      used as required density, if it is larger than ``maxlim`` then
      ``maxdens`` is used as required density.

   -  ``remeshingdensityratio`` parameter determines the allowed ratio
      between proposed density and actual density. The remeshing is
      forced, whenever the actual ratio is smaller than this value.
      Default value is equal to 0.80.

-  | ``EET_ZZEE``
   | [``regionskipmap #(ia)``]  ``normtype #(in)``
     ``requirederror #(rn)`` ``minelemsize #(rn)``

   -  ``regionskipmap`` parameter allows to skip some regions. The error
      is not evaluated in these regions and default mesh density is
      used. The size of this array should be equal to number of regions
      and nonzero entry indicates region to skip.

   -  ``normtype`` Allows select the type of norm used in evaluation of
      error. Default value is to use L2 norm (equal to 0), value equal
      to 1 uses the energy norm.

   -  ``requirederror`` parameter determines the required error to
      obtain (in percents/100).

   -  minelemsize parameter allows to set minimum limit on element size.

-  EET_CZZSI - combination of parameters for EET_SEI and EET_ZZEE; the
   in elastic regions are driven using EET_SEI, the elastic are driven
   by EET_ZZEE.

.. _materialinterfaces:

Material interfaces
-------------------

The material interfaces are used to represent and track the position of
various interfaces on fixed grids. Typical examples include free
surface, evolving interface between two materials, etc. Available
representations include:

======== ====== =============
MI       miflag Compatibility
======== ====== =============
LEPlic   0      2D triangular
LevelSet 1      2D triangular
======== ====== =============

-  | LEPlic- representation based on Volume-Of-Fluid approach; the
     initial distribution of VOF fractions should be specified for each
     element (see element manual)
   | [``refvol #(rn)``]

   -  parameter ``refvol`` allows to set initial volume of reference
      fluid, then the reference volume is computed in each step and
      printed, so the accuracy and mass conservation can be monitored.

-  | LevelSet- level set based representation
   | ``levelset #(ra)`` OR ``refmatpolyx #(ra)`` ``refmatpolyy #(ra)``
   | [``lsra #(in)``] [``rdt #(rn)``] [``rerr #(rn)``]

   -  ``levelset`` allows to specify the initial level set values for
      all nodes directly. The size should be equal to total number of
      nodes within the domain.

   -  Parameters ``refmatpolyx`` and ``refmatpolyy`` allow to initialize
      level set by specifying interface geometry as 2d polygon. Then
      polygon describes the initial zero level set, and level set values
      are then defined as signed distance from this polygon. Positive
      values are on the left side when walking along polygon. The
      parameter ``refmatpolyx`` specifies the x-coordinates of polygon
      vertices, parameter ``refmatpolyy`` y-corrdinates. Please note,
      that level set must be initialized, either using ``levelset``
      parameter or using ``refmatpolyx`` and ``refmatpolyy``.

   -  Parameter ``lsra`` allows to select level set reinitialization
      algorithm. Currently supported values are 0 (no
      re-initialization), 1 (re-initializes the level set representation
      by solving
      :math:`d_{\tau} = S(\phi)(1-\vert\boldsymbol{\nabla}d\vert)` to
      steady state, default), 2 (uses fast marching method to build
      signed distance level set representation).

   -  Parameters ``rdt`` ``rerr`` are used to control reinitialization
      algorithm for ``lsra`` = 0. ``rdt`` allows to change time step of
      integration algorithm and parameter ``rerr`` allows to change
      default error limit used to detect steady state.

.. _meshpackages:

Mesh generator interfaces
-------------------------

The mesh generator interface is responsible to provide a link to
specific mesh generator. The supported values of ``meshpackage``
parameter are

-  MPT_T3D: ``meshpackage`` = 0. T3d mesh interface. Default. Supports
   both 1d, 2d (triangles) and 3d (tetrahedras) meshes. Reliable.

-  MPT_TARGE2: ``meshpackage`` = 1. Interface to Targe2 2D mesh
   generator.

-  MPT_SUBDIVISION: ``meshpackage``\ =3. Built-in subdivision algorithm.
   Supports triangular 2D and tetrahedral 3D meshes. Can operate in
   parallel mode.

.. _InitModulesSec:

Initialization modules
----------------------

| Initialization modules allow to initialize the state variables using
  data previously computed by external software. The number of
  initialization module records is specified in analysis record using
  ``ninitmodules`` parameter (see the initial part of section
  :ref:`AnalysisRecord`. The general format is the following:

| ``EntType``   ``initfile #(string)``
  The file name following the keyword
  “initfile” specifies the path to the file that contains the
  initialization data and should be given without quotes.

Currently, the only supported initialization module is

-  Gauss point initialization module

   | ``GPInitModule`` ``initfile #(string)``
   |
   | -  Each Gauss point is represented by one line in the initialization file.
   | -  The Gauss points should be given in a specific order, based on the
        element number and the Gauss point number, in agreement with the
        mesh specified in later sections.

   | -  Each line referring to a Gauss point should contain the following data:
      |  ``elnum #(in)`` ``gpnum #(in)`` ``coords #(ra)``
      | ``ng #(in)`` ``var_1_id #(in)`` ``values_1 #(ra)`` ``...`` ``var_ng_id #(in)`` ``values_ng #(ra)``

       -  ``elnum`` is the element number

       -  ``gpnum`` is the Gauss point number

       -  ``coords`` are the coordinates of the Gauss point

       -  ``ng`` is the number of groups of variables that will follow

       -  ``var_1_id`` is the identification number of variable group number
           1 (according to the definitions in internalstatetype.h)

       -  ``values_1`` are the values of variables in group number 1

       -  ``var_ng_id`` is the identification number of variable group number ng

       -  ``values_ng``  are the values of variables in group number ng

   | -  Example:
      |   ``37 4 3 0.02 0.04 0.05 3 52 1 0.23 62 1 0.049 1 6 0 -2.08e+07 0 0 0 0``
          means that Gauss point number 4 of element number 37 has
         coordinates :math:`x=0.02`, :math:`y=0.04` and :math:`z=0.05`
         and the initial values are specified for 3 groups of variables;
         the first group (variable ID 52) is of type IST_DamageScalar
         (see internalstatetype.h) and contains 1 variable (since it is a
         scalar) with value 0.23;
         the second group (ID 62) is of type IST_CumPlasticStrain and
         contains 1 variable with value 0.049;
         the third group is of type IST_StressTensor and contains 6
         variables (stress components :math:`\sigma_x`, :math:`\sigma_y`,
         etc.) with values 0, -2.08e+07, 0, 0, 0, 0

.. _ExportModulesSec:

Export modules
--------------

| Export modules allow to export computed data into external software
  for post-processing. The number of export module records is specified
  in analysis record using ``nmodules`` parameter (see the initial part
  of section AnalysisRecord_. The general format is the
  following:
| ``EntType`` [``tstep_all``] [``tstep_step #(in)``] [``tsteps_out #(rl)``]
  [``subtsteps_out #(in)``] [``domain_all``]
  [``domain_mask #(in)``] [``regionsets #(ia)``]
  [``timeScale #(rn)``]

| To select all solution steps, in which output will be performed, use
   ``tstep_all``. To select each ``tstep_step``-nth step, use
   ``tstep_step`` parameter. In order to select only specific solution
   steps, the ``tsteps_out`` list can be specified, supplying solution step
   number list in which output will be done. To select output for all
   domain of the problem the ``domain_all`` keyword can be used. To select
   only specific domains, ``domain_mask`` array can be used, where the
   values of the array specify the domain numbers to be exported. If the
   parameter ``subtsteps_out`` = 1, it turns on the export of intermediate
   results, for example during the substepping or individual equilibrium
   iterations. This option requires support from the solver.

The export is done on region basis, on each region, the nodal
   recovery is performed independently and results are exported in a
   separate piece. This allows to take into account for
   discntinuities, or to export variables defined only by particular
   material model. The region volumes are defined using sets
   containing individual elements. By default the one region is
   created, containing all element in the problem domain. The
   optional parameter ``regionsets`` allows to use user-defined. The
   individual values refer to numbers (ids) of domain sets. Note,
   that regions are determined solely using elements.

   vtkxml tstep_all cellvars 1 46 vars 1 1 primvars 1 1 stype 2
   regionsets 2 1 2

Optional parameter ``timeScale`` scales time in output. In transport problem, basic
   units are seconds. Setting timeScale = 2.777777e-4 (=1/3600.)
   converts all time data in vtkXML from seconds to hours.


Currently, the supported export modules are following

-  VTK export, **DEPRECATED - Use VTKXML or vtkhdf5**

   ``vtk`` [``vars #(ia)``] [``primvars #(ia)``] [``cellvars #(ia)``]
   [``stype #(in)``] [``regionstoskip #(ia)``]

   ``vtkxml`` [``vars #(ia)``] [``primvars #(ia)``] [``cellvars #(ia)``]
   [``ipvars #(ia)``] [``stype #(in)``] [``setmembership``]

   | <``ver 1.6``> 
   | ``vtkhdf5`` [``vars #(ia)``] [``primvars #(ia)``] [``cellvars #(ia)``]
     [``ipvars #(ia)``] [``stype #(in)``] 

   -  The vtk module is obsolete, use vtkxml or vtkhdf5 instead. Vtkxml allows to
      export results recovered on region by region basis and has more
      features. The vtkxml export module exports result into vtu files, one file for each solution step (vtk ustructured grids).
      The vtkhdf5 export module exports results into vtk hdf5 file (single hdf file fro all time steps). 
      Both output formats can be visualized using Paraview.

      The vtkhdf5 export module requires HDF5 library and oofem has to be configured with USE_HDF5=ON. 
      The HDF5 library version 1.14 (or later) is recommended.


   -  The array ``vars`` contains identifiers for those internal
      variables which are to be exported. These variables will be
      smoothed and transfered to nodes. The id values are defined by
      InternalStateType enumeration, which is defined in include file
      “src/core/internalstatetype.h”.

   -  The array ``primvars`` contains identifiers of primary variables
      to be exported. The possible values correspond to the values of
      enumerated type UnknownType, which is again defined in
      “src/core/unknowntype.h”. Please note, that the values
      corresponding to enumerated type values start from zero, if not
      specified directly and that not all values are supported by
      particular material model or analysis type.

   -  The array ``cellvars`` contains identifiers of constant variables
      defined on an element (cell), e.g. a material number. Identifier
      numbers are specified in “src/core/internalstatetype.h”.

   -  The array ``ipvars`` contains identifiers for those internal
      variables which are to be exported. These variables will be
      directly exported (no smoothing) as point dataset, where each
      point corresponds to individual integration point. A separate vtu
      file for these raw, point data will be created. The id values are
      defined by InternalStateType enumeration, which is defined in
      include file “src/core/internalstatetype.h”.

   -  The parameter ``stype`` allows to select smoothing procedure for
      internal variables, which is used to compute nodal values from
      values in integration points. The supported values are :math:`0`
      for simple nodal averageing (generally supported only by
      triangular and tetrahedral elements), :math:`1` for Zienkiewicz
      Zhu recovery (default), and :math:`2` for Superconvergent Patch
      Recovery (SPR, based on least square fitting).
   -  The paramneter ``setmembership`` will trigger export of set
      membership information. The set membership information (at present for DofManagers and cells) is
      exported as a point/cell variable. As vtk format does not support sparse sets, 
      the set membership is encoded into datasets containing byte values (UINT8, named VertexSetMembership and CellSetMembership), 
      where i-th bit is set to 1 if the given component is member of set, zero otherwise. 

   
-  VTK pfem (particle FEM) export. Exports particle positions to vtk as a point dataset.

   ``vtkpfem`` [``vars #(ia)``] [``primvars #(ia)``] [``cellvars #(ia)``]
   [``ipvars #(ia)``] [``stype #(in)``] 

-  VTK memory export. This module is not producing any output, but prepares necessary data structures to suport vtk export or vtk visualization. It is used by Python interface to access vtk datasets.   
  
   ``vtkmemory`` [``vars #(ia)``] [``primvars #(ia)``] [``cellvars #(ia)``] [``ipvars #(ia)``] 

-  VTK xfem export module. Exports xfem related data. The data exported are determined
   by XfemManager vtkExportFields parameter (see ``exportfields`` keyword).

   ``vtkxmlxfem``   

-  Homogenization of IP quantities in the global coordinate system (such
   as stress, strain, damage, heat flow). Corresponding IP quantities
   are summed and averaged over the volume. It is possible to select
   region sets from which the averaging occurs. The averaging works for
   all domains with an extension to trusses. A truss is considered as a
   volume element with oriented stress and strain components along the
   truss axis. The transformation to global components occurs before
   averaging.

   ``hom`` ``ists #(ia)``   [``scale #(rn)``] [``regionSets #(ia)``]
   [``strain_energy``]

   -  An integer array ``ists`` specifies internal state types for
      export which are defined in internalstatetype.h file.

   -  The parameter ``scale`` multiplies all averaged IP quantities.
      ``scale``\ =1 by default.

   -  An integer array ``regionSets`` specifies region sets for
      averaging. All domain is averaged by default.

   -  ``strain_energy`` calculates strain energy over selected elements
      (defined by sets) by

      .. math:: W^*=\int_V \int \sigma \mathrm{d} (\varepsilon-\varepsilon_{eig}) \mathrm{d} V

      where :math:`\sigma` is the stress tensor, :math:`\varepsilon`
      stands for the strain tensor and :math:`\varepsilon_{eig}` is
      eigenstrain tensor (originates from temperature load or prescribed
      eigenstrain). Strain energy increment and total strain energy is
      reported in each step. The integration uses mid-point rule for
      stress and yields exact results for linear elastic materials.

-  Gauss point export is useful if one needs to plot a certain variable
   (such as damage) as a function of a spatial coordinate using tools
   like gnuplot. It generates files with data organized in columns, each
   row representing one Gauss point. In this way, one can plot e.g. the
   damage distribution along a one-dimensional bar.
   ``gpexportmodule`` [``vars #(ia)``] [``ncoords #(in)``]

   -  The array ``vars`` contains identifiers for those internal
      variables which are to be exported. The id values are defined by
      InternalStateType enumeration, which is defined in include file
      “src/core/internalstatetype.h”.

   -  Parameter ``ncoords`` specifies the number of spatial coordinates
      to be exported at each Gauss point. Depending on the spatial
      dimension of the domain, the points can have one, two or three
      coordinates. If ``ncoords`` is set to -1, only those coordinates
      that are actually used are exported. If ``ncoords`` is set to 0,
      no coordinates are exported. If ``ncoords`` is set to a positive
      integer, exactly ``ncoords`` coordinates are exported. If
      ``ncoords`` exceeds the actual number of coordinates, the actual
      coordinates are supplemented by zeros. For instance, if we deal
      with a 2D problem, the actual number of coordinates is 2. For
      ``ncoords``\ =3, the two actual coordinates followed by 0 will be
      exported. For ``ncoords``\ =1, only the first coordinate will be
      exported.

   The Gauss point export module creates a file with extension “gp”
   after each step for which the output is performed. This file contains
   a header with lines starting by the symbol #, followed by the actual
   data section. Each data line corresponds to one Gauss point and
   contains the following data:

   #. element number,

   #. material number,

   #. Gauss point number,

   #. contributing volume around Gauss point,

   #. Gauss point global coordinates (written as a real array of length
      ``ncoords``),

   #. internal variables according to the specification in ``vars``
      (each written as a real array of the corresponding length).

   | Example:
     ``GPExportModule 1 tstep_step 100 domain_all ncoords 2 vars 5 4 13 31 64 65``
      
     means that the \*.gp file will be written after each 100 steps and
     will contain for each of the Gauss points in the entire domain its
     2 coordinates and also internal variables of type 4, 13, 31, 64 and
     65, which are the strain tensor, damage tensor, maximum equivalent
     strain level, stress work density and dissipated work density. Of
     course, the material model must be able to deliver such variables.
     The size of the strain tensor depends on the spatial dimension, and
     the size of the damage tensor depends on the spatial dimension and
     type of model (e.g., for a simple isotropic damage model it will
     have just 1 component while for an anisotropic damage model it may
     have more). The other variables in this example are scalars, but
     they will be written as arrays of length 1, so the actual value
     will always be preceded by “1” as the length of the array. Since
     certain internal variables have the meaning of densities (per unit
     volume or area, again depending on the spatial dimension), it is
     useful to have access to the contributing volume of the Gauss
     point. The product of this contributing volume and the density
     gives an additive contribution to the total value of the
     corresponding variable. This can be exploited e.g. to evaluate the
     total dissipated energy over the entire domain.

- Solution status monitor <``ver 1.6``>
   | The Solution status monitor creates a report on problem solution. Its content is configurable, see below. 
   | The syntax of solutionstatus export module is following:
   | ``solutionstatus`` [``fmt #(s)``] 
   | 
   | The fmt parameter allows to control what will be reported for each solution step. The fmt string contains data codes separated by colon. 
   | The following data codes are supported:

+-----------+--------------------------------------------------------------------------+
| Data code | Description                                                              |
+===========+==========================================================================+
|m          | Meta step number                                                         |
+-----------+--------------------------------------------------------------------------+
|s          | Solution step number                                                     |
+-----------+--------------------------------------------------------------------------+
|a          | Attempt number                                                           |
+-----------+--------------------------------------------------------------------------+
|nite       | Number of iterations                                                     |
+-----------+--------------------------------------------------------------------------+
|t          | Solution step target time                                                |
+-----------+--------------------------------------------------------------------------+
|dt         | Solution step time increment                                             |
+-----------+--------------------------------------------------------------------------+
|st         | Time spend on solving solution step [s]                                  |
+-----------+--------------------------------------------------------------------------+
|cr         | | Convergence reason status                                              |
|           | | - *Converged*                                                          |
|           | | - *Diverged_I* (Convergence not reached with maximum iteration limit)  |
|           | | - *Diverged_T* (Diverged solution)                                     |
|           | | - *Failed*                                                             |
+-----------+--------------------------------------------------------------------------+
|\-         | Placeholder, will print ‘-‘                                              |
+-----------+--------------------------------------------------------------------------+

   | The default fmt string is "m:s:a:nite:t:dt:st:cr"
   | Example of input record:
   | ``solutionstatus`` ``tstep_all`` ``fmt m:s:a:nite:t:dt:st:cr``

- Error checking module
   Error checking rule used for regressions tests. It compares the user-selected computed results with the reference results.
   The individual rules are defined in a separate section with dedicated syntax, that can be part of the input file itself, or defined in external file.
   Their syntax follows the syntax of Extractor package, namely its `Value records section <https://www.oofem.org/resources/doc/extractorInput/html/node2.html>`_.
   
   The syntax of the error checking module is following:

   ``errorcheck`` [``filename #(s)``] [``extract``]
   
   By default the module performs the regression tests and issues the error when one or more tests fail.
   The module can also be used to extract user-selected quantities rather than testing them against expected values. This can be done by adding ``extract`` keyword to the input record of the module.
   The ``filename`` parameter allows to specify filename (with path) in which the rules are defined, insted of input file itself (the default).

.. _VariablesSec:

Variables (symbolic mpm module)
-------------------------------
The symbolic mpm module allows to define problem(s) by defining problem weak form using variables, terms and integrals. 
In this section we describe how to define variable or field apperaing in weak form. The general syntax is following:

``Variable`` ``name #(s)`` ``interpolation #(s)`` ``type #(in)`` ``quantity #(in)`` ``size #(in)`` ``dofs #(ia)``

where the parameters have following meaning:

- ``name`` is the string containing the name of the variable
- ``interpolation`` string, defining the interpolation of the varaible. The supported values are:
    
  - ``feiconst`` - constant interpolation
  - ``feilin`` - linear interpolation
  - ``feiquad`` - quadratic interpolation
- ``type`` defines the rank of the variable. The supported values are:
  
  - 0 - for scalar variable
  - 1 - for vector variable
- ``quantity`` attribute defines the physical meaning of variable. Supported values include
  
   - 0 - for displacement field
   - 1 - for velocity field
   - 2 - for temperature field
   - 3 - for pressure field
- ``size`` attribute determines the size (dimension) of variable.
- ``dofs`` array of integers, defining the physical meaning of variable DOFs. The size of the array should be equal to the size of the variable. 
    The supported values are defined in src/core/dofiditem.h file.

.. _TermsSec:

Terms (symbolic mpm module)
---------------------------
The integrals in the weak form integrate terms. The individual term record have the following generic syntax:

``TermType``  ``variable #(s)``  ``testvariable #(s)`` ``mmode #(in)``

The Supported TermType keywords are documented below.
The parameters have following meaning:

  - ``variable`` is the name of the unknown variable (field) of the term
  - ``testvariable`` is the name of the test variable (field) of the term
  - ``mmode`` allows to define material mode used to evaluate the term

Note that sopecific terms can introduce additional parameters to define the term.

Supported TermTypes
^^^^^^^^^^^^^^^^^^^^
+--------------+-----------------------------------------------------------------------------------------------------------------------------+
| Keyword      | Description & parameters                                                                                                    |
+==============+=============================================================================================================================+
|| BTSigmaTerm || :math:`\int_\Omega \nabla^s \mathbf{w}\ \mathbf{\sigma}(\nabla^s \mathbf{u})`,                                             |
||             || where :math:`\mathbf{\sigma}` is (nonlinar) operator evaluated by constitutive model.                                      |
||             || Supported material modes (``mmode``): _3dMat, _3dUP, _2dUP, _PlaneStress                                                   |
||             || Optional parameters                                                                                                        |
||             || * ``lhsmatmode #(in)``                                                                                                     |
+--------------+-----------------------------------------------------------------------------------------------------------------------------+
|| BTamNTerm   || :math:`\int_\Omega \nabla^s \mathbf{w}\ a\mathbf{m}\ p`,                                                                   |
||             || where ``a`` is material parameter, defined by contitutive model with meaning defined by ``atype`` parameter                |
||             || and :math:`\mathbf{m}^T=[1,1,1,0,0,0]^T`. Note that :math:`\mathbf{\sigma}_{D}=\mathbf{\sigma}-\mathbf{m}p`                |
+--------------+-----------------------------------------------------------------------------------------------------------------------------+
| NTamTBTerm   | This is transposed version of BTamNTerm :math:`=\int_\Omega w\ \alpha\mathbf{m}\ \nabla^s \mathbf{u}`                       |
+--------------+-----------------------------------------------------------------------------------------------------------------------------+
| NTcN         | :math:`\int_\Omega q c p`, where ``c`` is constant defined by material model with meaning determined by `ctype` parameter.  |
+--------------+-----------------------------------------------------------------------------------------------------------------------------+
| NTfTerm      | :math:`\int_\Omega \mathbf{w}\cdot\bar{\mathbf{t}}`, where ``t`` is given flux vector, defined by ``flux #(ra)`` parameter. |
+--------------+-----------------------------------------------------------------------------------------------------------------------------+

Notation
^^^^^^^^
+----------------------+-------------------------+
| Symbol               | Description             |
+======================+=========================+
| q,p                  | scalar valued functions |
+----------------------+-------------------------+
| :math:`\mathbf{w,u}` | vector valued functions |
+----------------------+-------------------------+





.. _IntergarlsSec:

Integrals (symbolic mpm module)
--------------------------------
The Integrals in the weak form integrate individual terms. The individual integral record have the following generic syntax:

``Integral #(in)``  ``domain #(in)`` ``set #(in)`` ``term #(in)``

The parameters have following meaning:
  - The record keyword is always ``Integral``, followed by its number,
  - ``domain`` determines the domain in which integral is defined (by providing domain number),
  - ``set`` is the set number, which defines integral domain using set of elements over which the integration is performed. Note that boundary integrals require the presence of boundary elements. 
  - ``term`` defines term to integrate by providing its number.



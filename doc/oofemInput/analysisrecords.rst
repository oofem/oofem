
.. _AnalysisRecord:

Analysis record
===============

This record describes the type of analysis, which should be performed.
The analysis record can be splitted into optional meta-step input
records (see below). Then certain attributes originally in analysis
record can be specified independently for each meta-step. This is marked
by adding “M” superscript to keyword. Then the attribute format is
``Keyword^M`` #(type).

The general format of this record can be specified using

-  | “standard-syntax”
   | ``nsteps #(in)`` [``renumber #(in)``]
     [``profileopt #(in)``] ``attributes #(string)``
     [``ninitmodules #(in)``] [``nmodules #(in)``]
     [``nxfemman #(in)``]

-  | “meta step-syntax”
   | ``nmsteps #(in)`` [``ninitmodules #(in)``]
     [``nmodules #(in)``] [``nxfemman #(in)``]
   | immediately followed by ``nmsteps`` meta step records with the
     following syntax:
   | ``nsteps #(in)`` ``attributes #(string)``
   | The ``nmsteps`` parameter determines the number of “metasteps”. The
     meta step represent sequence of solution steps with common
     attributes. There is expected to be ``nmsteps`` subsequent metastep
     records. The meaning of meta step record parameters (or analysis
     record parameters in “standard syntax”) is following:

   -  ``nsteps`` - determines number of subsequent solution steps within
      metastep.

   -  ``renumber`` - Turns out renumbering after each time step.
      Necessary when Dirichlet boundary conditions change during
      simulation. Can also be turned out by the executeable flag
      ``-rn``.

   -  ``profileopt`` - Nonzero value turns on the equation renumbering
      to optimize the profile of characteristic matrix (uses Sloan
      algorithm). By default, profile optimization is not performed. It
      will not work in parallel mode.

   -  ``attributes`` - contains the metastep related attributes of
      analysis (and solver), which are valid for corresponding solution
      steps within meta step. If used in standard syntax, the attributes
      are valid for all solution step.

   -  ``ninitmodules`` - number of initialization module records for
      given problem. The initialization modules are specified after meta
      step section (or after analysis record, if no metasteps are
      present). Initialization modules allow to initialize the state
      variables by values previously computed by external software. The
      available initialization modules are described in section
      :ref:`InitModulesSec`.

   -  ``nmodules`` - number of export module records for given problem.
      The export modules are specified after initialization modules.
      Export modules allow to export computed data into external
      software for postprocessing. The available export modules are
      described in section :ref:`ExportModulesSec`.

   -  ``nxfemman`` - 1 implies that an XFEM manager is created, 0
      implies that no XFEM manager is created. The XFEM manager stores a
      list of enrichment items. The syntax of the XFEM manager record
      and related records is described in section
      :ref:`XFEMManagerRecords`.

   -  ``eetype`` - optional error estimator type used for the problem.
      Used for adaptive analysis, but can also be used to compute and
      write error estimates to the output files. See adaptive
      engineering models for details.

Not all of analysis types support the metastep syntax, and if not
mentioned, the standard-syntax is expected. Currently, supported
analysis types are

-  Linear static analysis, see section :ref:`LinearStatic`,

-  Eigen value dynamic, see section :ref:`EigenValueDynamic`,

-  Direct explicit nonlinear dynamics, see section
   :ref:`NlDEIDynamic`,

-  Direct explicit (linear) dynamics, see section
   :ref:`DEIDynamic`,

-  Implicit linear dynamic, see section :ref:`DIIDynamic`,

-  Incremental **linear** static problem, see section
   :ref:`IncrementalLinearStatic`,

-  Non-linear static analysis, see section :ref:`NonLinearStatic`.

-  Dymmy problem, see section :ref:`DummyEngngModel` 

Structural Problems
-------------------

.. _StaticStructural:

StaticStructural
~~~~~~~~~~~~~~~~

``StaticStructural`` ``nsteps #(in)`` [``deltat #(...)``] [``prescribedtimes #(...)``] [``stiffmode #(...)``] [``nonlocalext #(...)``] [``sparselinsolverparams #(...)``]

Static structural analysis. Can be used to solve linear and nonlinear
static structural problems, supporting changes in boundary conditions
(applied load and supports). The problem can be solved under direct load
or displacement control, indirect control, or by their arbitrary
combination. Note, that the individual solution steps are used to
describe the history of applied incremental loading. The load cases are
not supported, for each load case the new analysis has to be performed.
To analyze linear static problem with multiple load combinations, please
use LinearStatic solver.

By default all material nonlinearities will be taken into account,
geometrical not. To include geometrically nonlinear effect one must
specify level of non-linearity in element records.

The ``sparselinsolverparams`` parameter describes the sparse linear
solver attributes and is explained in section
:ref:`sparselinsolver`. The optional parameter ``deltat`` defines
the length of time step (equal to 1.0 by default). The times
corresponding to individual solution times can be specified using
optional parameter ``prescribedtimes``, allowing to input array of
discrete solution times, the number of solution steps is then equal to
the size of this array. .

.. _LinearStatic:

Linear static analysis
~~~~~~~~~~~~~~~~~~~~~~

``LinearStatics`` ``nsteps #(in)`` [``sparselinsolverparams #(...)``] [``sparselinsolverparams #(...)``]

Linear static analysis. Parameter ``nsteps`` indicates the number of
loading cases. Problem supports multiple load cases, where number of
load cases correspods to number of solution steps, individual load
vectors are formed in individual time-steps. However, the static system
is assumed to be the same for all load cases. For each load case an
auxiliary time-step is generated with time equal to load case number.

The ``sparselinsolverparams`` parameter describes the sparse linear
solver attributes and is explained in section
:ref:`sparselinsolver`.

.. _LinearStability:

LinearStability
~~~~~~~~~~~~~~~

``LinearStability`` ``nroot #(in)`` ``rtolv #(rn)`` [``eigensolverparams #(...)``]

Solves linear stability problem. Only
first ``nroot`` smallest eigenvalues and corresponding eigenvectors will
be computed. Relative convergence tolerance is specified using ``rtolv``
parameter.

The ``eigensolverparams`` parameter describes the sparse linear solver
attributes and is explained in section :ref:`eigensolverssection`.

.. _EigenValueDynamic:

EigenValueDynamic
~~~~~~~~~~~~~~~~~

``EigenValueDynamic`` ``nroot #(in)`` ``rtolv #(rn)`` [``eigensolverparams #(...)``]

Represents the eigen value dynamic analysis. Only ``nroot`` smallest
eigenvalues and corresponding eigenvectors will be computed. Relative
convergence criteria is governed using ``rtolv`` parameter.

The ``eigensolverparams`` parameter describes the sparse linear solver
attributes and is explained in section :ref:`eigensolverssection`.

.. _NlDEIDynamic:

NlDEIDynamic
~~~~~~~~~~~~

``NlDEIDynamic`` ``nsteps #(in)`` ``dumpcoef #(rn)`` [``deltaT #(rn)``]

Represents the direct explicit nonlinear dynamic integration. The
central difference method with diagonal mass matrix is used, damping
matrix is assumed to be proportional to mass matrix,
:math:`\boldsymbol{C}
= \mathrm{dumpcoef} * \boldsymbol{M}`, where :math:`\boldsymbol{M}` is
diagonal mass matrix. Parameter ``nsteps`` specifies how many time steps
will be analyzed. ``deltaT`` is time step length used for integration,
which may be reduced by program in order to satisfy solution stability
conditions. Parameter ``reduct`` is a scaling factor (smaller than 1),
which is multiplied with the determined step length adjusted by the
program. If ``deltaT`` is reduced internally, then ``nsteps`` is
adjusted so that the total analysis time remains the same.

| The parallel version has the following additional syntax:
| &\ :math:`\langle`\ [``nonlocalext``]\ :math:`\rangle`\ &

.. _DEIDynamic:

DEIDynamic
~~~~~~~~~~

``DEIDynamic`` ``nsteps #(in)`` ``dumpcoef #(rn)`` [``deltaT #(rn)``]

Represent the **linear** explicit integration scheme for dynamic problem
solution. The central difference method with diagonal mass matrix is
used, damping matrix is assumed to be proportional to mass matrix,
:math:`\boldsymbol{C} = \mathrm{dumpcoef} * \boldsymbol{M}`, where
:math:`\boldsymbol{M}` is diagonal mass matrix. ``deltaT`` is time step
length used for integration, which may be reduced by program in order to
satisfy solution stability conditions. Parameter ``nsteps`` specifies
how many time steps will be analyzed.

.. _DIIDynamic:

DIIDynamic
~~~~~~~~~~

``DIIDynamic`` ``nsteps #(in)`` ``deltaT #(rn)`` [``ddtscheme #(in)``] [``gamma #(rn)``] [``beta #(rn)``] [``eta #(rn)``] [``delta #(rn)``] [``theta #(rn)``]

Represents direct implicit integration of linear dynamic problems. Solution procedure described in Solution procedure described in 
K. Subbaraj and M. A. Dokainish, A SURVEY OF DIRECT TIME-INTEGRATION METHODS IN COMPUTATIONAL STRUCTURAL DYNAMICS - II. IMPLICIT METHODS,
Computers & Structures Vol. 32. No. 6. pp. 1387-1401, 1989.

Parameter ``ddtscheme`` determines integration scheme, as defined in src/oofemlib/timediscretizationtype.h (TD_ThreePointBackward=0 (default), TD_TwoPointBackward =  1,
TD_Newmark =  2, TD_Wilson =  3, TD_Explicit  =  4).

Parameters ``beta`` and ``gamma`` determine the stability and acuracy of the integration algorithm, both have zero values as default. For ``gamma=0.5`` and ``beta = l/6``, the linear acceleration method is obtained. Unconditional stability is obtained, when :math:`2\beta \ge \gamma \ge 1/2`. 
The dafault values are ``beta=0.25`` and ``gamma=0.5``. The Wilson-theta metod requires additional ``theta`` parameter with default value equal to 1.37.
The damping is assumed to be modeled as Rayleigh damping :math:`\boldsymbol{C} = \eta \boldsymbol{M} + \delta \boldsymbol{K}`.

.. _IncrementalLinearStatic:

IncrementalLinearStatic
~~~~~~~~~~~~~~~~~~~~~~~

``IncrementalLinearStatic`` ``endOfTimeOfInterest #(rn)`` ``prescribedTimes #(ra)``

Represents incremental **linear** static problem. The problem is solved
as series of linear solutions and is intended to be used for solving
linear creep problems or incremental perfect plasticity.

Supports the changes of static scheme (applying, removing and changing
boundary conditions) during the analysis.

Response is computed in times defined by ``prescribedTimes`` array.
These times should include times, when generally the boundary conditions
are changing, and in other times of interest. (For linear creep
analysis, the values should be uniformly distributed on log-time scale,
if no change in loading or boundary conditions). The time at the end of
interested is specified using ``endOfTimeOfInterest`` parameter.

.. _NonLinearStatic:

NonLinearStatic
~~~~~~~~~~~~~~~

| **NonLinearStatic**
| Non-linear static analysis. The problem can be solved under direct
  load or displacement control, indirect control, or by their arbitrary
  combination. By default all material nonlinearities will be included,
  geometrical not. To include geometrically nonlinear effect one must
  specify level of non-linearity in element records. There are two
  different ways, how to specify the parameters - the extended and
  standard syntax.

Extended syntax
^^^^^^^^^^^^^^^

The extended syntax uses the “metastep” concept and has the following
format:

``NonLinearStatic`` [``nmsteps #(in)``] ``nsteps #(in)`` [``contextOutputStep #(in)``] [``sparselinsolverparams #(string)``] [``nonlinform #(in)``] <[``nonlocstiff #(in)``]>
<[``nonlocalext``]> <[``loadbalancing``]>

This record
is immediately followed by metastep records with the format described
below. The analysis parameters have following meaning

-  ``nmsteps`` - determines the number of “metasteps”, default is 1.

-  ``nsteps`` - determines number of solution steps.

-  ``contextOutputStep`` - causes the context file to be created for
   every contextOutputStep-th step and when needed. Useful for
   postprocessing.

-  The ``sparselinsolverparams`` parameter describes the sparse linear
   solver attributes and is explained in section
   :ref:`sparselinsolver`.

-  ``nonlinform`` - formulation of non-linear problem. If == 1
   (default), total Lagrangian formulation in undeformed original shape
   is used (first-order theory). If == 2, the equlibrated displacements
   are added to original ones and updated in each time step
   (second-order theory).

-

| The metastep record has the following general syntax:
| ``nsteps #(in)`` [``controlmode #(in)``] [``deltat #(rn)``]
  [``stiffmode #(in)``] [``refloadmode #(in)``]
  ``solverParams #()`` [``sparselinsolverparams #(string)``]
  [``donotfixload #()``]

where

-  ``controlmode`` - determines the type of solution control used for
   corresponding meta step. if == 0 then indirect control will be used
   to control solution process (arc-length method, default). if == 1
   then direct displacement or load control will be used (Newton-Raphson
   solver). In the later mode, one can apply the prescribed load
   increments as well as control displacements.

-  ``deltaT`` - is time step length. If not specified, it is set equal
   to 1,0. Each solution step has associated the corresponding intrinsic
   time, at which the loading is generated. The ``deltaT`` determines
   the spacing between solution steps on time scale.

-  ``stiffMode`` - If == 0 (default) then tangent stiffness will be used
   at new step beginning and whenever numerical method will ask for
   stiffness update. If == 1 the use of secant tangent will be forced.
   The secant stiffness will be used at new step beginning and whenever
   numerical method will ask for stiffness update. If == 2 then original
   elastic stiffness will be used during the whole solution process.

-  The ``refloadmode`` parameter determines how the reference force load
   vector is obtained from given totalLoadVector and initialLoadVector.
   The initialLoadVector describes the part of loading which does not
   scale. Works only for force loading, other non-force components
   (temperature, prescribed displacements should always given in total
   values). If ``refloadmode`` is 0 (rlm_total, default) then the
   reference incremental load vector is defined as totalLoadVector
   assembled at given time. If ``refloadmode`` is 1 (rlm_inceremental)
   then the reference load vector is obtained as incremental load vector
   at given time.

-  ``solverParams`` - parameters of solver. The solver type is
   determined using ``controlmode``.

-  The ``sparselinsolverparams`` parameter describes the sparse linear
   solver attributes and is explained in section
   :ref:`sparselinsolver`.

-  By default, reached load at the end of metastep will be maintained in
   subsequent steps as fixed, non scaling load and load level will be
   reset to zero. This can be changed using keyword ``donotfixload``,
   which if present, causes the loading to continue, not resetting the
   load level. For the indirect control the reached loading will not be
   fixed, however, the new reference loading vector will be assembled
   for the new metastep.

| The direct control corresponds to ``controlmode``\ =1 and the
  Newton-Raphson solver is used. Under the direct control, the total
  load vector assembled for specific solution step represents the load
  level, where equilibrium is searched. The implementation supports also
  displacement control - it is possible to prescribe one or more
  displacements by applying “quasi prescribed” boundary
  condition(s). The load level then represents the time, where the
  equilibrium has been found. The Newton-Raphson solver parameters
  (``solverParams``) for load-control are:
| ``maxiter #(in)`` [``minsteplength #(in)``]
  [``minIter #(in)``] [``manrmsteps #(in)``] [``ddm #(ia)``]
  [``ddv #(ra)``] [``ddltf #(in)``] [``linesearch #(in)``]
  [``lsearchamp #(rn)``] [``lsearchmaxeta #(rn)``]
  [``lsearchtol #(rn)``] [``nccdg #(in)`` ``ccdg1 #(ia)``
  ``ccdgN #(ia)`` ] ``rtolv #(rn)`` [``rtolf #(rn)``]
  [``rtold #(tn)``] [``initialGuess #(rn)``]  where

-  ``maxiter`` determines the maximum number of iterations allowed to
   reach equilibrium. If equilibrium is not reached, the step length
   (corresponding to time) is reduced.

-  ``minsteplength`` parameter is the minimum step length allowed.

-  ``minIter`` - minimum number of iterations which always proceed
   during the iterative solution.

-  If ``manrmsteps`` parameter is nonzero, then the modified N-R scheme
   is used, with the stiffness updated after ``manrmsteps`` steps.

-  ``ddm`` is array specifying the degrees of freedom, which
   displacements are controlled. Let the number of these DOFs is N. The
   format of ``ddm`` array is 2*N dofman1 idof1 dofman2 idof2 ...
   dofmanN idofN, where the dofmani is the number of i-th dof manager
   and idofi is the corresponding DOF number.

-  ``ddv`` is array of relative weights of controlled displacements, the
   size should be equal to N. The actual value of prescribed dofs is
   defined as a product of its weight and the value of load time
   function specified using ``ddltf`` parameter (see below).

-  ``ddltf`` number of load time function, which is used to evaluate the
   actual displacements of controlled dofs.

-  ``linesearch`` nonzero value turns on line search algorithm. The
   ``lsearchtol`` defines tolerance (default value is 0.8),
   amplification factor can be specified using ``lsearchamp`` parameter
   (should be in interval :math:`(1,10)`), and parameter
   ``lsearchmaxeta`` defines maximum limit on the length of iterative
   step (allowed range is :math:`(1.5,15)`).

-  ``nccdg`` allows to define one or more DOF groups, that are used for
   evaluation of convergence criteria. Each DOF is checked if it is a
   member of particular group and in this case its contribution is taken
   into account when evaluating the convergence criteria for that group.
   By default, if ``nccdg`` is not specified, one group containing all
   DOF types is created. The value of ``nccdg`` parameter defines the
   number of DOF type groups. For each group, the corresponding DOF
   types need to be specified using ``ccdg#`` parameter, where ’#’
   should be replaced by group number (numbering starts from 1). This
   array contains the DofIDItem values, that identify the physical
   meaning of DOFs in the group. The values and their physical meaning
   is defined by DofIDItem enum type (see src/oofemlib/dofiditem.h for
   reference).

-  ``rtolv`` determines relative convergence norm (both for displacement
   iterative change vector and for residual unbalanced force vector).
   Optionally, the ``rtolf`` and ``rtold`` parameters can be used to
   define independent relative convergence crteria for unbalanced forces
   and displacement iterative change. If the default convergence
   criteria is used, the parameters ``rtolv``,\ ``rtolf``, and ``rtold``
   are real values. If the convergence criteria DOF groups are used (see
   bellow the description of ``nccdg`` parameter) then they should be
   specified as real valued arrays of ``nccdg`` size, and individual
   values define relative convergence criteria for each individual dof
   group.

-  ``initialGuess`` is an optional parameter with default vaue 0, for
   which the first iteration of each step starts from the previously
   converged state and applies the prescribed displacement increments.
   This can lead to very high strains in elements connected to the nodes
   with changing prescribed displacements and the state can be far from
   equilibrium, which may results into slow convergence and strain
   localization near the boundary. If ``initialGuess`` is set to 1, the
   contribution of the prescribed displacement increments to the
   internal nodal forces is linearized and moved to the right-hand side,
   which often results into an initial solution closer to equilibrium.
   For instance, if the step is actually elastic, equilibrium is fully
   restored after the second iteration, while the default method may
   require more iterations.

| The indirect solver corresponds to ``controlmode``\ =0 and the CALM
  solver is used. The value of reference load vector is determined by
  ``refloadmode`` parameter mentioned above at the first step of each
  metastep. However, the user must ensure that the same value of
  reference load vector could be obtained for all solution steps of
  particular metastep (this is necessary for restart and adaptivity to
  work). The corresponding meta step solver parameters
  (``solverParams``) are:
| ``Psi #(rn)`` ``MaxIter #(in)`` ``stepLength #(rn)``
  [``minStepLength #(in)``] [``initialStepLength #(rn)``]
  [``forcedInitialStepLength #(rn)``] [``reqIterations #(in)``]
  [``maxrestarts #(in)``] [``minIter #(in)``]
  [``manrmsteps #(in)``] [``hpcmode #(in)``] [``hpc #(ia)``]
  [``hpcw #(ra)``] [``linesearch #(in)``] [``lsearchamp #(rn)``]
  [``lsearchmaxeta #(rn)``] [``lsearchtol #(rn)``] [``nccdg #(in)``
  ``ccdg1 #(ia)`` ... ``ccdgN #(ia)``] ``rtolv #(rn)``
  [``rtolf #(rn)``] [``rtold #(rn)``] [``pert #(ia)``]
  [``pertw #(ra)``] [``rpa #(rn)``] [``rseed #(in)``] where

-  ``Psi`` - CALM :math:`\Psi` control parameter. For :math:`\Psi` = 0
   displacement control is applied. For nonzero values the load control
   applies together with displacement control (ALM). For large
   :math:`\Psi` load control apply.

-  ``MaxIter`` - determines the maximum number of iteration allowed to
   reach equilibrium state. If this limit is reached, restart follows
   with smaller step length.

-  ``stepLength`` - determines the maximum value of arc-length (step
   length).

-  ``minStepLength`` - minimum step length. The step length will never
   be smaller. If convergence problems are encountered and step length
   cannot be decreased, computation terminates.

-  ``initialsteplength`` - determines the initial step length (the
   arc-length). If not provided, the maximum step length (determined by
   ``stepLength`` parameter) will be used as the value of initial step
   length.

-  ``forcedInitialStepLength`` - When simulation is restarted, the last
   predicted step length is used. Use ``forcedInitialStepLength``
   parameter to override the value of step length. This parameter will
   also override the value of initial step length set by
   ``initialsteplength`` parameter.

-  ``reqIterations`` - approximate number of iterations controlled by
   changing the step length.

-  ``maxrestarts`` - maximum number of restarting computation when
   convergence not reached up to ``MaxIter``.

-  ``minIter`` - minimum number of iterations which always proceed
   during the iterative solution. ``reqIterations`` are set to be the
   same, ``MaxIter`` are increased if lower.

-  ``manrmsteps`` - Forces the use of accelerated Newton Raphson method,
   where stiffness is updated after ``manrmsteps`` steps. By default,
   the modified NR method is used (no stiffness update).

-  ``hpcmode`` Parameter determining the alm mode. Possible values are:
   0 - (default) full ALM with quadratic constrain and all dofs, 1 -
   (default, if ``hpc`` parameter used) full ALM with quadratic
   constrain, taking into account only selected dofs (see ``hpc``
   param), 2 - linearized constrain in displacements only, taking into
   account only selected dofs with given weight (see ``hpc`` and
   ``hpcw`` parameters).

-  ``hpc`` - Special parameter for Hyper-plane control, when only
   selected DOFs are taken account in ALM step length condition.
   Important mainly for material nonlinear problems with strong
   localization. This array selects the degrees of freedom, which
   displacements are controlled. Let the number of these DOFs be N. The
   format of ``ddm`` array is 2*N dofman1 idof1 dofman2 idof2 ...
   dofmanN idofN, where the dofmani is the number of i-th dof manager
   and idofi is the corresponding DOF number.

-  ``hpcw`` Array of DOF weights in linear constraint. The dof ordering
   is determined by ``hpc`` parameter, the size of the array should be
   N.

-  ``linesearch`` nonzero value turns on line search algorithm. The
   ``lsearchtol`` defines tolerance, amplification factor can be
   specified using ``lsearchamp`` parameter (should be in interval
   :math:`(1,10)`), and parameter ``lsearchmaxeta`` defines maximum
   limit on the length of iterative step (allowed range is
   :math:`(1.5,15)`).

-  ``nccdg`` allows to define one or more DOF groups, that are used for
   evaluation of convergence criteria. Each DOF is checked if it is a
   member of particular group and in this case its contribution is taken
   into account when evaluating the convergence criteria for that group.
   By default, if ``nccdg`` is not specified, one group containing all
   DOF types is created. The value of ``nccdg`` parameter defines the
   number of DOF type groups. For each group, the corresponding DOF
   types need to be specified using ``ccdg#`` parameter, where ’#’
   should be replaced by group number (numbering starts from 1). This
   array contains the DofIDItem values, that identify the physical
   meaning of DOFs in the group. The values and their physical meaning
   is defined by DofIDItem enum type (see src/oofemlib/dofiditem.h for
   reference).

-  ``rtolv`` determines relative convergence norm (both for displacement
   iterative change vector and for residual unbalanced force vector).
   Optionally, the ``rtolf`` and ``rtold`` parameters can be used to
   define independent relative convergence crteria for unbalanced forces
   and displacement iterative change. If the default convergence
   criteria is used, the parameters ``rtolv``,\ ``rtolf``, and ``rtold``
   are real values. If the convergence criteria DOF groups are used (see
   bellow the description of ``nccdg`` parameter) then they should be
   specified as real valued arrays of ``nccdg`` size, and individual
   values define relative convergence criteria for each individual dof
   group.

-  ``pert`` Array specifying DOFs that should be perturbed after the
   first iteration of each step. Let the number of these DOFs be M. The
   format of ``ddm`` array is 2*M dofman1 idof1 dofman2 idof2 ...
   dofmanN idofN, where the dofmani is the number of i-th dof manager
   and idofi is the corresponding DOF number.

-  ``pertw`` Array of DOF perturbations. The dof ordering is determined
   by ``pert`` parameter, the size of the array should be M.

-  ``rpa`` Amplitude of random perturbation that is applied to each DOF.

-  ``rseed`` Seed for the random generator that generates random
   perturbations.

Standard syntax
^^^^^^^^^^^^^^^

| In this case, all parameters (for analysis as well as for the solver)
  are supplied in analysis record. The default meta step is created for
  all solution steps required. Then the meta step attributes are
  specified within analysis record. The format of analysis record is
  then following
| ``NonLinearStatic`` ``nsteps #(in)`` [``nonlocstiff #(in)``]
  [``contextOutputStep #(in)``] [``controlmode #(in)``]
  [``deltat`` #(rn)``] ``rtolv #(rn)`` [``stiffmode #(in)``]
  ``lstype #(in)`` ``smtype #(in)``
  ``solverParams #()`` [``nonlinform #(in)``]
  <[``nonlocstiff #(in)``]> <[``nonlocalext``]> <[``loadbalancing``]

The meaning of parameters is the same as for extended syntax.

Parameter ``lstype`` allows to select the solver for the linear system
of equations. Parameter ``smtype`` allows to select the sparse matrix
storage scheme. The scheme should be compatible with the solver type.
See section :ref:`sparselinsolver` for further details.

.. _AdaptiveLinearStatic:

Adaptive linear static
~~~~~~~~~~~~~~~~~~~~~~

``Adaptlinearstatic`` ``nsteps #(in)`` [``sparselinsolverparams #(...)``]
[``meshpackage #(in)``] ``errorestimatorparams #(...)``

Adaptive
linear static analysis. Multiple loading cases are not supported. Due to
linearity of a problem, the complete reanalysis from the beginning is
done after adaptive remeshing. After first step the error is estimated,
information about required density is generated (using mesher interface)
and solution terminates. If the error criteria is not satisfied, then
the new mesh and corresponding input file is generated and new analysis
should be performed until the error is acceptable. Currently, the
available error estimator for linear problems is Zienkiewicz-Zhu. Please
note, that adaptive framework requires specific functionality provided
by elements and material models. For details, see element and material
model manuals.

-  Parameter ``nsteps`` indicates the number of loading cases. Should be
   set to 1.

-  The ``sparselinsolverparams`` parameter describes the sparse linear
   solver attributes and is explained in section
   :ref:`sparselinsolver`.

-  The ``meshpackage`` parameter selects the mesh package interface,
   which is used to generate information about required mesh density for
   new remeshing. The supported interfaces are explained in section
   :ref:`meshpackages`. By default, the T3d interface is used.

-  The ``errorerestimatorparams`` parameter contains the parameters of
   Zienkiewicz Zhu Error Estimator. These are described in section
   :ref:`errorestimators`.

Adaptive nonlinear static
~~~~~~~~~~~~~~~~~~~~~~~~~

``Adaptnlinearstatic`` ``Nonlinearstaticparams #()`` [``equilmc #(in)``]
[``meshpackage #(in)``] [``eetype #(in)``] ``errorestimatorparams #(...)``

Represents Adaptive Non-LinearStatic problem. Solution is performed as a
series of increments (loading or displacement). The error is estimated
at the end of each load increment (after equilibrium is reached), and
based on reached error, the computation continues, or the new mesh
densities are generated and solution stops. Then the new discretization
should be generated. The truly adaptive approach is supported, so the
computation can be restarted from the last step (see section
:ref:`running-the-code`), solution is mapped to new mesh (separate
solution step) and new load increment is applied. Of course, one can
start the analysis from the very beginning using new mesh. Currently,
the available estimators/indicators include only linear Zienkiewicz-Zhu
estimator and scalar error indicator. Please note, that adaptive
framework requires specific functionality provided by elements and
material models. For details, see element and material model manuals.

-  Set of parameters ``Nonlinearstaticparams`` are related to nonlinear
   analysis. They are described in section :ref:`NonLinearStatic`.

-  Parameter ``equilmc`` determines, whether after mapping of primary
   and internal variables to new mesh the equilibrium is restored or not
   before new load increment is applied. The possible values are: 0
   (default), when no equilibrium is restored, and 1 forcing the
   equilibrium to be restored before applying new step.

-  The ``meshpackage`` parameter selects the mesh package interface,
   which is used to generate information about required mesh density for
   new remeshing. The supported interfaces are explained in section
   :ref:`meshpackages`. By default, the T3d interface is used.

-  Parameter ``eetype`` determines the type of error estimator/indicator
   to be used. The parameters ``errorestimatorparams`` represent set of
   parameters corresponding to selected error estimator. For
   description, follow to section :ref:`errorestimators`.

.. _FreeWarping:

Free warping analysis
~~~~~~~~~~~~~~~~~~~~~

``FreeWarping`` ``nsteps #(in)``

Free warping analysis computes the deplanation function of cross section
with arbitrary shape. It is done by solving the Laplace’s equation with
automatically generated boundary conditions corresponding to the free
warping problem.

This type of analysis supports only ``TrWarp`` elements and
``WarpingCS`` cross sections. One external node must be defined for each
warping cross section. The coordinates of this node can be arbitrary but
this node must be defined with parametr ``DofIDMask 1 24`` and one
boundary condition which represents relative twist acting on
corresponding warping cross section. No additional loads make sence in
free warping analysis.

Parameter ``nsteps`` indicates the number of loading cases. Series of
loading cases is maintained as sequence of time-steps. For each load
case an auxiliary time-step is generated with time equal to load case
number. Load vectors for each load case are formed as load vectors at
this auxiliary time.

Transport Problems
------------------

.. _StationaryTransport:

Stationary transport problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``StationaryProblem`` ``nsteps #(in)`` [``sparselinsolverparams #(...)``] [``exportfields #(ia)``]

Stationary transport problem. Series of loading cases is maintained as
sequence of time-steps. For each load case an auxiliary time-step is
generated with time equal to load case number. Load vectors for each
load case are formed as load vectors at this auxiliary time. The
``sparselinsolverparams`` parameter describes the sparse linear solver
attributes and is explained in section :ref:`sparselinsolver`.

If the present problem is used within the context of staggered-like
analysis, the temperature field obtained by the solution can be exported
and made available to any subsequent analyses. For example, temperature
field obtained by present analysis can be taken into account in
subsequent mechanical analysis. To allow this, the temperature must be
“exported”. This can be done by adding array ``exportfields``. This
array contains the field identifiers, which tell the problem to register
its primary unknowns under given identifiers. See file field.h. Then the
subsequent analyses can get access to exported fields and take them into
account, if they support such feature.

.. _TransientTransport:

Transient transport problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``TransientTransport`` ``nsteps #(in)`` ``deltaT #(rn)`` *or* ``dTfunction #(in)`` *or* ``prescribedtimes #(ra)``
``alpha #(rn)`` [``initT #(rn)``] [``lumped``]
[``keeptangent``] [``exportfields #(ia)``]

**Nonlinear** implicit integration scheme for transient transport
problems. The generalized midpoint rule (sometimes called
:math:`\alpha`-method) is used for time discretization, with alpha
parameter, which has limits :math:`0\le\alpha\le1`. For :math:`\alpha=0`
explicit Euler forward method is obtained, for :math:`\alpha=0.5`
implicit trapezoidal rule is recovered, which is unconditionally stable,
second-order accurate in :math:`\Delta t`, and :math:`\alpha=1.0` yields
implicit Euler backward method, which is unconditionally stable, and
first-order accurate in :math:`\Delta t`. ``deltaT`` is time step length
used for integration, ``nsteps`` parameter specifies number of time
steps to be solved. It is possible to define ``dTfunction`` with a
number referring to corresponding time function, see
section :ref:`TimeFunctionsRecords`. Variable time step is
advantageous when calculating large time intervals.

The ``initT`` sets the initial time for integration, 0. by default. If
``lumped`` is set, then the stabilization of numerical algorithm using
lumped capacity matrix will be used, reducing the initial oscillations.
See section :ref:`StationaryTransport` for an explanation on
``exportfields``.

This transport problem supports sets and changes in number of equations.
It is possible to impose/remove Dirichlet boundary conditions during
solution.

.. _LinearTransientTransport:

Transient transport problem - linear case - obsolete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NonStationaryProblem`` ``nsteps #(in)`` ``deltaT #(rn)`` |
``deltaTfunction #(in)`` ``alpha #(rn)`` [``initT #(rn)``]
[``lumpedcapa``] [``sparselinsolverparams #(..)``]
[``exportfields #(ia)``] [``changingProblemSize``]

**Linear** implicit integration scheme for transient transport problems.
The generalized midpoint rule (sometimes called :math:`\alpha`-method)
is used for time discretization, with alpha parameter, which has limits
:math:`0\le\alpha\le1`. For :math:`\alpha=0` explicit Euler forward
method is obtained, for :math:`\alpha=0.5` implicit trapezoidal rule is
recovered, which is unconditionally stable, second-order accurate in
:math:`\Delta t`, and :math:`\alpha=1.0` yields implicit Euler backward
method, which is unconditionally stable, and first-order accurate in
:math:`\Delta t`. ``deltaT`` is time step length used for integration,
``nsteps`` parameter specifies number of time steps to be solved. It is
possible to define ``deltaTfunction`` with a number referring to
corresponding time function, see
section :ref:`TimeFunctionsRecords`. Variable time step is
advantageous when calculating large time intervals. It is strongly
suggested to use nonlinear transport solver due to stability reasons,
see section :ref:`TransientTransport`.

The ``initT`` sets the initial time for integration, 0 by default. If
``lumpedcapa`` is set, then the stabilization of numerical algorithm
using lumped capacity matrix will be used, reducing the initial
oscillations. See section :ref:`StationaryTransport` for an
explanation on ``exportfields``.

This linear transport problem supports changes in number of equations.
It is possible to impose/remove Dirichlet boundary conditions during
solution. This feature is enabled with ``changingProblemSize``, which
ensures storing solution values on nodes (DoFs) directly. If the problem
does not grow/decrease during solution, it is more efficient to use
conventional solution strategy and the parameter should not be
mentioned.

Note: This problem type **requires transport module** and it can be used
only when this module is included in your oofem configuration.

.. _TransientTransport2:

Transient transport problem - nonlinear case - obsolete
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NlTransientTransportProblem``  ``nsteps #(in)`` ``deltaT #(rn)`` |
``deltaTfunction #(in)`` ``alpha #(rn)`` [``initT #(rn)``]
[``lumpedcapa #()``] [``nsmax #(in)``] ``rtol #(rn)``
[``manrmsteps #(in)``] [``sparselinsolverparams #(...)``]
[``exportfields #(ia)``] [``changingProblemSize``]

Implicit integration scheme for transient transport problems. The
generalized midpoint rule (sometimes called :math:`\alpha`-method) is
used for time discretization, with alpha parameter, which has limits
:math:`0\le\alpha\le1`. For :math:`\alpha=0` explicit Euler forward
method is obtained, for :math:`\alpha=0.5` implicit trapezoidal rule is
recovered, which is unconditionally stable, second-order accurate in
:math:`\Delta t`, and :math:`\alpha=1.0` yields implicit Euler backward
method, which is unconditionally stable, and first-order accurate in
:math:`\Delta t`. See matlibmanual.pdf for solution algorithm.

``deltaT`` is time step length used for integration, ``nsteps``
parameter specifies number of time steps to be solved. For
``deltaTfunction`` and ``initT`` see
section :ref:`LinearTransientTransport`. Parameter ``maxiter``
determines the maximum number of iterations allowed to reach equilibrium
(default is 30). Norms of residual physical quantity (heat, mass)
described by solution vector and the change of solution vector are
determined in each iteration. The convergence is reached, when the norms
are less than the value given by ``rtol``. If ``manrmsteps`` parameter
is nonzero, then the modified N-R scheme is used, with the left-hand
side matrix updated after ``manrmsteps`` steps. ``nsmax`` maximum number
of iterations per time step, default is 30. If ``lumpedcapa`` is set,
then the stabilization of numerical algorithm using lumped capacity
matrix will be used, reducing the initial oscillations.

See the Section :ref:`StationaryTransport` for an explanation on
``exportfields``. The meaning of ``changingProblemSize`` is given in
Section :ref:`LinearTransientTransport`.

Note: This problem type **requires transport module** and it can be used
only when this module is included in your oofem configuration.

Fluid Dynamic Problems
----------------------

.. _cbsIncomp:

Transient incompressible flow - CBS Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``CBS`` ``nsteps #(in)`` ``deltaT #()`` [``theta1 #(in)``]
[``theta2 #(in)``] [``cmflag #(in)``] [``scaleflag #(in)``
``lscale #(in)`` ``uscale #(in)`` ``dscale #(in)``] [``lstype #(in)``]
[``smtype #(in)``]

Solves the transient incompressible flow using algorithm based on
Characteristics Based Split (CBS, for reference see O.C.Zienkiewics and
R.L.Taylor: The Finite Element Method, 3rd volume,
Butterworth-Heinemann, 2000). At present, only semi-implicit form of the
algorithm is available and energy equation, yielding the temperature
field, is not solved. Parameter ``nsteps`` determines number of solution
steps. Parameter ``deltaT`` is time step length used for integration.
This time step will be automatically adjusted to satisfy integration
stability limits
:math:`\Delta t \le {\frac{h}{\vert\boldsymbol{u}\vert}}` and
:math:`\Delta t \le {\frac{h^2}{2\nu}}`, if necessary. Parameters
``theta1`` and ``theta2`` are integration constants,
:math:`\theta_1, \theta_2 \in \langle{\frac12}, 1\rangle`. If ``cmflag``
is given a nonzero value, then consistent mass matrix will be used
instead of (default) lumped one.

The characteristic equations can be solved in non-dimensional form. To
enable this, the ``scaleflag`` should have a nonzero value, and the
following parameters should be provided: ``lscale``, ``uscale``, and
``dscale`` representing typical length, velocity, and density scales.

Parameter ``lstype`` allows to select the solver for the linear system
of equations. Parameter ``smtype`` allows to select the sparse matrix
storage scheme. The scheme should be compatible with the solver type.
See section :ref:`sparselinsolver` for further details.

.. _supgIncomp:

Transient incompressible flow SUPG/PSPG Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SUPG`` ``nsteps #(in)`` ``deltaT #(rn)`` ``rtolv #(rn)``
[``atolv #(rn)``] [``stopmaxiter #(in)``] [``alpha #(rn)``]
[``cmflag #(in)``] [``deltatltf #(in)``] [``miflag #(in)``]
[``scaleflag #(in)`` ``lscale #(in)`` ``uscale #(in)``
``dscale #(in)``] [``lstype #(in)``] [``smtype #(in)``]

Solves the transient incompressible flow using stabilized formulation
based on SUPG and PSPG stabilization terms. The stabilization provides
stability and accuracy in the solution of advection-dominated problems
and permits usage of equal-order interpolation functions for velocity
and pressure. Furthermore, stabilized formulation significantly improves
convergence rate in iterative solution of large nonlinear systems of
equations.

By changing the value :math:`\alpha`, different methods from
“Generalized mid-point family” can be chosen, i.e., Forward Euler
(:math:`\alpha=0`), Midpoint rule (:math:`\alpha=0.5`), Galerkin
(:math:`\alpha=2/3`), and Backward Euler (:math:`\alpha=1`). Except the
first one, all the methods are implicit and require matrix inversion for
solution. Some results form an energy method analysis suggest
unconditional stability for :math:`\alpha\ge 0.5` for the generalized
mid-point family. As far as accuracy is concerned, the midpoint rule is
to be generally preferred.

Parameter ``nsteps`` determines number of solution steps. Parameter
``deltaT`` is time step length used for integration. Alternatively, the
load time function can be used to determine time step length for
particular solution step. The load time function number is determined by
parameter ``deltatltf`` and its value evaluated for solution step number
should yield the step length.

Parameters ``rtolv`` and ``atolv`` allow to specify relative and
absolute errors norms for residual vector. The equilibrium iteration
process will stopped when both error limits are satisfied or when the
number of iteration exceeds the value given by parameter
``stopmaxiter``.

If ``cmflag`` is given a nonzero value, then consistent mass matrix will
be used instead of (default) lumped one.

The algorithm allows to solve the flow of two immiscible fluids in fixed
spatial domain (currently only in 2d). This can be also used for solving
free surface problems, where one of the fluids should represent air. To
enable multi-fluid analysis, user should set parameter ``miflag``. The
supported values are described in
section :ref:`materialinterfaces`. Please note, that the initial
distribution of reference fluid volume should be provided as well as
constitutive models for both fluids.

The characteristic equations can be solved in non-dimensional form. To
enable this, the ``scaleflag`` should have a nonzero value, and the
following parameters should be provided: ``lscale``, ``uscale``, and
``dscale`` representing typical length, velocity, and density scales.

Parameter ``lstype`` allows to select the solver for the linear system
of equations. Parameter ``smtype`` allows to select the sparse matrix
storage scheme. Please note that the present algorithm leads to a
non-symmetric matrix. The scheme should be compatible with the solver
type. See section :ref:`sparselinsolver` for further details.

.. _pfemIncomp:

Transient incompressible flow (PFEM Algorithm)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``PFEM`` ``nsteps #(in)``  ``deltaT #(rn)`` ``material #(in)``
``cs #(in)`` ``pressure #(in)`` [``mindeltat #(rn)``]
[``maxiter #(in)``] [``rtolv #(rn)``] [``rtolp #(rn)``]
[``alphashapecoef #(rn)``] [``removalratio #(rn)``]
[``scheme #(in)``] [``lstype #(in)``] [``smtype #(in)``]

Solves the transient incompressible flow using particle finite element
method based on the Lagrangian formulation of Navier-Stokes equations.

Mesh nodes are represented by
PFEMParticles (see :ref:`pfemparticles`), which can freely
move and even separate from the main domain. To integrate governing
equations in each solution step, a temporary mesh, built from particles,
is needed. The mesh is rebuilt from scratch in each solution step to
prevent large distortion of elements. Paramters ``cs`` and ``material``
assign types from cross section and material record to created elements.
Thus, the problem is defined without any elements in the input file.

Mesh is generated using Delaunay triangulation and Alpha shape technique
for the identification of the free surface. The parameter
``alphashapecoef`` should reflect initial distribution of PFEMParticles.
Value approximately equal to 1,5-multiple of shortest distance of two
neighboring particles has been found well. On the free surface the
zero-pressure boundary condition is enforced. This must be defined in
boundary condition record under the number defined by ``pressure``.

Parameter ``scheme`` controls whether the equation system for the
components of the auxiliary velocity is solved explicitly (0) or
implicitly (1). The last is the default option.

Parameter ``nsteps`` determines number of solution steps. Parameter
``deltaT`` is time step length used for integration. To ensure numerical
stability, step length is adapted upon mesh geometry and velocity of
paricular nodes. To avoid to short time length a minimal size can be
defined by ``mindeltat``. Alternatively prescribing limit
``removalratio`` of the element edge length too close particles can be
removed from solution.

Optional parameters ``rtolv`` and ``rtolp`` allow to specify relative
norms for velocity and pressure difference of two subsequent iteration
step. Default values are 1.e-8. By default maximal 50 iterations are
performed, if not specified by ``maxiter``.

Parameter ``lstype`` allows to select the solver for the linear system
of equations. Parameter ``smtype`` allows to select the sparse matrix
storage scheme. Please note that the present algorithm leads to a
non-symmetric matrix. The scheme should be compatible with the solver
type. See section :ref:`sparselinsolver` for further details.

Coupled Problems
----------------

.. _staggeredproblem:

Staggered Problem
~~~~~~~~~~~~~~~~~

``StaggeredProblem`` (``nsteps #(in)`` ``deltaT #(rn))`` :math:`|`
``timeDefinedByProb #(in)`` ``prob1 #(s)`` ``prob2 #(s)``
[``stepMultiplier #(rn)``]

Represent so-called staggered analysis. This can be described as an
sequence of sub-problems, where the result of some sub-problem in the
sequence can depend on results of previous sub-problems in sequence.
Typical example is heat transfer analysis followed by mechanical
analysis taking into account the temperature field generated by the heat
transfer analysis. Similar analysis can be done when coupling moisture
transport with concrete drying strain.

The actual implementation supports only sequence of two sub-problems.
The sub-problems are described using sub-problem input files. The syntax
of sub-problem input file is the same as for standalone problem. The
only addition is that sub-problems should export their solution fields
so that they became available for subsequent sub-problems. See the
Section :ref:`StationaryTransport`.

The subproblem input files are described using ``prob1`` and ``prob2``
parameters, which are strings containing a path to sub-problem input
files, the ``prob1`` contains input file path of the first sub-problem,
which runs first for each solution step, the ``prob2`` contains input
file path of the second sub-problem.

There are two options how to control a time step sequence. The first
approach uses ``timeDefinedByProb`` which uses time sequence from the
corresponding subproblem. The subproblem may specify arbitrary loading
steps and allows high flexibility. The second approach uses the
staggered problem to take control over time. Therefore any sub-problem
time-stepping parameters are ignored (even if they are required by
sub-problem input syntax) and only staggered-problem parameters are
relevant. ``deltaT`` is than a time step length used for integration,
``nsteps`` parameter specifies number of time steps to be solved.
``stepMultiplier`` multiplies all times with a given constant. Default
is 1.

Note: This problem type **is included in transport module** and it can
be used only when this module is configured. Note: All material models
derived from StructuralMaterial base will take into account the external
registered temperature field, if provided.

.. _fluidstructureproblem:

FluidStructure Problem
~~~~~~~~~~~~~~~~~~~~~~

``FluidStructureProblem``  ``nsteps #(in)`` ``deltaT #(rn)`` ``prob1 #(s)``
``prob2 #(s)`` [``maxiter #(in)``] [``rtolv #(rn)``]
[``rtolp #(rn)``]

Represents a fluid-structure analysis based on StaggeredProblem but
providing iterative synchronization of sub-problems. The implementation
uses the the PFEM model :ref:`pfemIncomp` for the fluid part.
For the structural part a full dynamic analysis using implicit direct
integration DIIDynamic :ref:`DIIDynamic` is considered.

The coupling of both phases is based on the idea of enforcing
compatibility on the interface. Special fluid particle are attached to
every structural node on the interface that can be hit by the fluid.
These special particles have no degrees of freedom associated, so no
equations are solved on them. However, their movement is fully
determined by associated structural nodes. Their velocities governed by
the solid part affect the fluid equation naturally.

This iterative procedure is based on the so-called Dirichlet-Neumann
approach. Dirichlet boundary conditions are the prescribed velocities on
the fluid side of the interface, whereas applied forces on the
structural side represent the Neumann boundary conditions.

The convergence criterion is based on the difference of the pressure and
velocity values on the interface from the subsequent iterative steps.
Once they are smaller than prescribed tolerance, the iteration is
terminated and solution can proceed to the next step.

The subproblem input files are described using ``prob1`` and ``prob2``
parameters, which are strings containing a path to sub-problem input
files, the ``prob1`` contains input file path of the first sub-problem,
which runs first for each solution step, the ``prob2`` contains input
file path of the second sub-problem. The time step sequence is
controlled by the number of steps ``nsteps`` and the time step length
``deltaT``.

Optional parameters ``rtolv`` and ``rtolp`` allow to specify relative
norms for velocity and pressure differnce of two subsequent iteration
step. Default values are 1.e-3. By default maximal 50 iterations are
performed, if not specified by ``maxiter``.

Note: This problem type **is included in PFEM module** and it can be
used only when this module is configured.

.. _DummyEngngModel:

DummyEngngModel
~~~~~~~~~~~~~~~~

``Dummy``  ``nnmodules #(in)`` 

Represents a dummy model, whch is not capable to perform any analysis. 
Its intended use is to invoke the configured export modules, 
so that the problem geometry can be exported without requiring to actually solve the problem.

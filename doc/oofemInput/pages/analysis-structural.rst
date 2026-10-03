Structural problems
===================

.. _StaticStructural:

StaticStructural
----------------

Static structural analysis.  Solves linear and nonlinear static structural
problems, supporting changes in boundary conditions (applied load and
supports).  The problem can be solved under direct load or displacement
control, indirect control, or any combination of the two.

The individual solution steps describe the history of the applied incremental
loading.  Load cases are *not* supported: a new analysis has to be run for each
load case.  To analyse a linear static problem with several load combinations,
use :ref:`LinearStatic` instead.

By default all material nonlinearities are taken into account, geometrical ones
are not.  To include geometrically nonlinear effects, specify the level of
non-linearity in the element records.

.. record::

   :descitem:`StaticStructural` :elemparam:`nsteps{in}`
   :optelemparam:`deltat{rn}` :optelemparam:`prescribedtimes{ra}`
   :optelemparam:`solvertype{s}` :optelemparam:`stiffmode{in}`
   :optelemparam:`initialguess{in}` :optelemparam:`smtype{in}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`recomputeaftercrackprop{}`
   ``<``\ :optelemparam:`nonlocalext{}`\ ``>``

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.  Ignored when :param:`prescribedtimes` is given.

:optparam:`deltat{rn}`
    Length of the time step.  Default ``1.0``.  Ignored when
    :param:`prescribedtimes` is given.

:optparam:`prescribedtimes{ra}`
    Array of discrete solution times.  When present, the number of solution
    steps equals the size of this array.

:optparam:`solvertype{s}`
    Name of the non-linear solver.  Default ``nrsolver``, the Newton-Raphson
    solver whose attributes are documented under :ref:`NonLinearStatic`.

:optparam:`stiffmode{in}`
    Which stiffness the solver requests: ``0`` tangent (default), ``1`` secant,
    ``2`` the original elastic stiffness throughout.

:optparam:`initialguess{in}`
    How each step starts; see :ref:`common-analysis-parameters`.  For
    ``StaticStructural`` the default is ``1`` (approximated tangent problem),
    which is the better choice when prescribed displacements change.

:optparam:`smtype{in}`
    Sparse matrix storage scheme.  Default ``0``, the symmetric skyline.  See
    :ref:`sparselinsolver`.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`recomputeaftercrackprop{}`
    A flag; when present, the solution step is recomputed after crack
    propagation.

:optparam:`nonlocalext{}`
    Parallel build only.  Turns on the nonlocal constitutive extension, which
    maintains a band of remote elements taking part in the computation of
    nonlocal variables — see :numref:`fig-nodecut-nlm`.  In parallel mode the
    analysis record must be identical in every partition input file, and the
    parallel version requires the PETSc module.

Adaptive time stepping is available through the meta-step attributes described
in :ref:`timestepping`.

Contact analysis
~~~~~~~~~~~~~~~~

In addition to the standard boundary condition types, ``StaticStructural``
supports **contact analysis** through a penalty-based contact formulation.
Contact interactions are defined with ``structuralpenaltycontactbc`` boundary
conditions (see :ref:`structural-penalty-contact`) together with
``StructuralContactElement_*`` elements and ``StructuralFEContactSurface``
records (see :ref:`ContactSurfaceRecords`).  This makes frictionless — and,
experimentally, frictional — contact between deformable bodies possible within a
static framework.

.. _LinearStatic:

Linear static analysis
----------------------

Linear static analysis.  The problem supports multiple load cases: the number of
load cases equals the number of solution steps, and the individual load vectors
are formed in the individual time steps.  The static system, however, is assumed
to be the same for all load cases.  For each load case an auxiliary time step is
generated whose time equals the load case number.

.. record::

   :descitem:`LinearStatic` :elemparam:`nsteps{in}`
   :optelemparam:`sparselinsolverparams{...}`

**Parameters**

:param:`nsteps{in}`
    Number of load cases.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

.. _Pdelta:

Pdelta
------

Performs P-delta analysis, a form of geometrically nonlinear analysis in which
the effect of axial (normal) forces on the deformed configuration is taken into
account.  It solves the nonlinear equation

.. math::

   \left(\boldsymbol{K}+\boldsymbol{K}_s(\boldsymbol{r})\right)\boldsymbol{r}
   = \boldsymbol{f},

where :math:`\boldsymbol{K}` is the stiffness matrix,
:math:`\boldsymbol{K}_s` the initial stress matrix, :math:`\boldsymbol{r}` the
displacement vector and :math:`\boldsymbol{f}` the load vector.

.. record::

   :descitem:`pdelta` :elemparam:`nsteps{in}` :elemparam:`rtolv{rn}`
   :optelemparam:`lumped{in}` :optelemparam:`maxiter{in}`
   :optelemparam:`sparselinsolverparams{...}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps; set it to ``1``.

:param:`rtolv{rn}`
    Relative convergence tolerance.

:optparam:`lumped{in}`
    Nonzero (the default) uses the lumped initial stress matrix; set it to zero
    to use the consistent initial stress matrix.

:optparam:`maxiter{in}`
    Maximum number of iterations.  Default ``50``.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

.. _EigenValueDynamic:

EigenValueDynamic
-----------------

Eigenvalue dynamic analysis.  Only the :param:`nroot` smallest eigenvalues and
their corresponding eigenvectors are computed.

.. record::

   :descitem:`EigenValueDynamic` :elemparam:`nroot{in}`
   :elemparam:`rtolv{rn}` :optelemparam:`eigensolverparams{...}`

**Parameters**

:param:`nroot{in}`
    Number of eigenvalues to compute.

:param:`rtolv{rn}`
    Relative convergence criterion.

:optparam:`eigensolverparams{...}`
    Attributes of the eigenvalue solver; see :ref:`eigensolverssection`.

.. _NlDEIDynamic:

NlDEIDynamic
------------

Direct explicit **nonlinear** dynamic integration.  The central difference
method with a diagonal mass matrix is used, and the damping matrix is assumed
proportional to the mass matrix,
:math:`\boldsymbol{C} = \mathrm{dumpcoef}\cdot\boldsymbol{M}`, where
:math:`\boldsymbol{M}` is the diagonal mass matrix.

.. record::

   :descitem:`NlDEIDynamic` :elemparam:`nsteps{in}`
   :elemparam:`dumpcoef{rn}` :optelemparam:`deltaT{rn}`
   :optelemparam:`reduct{rn}` :optelemparam:`drflag{in}`
   :optelemparam:`tau{rn}` :optelemparam:`py{rn}`
   ``<``\ :optelemparam:`nonlocalext{}`\ ``>``

**Parameters**

:param:`nsteps{in}`
    Number of time steps to analyse.  If :param:`deltaT` is reduced
    internally, :param:`nsteps` is adjusted so that the total analysis time
    stays the same.

:param:`dumpcoef{rn}`
    Mass-proportional damping coefficient.

:optparam:`deltaT{rn}`
    Time step length used for integration.  The program may reduce it to
    satisfy the stability conditions of the solution.

:optparam:`reduct{rn}`
    Scaling factor, smaller than 1, multiplying the step length the program
    determined.

:optparam:`drflag{in}`
    A nonzero value turns on dynamic relaxation, used to find a static
    solution by damping the transient response.  When it is nonzero,
    :param:`tau` and :param:`py` become **mandatory**.

:optparam:`tau{rn}`
    Dynamic relaxation: the :math:`\tau` parameter.  Mandatory when
    :param:`drflag` is nonzero.

:optparam:`py{rn}`
    Dynamic relaxation: the estimate of the loading magnitude.  Mandatory when
    :param:`drflag` is nonzero.

:optparam:`nonlocalext{}`
    Parallel build only.  Turns on the nonlocal constitutive extension; see
    :param:`nonlocalext` under :ref:`StaticStructural`.

.. _DEIDynamic:

DEIDynamic
----------

Direct explicit **linear** integration scheme for dynamic problems.  The central
difference method with a diagonal mass matrix is used, and the damping matrix is
assumed proportional to the mass matrix,
:math:`\boldsymbol{C} = \mathrm{dumpcoef}\cdot\boldsymbol{M}`.

.. record::

   :descitem:`DEIDynamic` :elemparam:`nsteps{in}`
   :elemparam:`dumpcoef{rn}` :optelemparam:`deltaT{rn}`

**Parameters**

:param:`nsteps{in}`
    Number of time steps to analyse.

:param:`dumpcoef{rn}`
    Mass-proportional damping coefficient.

:optparam:`deltaT{rn}`
    Time step length used for integration.  The program may reduce it to
    satisfy the stability conditions of the solution.

.. _DIIDynamic:

DIIDynamic
----------

Direct implicit integration of linear dynamic problems.  The solution procedure
is described in K. Subbaraj and M. A. Dokainish, *A survey of direct
time-integration methods in computational structural dynamics — II. Implicit
methods*, Computers & Structures **32**\ (6), 1387–1401, 1989.

Damping is modelled as Rayleigh damping,
:math:`\boldsymbol{C} = \eta\,\boldsymbol{M} + \delta\,\boldsymbol{K}`.

.. record::

   :descitem:`DIIDynamic` :elemparam:`nsteps{in}` :elemparam:`deltaT{rn}`
   :optelemparam:`ddtscheme{in}` :optelemparam:`gamma{rn}`
   :optelemparam:`beta{rn}` :optelemparam:`eta{rn}`
   :optelemparam:`delta{rn}` :optelemparam:`theta{rn}`

**Parameters**

:param:`nsteps{in}`
    Number of time steps to analyse.

:param:`deltaT{rn}`
    Time step length used for integration.

:optparam:`ddtscheme{in}`
    Integration scheme, as defined in ``src/core/timediscretizationtype.h``:
    ``0`` ``TD_ThreePointBackward`` (default), ``1`` ``TD_TwoPointBackward``,
    ``2`` ``TD_Newmark``, ``3`` ``TD_Wilson``, ``4`` ``TD_Explicit``.

:optparam:`gamma{rn}`, :optparam:`beta{rn}`
    Newmark parameters governing the stability and accuracy of the integration
    algorithm.  Defaults are :math:`\gamma=0.5` and :math:`\beta=0.25`.  With
    :math:`\gamma=0.5` and :math:`\beta=1/6` the linear acceleration method is
    obtained.  Unconditional stability requires
    :math:`2\beta \ge \gamma \ge 1/2`.

:optparam:`theta{rn}`
    Required by the Wilson-:math:`\theta` method.  Default ``1.37``.

:optparam:`eta{rn}`, :optparam:`delta{rn}`
    Mass- and stiffness-proportional Rayleigh damping coefficients
    :math:`\eta` and :math:`\delta`.

.. _IncrementalLinearStatic:

IncrementalLinearStatic
-----------------------

Incremental **linear** static problem.  The problem is solved as a series of
linear solutions and is intended for linear creep problems and incremental
perfect plasticity.  Changes of the static scheme — applying, removing and
changing boundary conditions — are supported during the analysis.

.. note::

   The record keyword is ``incrlinearstatic``, not the class name
   ``IncrementalLinearStatic``.

.. record::

   :descitem:`incrlinearstatic` :elemparam:`endOfTimeOfInterest{rn}`
   :elemparam:`prescribedTimes{ra}` :optelemparam:`deltat{rn}`

**Parameters**

:param:`endOfTimeOfInterest{rn}`
    Time at the end of the interval of interest.

:param:`prescribedTimes{ra}`
    Times at which the response is computed.  These should include the times at
    which boundary conditions change, plus any other times of interest.  For
    linear creep analysis, and in the absence of changes in loading or boundary
    conditions, the values should be distributed uniformly on a logarithmic
    time scale.

:optparam:`deltat{rn}`
    Length of the time step, used when :param:`prescribedTimes` is not given.

.. _NonLinearStatic:

NonLinearStatic
---------------

Non-linear static analysis.  The problem can be solved under direct load or
displacement control, indirect control, or any combination of the two.  By
default all material nonlinearities are included, geometrical ones are not; to
include geometrically nonlinear effects, specify the level of non-linearity in
the element records.

There are two ways of supplying the parameters: the *extended* syntax, which
uses meta-steps, and the *standard* syntax, which puts everything in the
analysis record.

Extended syntax
~~~~~~~~~~~~~~~

.. record::

   :descitem:`NonLinearStatic` :elemparam:`nsteps{in}`
   :optelemparam:`nmsteps{in}` :optelemparam:`contextOutputStep{in}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`nonlinform{in}`
   ``<``\ :optelemparam:`nonlocstiff{in}`\ ``>``
   ``<``\ :optelemparam:`nonlocalext{}`\ ``>``
   ``<``\ :optelemparam:`loadbalancing{}`\ ``>``

This record is immediately followed by the meta-step records described below.

**Parameters**

:optparam:`nmsteps{in}`
    Number of meta-steps.  Default ``1``.

:param:`nsteps{in}`
    Number of solution steps.

:optparam:`contextOutputStep{in}`
    Create a context file every :param:`contextOutputStep`-th step, and
    whenever one is needed.  Useful for post-processing.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`nonlinform{in}`
    Formulation of the non-linear problem.  ``1`` (default) uses the total
    Lagrangian formulation on the undeformed original shape (first-order
    theory).  ``2`` adds the equilibrated displacements to the original ones
    and updates them in each time step (second-order theory).

:optparam:`nonlocstiff{in}`
    Parallel build only.  Whether the tangent stiffness extension for nonlocal
    models is active: ``0`` (default) inactive, ``1`` active.

:optparam:`nonlocalext{}`
    Parallel build only.  Turns on the nonlocal constitutive extension, which
    maintains a band of remote elements taking part in the computation of
    nonlocal variables — see :numref:`fig-nodecut-nlm`.

:optparam:`loadbalancing{}`
    Parallel build only.  Turns on dynamic load balancing; its attributes are
    described in :ref:`dynamicloadbalancing`.

Meta-step record
~~~~~~~~~~~~~~~~

.. record::

   :elemparam:`nsteps{in}` :elemparam:`solverParams{...}`
   :optelemparam:`controlmode{in}` :optelemparam:`deltat{rn}`
   :optelemparam:`stiffmode{in}` :optelemparam:`refloadmode{in}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`donotfixload{}`

**Parameters**

:optparam:`controlmode{in}`
    Type of solution control used for the meta-step.  ``0`` (default) uses
    indirect control — the arc-length method, i.e. the CALM solver.  ``1`` uses
    direct displacement or load control, i.e. the Newton-Raphson solver; in
    that mode prescribed load increments and controlled displacements can both
    be applied.

:optparam:`deltat{rn}`
    Time step length.  Default ``1.0``.  Each solution step has an associated
    intrinsic time at which the loading is generated, and :param:`deltat`
    determines the spacing of the solution steps on the time scale.

:optparam:`stiffmode{in}`
    ``0`` (default) uses the tangent stiffness at the beginning of a new step
    and whenever the numerical method asks for a stiffness update.  ``1``
    forces the use of the secant stiffness in the same places.  ``2`` uses the
    original elastic stiffness throughout the solution.

:optparam:`refloadmode{in}`
    How the reference force load vector is obtained from the given total and
    initial load vectors.  The initial load vector is the part of the loading
    that does not scale.  This works only for force loading; other,
    non-force components such as temperature and prescribed displacements
    should always be given as total values.  ``0`` (``rlm_total``, default)
    defines the reference incremental load vector as the total load vector
    assembled at the given time; ``1`` (``rlm_incremental``) obtains the
    reference load vector as the incremental load vector at the given time.

:param:`solverParams{...}`
    Parameters of the solver.  Which solver they belong to is determined by
    :param:`controlmode`: see :ref:`nrsolver-params` for direct control and
    :ref:`calm-params` for indirect control.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`donotfixload{}`
    By default the load reached at the end of a meta-step is maintained in the
    subsequent steps as a fixed, non-scaling load, and the load level is reset
    to zero.  When this keyword is present the loading continues instead and
    the load level is not reset.  Under indirect control the load reached is
    not fixed in either case; a new reference loading vector is assembled for
    the new meta-step.

.. _nrsolver-params:

Solver parameters for direct control (Newton-Raphson)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Direct control corresponds to :param:`controlmode`\ ``=1``, for which the
Newton-Raphson solver is used.  The total load vector assembled for a solution
step is the load level at which equilibrium is sought.  The implementation also
supports displacement control: one or more displacements can be prescribed by
applying "quasi prescribed" boundary conditions, and the load level then
represents the time at which equilibrium was found.

.. note::

   The problem does not support changes of the static system.  Direct
   displacement control can, however, be applied without a boundary condition
   being present (see the ``nrsolver`` documentation), so direct displacement
   control can be combined with direct load control or with indirect control.

.. record::

   :elemparam:`maxiter{in}` :elemparam:`rtolv{rn}`
   :optelemparam:`minsteplength{rn}` :optelemparam:`minIter{in}`
   :optelemparam:`manrmsteps{in}` :optelemparam:`ddm{ia}`
   :optelemparam:`ddv{ra}` :optelemparam:`ddltf{in}`
   :optelemparam:`linesearch{in}` :optelemparam:`lsearchamp{rn}`
   :optelemparam:`lsearchmaxeta{rn}` :optelemparam:`lsearchtol{rn}`
   :optelemparam:`nccdg{in}` :optelemparam:`ccdg1{ia}` ``...``
   :optelemparam:`ccdgN{ia}` :optelemparam:`rtolf{rn}`
   :optelemparam:`rtold{rn}` :optelemparam:`initialGuess{in}`

**Parameters**

:param:`maxiter{in}`
    Maximum number of iterations allowed to reach equilibrium.  If equilibrium
    is not reached, the step length — corresponding to time — is reduced.

:optparam:`minsteplength{rn}`
    Minimum step length allowed.

:optparam:`minIter{in}`
    Minimum number of iterations always performed during the iterative
    solution.

:optparam:`manrmsteps{in}`
    When nonzero, the modified Newton-Raphson scheme is used, with the
    stiffness updated after :param:`manrmsteps` steps.

:optparam:`ddm{ia}`
    Degrees of freedom whose displacements are controlled.  With *N* such DOFs,
    the array has the form ``2*N dofman1 idof1 dofman2 idof2 ... dofmanN
    idofN``, where *dofmani* is the number of the *i*-th dof manager and
    *idofi* the corresponding DOF number.

:optparam:`ddv{ra}`
    Relative weights of the controlled displacements; its size should be *N*.
    The actual value of a prescribed DOF is the product of its weight and the
    value of the time function given by :param:`ddltf`.

:optparam:`ddltf{in}`
    Number of the time function used to evaluate the actual displacements of
    the controlled DOFs.

:optparam:`linesearch{in}`
    A nonzero value turns on the line search algorithm.

:optparam:`lsearchtol{rn}`
    Line search tolerance.  Default ``0.8``.

:optparam:`lsearchamp{rn}`
    Line search amplification factor, which should lie in :math:`(1,10)`.

:optparam:`lsearchmaxeta{rn}`
    Maximum limit on the length of the iterative step; the allowed range is
    :math:`(1.5,15)`.

:optparam:`nccdg{in}`
    Number of DOF groups used for evaluating the convergence criteria.  Each
    DOF is tested for membership of a group, and if it is a member its
    contribution is taken into account when the criteria for that group are
    evaluated.  When :param:`nccdg` is omitted, one group containing all DOF
    types is created.

:optparam:`ccdg1{ia}` ... :optparam:`ccdgN{ia}`
    The DOF types of each convergence group, where the ``#`` in ``ccdg#`` is
    the group number, counted from 1.  The array holds ``DofIDItem`` values
    identifying the physical meaning of the DOFs in the group; see
    ``src/core/dofiditem.h``.

:param:`rtolv{rn}`
    Relative convergence norm, for both the iterative displacement change
    vector and the residual unbalanced force vector.

:optparam:`rtolf{rn}`, :optparam:`rtold{rn}`
    Independent relative convergence criteria for the unbalanced forces and
    for the iterative displacement change respectively.  With the default
    convergence criteria, :param:`rtolv`, :param:`rtolf` and :param:`rtold` are
    real values; when convergence criteria DOF groups are used, they must be
    real arrays of size :param:`nccdg`, each value giving the relative
    criterion for one DOF group.

:optparam:`initialGuess{in}`
    Default ``0``, for which the first iteration of each step starts from the
    previously converged state and applies the prescribed displacement
    increments.  That can produce very high strains in the elements connected
    to nodes with changing prescribed displacements, leaving the state far from
    equilibrium, which may lead to slow convergence and to strain localisation
    near the boundary.  With ``1``, the contribution of the prescribed
    displacement increments to the internal nodal forces is linearised and
    moved to the right-hand side, which often gives an initial solution closer
    to equilibrium.  If the step is in fact elastic, equilibrium is then fully
    restored after the second iteration, whereas the default method may need
    more iterations.

.. _calm-params:

Solver parameters for indirect control (CALM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Indirect control corresponds to :param:`controlmode`\ ``=0``, for which the CALM
solver is used.  The value of the reference load vector is determined by
:param:`refloadmode` at the first step of each meta-step.  You must ensure that
the same reference load vector can be obtained for every solution step of the
meta-step; this is necessary for restart and adaptivity to work.

.. record::

   :elemparam:`Psi{rn}` :elemparam:`MaxIter{in}`
   :elemparam:`stepLength{rn}` :elemparam:`rtolv{rn}`
   :optelemparam:`minStepLength{rn}` :optelemparam:`initialStepLength{rn}`
   :optelemparam:`forcedInitialStepLength{rn}`
   :optelemparam:`reqIterations{in}` :optelemparam:`maxrestarts{in}`
   :optelemparam:`minIter{in}` :optelemparam:`manrmsteps{in}`
   :optelemparam:`hpcmode{in}` :optelemparam:`hpc{ia}`
   :optelemparam:`hpcw{ra}` :optelemparam:`linesearch{in}`
   :optelemparam:`lsearchamp{rn}` :optelemparam:`lsearchmaxeta{rn}`
   :optelemparam:`lsearchtol{rn}` :optelemparam:`nccdg{in}`
   :optelemparam:`ccdg1{ia}` ``...`` :optelemparam:`ccdgN{ia}`
   :optelemparam:`rtolf{rn}` :optelemparam:`rtold{rn}`
   :optelemparam:`pert{ia}` :optelemparam:`pertw{ra}`
   :optelemparam:`rpa{rn}` :optelemparam:`rseed{in}`

**Parameters**

:param:`Psi{rn}`
    CALM :math:`\Psi` control parameter.  For :math:`\Psi = 0` displacement
    control is applied.  For nonzero values load control applies together with
    displacement control (ALM), and for large :math:`\Psi` load control
    dominates.

:param:`MaxIter{in}`
    Maximum number of iterations allowed to reach the equilibrium state.  When
    the limit is reached, a restart with a smaller step length follows.

:param:`stepLength{rn}`
    Maximum value of the arc length (step length).

:optparam:`minStepLength{rn}`
    Minimum step length; the step length never falls below it.  If convergence
    problems are encountered and the step length cannot be decreased further,
    the computation terminates.

:optparam:`initialStepLength{rn}`
    Initial arc length.  If not given, the maximum step length from
    :param:`stepLength` is used.

:optparam:`forcedInitialStepLength{rn}`
    When a simulation is restarted the last predicted step length is used;
    this parameter overrides it, and also overrides
    :param:`initialStepLength`.

:optparam:`reqIterations{in}`
    Approximate number of iterations, maintained by adjusting the step length.

:optparam:`maxrestarts{in}`
    Maximum number of restarts when convergence is not reached within
    :param:`MaxIter`.

:optparam:`minIter{in}`
    Minimum number of iterations always performed during the iterative
    solution.  :param:`reqIterations` is set to the same value, and
    :param:`MaxIter` is raised if it is lower.

:optparam:`manrmsteps{in}`
    Forces the use of the accelerated Newton-Raphson method, where the
    stiffness is updated after :param:`manrmsteps` steps.  By default the
    modified NR method is used, with no stiffness update.

:optparam:`hpcmode{in}`
    ALM mode.  ``0`` (default) full ALM with a quadratic constraint on all
    DOFs; ``1`` (the default once :param:`hpc` is used) full ALM with a
    quadratic constraint taking only the selected DOFs into account; ``2``
    constraint linearised in displacements only, taking only the selected DOFs
    with the given weights into account.

:optparam:`hpc{ia}`
    Hyper-plane control: the degrees of freedom taken into account in the ALM
    step length condition.  Important mainly for materially nonlinear problems
    with strong localisation.  With *N* such DOFs, the array has the form
    ``2*N dofman1 idof1 ... dofmanN idofN``.

:optparam:`hpcw{ra}`
    DOF weights in the linear constraint.  The DOF ordering follows
    :param:`hpc` and the array size should be *N*.

:optparam:`linesearch{in}`, :optparam:`lsearchtol{rn}`, :optparam:`lsearchamp{rn}`, :optparam:`lsearchmaxeta{rn}`
    Line search control; see the identically named parameters under
    :ref:`nrsolver-params`.

:optparam:`nccdg{in}`, :optparam:`ccdg1{ia}` ... :optparam:`ccdgN{ia}`
    Convergence criteria DOF groups; see :ref:`nrsolver-params`.

:param:`rtolv{rn}`
    Relative convergence norm, for both the iterative displacement change
    vector and the residual unbalanced force vector.

:optparam:`rtolf{rn}`, :optparam:`rtold{rn}`
    Independent relative convergence criteria for the unbalanced forces and
    for the iterative displacement change; see :ref:`nrsolver-params`.

:optparam:`pert{ia}`
    DOFs to be perturbed after the first iteration of each step.  With *M* such
    DOFs, the array has the form ``2*M dofman1 idof1 ... dofmanM idofM``.

:optparam:`pertw{ra}`
    DOF perturbations.  The DOF ordering follows :param:`pert` and the array
    size should be *M*.

:optparam:`rpa{rn}`
    Amplitude of the random perturbation applied to each DOF.

:optparam:`rseed{in}`
    Seed for the random generator producing the random perturbations.

Standard syntax
~~~~~~~~~~~~~~~

All parameters, for the analysis as well as for the solver, are supplied in the
analysis record.  A default meta-step is created for all the solution steps
required, and the meta-step attributes are then specified within the analysis
record.

.. record::

   :descitem:`NonLinearStatic` :elemparam:`nsteps{in}`
   :elemparam:`rtolv{rn}` :elemparam:`lstype{in}` :elemparam:`smtype{in}`
   :elemparam:`solverParams{...}` :optelemparam:`contextOutputStep{in}`
   :optelemparam:`controlmode{in}` :optelemparam:`deltat{rn}`
   :optelemparam:`stiffmode{in}` :optelemparam:`nonlinform{in}`
   ``<``\ :optelemparam:`nonlocstiff{in}`\ ``>``
   ``<``\ :optelemparam:`nonlocalext{}`\ ``>``
   ``<``\ :optelemparam:`loadbalancing{}`\ ``>``

The meaning of the parameters is the same as for the extended syntax.
:param:`lstype` selects the solver for the linear system of equations and
:param:`smtype` the sparse matrix storage scheme; the scheme must be compatible
with the solver type.  See :ref:`sparselinsolver` for details.

.. _AdaptiveLinearStatic:

Adaptive linear static
----------------------

Adaptive linear static analysis.  Multiple loading cases are not supported.
Because the problem is linear, a complete reanalysis from the beginning is done
after adaptive remeshing: after the first step the error is estimated,
information about the required density is generated through the mesher
interface, and the solution terminates.  If the error criterion is not
satisfied, a new mesh and a corresponding input file are generated and a new
analysis should be run, until the error is acceptable.  The error estimator
currently available for linear problems is Zienkiewicz-Zhu.

.. note::

   The adaptive framework needs specific functionality from the elements and
   material models.  See the Element Library Manual and the Material Library
   Manual.

.. record::

   :descitem:`Adaptlinearstatic` :elemparam:`nsteps{in}`
   :elemparam:`errorestimatorparams{...}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`meshpackage{in}`

**Parameters**

:param:`nsteps{in}`
    Number of loading cases; set it to ``1``.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`meshpackage{in}`
    Mesh package interface used to generate the required mesh density for
    remeshing.  The supported interfaces are described in :ref:`meshpackages`.
    The T3d interface is used by default.

:param:`errorestimatorparams{...}`
    Parameters of the Zienkiewicz-Zhu error estimator; see
    :ref:`errorestimators`.

.. _AdaptiveNonLinearStatic:

Adaptive nonlinear static
-------------------------

Adaptive non-linear static problem.  The solution proceeds as a series of
loading or displacement increments.  The error is estimated at the end of each
load increment, once equilibrium has been reached; depending on the error
reached, the computation either continues or generates new mesh densities and
stops, in which case a new discretization should be generated.

The truly adaptive approach is supported: the computation can be restarted from
the last step (see :ref:`running-the-code`), the solution is mapped onto the new
mesh in a separate solution step, and a new load increment is applied.  You may
of course also start the analysis from the beginning with the new mesh.  The
estimators and indicators currently available are the linear Zienkiewicz-Zhu
estimator and the scalar error indicator.

.. note::

   The adaptive framework needs specific functionality from the elements and
   material models.  See the Element Library Manual and the Material Library
   Manual.

.. record::

   :descitem:`Adaptnlinearstatic` :elemparam:`Nonlinearstaticparams{...}`
   :elemparam:`errorestimatorparams{...}` :optelemparam:`equilmc{in}`
   :optelemparam:`meshpackage{in}` :optelemparam:`eetype{in}`

**Parameters**

:param:`Nonlinearstaticparams{...}`
    The parameters related to the nonlinear analysis itself, described in
    :ref:`NonLinearStatic`.

:optparam:`equilmc{in}`
    Whether equilibrium is restored after the primary and internal variables
    have been mapped onto the new mesh, before the new load increment is
    applied.  ``0`` (default) does not restore equilibrium; ``1`` forces it to
    be restored before the new step is applied.

:optparam:`meshpackage{in}`
    Mesh package interface used to generate the required mesh density for
    remeshing.  The supported interfaces are described in :ref:`meshpackages`.
    The T3d interface is used by default.

:optparam:`eetype{in}`
    Type of error estimator or indicator to be used; see
    :ref:`errorestimators`.

:param:`errorestimatorparams{...}`
    The set of parameters belonging to the selected error estimator; see
    :ref:`errorestimators`.

.. _FreeWarping:

Free warping analysis
---------------------

Free warping analysis computes the deplanation function of a cross section of
arbitrary shape, by solving Laplace's equation with automatically generated
boundary conditions corresponding to the free warping problem.

This analysis supports only ``TrWarp`` elements and ``WarpingCS`` cross
sections.  One external node must be defined for each warping cross section.
The coordinates of this node are arbitrary, but it must be defined with
``DofIDMask 1 24`` and with one boundary condition representing the relative
twist acting on the corresponding warping cross section.  No additional loads
are meaningful in free warping analysis.

.. record::

   :descitem:`FreeWarping` :elemparam:`nsteps{in}`

**Parameters**

:param:`nsteps{in}`
    Number of loading cases.  The series of loading cases is maintained as a
    sequence of time steps: for each load case an auxiliary time step is
    generated whose time equals the load case number, and the load vector for
    each load case is formed at that auxiliary time.

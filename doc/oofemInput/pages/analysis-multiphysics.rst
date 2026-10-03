Multi-physics problems
======================

OOFEM solves a multi-physics problem in one of two ways, and this section covers
both.

A *partitioned* (staggered) problem is a sequence of sub-problems, each a
complete analysis with its own input file, run one after another in every time
step; the fields one sub-problem computes are exported and picked up by the
next.  :ref:`staggeredproblem`, :ref:`fluidstructureproblem` and
:ref:`AdditiveManufacturingProblem` work this way.

A *monolithic* problem is stated as a single weak form, assembled and solved as
one system.  :ref:`StationaryMPMSProblem` and :ref:`NonStationaryMPMSProblem`
work this way, with the weak form written symbolically in the input file out of
variables, terms and integrals.

.. _staggeredproblem:

Staggered problem
-----------------

.. note::

   This problem type is part of the **transport module** and can be used only
   when that module is configured.

Represents a so-called staggered analysis: a sequence of sub-problems in which
the result of one sub-problem may depend on the results of those before it.  The
typical example is a heat transfer analysis followed by a mechanical analysis
that takes the computed temperature field into account.  A similar analysis
couples moisture transport with the drying strain of concrete.

The current implementation supports a sequence of two sub-problems.  The
sub-problems are described by their own input files, whose syntax is the same as
for a standalone problem; the only addition is that a sub-problem must *export*
its solution fields so that they become available to the subsequent
sub-problem — see :ref:`exportfields`.

.. note::

   All material models derived from the ``StructuralMaterial`` base class take
   an externally registered temperature field into account when one is
   provided.

.. record::

   :descitem:`StaggeredProblem` :elemparam:`prob1{s}` :elemparam:`prob2{s}`
   ((:elemparam:`nsteps{in}` :elemparam:`deltaT{rn}`) |
   :elemparam:`timeDefinedByProb{in}`)
   :optelemparam:`stepMultiplier{rn}`

**Parameters**

:param:`prob1{s}`, :param:`prob2{s}`
    Paths to the sub-problem input files.  :param:`prob1` runs first in each
    solution step, :param:`prob2` second.

:param:`timeDefinedByProb{in}`
    Take the time sequence from the given sub-problem.  The sub-problem may
    then specify arbitrary loading steps, which allows a great deal of
    flexibility.  Use either this parameter or the :param:`nsteps` /
    :param:`deltaT` pair below.

:param:`nsteps{in}`, :param:`deltaT{rn}`
    Let the staggered problem control the time instead.  Any time-stepping
    parameters of the sub-problems are then ignored — even where the
    sub-problem input syntax requires them — and only these parameters matter.
    :param:`nsteps` is the number of time steps to solve and :param:`deltaT`
    the time step length used for integration.

:optparam:`stepMultiplier{rn}`
    Multiplies all times by this constant.  Default ``1``.

.. _fluidstructureproblem:

FluidStructure problem
----------------------

.. note::

   This problem type is part of the **PFEM module** and can be used only when
   that module is configured.

.. warning::

   The record keyword is ``fluidstuctureproblem`` — the spelling in the OOFEM
   sources is missing an "r", and the input file must match it.

Represents a fluid-structure analysis based on the staggered problem, but with
iterative synchronization of the sub-problems.  The implementation uses the PFEM
model (:ref:`pfemIncomp`) for the fluid part, and for the structural part a full
dynamic analysis by implicit direct integration (:ref:`DIIDynamic`).

The coupling of the two phases enforces compatibility on the interface.  Special
fluid particles are attached to every structural node on the interface that the
fluid can reach.  These particles have no degrees of freedom of their own, so no
equations are solved on them; their movement is determined entirely by the
associated structural nodes, and their velocities — governed by the solid
part — affect the fluid equations naturally.

The iterative procedure follows the Dirichlet-Neumann approach: prescribed
velocities on the fluid side of the interface are the Dirichlet conditions,
while applied forces on the structural side are the Neumann conditions.  The
convergence criterion is based on the difference of the pressure and velocity
values on the interface between successive iterations; once these are smaller
than the prescribed tolerance the iteration stops and the solution proceeds to
the next step.

.. record::

   :descitem:`fluidstuctureproblem` :elemparam:`nsteps{in}`
   :elemparam:`deltaT{rn}` :elemparam:`prob1{s}` :elemparam:`prob2{s}`
   :optelemparam:`maxiter{in}` :optelemparam:`rtolv{rn}`
   :optelemparam:`rtolp{rn}`

**Parameters**

:param:`nsteps{in}`, :param:`deltaT{rn}`
    Number of time steps and the time step length.

:param:`prob1{s}`, :param:`prob2{s}`
    Paths to the sub-problem input files.  :param:`prob1` runs first in each
    solution step, :param:`prob2` second.

:optparam:`maxiter{in}`
    Maximum number of iterations.  Default ``50``.

:optparam:`rtolv{rn}`, :optparam:`rtolp{rn}`
    Relative norms for the velocity and pressure difference between two
    successive iterations.  Defaults ``1.e-3``.

.. _AdditiveManufacturingProblem:

Additive manufacturing problem
------------------------------

Represents an additive manufacturing analysis, performed as a sequence of
sub-problems in each time step over a time-evolving domain defined by a G-code
file.

The time-evolving computational domain is defined by the G-code file
(:param:`gcode`), which is interpreted by a built-in G-code processor.  The
processor emulates the machine G-code precisely in a virtual domain represented
by a structured grid, with user-defined resolution in space
(:param:`stepx`, :param:`stepy`, :param:`stepz`) and in time.  It tracks the
evolution of the Representative Volume of Deposited Material (RVOM) in the
individual elements of the grid, taking into account the realistic extruder
path, speed and deposition profile defined by the G-code.

.. note::

   G-code units are millimetres and mm/s, but they are converted internally to
   the basic units of metres and m/s.

Initially the RVOM of every element is zero, so all elements are empty and
defined as inactive.  In each step the updated computational domain is
identified as the set of elements whose RVOM exceeds the threshold given by
:param:`minvof`.  Active nodes are then the nodes shared by at least one active
element.  The boundary conditions — ambient temperature and prescribed bottom
temperature — are generated by the G-code processor.

See the ``tests/am`` directory for example input files.

.. record::

   :descitem:`AdditiveManufacturingProblem` :elemparam:`gcode{s}`
   :elemparam:`prob1{s}` :elemparam:`prob2{s}` :elemparam:`stepx{rn}`
   :elemparam:`stepy{rn}` :elemparam:`stepz{rn}`
   (:elemparam:`nsteps{in}` :elemparam:`deltaT{rn}` |
   :elemparam:`timeDefinedByProb{in}` | :elemparam:`dtf{in}` |
   :elemparam:`prescribedtimes{ra}`)
   :optelemparam:`prob3{s}` :optelemparam:`coupling{ia}`
   :optelemparam:`minvof{rn}` :optelemparam:`skipsm{in}`
   :optelemparam:`maxcommands{in}` :optelemparam:`stepMultiplier{rn}`
   :optelemparam:`adaptivesteplength{}` :optelemparam:`minsteplength{rn}`
   :optelemparam:`maxsteplength{rn}` :optelemparam:`reqIterations{in}`
   :optelemparam:`endoftimeofinterest{rn}`
   :optelemparam:`adaptivestepsince{rn}` :optelemparam:`layerheight{rn}`
   :optelemparam:`extrusionwidth{rn}`
   :optelemparam:`chambertemperature{rn}`
   :optelemparam:`depositiontemperature{rn}`
   :optelemparam:`heatbedtemperature{rn}`
   :optelemparam:`heattransferfilmcoefficient{rn}`
   :optelemparam:`depositedmaterialheatpower{rn}`

**Parameters**

:param:`gcode{s}`
    Path to the G-code file defining the time-evolving domain.

:param:`prob1{s}`, :param:`prob2{s}`
    Paths to the sub-problem input files, run in a staggered manner for each
    solution step.  :param:`prob1` is typically a transient heat analysis,
    followed by the mechanical analysis in :param:`prob2`.

:optparam:`prob3{s}`
    Path to an optional third sub-problem input file.

:optparam:`coupling{ia}`
    Numbers of the coupled sub-problems, used to resolve which sub-problem
    supplies a field to which.

:param:`stepx{rn}`, :param:`stepy{rn}`, :param:`stepz{rn}`
    Spatial resolution of the structured grid representing the virtual domain.

:param:`nsteps{in}`, :param:`deltaT{rn}`
    Number of time steps and the time step length.

:param:`timeDefinedByProb{in}`
    Take the time sequence from the given sub-problem instead.

:param:`dtf{in}`
    Number of a time function determining the time step length instead.

:param:`prescribedtimes{ra}`
    Explicit array of solution times instead.

:optparam:`minvof{rn}`
    Threshold RVOM value above which an element counts as active and becomes
    part of the computational domain.

:optparam:`skipsm{in}`
    Skip the structural (mechanical) sub-problem, running the thermal analysis
    alone.

:optparam:`maxcommands{in}`
    Maximum number of G-code commands to process, which limits how far into the
    print the emulation runs.

:optparam:`stepMultiplier{rn}`
    Multiplies all times by this constant.  Default ``1``.

Adaptive step length
~~~~~~~~~~~~~~~~~~~~

:optparam:`adaptivesteplength{}`
    A flag; when present, the time step length is adapted during the analysis.
    The remaining parameters in this group are read only when it is given.

:optparam:`minsteplength{rn}`
    Minimum time step length.  Default ``0.``.

:optparam:`maxsteplength{rn}`
    Maximum time step length.

:optparam:`reqIterations{in}`
    Approximate number of iterations, maintained by adjusting the step length.
    Default ``1``.

:optparam:`endoftimeofinterest{rn}`
    Time at which the analysis stops.  Default ``1.e32``.

:optparam:`adaptivestepsince{rn}`
    Time from which the step length starts being adapted.  Default ``0.``.

Printer process parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

:optparam:`layerheight{rn}`
    Layer height, in millimetres.

:optparam:`extrusionwidth{rn}`
    Extrusion width, in millimetres.

:optparam:`chambertemperature{rn}`
    Ambient temperature in the printer chamber.

:optparam:`depositiontemperature{rn}`
    Temperature of the deposited material.

:optparam:`heatbedtemperature{rn}`
    Prescribed temperature of the heat bed.

:optparam:`heattransferfilmcoefficient{rn}`
    Heat transfer film coefficient between the deposited material and the air
    surrounding it in the printer chamber.

:optparam:`depositedmaterialheatpower{rn}`
    Heat power of the deposited material, defined as the product of the
    material's specific heat and its density.

.. _StationaryMPMSProblem:

Stationary symbolic MPM problem
-------------------------------

The stationary symbolic multiphysics problem lets the problem be formulated
symbolically, by giving its weak form in terms of variables, terms and
integrals.

.. warning::

   The record keyword is ``mpmsymbolicstationaryproblem``, not the class name
   ``StationaryMPMSProblem``.

The problem is defined by its weak form as a set of integrals of specific,
predefined terms, which are themselves functions of the defined variables.  The
resulting nonlinear system is solved with the Newton-Raphson solver
(``NRSolver``).

.. record::

   :descitem:`mpmsymbolicstationaryproblem` :elemparam:`nsteps{in}`
   :elemparam:`nvariables{in}` :elemparam:`nterms{in}`
   :elemparam:`nintegrals{in}` :elemparam:`lhsterms{ia}`
   :elemparam:`rhsterms{ia}` :optelemparam:`smtype{in}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.

:param:`nvariables{in}`, :param:`nterms{in}`, :param:`nintegrals{in}`
    Numbers of Variable, Term and Integral records that follow.  See
    :ref:`VariablesSec`, :ref:`TermsSec` and :ref:`IntegralsSec`.

:param:`lhsterms{ia}`
    Numbers of the integrals contributing to the left-hand side of the
    nonlinear problem.

:param:`rhsterms{ia}`
    Numbers of the integrals contributing to the right-hand side.

:optparam:`smtype{in}`
    Sparse matrix storage scheme.  Default ``0``, the symmetric skyline.  See
    :ref:`sparselinsolver`.

**Example**

See ``tests/regression/mpm/mpms06.in`` for a complete input file, and
`the commented up-formulation example
<https://oofem.github.io/blog/mpm-incompressible-elasticity-up-formulation/>`_
on the OOFEM blog.

.. _NonStationaryMPMSProblem:

Non-stationary symbolic MPM problem
-----------------------------------

The transient counterpart of the stationary symbolic MPM problem.  Time is
discretized with the generalized midpoint rule using the parameter
:param:`alpha`, as in :ref:`TransientTransport`.

Unlike the stationary problem, this model has three formulations, selected by
:param:`ptype`.  Only ``symbolic`` reads the integral lists; the other two are
built-in formulations retained for compatibility.

.. note::

   The keyword is ``mpmnonstationaryproblem``.  The alias ``mpmproblem`` is
   also accepted, for compatibility with older input files.

.. record::

   :descitem:`mpmnonstationaryproblem` :elemparam:`nsteps{in}`
   :elemparam:`alpha{rn}`
   (:elemparam:`deltat{rn}` | :elemparam:`deltatfunction{in}` |
   :elemparam:`prescribedtimes{ra}`)
   :optelemparam:`ptype{s}` :optelemparam:`nvariables{in}`
   :optelemparam:`nterms{in}` :optelemparam:`nintegrals{in}`
   :optelemparam:`lhsterms{ia}` :optelemparam:`lhsdotterms{ia}`
   :optelemparam:`rhsterms{ia}` :optelemparam:`initt{rn}`
   :optelemparam:`keeptangent{}` :optelemparam:`exportfields{ia}`
   :optelemparam:`smtype{in}` :optelemparam:`debug{in}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.

:param:`alpha{rn}`
    The :math:`\alpha` parameter of the generalized midpoint rule.

:param:`deltat{rn}`
    Time step length.  Exactly one of :param:`deltat`,
    :param:`deltatfunction` or :param:`prescribedtimes` must be given; the
    analysis stops with "Time step not defined" if none is.

:param:`deltatfunction{in}`
    Number of a time function returning the time step length.

:param:`prescribedtimes{ra}`
    Explicit array of solution times.

:optparam:`ptype{s}`
    Problem formulation: ``up`` (the default, for compatibility) for the
    displacement-pressure formulation, ``tm`` for the thermo-mechanical one, or
    ``symbolic`` for the symbolic weak form.  Any other value is rejected.

:optparam:`lhsterms{ia}`, :optparam:`lhsdotterms{ia}`, :optparam:`rhsterms{ia}`
    Numbers of the integrals contributing to the left-hand side, to the
    left-hand side time-derivative term, and to the right-hand side.  All three
    are **mandatory** when :param:`ptype` is ``symbolic``, and unused
    otherwise.

:optparam:`nvariables{in}`, :optparam:`nterms{in}`, :optparam:`nintegrals{in}`
    Numbers of Variable, Term and Integral records that follow; needed for the
    ``symbolic`` formulation.  See :ref:`VariablesSec`, :ref:`TermsSec` and
    :ref:`IntegralsSec`.

:optparam:`initt{rn}`
    Initial time for the integration.

:optparam:`keeptangent{}`
    A flag; when present, the tangent is assembled once and reused in every
    step.

:optparam:`exportfields{ia}`
    Fields to export to a subsequent analysis; see :ref:`exportfields`.
    Displacement and pressure fields are supported.

:optparam:`smtype{in}`
    Sparse matrix storage scheme.  Default ``0``, the symmetric skyline.

:optparam:`debug{in}`
    A nonzero value turns on additional diagnostic output.

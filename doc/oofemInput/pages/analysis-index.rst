.. _AnalysisRecord:

Analysis record
===============

The analysis record describes the type of analysis to be performed.  Its record
keyword is the name of an *engineering model*; the pages listed at the bottom of
this page document one group of models each.

The analysis record may be followed by optional *meta-step* records.  Certain
attributes that would otherwise appear in the analysis record can then be given
independently for each meta-step.  Such attributes are marked in this manual by
a superscript "M", as in ``keyword``\ :sup:`M`.

.. note::

   Record keywords are matched case-insensitively (see
   :ref:`anatomy-of-an-input-file`), but the spelling must be exact.  This
   manual gives the keyword each engineering model actually registers, which is
   not always its class name — ``incrlinearstatic`` rather than
   ``IncrementalLinearStatic``, for instance.

Standard syntax
---------------

All parameters, for the analysis as well as for the solver, are given in the
analysis record itself.  A single default meta-step is created covering all
solution steps.

.. record::

   :descitem:`EngngModelType` :elemparam:`nsteps{in}`
   :optelemparam:`renumber{in}` :optelemparam:`profileopt{in}`
   :elemparam:`attributes{s}` :optelemparam:`ninitmodules{in}`
   :optelemparam:`nmodules{in}` :optelemparam:`nxfemman{in}`
   :optelemparam:`nvariables{in}` :optelemparam:`nterms{in}`
   :optelemparam:`nintegrals{in}`

Meta-step syntax
----------------

Not every analysis type supports meta-steps; where a page does not say
otherwise, only the standard syntax is available.

.. record::

   :descitem:`EngngModelType` :elemparam:`nmsteps{in}`
   :optelemparam:`ninitmodules{in}` :optelemparam:`nmodules{in}`
   :optelemparam:`nxfemman{in}` :optelemparam:`nvariables{in}`
   :optelemparam:`nterms{in}` :optelemparam:`nintegrals{in}`

immediately followed by :param:`nmsteps` meta-step records, each of the form

.. record:: Meta-step record

   :elemparam:`nsteps{in}` :elemparam:`attributes{s}`

A meta-step represents a sequence of solution steps sharing common attributes.
Exactly :param:`nmsteps` meta-step records are expected.

Parameters
----------

These parameters apply to the analysis record in standard syntax, and to the
meta-step records in meta-step syntax.

:param:`nsteps{in}`
    Number of solution steps.  Within a meta-step, the number of solution steps
    belonging to that meta-step.

:param:`nmsteps{in}`
    Number of meta-steps.

:optparam:`renumber{in}`
    Turns on renumbering after each time step.  Necessary when Dirichlet
    boundary conditions change during the simulation.  Equivalent to the
    ``-rn`` command line flag.

:optparam:`profileopt{in}`
    A nonzero value turns on equation renumbering to optimise the profile of
    the characteristic matrix, using Sloan's algorithm.  Profile optimisation is
    off by default, and does not work in parallel mode.

:param:`attributes{s}`
    The analysis and solver attributes valid for the solution steps of the
    corresponding meta-step.  In standard syntax these attributes are valid for
    all solution steps.

:optparam:`ninitmodules{in}`
    Number of initialization module records.  Initialization modules let the
    state variables be initialised from values computed previously by external
    software.  Their records follow the meta-step section, or the analysis
    record if there are no meta-steps.  See :ref:`InitModulesSec`.

:optparam:`nmodules{in}`
    Number of export module records.  Export modules write computed data for
    post-processing by external software.  Their records follow the
    initialization module records.  See :ref:`ExportModulesSec`.

:optparam:`nvariables{in}`, :optparam:`nterms{in}`, :optparam:`nintegrals{in}`
    Number of Variables, Terms and Integrals appearing in the weak form, for
    the symbolic computations provided by the ``mpm`` module.  The corresponding
    Variable, Term and Integral records follow the export module records.  See
    :ref:`VariablesSec`, :ref:`TermsSec` and :ref:`IntegralsSec`.

:optparam:`nxfemman{in}`
    ``1`` creates an XFEM manager, ``0`` (the default) does not.  The XFEM
    manager holds the list of enrichment items; see :ref:`XFEMManagerRecords`.

:optparam:`eetype{in}`
    Error estimator type for the problem.  Used for adaptive analysis, but it
    can also be used to compute error estimates and write them to the output
    files.  See :ref:`errorestimators` and the adaptive engineering models.

.. _timestepping:

Time stepping and adaptive time step
------------------------------------

The attributes below are read by the meta-step record, and therefore apply to
every engineering model that uses meta-steps — in standard syntax, to the
single default meta-step created from the analysis record.  A model may
additionally read :param:`deltat` or :param:`prescribedtimes` itself, in which
case its own section says so.

:optparam:`deltat{rn}`
    Length of the time step.  Default ``1.0``.

:optparam:`finalt{rn}`
    Final time of the meta-step.  When given and nonzero, it determines how far
    the meta-step runs and :param:`nsteps` is not read; it must not be
    negative.  When omitted, the final time is derived from :param:`nsteps`
    together with :param:`deltat`, :param:`prescribedtimes` or
    :param:`dtfunction`.

:optparam:`prescribedtimes{ra}`
    Explicit array of solution times.  The number of solution steps is then the
    size of this array.

:optparam:`dtfunction{in}`
    Number of a time function (see :ref:`TimeFunctionsRecords`) that returns
    the length of the time step for a given solution step.  A variable time
    step is useful when a long time interval has to be covered.  Takes
    precedence over :param:`prescribedtimes`.

:optparam:`treductiontype{s}`
    Name of the time-step reduction strategy, which controls what happens when
    a step fails to converge.  Default ``noreduction``.

    ``noreduction``
        No reduction: a step that fails to converge terminates the analysis.
        This is the default, so adaptive time stepping is **off** unless
        :param:`treductiontype` is given.

    ``simplereduction``
        Reduce the step by a fixed factor on failure, and adjust it according
        to the number of iterations the previous step needed.  It reads the
        additional attributes listed below from the same record.

Attributes of the ``simplereduction`` strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are read from the same record as :param:`treductiontype`.  Note that
:param:`dtmax` and :param:`dtmin` are **mandatory** once ``simplereduction`` is
selected.

:elemparam:`dtmax{rn}`
    Largest permitted time step.  The step is never enlarged beyond this value.

:elemparam:`dtmin{rn}`
    Smallest permitted time step.  The step is never reduced below this value.

:optparam:`ncrf{rn}`
    No-convergence reduction factor.  When a step fails to converge, the time
    step is multiplied by this factor and the step is attempted again.
    Default ``0.5``.

:optparam:`nmaxrestarts{in}`
    Maximum number of times a step may be restarted with a reduced time step
    before the analysis gives up.  Default ``5``.

The number of iterations the solver needed in the previous step decides whether
the next time step is enlarged or reduced.  Within the interval the step is
kept; below it the step is doubled, above it the step is halved.

:optparam:`minrequirediter{in}`
    If convergence was reached in fewer iterations than this, the step was
    "too easy" and the time step is multiplied by ``2.0``.  Default ``3``.

:optparam:`maxrequirediter{in}`
    If convergence needed more iterations than this, the step was too large and
    the time step is multiplied by ``0.5``.  Default ``8``.

.. _common-analysis-parameters:

Parameters common to all engineering models
-------------------------------------------

The following optional attributes are read by the ``EngngModel`` base class, so
every analysis type accepts them even where its own section does not repeat
them.

:optparam:`contextoutputstep{in}`
    Create a context file every :param:`contextoutputstep`-th step, and
    whenever one is otherwise needed.  Useful for restarts and for
    post-processing.

:optparam:`suppress_output`
    A flag; when present, writing to the ``.out`` file is suppressed.  Useful
    when the results are taken from an export module instead and the text
    output would only cost time and disk space.

:optparam:`initialguess{in}`
    How the first iteration of each step is started.  ``0`` (``IG_None``)
    continues from the previously converged state; ``1`` (``IG_Tangent``)
    solves an approximated tangent problem from the last iteration, which helps
    when Dirichlet boundary conditions change.  The default differs between
    engineering models and is stated in each section.

:optparam:`lstype{in}`, :optparam:`smtype{in}`
    Linear solver type and sparse matrix storage scheme.  These are the
    individual members of the ``sparselinsolverparams`` group; see
    :ref:`sparselinsolver` for the permitted values and combinations.

:optparam:`nonlinform{in}`
    Formulation of a non-linear problem.  ``1`` (default) uses the total
    Lagrangian formulation on the undeformed original shape (first-order
    theory).  ``2`` adds the equilibrated displacements to the original ones
    and updates them in each time step (second-order theory).

:optparam:`parallelflag{in}`
    Nonzero forces the model to run in parallel mode even without the ``-p``
    command line flag.  Parallel build only.

:optparam:`lbflag{in}`, :optparam:`forcelb1{in}`
    Dynamic load balancing.  See :ref:`dynamicloadbalancing` for these and the
    remaining load balancing attributes.

Supported analysis types
------------------------

.. toctree::
    :maxdepth: 2

    analysis-structural
    analysis-transport
    analysis-fluid
    analysis-multiphysics
    analysis-other

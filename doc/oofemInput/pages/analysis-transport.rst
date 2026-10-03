Transport problems
==================

.. note::

   These problem types require the **transport module** and can be used only
   when that module is included in your OOFEM configuration.

.. _StationaryTransport:

Stationary transport problem
----------------------------

Stationary transport problem.  A series of loading cases is maintained as a
sequence of time steps: for each load case an auxiliary time step is generated
whose time equals the load case number, and the load vector for each load case
is formed at that auxiliary time.

.. record::

   :descitem:`StationaryProblem` :elemparam:`nsteps{in}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`exportfields{ia}` :optelemparam:`keeptangent{}`

**Parameters**

:param:`nsteps{in}`
    Number of load cases.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`exportfields{ia}`
    See :ref:`exportfields` below.

:optparam:`keeptangent{}`
    A flag; when present, the tangent is assembled once and reused in every
    step, which saves time when it does not in fact change.

.. _exportfields:

Exporting fields to a subsequent analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the problem is used within a staggered analysis, the field it computes can
be exported and made available to subsequent analyses.  A temperature field
computed here, for instance, can then be taken into account by a subsequent
mechanical analysis.

To make that possible the field must be "exported", which is what the
:param:`exportfields` array does: it holds the field identifiers under which
the problem registers its primary unknowns.  The identifiers are those of the
``FieldType`` enumeration in ``src/core/field.h``.  Subsequent analyses can then
access the exported fields, if they support the feature.

.. _TransientTransport:

Transient transport problem
---------------------------

**Nonlinear** implicit integration scheme for transient transport problems.
This is the recommended transient transport solver.

The generalized midpoint rule, sometimes called the :math:`\alpha`-method, is
used for time discretization, with :math:`0\le\alpha\le1`:

* :math:`\alpha=0` gives the explicit Euler forward method;
* :math:`\alpha=0.5` recovers the implicit trapezoidal rule, which is
  unconditionally stable and second-order accurate in :math:`\Delta t`;
* :math:`\alpha=1.0` gives the implicit Euler backward method, which is
  unconditionally stable and first-order accurate in :math:`\Delta t`.

This problem supports sets and changes in the number of equations, so Dirichlet
boundary conditions may be imposed and removed during the solution.

.. record::

   :descitem:`TransientTransport` :elemparam:`nsteps{in}`
   :elemparam:`alpha{rn}`
   (:elemparam:`deltaT{rn}` | :elemparam:`dTfunction{in}` |
   :elemparam:`prescribedtimes{ra}`)
   :optelemparam:`initT{rn}` :optelemparam:`lumped{}`
   :optelemparam:`keeptangent{}` :optelemparam:`exportfields{ia}`

**Parameters**

:param:`nsteps{in}`
    Number of time steps to solve.

:param:`alpha{rn}`
    The :math:`\alpha` parameter of the generalized midpoint rule, as above.

:param:`deltaT{rn}`
    Time step length used for integration.  Give exactly one of
    :param:`deltaT`, :param:`dTfunction` or :param:`prescribedtimes`.

:param:`dTfunction{in}`
    Number of a time function (see :ref:`TimeFunctionsRecords`) returning the
    length of the time step.  A variable time step is advantageous when a large
    time interval has to be covered.

:param:`prescribedtimes{ra}`
    Explicit array of solution times.

:optparam:`initT{rn}`
    Initial time for the integration.  Default ``0.``.

:optparam:`lumped{}`
    A flag; when present, the numerical algorithm is stabilised using a lumped
    capacity matrix, which reduces the initial oscillations.

:optparam:`keeptangent{}`
    A flag; when present, the tangent is assembled once and reused in every
    step.

:optparam:`exportfields{ia}`
    See :ref:`exportfields`.

.. _LinearTransientTransport:

Transient transport problem — linear case (obsolete)
----------------------------------------------------

.. warning::

   Obsolete.  Use :ref:`TransientTransport` instead; the nonlinear transport
   solver is strongly preferred for reasons of stability.

**Linear** implicit integration scheme for transient transport problems, using
the same generalized midpoint rule as :ref:`TransientTransport`.

.. record::

   :descitem:`NonStationaryProblem` :elemparam:`nsteps{in}`
   :elemparam:`alpha{rn}`
   (:elemparam:`deltaT{rn}` | :elemparam:`deltaTfunction{in}` |
   :elemparam:`prescribedtimes{ra}`)
   :optelemparam:`initT{rn}` :optelemparam:`lumpedcapa{}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`exportfields{ia}`
   :optelemparam:`changingProblemSize{}`

**Parameters**

:param:`nsteps{in}`
    Number of time steps to solve.

:param:`alpha{rn}`
    The :math:`\alpha` parameter of the generalized midpoint rule.

:param:`deltaT{rn}`
    Time step length used for integration.

:param:`deltaTfunction{in}`
    Number of a time function (see :ref:`TimeFunctionsRecords`) returning the
    length of the time step.

:param:`prescribedtimes{ra}`
    Explicit array of solution times.

:optparam:`initT{rn}`
    Initial time for the integration.  Default ``0``.

:optparam:`lumpedcapa{}`
    A flag; when present, the numerical algorithm is stabilised using a lumped
    capacity matrix, which reduces the initial oscillations.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`exportfields{ia}`
    See :ref:`exportfields`.

:optparam:`changingProblemSize{}`
    A flag enabling changes in the number of equations, so that Dirichlet
    boundary conditions can be imposed and removed during the solution.  It
    works by storing the solution values directly on the nodes (DOFs).  If the
    problem does not grow or shrink during the solution it is more efficient to
    use the conventional strategy and omit this flag.

.. _TransientTransport2:

Transient transport problem — nonlinear case (obsolete)
-------------------------------------------------------

.. warning::

   Obsolete.  Use :ref:`TransientTransport` instead.

Implicit integration scheme for transient transport problems, using the same
generalized midpoint rule as :ref:`TransientTransport`.  See the Material
Library Manual for the solution algorithm.

.. record::

   :descitem:`NlTransientTransportProblem` :elemparam:`nsteps{in}`
   :elemparam:`alpha{rn}` :elemparam:`rtol{rn}`
   (:elemparam:`deltaT{rn}` | :elemparam:`deltaTfunction{in}`)
   :optelemparam:`initT{rn}` :optelemparam:`lumpedcapa{}`
   :optelemparam:`nsmax{in}` :optelemparam:`manrmsteps{in}`
   :optelemparam:`sparselinsolverparams{...}`
   :optelemparam:`exportfields{ia}`
   :optelemparam:`changingProblemSize{}`

**Parameters**

:param:`nsteps{in}`
    Number of time steps to solve.

:param:`alpha{rn}`
    The :math:`\alpha` parameter of the generalized midpoint rule.

:param:`rtol{rn}`
    Relative convergence tolerance.  The norms of the residual physical
    quantity — heat or mass — described by the solution vector, and of the
    change of the solution vector, are evaluated in each iteration;
    convergence is reached when both are smaller than :param:`rtol`.

:param:`deltaT{rn}`
    Time step length used for integration.

:param:`deltaTfunction{in}`
    Number of a time function returning the length of the time step; see
    :ref:`LinearTransientTransport`.

:optparam:`initT{rn}`
    Initial time for the integration; see :ref:`LinearTransientTransport`.

:optparam:`lumpedcapa{}`
    A flag; when present, the numerical algorithm is stabilised using a lumped
    capacity matrix, which reduces the initial oscillations.

:optparam:`nsmax{in}`
    Maximum number of iterations per time step.  Default ``30``.

:optparam:`manrmsteps{in}`
    Default ``0``, for which the modified Newton-Raphson scheme is used and the
    left-hand side matrix is never updated.  A nonzero value switches to the
    accelerated Newton-Raphson scheme, updating the left-hand side matrix after
    :param:`manrmsteps` steps.

:optparam:`sparselinsolverparams{...}`
    Attributes of the sparse linear solver; see :ref:`sparselinsolver`.

:optparam:`exportfields{ia}`
    See :ref:`exportfields`.

:optparam:`changingProblemSize{}`
    See :ref:`LinearTransientTransport`.

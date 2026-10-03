Fluid dynamic problems
======================

.. _nondimensional-scaling:

Non-dimensional form
--------------------

The characteristic equations of the CBS and SUPG/PSPG algorithms can be solved
in non-dimensional form.  To enable this, give :param:`scaleflag` a nonzero
value and supply the three scales:

:optparam:`scaleflag{in}`
    A nonzero value turns on the non-dimensional form.

:param:`lscale{rn}`, :param:`uscale{rn}`, :param:`dscale{rn}`
    Typical length, velocity and density scales.  Required when
    :param:`scaleflag` is nonzero.

.. _cbsIncomp:

Transient incompressible flow — CBS algorithm
---------------------------------------------

Solves transient incompressible flow with an algorithm based on the
Characteristics Based Split (CBS); for reference see O. C. Zienkiewicz and
R. L. Taylor, *The Finite Element Method*, volume 3, Butterworth-Heinemann,
2000.

At present only the semi-implicit form of the algorithm is available, and the
energy equation — and with it the temperature field — is not solved.

.. record::

   :descitem:`CBS` :elemparam:`nsteps{in}` :elemparam:`deltaT{rn}`
   :optelemparam:`theta1{rn}` :optelemparam:`theta2{rn}`
   :optelemparam:`cmflag{in}` :optelemparam:`scaleflag{in}`
   :optelemparam:`lscale{rn}` :optelemparam:`uscale{rn}`
   :optelemparam:`dscale{rn}` :optelemparam:`lstype{in}`
   :optelemparam:`smtype{in}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.

:param:`deltaT{rn}`
    Time step length used for integration.  It is adjusted automatically, if
    necessary, to satisfy the stability limits
    :math:`\Delta t \le h/\vert\boldsymbol{u}\vert` and
    :math:`\Delta t \le h^2/(2\nu)`.

:optparam:`theta1{rn}`, :optparam:`theta2{rn}`
    Integration constants,
    :math:`\theta_1, \theta_2 \in \langle\tfrac12, 1\rangle`.

:optparam:`cmflag{in}`
    A nonzero value uses the consistent mass matrix instead of the default
    lumped one.

:optparam:`scaleflag{in}`, :optparam:`lscale{rn}`, :optparam:`uscale{rn}`, :optparam:`dscale{rn}`
    See :ref:`nondimensional-scaling`.

:optparam:`lstype{in}`, :optparam:`smtype{in}`
    Linear solver type and sparse matrix storage scheme; the scheme must be
    compatible with the solver type.  See :ref:`sparselinsolver`.

.. _supgIncomp:

Transient incompressible flow — SUPG/PSPG algorithm
---------------------------------------------------

Solves transient incompressible flow using a stabilized formulation based on
SUPG and PSPG stabilization terms.  The stabilization provides stability and
accuracy for advection-dominated problems and permits equal-order interpolation
for velocity and pressure.  It also improves the convergence rate substantially
when large nonlinear systems are solved iteratively.

Varying :math:`\alpha` selects different methods from the generalized midpoint
family: forward Euler (:math:`\alpha=0`), the midpoint rule
(:math:`\alpha=0.5`), Galerkin (:math:`\alpha=2/3`) and backward Euler
(:math:`\alpha=1`).  Except the first one, all are implicit and require a matrix
inversion.  Energy-method analysis suggests unconditional stability for
:math:`\alpha\ge 0.5`; as far as accuracy goes, the midpoint rule is generally
to be preferred.

.. record::

   :descitem:`SUPG` :elemparam:`nsteps{in}` :elemparam:`deltaT{rn}`
   :elemparam:`rtolv{rn}` :optelemparam:`atolv{rn}`
   :optelemparam:`stopmaxiter{in}` :optelemparam:`alpha{rn}`
   :optelemparam:`cmflag{in}` :optelemparam:`deltatltf{in}`
   :optelemparam:`miflag{in}` :optelemparam:`scaleflag{in}`
   :optelemparam:`lscale{rn}` :optelemparam:`uscale{rn}`
   :optelemparam:`dscale{rn}` :optelemparam:`lstype{in}`
   :optelemparam:`smtype{in}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.

:param:`deltaT{rn}`
    Time step length used for integration.

:optparam:`deltatltf{in}`
    Number of a time function used instead of :param:`deltaT` to determine the
    time step length; its value evaluated for the solution step number yields
    the step length.

:optparam:`alpha{rn}`
    The :math:`\alpha` parameter of the generalized midpoint family, as above.

:param:`rtolv{rn}`, :optparam:`atolv{rn}`
    Relative and absolute error norms for the residual vector.  The equilibrium
    iteration stops when both limits are satisfied, or when the number of
    iterations exceeds :param:`stopmaxiter`.

:optparam:`stopmaxiter{in}`
    Maximum number of equilibrium iterations.

:optparam:`cmflag{in}`
    A nonzero value uses the consistent mass matrix instead of the default
    lumped one.

:optparam:`miflag{in}`
    Enables multi-fluid analysis, letting the algorithm solve the flow of two
    immiscible fluids in a fixed spatial domain (currently in 2D only).  This
    can also be used for free surface problems, with one of the fluids
    representing air.  The supported values are described in
    :ref:`materialinterfaces`.  The initial distribution of the reference fluid
    volume must be provided, as well as constitutive models for both fluids.

:optparam:`scaleflag{in}`, :optparam:`lscale{rn}`, :optparam:`uscale{rn}`, :optparam:`dscale{rn}`
    See :ref:`nondimensional-scaling`.

:optparam:`lstype{in}`, :optparam:`smtype{in}`
    Linear solver type and sparse matrix storage scheme.  Note that this
    algorithm leads to a non-symmetric matrix, so the scheme must be
    compatible with the solver type.  See :ref:`sparselinsolver`.

.. _pfemIncomp:

Transient incompressible flow — PFEM algorithm
----------------------------------------------

.. note::

   This problem type is part of the **PFEM module** and can be used only when
   that module is configured.

Solves transient incompressible flow using the particle finite element method,
based on the Lagrangian formulation of the Navier-Stokes equations.

Mesh nodes are represented by PFEM particles (see :ref:`pfemparticles`), which
can move freely and even separate from the main domain.  A temporary mesh built
from the particles is needed to integrate the governing equations in each
solution step; it is rebuilt from scratch in every step to prevent large element
distortion.  The problem is therefore defined without any elements in the input
file, and :param:`cs` and :param:`material` assign the cross section and
material to the elements that are created.

The mesh is generated by Delaunay triangulation, with the alpha-shape technique
identifying the free surface.  On the free surface a zero-pressure boundary
condition is enforced, which must be defined in the boundary condition record
whose number is given by :param:`pressure`.

.. record::

   :descitem:`PFEM` :elemparam:`nsteps{in}` :elemparam:`deltaT{rn}`
   :elemparam:`material{in}` :elemparam:`cs{in}`
   :elemparam:`pressure{in}` :optelemparam:`mindeltat{rn}`
   :optelemparam:`maxiter{in}` :optelemparam:`rtolv{rn}`
   :optelemparam:`rtolp{rn}` :optelemparam:`alphashapecoef{rn}`
   :optelemparam:`removalratio{rn}` :optelemparam:`scheme{in}`
   :optelemparam:`lstype{in}` :optelemparam:`smtype{in}`

**Parameters**

:param:`nsteps{in}`
    Number of solution steps.

:param:`deltaT{rn}`
    Time step length used for integration.  To ensure numerical stability the
    step length is adapted to the mesh geometry and to the velocity of the
    individual nodes.

:param:`material{in}`, :param:`cs{in}`
    Numbers of the material and cross section records assigned to the elements
    created from the particles.

:param:`pressure{in}`
    Number of the boundary condition record enforcing zero pressure on the free
    surface.

:optparam:`mindeltat{rn}`
    Minimum time step length, preventing the step from becoming too short.

:optparam:`removalratio{rn}`
    Limit, relative to the element edge length, below which particles that are
    too close together are removed from the solution.

:optparam:`maxiter{in}`
    Maximum number of iterations.  Default ``50``.

:optparam:`rtolv{rn}`, :optparam:`rtolp{rn}`
    Relative norms for the velocity and pressure difference between two
    successive iterations.  Defaults ``1.e-8``.

:optparam:`alphashapecoef{rn}`
    Alpha-shape coefficient, which should reflect the initial distribution of
    the PFEM particles.  A value of about 1.5 times the shortest distance
    between two neighbouring particles has been found to work well.

:optparam:`scheme{in}`
    Whether the equation system for the components of the auxiliary velocity is
    solved explicitly (``0``) or implicitly (``1``, the default).

:optparam:`lstype{in}`, :optparam:`smtype{in}`
    Linear solver type and sparse matrix storage scheme.  Note that this
    algorithm leads to a non-symmetric matrix, so the scheme must be
    compatible with the solver type.  See :ref:`sparselinsolver`.

Nonlinear solver parameters
===========================

.. _petscsnessolver:

PETSc SNES
----------

The ``petscsnes`` adapter delegates nonlinear globalization to PETSc SNES
while retaining OOFEM's residual and tangent assembly.  It is selected in a
``StaticStructural`` engineering-model or meta-step record.

.. record:: PETSc SNES solver options

   :elemparam:`solvertype{s}` :elemparam:`smtype{in}`
   :optelemparam:`snesmaxiter{in}` :optelemparam:`snesmaxfunc{in}`
   :optelemparam:`snesatol{rn}` :optelemparam:`snesrtol{rn}`
   :optelemparam:`snesstol{rn}` :optelemparam:`snestype{s}`
   :optelemparam:`snesoptionsprefix{s}`
   :optelemparam:`snesmonitor{in}` :optelemparam:`soldepextforces{}`
   :optelemparam:`snesfieldscaling{in}`
   :optelemparam:`snesminscalesquared{rn}`
   :optelemparam:`forcescaledofs{ia}` :optelemparam:`forcescale{ra}`
   :optelemparam:`snesbounddofs{ia}`
   :optelemparam:`sneslowerbounds{ra}`
   :optelemparam:`snesupperbounds{ra}`
   :optelemparam:`snesirreversibledofs{ia}`
   :optelemparam:`snestrustregion{in}`
   :optelemparam:`snestrdelta0{rn}`
   :optelemparam:`snestrdeltamin{rn}`
   :optelemparam:`snestrdeltamax{rn}`
   :optelemparam:`snestrmaxtrials{in}`

:param:`solvertype{s}`
    Set to ``petscsnes``.

:param:`smtype{in}`
    Set to ``7``, the PETSc sparse matrix storage.

This solver requires an OOFEM build with PETSc.  It currently supports only an
OOFEM model running in serial; shared-memory assembly may still be used.  The
adapter fixes the load factor at one and does not provide arc-length or other
indirect load control.  Each attempted solve invokes ``SNESSolve`` once;
cutbacks invoke it again for the reduced step.  A failed SNES solve is reported
as nonconvergence so that the engineering model can apply its usual cutback and
restart procedure.

The residual supplied to PETSc is the internal force minus the external force,
and its Jacobian is the OOFEM tangent.  SNES uses a separate PETSc matrix so
that residual row scaling does not modify the physical OOFEM tangent.

Basic options
~~~~~~~~~~~~~

:optparam:`snesmaxiter{in}`
    Maximum number of nonlinear iterations.  The default is ``50`` and the
    value must be positive.

:optparam:`snesmaxfunc{in}`
    Maximum number of residual evaluations.  The default is ``10000`` and
    the value must be positive.

:optparam:`snesatol{rn}`
    Absolute residual tolerance.  The default is ``1.e-50`` and the value
    must be non-negative.

:optparam:`snesrtol{rn}`
    Relative residual tolerance.  The default is ``1.e-8`` and the value
    must be non-negative.

:optparam:`snesstol{rn}`
    Relative step tolerance.  The default is ``1.e-12`` and the value must
    be non-negative.

:optparam:`snestype{s}`
    PETSc SNES type name.  The default is ``newtonls`` and the string must
    not be empty.  Other PETSc types can also be selected through PETSc
    run-time options.

:optparam:`snesmonitor{in}`
    If nonzero, print the OOFEM SNES iteration monitor.  The default is ``1``.
    The output also follows OOFEM's iteration-log and problem-scale filtering.
    PETSc's own monitors can independently be enabled at run time.

:optparam:`soldepextforces{}`
    Reassemble ``ExternalRhs`` at every SNES trial solution and after state
    restoration.  By default the external-force vectors passed to the solver
    remain fixed.  This option disables support for the compact custom equation
    numbering used by staggered solves.

Residual scaling
~~~~~~~~~~~~~~~~

:optparam:`snesfieldscaling{in}`
    If nonzero, scale the residual and corresponding Jacobian rows by DOF
    group.  The default is ``1``.  The scale is constructed at the first
    admissible residual evaluation and remains fixed during that nonlinear
    solve.  Consequently, the SNES residual tolerances apply to the scaled
    norm.

:optparam:`snesminscalesquared{rn}`
    Positive threshold below which a DOF group is left unscaled.  The default
    is ``1.e-6`` and the value must be finite and positive.

:optparam:`forcescaledofs{ia}`
    DOF identifiers for which characteristic forces are supplied.  Specify
    together with :optparam:`forcescale{ra}`; both arrays must have the same
    length and a DOF identifier may occur only once.

:optparam:`forcescale{ra}`
    Positive, finite characteristic force corresponding to every entry in
    :optparam:`forcescaledofs{ia}`.

For a DOF group :math:`d`, the squared characteristic force is formed from
the squared external-force norm, the sum of squared element-force
contributions, and, when specified, :math:`n_d f_d^2`, where :math:`n_d` is
the number of active equations in the group and :math:`f_d` is its
``forcescale`` value.  If this sum is at least
``snesminscalesquared``, the residual and tangent rows in the group are
multiplied by the inverse square root of the sum.  Otherwise they are left
unscaled.  ``forcescale`` and ``snesminscalesquared`` have no effect when
``snesfieldscaling`` is zero.

Variable bounds
~~~~~~~~~~~~~~~

:optparam:`snesbounddofs{ia}`
    DOF identifiers to be bounded.  Specify together with
    :optparam:`sneslowerbounds{ra}` and :optparam:`snesupperbounds{ra}`; all
    three arrays must have equal lengths and a DOF identifier may occur only
    once.  A bound is applied to every active equation having the corresponding
    DOF identifier.

:optparam:`sneslowerbounds{ra}`
    Finite lower bound corresponding to every entry in
    :optparam:`snesbounddofs{ia}`.

:optparam:`snesupperbounds{ra}`
    Finite upper bound corresponding to every entry in
    :optparam:`snesbounddofs{ia}`.  Each lower bound must not exceed its upper
    bound.

Bounds require the SNES type ``vinewtonrsls`` or ``vinewtonssls``.  This
requirement is checked after PETSc run-time options have been applied.

:optparam:`snesirreversibledofs{ia}`
    Bounded DOF identifiers whose value may not decrease during a step.  The
    default is an empty list.  Every identifier must also occur in
    ``snesbounddofs``.  Its effective lower bound is the larger of the
    specified lower bound and the accepted value at the beginning of the step.
    That reference value is retained across cutbacks and staggered sweeps.

PETSc's VI active-bound tolerance defaults to ``1.e-8`` and can be changed
with ``-snes_vi_zero_tolerance`` or its prefixed equivalent.

Bounded trust region
~~~~~~~~~~~~~~~~~~~~

:optparam:`snestrustregion{in}`
    Enable OOFEM's bounded reduced-space trust-region line search when set to
    ``1``.  The default is ``0``.  This is distinct from selecting PETSc's
    ``newtontr`` SNES type.  It requires at least one variable bound, and the
    selected SNES type must resolve to ``vinewtonrsls`` after run-time
    options have been processed.

:optparam:`snestrdelta0{rn}`
    Initial trust-region radius.  The default is ``1``.

:optparam:`snestrdeltamin{rn}`
    Minimum trust-region radius.  The default is ``1.e-8`` and it must be
    positive.

:optparam:`snestrdeltamax{rn}`
    Maximum trust-region radius.  The default is ``1.e8``.

:optparam:`snestrmaxtrials{in}`
    Maximum number of trial points in one trust-region step.  The default is
    ``20`` and the value must be positive.

The three radii must be finite and satisfy

.. math::

   0 < \mathtt{snestrdeltamin}
     \leq \mathtt{snestrdelta0}
     \leq \mathtt{snestrdeltamax}.

The trust-region metric is balanced by DOF group and is frozen from the first
full Newton direction of the nonlinear solve.

PETSc run-time options
~~~~~~~~~~~~~~~~~~~~~~

:optparam:`snesoptionsprefix{s}`
    PETSc options prefix assigned to this SNES object.  The default is an empty
    string.  It must not start with ``-``, which PETSc adds automatically; a
    trailing underscore is customary.  For example,
    ``snesoptionsprefix "nl_"`` selects options such as
    ``-nl_snes_type``, ``-nl_ksp_type``, ``-nl_pc_type``,
    ``-nl_snes_linesearch_type``, and ``-nl_snes_vi_zero_tolerance``.

PETSc options are processed after the OOFEM input fields and therefore
override the input SNES type, tolerances, iteration limits, and PETSc KSP/PC
choices.  When ``snestrustregion 1`` is used, its shell line search is
installed after PETSc options are processed and therefore replaces any
run-time line-search type.

The tangent predictor used by ``StaticStructural`` has a separate,
unprefixed PETSc KSP.  Without ``snesoptionsprefix``, global ``-ksp_*``
and ``-pc_*`` options affect both that predictor and SNES's inner linear
solve.  With a prefix, unprefixed KSP/PC options affect the predictor, while
prefixed KSP/PC options affect SNES.

Each meta-step record is self-contained.  Omitted SNES fields are reset to the
defaults above instead of inheriting values from the preceding meta-step.

Example
~~~~~~~

A minimal ``StaticStructural`` analysis record using the default tangent
predictor and default ``newtonls`` SNES type is

.. code-block:: text

   StaticStructural nsteps 1 smtype 7 solvertype "petscsnes" snesrtol 1.e-10 snesatol 1.e-12

PETSc algorithms and preconditioners may then be selected at run time, for
example:

.. code-block:: console

   oofem -f model.in -snes_type newtonls -ksp_type preonly -pc_type lu

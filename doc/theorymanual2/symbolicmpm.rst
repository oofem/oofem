Symbolic MPM module (Concept & Formulation)
===========================================================================

The symbolic MPM module allows to define problems using syntax describing their weak form.
The weak form is represented using

   * *Variables*, defining unknown fields and test fields appearing in the weak form. Variables define physical meaning of the field, its interpolation and degrees of freedom.
   * *Terms*, representing mathematical expressions appearing in the weak form. Terms are evaluated for specific variables. Beside the library of pre-defined terms, an arbitrary term can be entered directly in the input file as a symbolic expression (*SymbolicTerm*, see below).
   * *Integrals*, representing integrals of terms over specific domains. Integral domains are defined using sets.
   * *Materials*, providing the constitutive response. A material receives the generalized state of an integration point once per equilibrium iteration and the terms then query the cached response (see :ref:`mpm-material-interface`).

The characteristic discretized equations are assembled from individual integral contributions by looping over individual elements of each integral domain.
Each term can evaluate itself for known variables at given time and also evaluate its linearization with respect to unknown variables.
How these contributions are assembled is defined by the engineering model, representing  the problem type to be solved.

We will illustrate the concept of MPM Symbolic module on an example of formulation
of the elasticity problem able to represent the
incompressible behavior. This formulation can be written using mixed approach involving
pressure and displacement fields

.. math:: \nabla p + 2\mu\nabla\cdot\rm{dev}[\nabla^s u]+f = 0\;\rm{in}\ \Omega

.. math:: \frac{1}{K}p-\nabla\cdot u=0\;\rm{in}\ \Omega

.. math:: u=\bar u\;\rm{on}\ \Gamma_u

.. math:: \sigma n = \bar t\;\rm{on}\ \Gamma_t

The weak form of the above equations can be written as

.. math:: \int_\Omega\overbrace{2\mu\nabla^sw:\rm{dev}[\nabla^su]}^{T_1}\ d\Omega-\int_\Omega\overbrace{ \nabla\cdot w\ p}^{T_2}\ d\Omega=\int_{\Gamma_t}\overbrace{w\cdot \bar t}^{T_3}\ d\Gamma

.. math::

   -\int_\Omega\underbrace{q\ \nabla\cdot u}_{T_4 = T_2^T}\ d\Omega+\int_\Omega\underbrace{\frac{1}{K}q\ p}_{T_5}\ d\Omega = 0

where

   * :math:`u, p` are unknown displacement and pressure fields (*Variables*),
   * :math:`w, q` are the corresponding test fields (*Variables*).

In order to
obtain stable results, the approximations should satisfy Babuška-Brezzi
condition or some form of stabilization would be required. To satisfy
the B-B condition, the approximation order of displacement should be
higher than that of hydrostatic pressure. In the following we will
consider 2D case.

Variables
---------

In the input deck, the *Variables* (:math:`u,p,w,q`) can be set up using
following syntax

::

   Variable name "u"  interpolation "feilin"   type 1 quantity "Displacements" size 2 dofs 2 1 2
   Variable name "w"  interpolation "feilin"   type 1 quantity "Displacements" size 2 dofs 2 1 2
   Variable name "p"  interpolation "feiconst" type 0 quantity "Pressure"      size 1 dofs 1 11
   Variable name "dp" interpolation "feiconst" type 0 quantity "Pressure"      size 1 dofs 1 11

where *interpolation* determines the interpolation used for specific
field. Here we use linear approximation (*interpolation “feilin”*) for
displacement field and related test field and constant approximation
(*interpolation “feiconst”*) for pressure (so called Q1P0 element). The
:math:`u, w` fields are vector fields (*type 1*) with physical meaning
of displacement (*quantity "Displacements"*) and two degrees of freedom
(*size 2 dofs 2 1 2*). And :math:`p, dp` are scalar fields (*type 0*)
with physical meaning of pressure (*quantity "Pressure" size 1 dofs 1 11*).
The test field of the pressure, denoted :math:`q` in the weak form above,
is named ``dp`` here, following the shipped input decks.

The *quantity* attribute names the physical field the variable represents.
It is a ``FieldType``, written as the enumerator name with the ``FT_``
prefix dropped — ``"Displacements"``, ``"Pressure"``, ``"Pressure2"``,
``"Temperature"``, ``"Concentration1"``, ``"Concentration2"``. This is what
the assembly matches against the state layout declared by the material
(:ref:`mpm-material-interface`), so it has to say what the field *is*;
*dofs* says only on which degrees of freedom that field is numbered. The
two declarations are independent and are cross-checked at initialization: a
deck storing, say, temperature on ``D_w`` is legitimate and only produces a
warning. Numeric values of *quantity* are still accepted but should not be
used — they are an implementation detail, and the meaning of a given number
has changed in the past.

A variable that is a test (weighting) function should say so, with
*dualto* naming the unknown it weights

::

   Variable name "dp" interpolation "feiconst" type 0 quantity "Pressure" size 1 dofs 1 11 dualto "p"

Both records are otherwise
indistinguishable — they may even carry different interpolations, which is
what a non-symmetric (Petrov-Galerkin) formulation looks like — and without
it the assembly cannot tell which field to read nodal unknowns from when
pushing the state to the material. It is required only where a deck names a
test function as some term's *variable*, which happens for pure source
terms (those with a zero linearization), since they have no unknown to
depend on.

Pre-defined terms
-----------------

The weak form above consists of several terms to be evaluated

   * :math:`T_1`: This is represented by *BTSigmaTerm*, evaluated for
     :math:`w` test field and :math:`u` as unknown field, under plain strain
     assumptions (*mmode 7*) and considering deviatoric part only
     (*lhsmatmode 29*), with the stress read back for the residual
     (*rhsmatmode 30*)

     ::

        BTSigmaTerm 1 variable "u"  testvariable "w" mmode 7 lhsmatmode 29 rhsmatmode 30

   * :math:`T_2`: represented by *BTamNTerm*, evaluated for :math:`w` test
     field and :math:`p` as unknown field, under plain strain assumptions
     (*mmode 7*) and with scalar parameter equal to 1.0 (*atype 28*)

     ::

        BTamNTerm 2 variable "p" testvariable "w" mmode 7 atype 28

   * Similarly, we set up remaining three terms :math:`T_3, T_4=T_2^T` and
     :math:`T_5`. The *uvmt 1* attribute selects the value mode in which the
     unknown field is read (total values instead of the default rates)

     ::

        NTamTBTerm 3 variable "u" testvariable "dp" mmode 7 atype 28 uvmt 1
        NTcN 4 variable "p" testvariable "dp" mmode 7 ctype 27 uvmt 1
        NTfTerm 5 variable "u" testvariable "w" mmode 6 flux 2 0. 6.25

The *mmode* attribute takes a ``MaterialMode``, which may be given by name
(``mmode "_PlaneStrain"``) instead of by number. The library of available
terms is documented in the *Term documentation* chapter.
The terms are integrated over specific domains
(:math:`\Omega,\ \Gamma_t`), defined using corresponding sets (defined
bellow in the example):

::

   Integral 1 domain 1 set 1 term 1
   Integral 2 domain 1 set 1 term 2 factor -1.0
   Integral 3 domain 1 set 1 term 3 factor -1.0
   Integral 4 domain 1 set 1 term 4
   Integral 5 domain 1 set 2 term 5

User-defined terms: *SymbolicTerm*
----------------------------------

A term does not have to exist as a C++ class. *SymbolicTerm* takes the
mathematical expressions of the term itself, written in a small matrix
language, and compiles them into bytecode at initialization; the assembly
loops then execute the compiled form:

::

   SymbolicTerm <id> variable "<name>" testvariable "<name>" mmode <mode>
                lexpression "<lhs expression>" rexpression "<rhs expression>"

*lexpression* defines the linearization (the tangent/stiffness
contribution, a matrix) and *rexpression* the residual/internal force
contribution (a column matrix, converted to a vector). The operators of the
weak form are available as functors: ``Grad_s``, ``Grad``, ``Div`` and
``N`` build the interpolation operator matrices of a variable at the
integration point, ``MDer``, ``MVec`` and ``MProp`` read the material
response, and ``ru``/``rv`` read the nodal unknowns/rates of a field on the
cell.

The five terms of the Cook membrane formulation above read, in symbolic
form (compare with the pre-defined terms of the previous section)

::

   SymbolicTerm 1 variable "u" testvariable "w"  mmode 7 lexpression "Grad_s(w,gp).T*MDer(gp,ts, MatResponseMode::DeviatoricStiffness)*Grad_s(u,gp)" rexpression "Grad_s(w,gp).T*Sig_dev(u,gp,ts)"
   SymbolicTerm 2 variable "p" testvariable "w"  mmode 7 lexpression "Div(w,gp).T*N(p,gp)" rexpression "Div(w,gp).T*N(p,gp)*ru(p, cell, ts)"
   SymbolicTerm 3 variable "u" testvariable "dp" mmode 7 lexpression "N(dp,gp).T*Div(u,gp)" rexpression "N(dp,gp).T*Div(u,gp)*ru(u, cell, ts)"
   SymbolicTerm 4 variable "p" testvariable "dp" mmode 7 lexpression "N(dp,gp).T*2.4e-7*N(p,gp)" rexpression "N(dp,gp).T*2.4e-7*N(p,gp)*ru(p, cell, ts)"
   NTfTerm     5 variable "u" testvariable "w"  mmode 6 flux 2 0. 6.25

with the deck otherwise unchanged; the complete input file is
`tests/regression/mpm/mpms_cook2_u2p1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms_cook2_u2p1.in>`__.
Response modes are best given by name
(``MatResponseMode::DeviatoricStiffness``) rather than by their integer
value. The full description of the expression language, of the available
functors and of the evaluation mechanics is in ``doc/symbolic_term.md``.

.. _mpm-material-interface:

Material interface: state push and declared state layout
--------------------------------------------------------

The material interface used by the module is split into two halves. Once
per equilibrium iteration, before any term is evaluated, the solver hands
each integration point's material its *generalized state* — strains, field
gradients, field values — and the material performs all of its constitutive
work there, caching the results:

.. code-block:: cpp

   virtual void updateTempState(const FloatArray &stateVector, GaussPoint *gp, TimeStep *tStep);

The terms then only *read* that cache, selecting what they need by
``MatResponseMode``:

.. code-block:: cpp

   virtual void   giveCharacteristicVector(FloatArray &answer,  MatResponseMode, GaussPoint*, TimeStep*) const;
   virtual void   giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode, GaussPoint*, TimeStep*) const;
   virtual double giveCharacteristicValue (              MatResponseMode, GaussPoint*, TimeStep*) const;

This matters for the formulation because the tangent and the residual are
assembled in separate sweeps over the integration points: with the state
pushed once per iteration, the response a term reads no longer depends on
the order in which the sweeps happen to run, and a term never has to
assemble and pass the state itself.

What the state vector contains, and in which order, is declared by the
material rather than fixed by convention:

.. code-block:: cpp

   virtual StateVariableLayout giveStateVariableIDs(MaterialMode mmode) const;

Each entry of the returned layout is a ``{FieldType field, StateOperator op}``
pair, in packing order: *which* primary field, and *what is taken of it*.
For each entry the assembly finds the unknown variable whose *quantity* is
``field``, reads its nodal unknowns on the cell, applies the operator at the
integration point and appends the resulting block. Four operators are
available:

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - ``StateOperator``
     - Meaning
     - Block size
   * - ``SO_Value``
     - the field itself, :math:`Nu`
     - ``Variable::size``
   * - ``SO_Gradient``
     - gradient of a scalar field, :math:`\nabla u`
     - number of spatial dimensions
   * - ``SO_SymmetricGradient``
     - symmetric gradient of a vector field, :math:`\nabla^s u` (strain, Voigt form)
     - material mode dependent (6 in ``_3dMat`` and the ``UP`` modes, 4 in ``_PlaneStrain``, 3 in ``_PlaneStress``, 1 in ``_1dMat``)
   * - ``SO_Divergence``
     - divergence of a vector field, :math:`\nabla\cdot u`
     - 1

So a poromechanical material consuming strain, pressure gradient and
pressure declares

.. code-block:: cpp

   StateVariableLayout giveStateVariableIDs(MaterialMode mmode) const override {
       return { { FT_Displacements, SO_SymmetricGradient },
                { FT_Pressure,      SO_Gradient },
                { FT_Pressure,      SO_Value } };
   }

and receives :math:`[\,\varepsilon,\ \nabla p,\ p\,]` packed in that order.
A material returning an empty layout — the default — does not participate
in this protocol and keeps using its physics-specific entry points; a cell
that cannot supply every declared field (a thermo-mechanical material
driven for the thermal sub-problem alone, say) is skipped rather than given
a partial state. The same interface is available to materials implemented in
python, where ``giveStateVariableIDs(mmode)`` returns a sequence of
``(FieldType, StateOperator)`` pairs. The design rationale is described in
``doc/unified_material_interface.md``.

Example: Cook membrane
----------------------

The Cook’s membrane is a standard benchmark problem. It consists of a
tapered plate clamped on one of its sides with a transversal distributed
load applied to the opposite side. The plate is in plain strain and its
dimensions, as well as the material parameters and boundary conditions,
are shown in figure bellow

.. figure:: https://github.com/oofem/oofem.github.io/blob/master/assets/images/cookGeometry.png?raw=true
   :alt: Cook membrane geometry and boundary conditions

   Cook membrane geometry and boundary conditions

The complete OOFEM input deck for mesh consisting of 2x2 is shown below

::

   cook2_u1p0_2.out
   Demo of symbolic mpm problem; Cook membrane benchmark
   #
   mpmsymbolicstationaryproblem nsteps 1 nvariables 4 nterms 5 nintegrals 5 lhsterms 4 1 2 3 4 rhsterms 1 5 nmodules 1 profileopt 1 maxiter 10
   errorcheck
   Variable name "u" interpolation "feilin" type 1 quantity "Displacements" size 2 dofs 2 1 2 # displacement
   Variable name "w" interpolation "feilin" type 1 quantity "Displacements" size 2 dofs 2 1 2 # test function
   Variable name "p" interpolation "feiconst"  type 0 quantity "Pressure" size 1 dofs 1 11 # pressure
   Variable name "dp" interpolation "feiconst" type 0 quantity "Pressure" size 1 dofs 1 11 # test function
   BTSigmaTerm 1 variable "u"  testvariable "w" mmode 7 lhsmatmode 29 rhsmatmode 30
   BTamNTerm 2 variable "p" testvariable "w" mmode 7 atype 28
   NTamTBTerm 3 variable "u" testvariable "dp" mmode 7 atype 28 uvmt 1
   NTcN 4 variable "p" testvariable "dp" mmode 7 ctype 27 uvmt 1
   NTfTerm 5 variable "u" testvariable "w" mmode 6 flux 2 0. 6.25
   Integral 1 domain 1 set 1 term 1
   Integral 2 domain 1 set 1 term 2 factor -1.0
   Integral 3 domain 1 set 1 term 3 factor -1.0
   Integral 4 domain 1 set 1 term 4
   Integral 5 domain 1 set 2 term 5
   domain HeatTransfer
   outputmanager tstep_all dofman_all element_all
   ndofman 9 nelem 6 nbc 1 ncrosssect 1 nic 0 nltf 2 nmat 1 nset 3
   Node 1 coords 3 0.0 0.0 0.0
   Node 2 coords 3 24.0 22.0 0.0
   Node 3 coords 3 48.0 44.0 0.0
   Node 4 coords 3 0.0 22.0 0.0
   Node 5 coords 3 24.0 37.0 0.0
   Node 6 coords 3 48.0 52.0 0.0
   Node 7 coords 3 0.0 44.0 0.0
   Node 8 coords 3 24.0 52.0 0.0
   Node 9 coords 3 48.0 60.0 0.0
   q1 1 nodes 4 1 2 5 4 mat 1 crosssect 1
   q1 2 nodes 4 2 3 6 5 mat 1 crosssect 1
   q1 3 nodes 4 4 5 8 7 mat 1 crosssect 1
   q1 4 nodes 4 5 6 9 8 mat 1 crosssect 1
   l1 5 nodes 2 3 6 mat 1 crosssect 1
   l1 6 nodes 2 6 9 mat 1 crosssect 1
   simplecs 1 thick 5.0
   isole 1 d 1 e 250 n 0.49999 talpha 1.
   # clamped-displacement
   boundarycondition 1 loadtimefunction 1 set 3 values 2 0 0   dofs 2 1 2
   constantfunction 1 f(t) 1
   PiecewiseLinFunction 2 nPoints 4 t 4 -10. 0. 1. 5. f(t) 4 0. 0. 1.0 1.0
   set 1 elementranges  {(1 4)}
   set 2 elementranges  {(5 6)}
   set 3 elementedges 4 1 4 3 4

The complete input deck can be found in
`tests/regression/mpm/cook2_u1p0_2.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u1p0_2.in>`__
file; the variant driven by the lightweight ``test`` engineering model is
`tests/regression/mpm/cook2_u1p0.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u1p0.in>`__.

Switching to high-order interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To switch to high order element Q2P1 with quadratic interpolation of
displacement and linear interpolation of pressure we need just to use
different interpolations for variables. The updated variable input
records are following

::

   Variable name "u" interpolation "feiquad" type 1 quantity "Displacements" size 2 dofs 2 1 2 # displacement
   Variable name "w" interpolation "feiquad" type 1 quantity "Displacements" size 2 dofs 2 1 2 # test function
   Variable name "p" interpolation "feilin"  type 0 quantity "Pressure" size 1 dofs 1 11 # pressure
   Variable name "dp" interpolation "feilin" type 0 quantity "Pressure" size 1 dofs 1 11 # test function

The complete input deck can be found in
`tests/regression/mpm/cook2_u2p1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u2p1.in>`__
file, and its *SymbolicTerm* counterpart in
`tests/regression/mpm/mpms_cook2_u2p1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms_cook2_u2p1.in>`__.


To illustrate the convergence, sequence of uniform meshes of the plate
is considered, starting from a mesh consisting of a single elements and
proceeding by uniform refinement.

The figure below shows the vertical displacement of the plane tip
plotted against the number of element segments along each side. The
solution is compared to reference solution  [1]_.

.. figure:: https://github.com/oofem/oofem.github.io/blob/master/assets/images/cookConvergence.png?raw=true
   :alt: Cook membrane convergence graph


References
~~~~~~~~~~

.. [1]
   Ignacio Romero, Manfred Bischoff, Incompatible Bubbles: A
   non-conforming finite element formulation for linear elasticity,
   Computer Methods in Applied Mechanics and Engineering, Volume 196,
   Issues 9–12, 2007, Pages 1662-1672, ISSN 0045-7825,
   https://doi.org/10.1016/j.cma.2006.09.010.

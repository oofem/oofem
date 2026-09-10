Plasticity-based material models
================================


.. _DPmodel:

Drucker-Prager model - DruckerPrager
------------------------------------

The Drucker-Prager plasticity model (Contributed by Simon Rolshoven, LSC, FENAC, EPFL.) is an isotropic elasto-plastic model based
on a yield function

.. math::

   f \left(\bsig, \tauY\right) = F\left(\bsig\right) - \tauY

with the pressure-dependent equivalent stress

.. math::

   F\left(\bsig\right) = \alpha I_1 + \sqrt{J_2}

As usual, :math:`\bsig` is the stress tensor, :math:`\tauY` is the yield stress
under pure shear, and :math:`I_1` and :math:`J_2` are the first invariant and second
deviatoric invariant of the stress tensor.
The friction coefficient :math:`\alpha` is a positive parameter that
controls the influence of the pressure on the yield limit, important for
cohesive-frictional materials such as concrete, soils or other
geomaterials. Regarding Mohr-Coulomb plasticity, relation to cohesion, :math:`c`, and 
the angle of friction, :math:`\theta`, exists for the Drucker-Prager model

.. math::

   \alpha = \frac{2\sin\theta}{(3-\sin\theta)\sqrt{3}}\\
   \tauY = \frac{6c\cos\theta}{(3-\sin\theta)\sqrt{3}}

The flow rule is derived from the plastic potential

.. math::

   g\left(\bsig\right) = \alphaPsi I_1 + \sqrt{J_2}

where :math:`\alphaPsi` is the dilatancy coefficient. An associated
model with :math:`\alphaPsi=\alpha` would overestimate the dilatancy of
concrete, so the dilatancy coefficient is usually chosen smaller than the
friction coefficient.

The model is described by the equations

.. math::

   \bsig & = &\mbf{D} : \left(\eps - \epsp \right)
   \\
   \tauY & = & h(\kappa)
   \\
   \epspd = \dot{\lambda} \frac{\partial g}{\partial \bsig} & = &
   \dot{\lambda} \left( \alphaPsi \mbf{\delta} + \frac{\mbf{s}}{2\sqrt{J_2}} \right)
   \\
   \dot{\kappa} & = & \sqrt{\frac{2}{3}} \; \| \epspd \|

and

.. math::

   \dot{\lambda} \ge 0, \;\;\; f \left(\bsig, \tauY\right) \le 0, \;\;\; \dot{\lambda}\, f \left(\bsig, \tauY\right)  = 0

which represent the linear elastic law, hardening law, evolution laws
for plastic strain and hardening variable,  and the
loading-unloading conditions.

In the above, :math:`\mbf{D}` is the elastic stiffness
tensor, :math:`\eps` is the strain tensor, :math:`\epsp` is the plastic strain tensor,
:math:`\lambda` is the plastic multiplier, :math:`\mbf{\delta}` is the unit
second-order tensor, :math:`\mbf{s}` is the
deviatoric stress tensor, :math:`\kappa` is the hardening variable, and a
superior dot marks the derivative with respect to time.

The flow rule has the form given in Eq.~(``eq:flowRule``) at all
points of the conical yield surface with the exception of its vertex,
located on the hydrostatic axis.

For the present model, the evolution
of the hardening variable can be explicitly linked to the plastic
multiplier.
Substituting the flow rule
(``eq:flowRule``) into Eq.~(``eq:hardening``) and computing the norm
leads to

.. math::

   \dot\kappa = k \dot \lambda

with a constant parameter :math:`k = \sqrt{1/3 + 2 \alphaPsi^2}`, so the
hardening variable is proportional to the plastic multiplier.
For :math:`\alpha=\alphaPsi=0`, the associated :math:`J_2`-plasticity model
is recovered as a special case.

In the simplest case of linear hardening, the hardening function is a linear function of :math:`\kappa`, given by

.. math::
   :label: bilin-soft

   h(\kappa) = \tau_0 + H E \kappa

where :math:`\tau_0` is the initial yield stress, and :math:`H` is the hardening modulus normalized with the elastic modulus. Alternatively, an exponential hardening function

.. math::

   h(\kappa) = \tau_{\mathrm{limit}} + \left( \tau_0 - \tau_{\mathrm{limit}} \right) \mathrm{e}^{-\kappa/\kappac}

can be used for a more realistic description of hardening.

The stress-return algorithm is based on the Newton-iteration. In plasticity, this is commonly called Closest-Point-Projection (CPP), and it generally leads to quadratic convergence. The implemented algorithm is convergent in any stress case, but in the vicinity of the vertex region, quadratic convergence might be lost because of insufficient regularity of the yield function.

The algorithmic tangent stiffness matrix is implemented for both the regular case and the vertex region. Generally, the error decreases quadratically (of course only asymptotically). Again, in the vicinity of the vertex region, quadratic convergence might be lost due to insufficient regularity. Furthermore, the tangent stiffness matrix does not always exist for the vertex case. In these cases, the elastic stiffness is used instead. It is generally safer (but slower) to use the elastic stiffness if you encounter any convergence problems, especially if your problem is tension-dominated.

.. table:: DP material - summary
   :name: DP_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | DP material                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`DruckerPrager` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`tAlpha{rn}`    |
   |                    | :elemparam:`E{rn}` :elemparam:`n{rn}` :elemparam:`alpha{rn}` :elemparam:`alphaPsi{rn}`       |
   |                    | :elemparam:`ht{in}` :elemparam:`iys{rn}` :elemparam:`lys{rn}` :elemparam:`hm{rn}`            |
   |                    | :elemparam:`kc{rn}` :optelemparam:`yieldtol{rn}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`alpha` friction coefficient                                                        |
   |                    | - :param:`alphaPsi` dilatancy coefficient                                                    |
   |                    | - :param:`ht` hardening type, 1: linear hardening, 2: exponential hardening                  |
   |                    | - :param:`iys` initial yield stress in shear, :math:`\tau_0`                                 |
   |                    | - :param:`lys` limit yield stress for exponential hardening, :math:`\tau_{\mathrm{limit}}`   |
   |                    | - :param:`hm` hardening modulus normalized with E-modulus (!)                                |
   |                    | - :param:`kc` :math:`\kappa_c` for the exponential softening law                             |
   |                    | - :param:`yieldtol` tolerance of the error in the yield criterion, default value 1.e-14      |
   |                    | - :param:`newtonIter` maximum number of iterations in :math:`\lambda` search, default value  |
   |                    |   30                                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStrain, 3dRotContinuum                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/DruckerPrager_01.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/DruckerPrager_01.in>`_        |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _DPCutmodel:


Drucker-Prager model with tension cut-off and isotropic damage - DruckerPragerCut
---------------------------------------------------------------------------------


The Drucker-Prager plasticity model with tension cut-off is a multisurface model, appropriate for cohesive-frictional materials such as concrete loaded both in compression and tension. The plasticity model is formulated for isotropic hardening and enhanced by isotropic damage, which is driven by the cumulative plastic strain. The model can be used only in the small-strain context, with additive split of the strain tensor into the elastic and plastic parts.

The basic equations include the additive decomposition of strain into elastic and plastic parts,

.. math::

   \veps = \veps_{\rm{e}} + \veps_{\rm{p}},

the stress strain law 

.. math::

   \vsig = (1-\omega)\bar{\vsig}=(1-\omega)\mbf{D} :(\veps-\veps_{\rm{p}}),

the definition of the yield function in terms of the effective stress,

.. math::

   f \left(\bar{\vsig},\kappa \right) = \alpha I_1 + \sqrt{J_2} - \tau_Y,

the flow rule

.. math::

   g \left( \bar{\vsig} \right) = \alphaPsi I_1 + \sqrt{J_2},

the linear hardening law

.. math::

   \tau_Y(\kappa) = \tau_0 + H\kappa,

where :math:`\tau_0` represents the initial yield stress under pure shear,
the damage law

.. math::
   :label: damagelawDP

   \omega(\kappa) = \omega_c(1-\mbox{e}^{-a\kappa}),

where :math:`\omega_c` is critical damage and :math:`a` is a positive dimensionless parameter.
More detailed descriptioin of some parameters is in Section :ref:`DPmodel`.

The dilatancy coefficient :math:`\alphaPsi` controls flow associativeness; if :math:`\alphaPsi=\alpha`, an associate model is recovered, which
overestimates the dilatancy of concrete. Hence, the dilatancy coefficient is usually chosen smaller, 
:math:`\alphaPsi \leq \alpha`, and the non-associated model is formulated.

.. table:: Drucker Prager material with tension cut-off - summary
   :name: DP_table_cut

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Drucker Prager material with tension cut-off                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`DruckerPragerCut` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`tAlpha{rn}` |
   |                    | :elemparam:`E{rn}` :elemparam:`n{rn}` :elemparam:`tau0{rn}` :elemparam:`alpha{rn}`           |
   |                    | [:elemparam:`alphaPsi{rn}`] [:elemparam:`H{rn}`] [:elemparam:`omega_crit{rn}`]               |
   |                    | [:elemparam:`a{rn}`] [:elemparam:`yieldtol{rn}`] [:elemparam:`NewtonIter{in}`]               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`tau0` initial yield stress in shear :math:`\tau_0`                                 |
   |                    | - :param:`alpha` friction coefficient                                                        |
   |                    | - :param:`alphaPsi` dilatancy coefficient, equals to :param:`alpha` by default               |
   |                    | - :param:`H` hardening modulus (can be negative in the case of plastic softening), 0 by      |
   |                    |   default                                                                                    |
   |                    | - :param:`omega_crit` critical damage in damage law (:eq:`damagelawDP`), 0 by default        |
   |                    | - :param:`a` exponent in damage law (:eq:`damagelawDP`), 0 by default                        |
   |                    | - :param:`yieldtol` tolerance of the error in the yield criterion, default value 1.e-14      |
   |                    | - :param:`newtonIter` maximum number of iterations in :math:`\lambda` search, default value  |
   |                    |   30                                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 1dMat, 3dMat, PlaneStrain                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _sec:misplast:

Mises plasticity model with isotropic damage - MisesMat
-------------------------------------------------------

This model is appropriate for the description of plastic yielding in ductile material such as metals, and it can also cover the effects of void growth. The model uses the Mises yield condition (in terms of the second deviatoric invariant, :math:`J_2`), associated flow rule, linear isotropic hardening driven by the cumulative plastic strain, and isotropic damage, also driven by the cumulative plastic strain. The model can be used in the small-strain context, with additive split of the strain tensor into the elastic and plastic parts, or in the large-strain context, with multiplicative split of the deformation gradient and with yield condition formulated in terms of Kirchhoff stress (which is the true Cauchy stress multiplied by the Jacobian).

.. _small-strain-formulation:

Small-strain formulation:
~~~~~~~~~~~~~~~~~~~~~~~~~

The small-strain version of hardening Mises plasticity can be combined with isotropic damage.
The basic equations include the additive decomposition of strain into elastic and plastic parts,

.. math::
   :label: VMadditiveDecomposition

   \veps = \veps_{\rm{e}} + \veps_{\rm{p}},

the stress strain law 

.. math::
   :label: VMconstitutiveEquation

   \vsig = (1-\omega)\bar{\vsig}=(1-\omega)\mbf{D} :(\veps-\veps_{\rm{p}}),

the definition of the yield function in terms of the effective stress,

.. math::
   :label: VMvonMisesCondition

   f(\bar{\vs},\kappa) = \sqrt{\frac{3}{2}\bar{\vs}:\bar{\vs}}-\sigma_Y(\kappa) = \sqrt{3 J_2(\bar{\vsig})}-\sigma_Y(\kappa),

the incremental definition of cumulative plastic strain

.. math::
   :label: VMcumPlasStrain

   \dot{\kappa} = \| \epspd \|,

the linear hardening law (for :param:`htype` = 0)

.. math::
   :label: VMlinearHardeningLaw

   \sigma_Y(\kappa) = \sigma_0 + H\kappa,

the evolution law for the plastic strain

.. math::

   \epspd = \dot{\lambda}\frac{\partial f}{\partial \bar{\vs}},

the loading-unloading conditions

.. math::
   :label: VMkuhnTucker

   \dot{\lambda} > 0 \qquad f(\bar{\vs},\kappa)\leq 0 \qquad \dot{\lambda} f(\bar{\vs},\kappa) = 0,

and the damage law 

.. math::
   :label: damagelawmp

   \omega(\kappa) = \omega_c(1-\mbox{e}^{-a\kappa}),

In the equations above, :math:`\veps` is the strain tensor, :math:`\veps_{\rm{e}}` is the elastic strain tensor, :math:`\veps_{p}` is the plastic strain tensor, :math:`\mbf{D}` is the elastic stiffness tensor, :math:`\vsig` is the nominal stress tensor, :math:`\bar{\vsig}` is the effective stress tensor, :math:`\bar{š}` is the effective deviatoric stress tensor, :math:`\sigma_Y` is the magnitude of stress at yielding under uniaxial tension (or compression), :math:`\kappa` is the cumulated plastic strain, :math:`H` is the hardening modulus, :math:`\lambda` is the plastic multiplier, :math:`\omega` is the damage variable, :math:`\omega_c` is critical damage and :math:`a` is a positive dimensionless parameter.

Large-strain formulation
~~~~~~~~~~~~~~~~~~~~~~~~

is based on the introduction of an intermediate local configuration, with respect to which the elastic response is characterized. This concept leads to a multiplicative decomposition of deformation gradient into elastic and plastic parts:

.. math::

   \mbf{F} = \mbf{F}^e\mbf{F}^p.

The stress-evaluation algorithm can be based on the 
classical radial return mapping; see [SimoHughes]_ for more details. Damage is not yet implemented in the large-strain version of the model.

The model description and parameters are summarized in :numref:`misesMat_table`.

.. table:: Mises plasticity -- summary.
   :name: misesMat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Mises plasticity model with isotropic hardening                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MisesMat` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}`                   |
   |                    | :elemparam:`n{rn}` :elemparam:`sig0{rn}` :elemparam:`H{rn}` :optelemparam:`htype{in}`        |
   |                    | :optelemparam:`h_eps{ra}` :optelemparam:`h(eps){ra}`                                         |
   |                    | :elemparam:`omega_crit{rn}`:elemparam:`a{rn}`                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression) (Required if htype =  |
   |                    |   0, which is default)                                                                       |
   |                    | - :param:`H` hardening modulus (can be negative in the case of plastic softening) (Required  |
   |                    |   if htype = 0, which is default)                                                            |
   |                    | - :param:`htype` hardening type (Optional parameter. Default = 0)                            |
   |                    | - :param:`h_eps` array of plastic strains (Required if htype = 1)                            |
   |                    | - :param:`h(eps)` array of yield stresses (Required if htype = 1)                            |
   |                    | - :param:`omega_crit` critical damage in damage law (:eq:`damagelawmp`)                      |
   |                    | - :param:`a` exponent in damage law (:eq:`damagelawmp`)                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 1dMat, PlaneStrain, 3dMat, 3dMatF                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/Mises01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/Mises01.in>`_, |
   |                    | `sm/Mises02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/Mises02.in>`_, |
   |                    | `sm/bond_link_2.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_2.in>`_,            |
   |                    | `sm/libeam3dboundary.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/libeam3dboundary.in>`_        |
   +--------------------+----------------------------------------------------------------------------------------------+

VTKxml output can report Mises stress, which equals to :math:`\sqrt{3J_2}`. When no hardening/softening exists, Mises stress reaches values up to given uniaxial yield stress :param:`sig0`.

.. _MisesMatNl:

Mises plasticity model with isotropic damage, nonlocal - MisesMatNl, MisesMatGrad
---------------------------------------------------------------------------------

The small-strain version of the model is regularized by the over-nonlocal formulation with damage driven by a combination of local and nonlocal cumulated plastic strain,

.. math::
   :label: overKappa1

   \hat{\kappa} = (1-m)\kappa + m\bar{\kappa},

where :math:`m` is a dimensionless material parameter (typically :math:`m>1`) and :math:`\bar{\kappa}` is the nonlocal cumulated plastic strain, which is evaluated either using the integral approach,
or using the implicit gradient approach.

.. _IntegralNonlocalFormulation:

Integral nonlocal formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One possible regularization technique is based on the integral definition of nonlocal
cumulated plastic strain

.. math::
   :label: nonLocalPlStrain1

   \bar{\kappa}(x) = \int\limits_V \alpha(x,s)\kappa(s)\,{\rm d}s

The nonlocal weight function is usually defined as

.. math::

   \alpha(x,s) = \frac{\alpha_0(\Vert x-s\Vert )}{\int\limits_V\alpha_0(\Vert x-t\Vert )\,{\rm d}t}

where

.. math::
   :label: alpha0

   \alpha_0(r) = \begin{cases} \left(1-\frac{r^2}{R^2}\right)^2 &\text{if $r<R$}\\ 
   \\
   0 & \text{if $r \ge R$}
   \end{cases}

is a nonnegative function, for :math:`r<R` monotonically decreasing with increasing distance :math:`r=\Vert x-s\Vert`, and :math:`V` denotes the domain occupied by the investigated material body.
The key idea is that the damage evolution at a certain point depends not only on the cumulated plastic strain at that point, but also on points at distances smaller than the interaction radius :math:`R`, considered as a new material parameter.

.. _ImplicitGradientFormulation:

Implicit gradient formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gradient formulation can be conceived as the differential counterpart to the integral formulation. The nonlocal cumulated plastic strain is computed from a Helmholtz-type differential equation 

.. math::
   :label: implicitGradient

   \bar{\kappa} - l^2\nabla^2\bar{\kappa} = \kappa

with homogeneous Neumann boundary condition

.. math::
   :label: neumannBC

   \frac{\partial\bar{\kappa}}{\partial n} = 0.

In (:eq:`implicitGradient`), :math:`l` is the length scale parameter and :math:`\nabla` is the Laplace operator.

The model description and parameters are summarized in :numref:`misesMatNl_table` and :numref:`misesMatGrad_table`. Note that the internal length parameter :param:`r` has the meaning of the radius of interaction :math:`R` for the integral version (and thus has the dimension of length) but for the gradient version it has the meaning of the coefficient :math:`l^2` multiplying the Laplacean in :eq:`implicitGradient`, and thus has the dimension of length squared.

.. table:: Nonlocal integral Mises plasticity -- summary.
   :name: misesMatNl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal Mises plasticity with isotropic hardening                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MisesMatNl` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}`                 |
   |                    | :elemparam:`n{rn}` :elemparam:`sig0{rn}` :elemparam:`H{rn}`                                  |
   |                    | :elemparam:`omega_crit{rn}`:elemparam:`a{rn}`:elemparam:`r{rn}`:elemparam:`m{rn}`\           |
   |                    | [:elemparam:`wft{in}`][:elemparam:`scalingType{in}`]                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression)                       |
   |                    | - :param:`H` hardening modulus                                                               |
   |                    | - :param:`omega_crit` critical damage                                                        |
   |                    | - :param:`a` exponent in damage law                                                          |
   |                    | - :param:`r` nonlocal interaction radius :math:`R` from eq. (:eq:`alpha0`)                   |
   |                    | - :param:`m` over-nonlocal parameter                                                         |
   |                    | - :param:`wft` selects the type of nonlocal weight function (see Section                     |
   |                    |   :ref:`sec:nidm`): ; 1 - default, quartic spline (bell-shaped function with bounded         |
   |                    |   support) ; 2 - Gaussian function ; 3 - exponential function (Green function in 1D) ; 4 -   |
   |                    |   uniform averaging up to distance :math:`R` ; 5 - uniform averaging over one finite element |
   |                    |   ; 6 - special function obtained by reducing the 2D exponential function to 1D (by          |
   |                    |   numerical integration)                                                                     |
   |                    | - :param:`scalingType` selects the type of scaling of the weight function (e.g. near a       |
   |                    |   boundary; see Section :ref:`sec:nidm`): ; 1 - default, standard scaling with integral      |
   |                    |   of weight function in the denominator ; 2 - no scaling (the weight function normalized in  |
   |                    |   an infinite body is used even near a boundary) ; 3 - Borino scaling (local complement)     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 1dMat, PlaneStrain, 3dMat                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+

.. table:: Gradient-enhanced Mises plasticity -- summary.
   :name: misesMatGrad_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Gradient-enhanced Mises plasticity with isotropic damage                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`MisesMatGrad` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}`               |
   |                    | :elemparam:`n{rn}` :elemparam:`sig0{rn}` :elemparam:`H{rn}`                                  |
   |                    | :elemparam:`omega_crit{rn}`:elemparam:`a{rn}`:elemparam:`r{rn}`:elemparam:`m{rn}`            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression)                       |
   |                    | - :param:`H` hardening modulus                                                               |
   |                    | - :param:`omega_crit` critical damage                                                        |
   |                    | - :param:`a` exponent in damage law                                                          |
   |                    | - :param:`r` internal length scale parameter :math:`l^2` from eq.                            |
   |                    |   (:eq:`implicitGradient`)                                                                   |
   |                    | - :param:`m` over-nonlocal parameter                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 1dMat, PlaneStrain, 3dMat                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _sec:rankplast:

Rankine plasticity model with isotropic damage and its nonlocal formulations - RankMat, RankMatNl, RankMatGrad
--------------------------------------------------------------------------------------------------------------

This model has a very similar structure to the model described in Section :ref:`sec:misplast`,
but is based on the Rankine yield condition. It is available in the small-strain version
only, and so far exclusively for plane stress analysis. 
The basic equations :eq:`VMadditiveDecomposition`--:eq:`VMconstitutiveEquation` 
and :eq:`VMcumPlasStrain`--:eq:`damagelawmp` remain valid, 
and the yield function :eq:`VMvonMisesCondition` is redefined as

.. math::
   :label: RankineCondition

   f(\bar{\vsig},\kappa) = \max_{I}\bar{\sigma_I}-\sigma_Y(\kappa)

where :math:`\bar{\sigma_I}` are the principal values of the effective stress tensor :math:`\bar{\vsig}`.
The hardening law can have either the linear form :eq:`VMlinearHardeningLaw`, or the exponential form

.. math::
   :label: VMexpHardeningLaw

   \sigma_Y(\kappa) = \sigma_0 + \Delta\sigma_Y\left(1-\exp(-H\kappa/\Delta\sigma_Y)\right),

where :math:`H` is now the initial plastic modulus and :math:`\Delta\sigma_Y` is the value of yield stress
increment asymptotically approached as :math:`\kappa\rightarrow\infty`.
In damage law :eq:`damagelawmp`, parameter :math:`\omega_c` is always set to 1.
If the plastic hardening is linear, the user can specify either the exponent :math:`a` from :eq:`damagelawmp`,
or the dissipation per unit volume, :math:`g_f`, which represents the area under the
stress-strain diagram (and parameter :math:`a` is then determined automatically such that
the area under the diagram has the prescribed value).
For exponential plastic hardening, the evaluation of :math:`a` from :math:`g_f` is not properly
implemented and it is better to specify :math:`a` directly.

The model description and parameters are summarized in :numref:`rankineMat_table`--:numref:`rankineMatGrad_table`. Note that the default value of parameter :math:`m` is equal to 1 for the
integral model but for the gradient model it is equal to 2. 
Also note that the internal length parameter :param{r} has the meaning of the
radius of interaction :math:`R` for the integral version (and thus has the dimension
of length) but for the gradient version it has the meaning of the coefficient :math:`l^2`
multiplying the Laplacean in :eq:`implicitGradient`, and thus has the dimension of length squared.

For the gradient model, it is possible to specify parameter :param:`negligible\_damage`, which sets the minimum value of damage that is considered as nonzero.

The approximate solution of Helmholtz equation :math:`(``\ implicitGradient``)` can lead to very small but nonzero nonlocal kappa at some points that are actually elastic. If such small values are positive, they lead to a very small but nonzero damage. If this is interpreted as "loading", the tangent terms are activated, but damage will not actually grow at such points and the convergence rate is slowed down. It is better to consider such points as elastic. By default, :param:`negligible\_damage` is set to 0, but it is recommended to set it e.g.\ to 1.e-6.

.. table:: Rankine plasticity -- summary.
   :name: rankineMat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Rankine plasticity with isotropic hardening and damage                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RankMat` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}` |
   |                    | :elemparam:`plasthardtype{in}` :elemparam:`sig0{rn}` :elemparam:`H{rn}`                      |
   |                    | :elemparam:`delSigY{rn}` :elemparam:`yieldtol{rn}` (:elemparam:`gf{rn}` :math:`\|`           |
   |                    | :elemparam:`a{rn}`)                                                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`plasthardtype` type of plastic hardening (0=linear=default, 1=exponential)         |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression)                       |
   |                    | - :param:`H` initial hardening modulus (default value 0.)                                    |
   |                    | - :param:`delSigY` final increment of yield stress (default value 0., needed only if         |
   |                    |   plasthardtype=1)                                                                           |
   |                    | - :param:`yieldtol` relative tolerance in the yield condition                                |
   |                    | - :param:`gf` dissipation per unit volume                                                    |
   |                    | - :param:`a` exponent in damage law (:eq:`damagelawmp`)                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | PlaneStress                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+

.. table:: Nonlocal integral Rankine plasticity -- summary.
   :name: rankineMatNl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal Rankine plasticity with isotropic hardening and damage                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RankMatNl` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}`                  |
   |                    | :elemparam:`n{rn}` :elemparam:`plasthardtype{in}` :elemparam:`sig0{rn}` :elemparam:`H{rn}`   |
   |                    | :elemparam:`delSigY{rn}` :elemparam:`yieldtol{rn}` (:elemparam:`gf{rn}` :math:`\|`           |
   |                    | :elemparam:`a{rn}`):elemparam:`r{rn}`:elemparam:`m{rn}`\                                     |
   |                    | [:elemparam:`wft{in}`][:elemparam:`scalingType{in}`]                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`plasthardtype` type of plastic hardening (0=linear=default, 1=exponential)         |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression)                       |
   |                    | - :param:`H` initial hardening modulus (default value 0.)                                    |
   |                    | - :param:`delSigY` final increment of yield stress (default value 0.)                        |
   |                    | - :param:`yieldtol` relative tolerance in the yield condition                                |
   |                    | - :param:`gf` dissipation per unit volume                                                    |
   |                    | - :param:`a` exponent in damage law (:eq:`damagelawmp`)                                      |
   |                    | - :param:`r` internal length scale parameter :math:`l^2` from eq.                            |
   |                    |   (:eq:`implicitGradient`)                                                                   |
   |                    | - :param:`m` over-nonlocal parameter (default value 1.)                                      |
   |                    | - :param:`wft` selects the type of nonlocal weight function (see Section                     |
   |                    |   :ref:`sec:nidm`): ; 1 - default, quartic spline (bell-shaped function with bounded         |
   |                    |   support) ; 2 - Gaussian function ; 3 - exponential function (Green function in 1D) ; 4 -   |
   |                    |   uniform averaging up to distance :math:`R` ; 5 - uniform averaging over one finite element |
   |                    |   ; 6 - special function obtained by reducing the 2D exponential function to 1D (by          |
   |                    |   numerical integration)                                                                     |
   |                    | - :param:`scalingType` selects the type of scaling of the weight function (e.g. near a       |
   |                    |   boundary; see Section :ref:`sec:nidm`): ; 1 - default, standard scaling with integral      |
   |                    |   of weight function in the denominator ; 2 - no scaling (the weight function normalized in  |
   |                    |   an infinite body is used even near a boundary) ; 3 - Borino scaling (local complement)     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | PlaneStress                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+

.. table:: Gradient-enhanced Rankine plasticity -- summary.
   :name: rankineMatGrad_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Gradient-enhanced Rankine plasticity with isotropic hardening and damage                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RankMatGrad` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}`                |
   |                    | :elemparam:`n{rn}` :elemparam:`plasthardtype{in}` :elemparam:`sig0{rn}` :elemparam:`H{rn}`   |
   |                    | :elemparam:`delSigY{rn}` :elemparam:`yieldtol{rn}` (:elemparam:`gf{rn}` :math:`\|`           |
   |                    | :elemparam:`a{rn}`):elemparam:`r{rn}`:elemparam:`m{rn}`:elemparam:`negligible_damage{rn}`    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`plasthardtype` type of plastic hardening (0=linear=default, 1=exponential)         |
   |                    | - :param:`sig0` initial yield stress in uniaxial tension (compression)                       |
   |                    | - :param:`H` hardening modulus (default value 0.)                                            |
   |                    | - :param:`delSigY` final increment of yield stress (default value 0.)                        |
   |                    | - :param:`yieldtol` relative tolerance in the yield condition                                |
   |                    | - :param:`gf` dissipation per unit volume                                                    |
   |                    | - :param:`a` exponent in damage law (:eq:`damagelawmp`)                                      |
   |                    | - :param:`r` internal length scale parameter :math:`l` from eq. (:eq:`implicitGradient`)     |
   |                    | - :param:`m` over-nonlocal parameter (default value 2.)                                      |
   |                    | - :param:`negligible_damage` optional parameter (default value 0.)                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | PlaneStress                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+


Perfectly plastic material with Mises yield condition - Steel1
--------------------------------------------------------------


This is an older model, kept here for compatibility with previous versions.
It uses Mises plasticity condition with no hardening and under small strain only.
The model description and parameters are summarized
in Tab.~:samp:`ref{Steel1_table}`. All its features are included in the model
described in Section~:samp:`ref{sec:misplast}`.

.. table:: Perfectly plastic material with Mises condition -- summary
   :name: Steel1_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Perfectly plastic material with Mises condition                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Steel1` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`E{rn}`                |
   |                    | :elemparam:`n{rn}` :elemparam:`tAlpha{rn}` :elemparam:`Ry{rn}`                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`Ry` uniaxial yield stress                                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer 3dBeam,      |
   |                    | PlaneStressRot                                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/sm/steel1.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/steel1.in>`_        |
   +--------------------+----------------------------------------------------------------------------------------------+


Composite plasticity model for masonry - Masonry02
--------------------------------------------------


Masonry is a composite material made of bricks and mortar. Nonlinear behavior of both components should be considered to obtain a realistic model able to describe cracking, slip, and crushing of the material. The model is based on paper by Lourenco and Rots [Rots]_. It is formulated on the basis of softening plasticity for tension, shear, and compression (see fig.(:numref:`compyieldsurffig`)). Numerical implementation is based on modern algorithmic concepts such as implicit integration of the rate equations and consistent tangent stiffness matrices.

.. figure:: /figures/constmodel.png
   :width: 80%
   :alt: Composite yield surface model for masonry
   :name: compyieldsurffig

   Composite yield surface model for masonry

The approach used in this work is based on the idea of concentrating all the damage in the relatively weak joints and, if necessary, in potential tension cracks in the bricks. The joint interface constitutive model should include all important damage mechanisms. Here, the concept of interface elements is used. An interface element allows to incorporate discontinuities in the displacement field and its behavior is described in terms of a relation between the tractions and relative displacement across the interface. In the present work, these quantities will be denoted as :math:`\sig`, generalized stress, and :math:`\e`, generalized strain. For 2D configuration, :math:`\sig=\{\sigma, \tau\}^T` and :math:`\e=\{u_n, u_s\}^T`, where :math:`\sigma` and :math:`\tau` are the normal and shear components of the traction interface vector and :math:`n` and :math:`s` subscripts distinguish between normal and shear components of displacement vector.

.. figure:: /figures/mmodel.png
   :width: 70%
   :alt: Modeling strategy for masonry

   Modeling strategy for masonry

The elastic response is characterized in terms of the elastic constitutive matrix :math:`\mbf{D}` as

.. math::

   \sig = \mbf{D} \e

For a 2D configuration, :math:`\mbf{D} = \text{diag}\{k_n, k_s\}`. The terms of the elastic stiffness matrix can be obtained from the properties of both masonry and joints as

.. math::

   k_n = \frac{E_b E_m}{t_m (E_b - E_m)};\; k_s = \frac{G_b G_m}{t_m (G_b - G_m)}

where :math:`E_b` and :math:`E_m` are Young's moduli, :math:`G_b` and :math:`G_m` are shear moduli for brick and mortar, and :math:`t_m` is the thickness of the joint. One should note that there is no contact algorithm assumed between bricks, this means that the overlap of neighboring units will be visible. On the other hand, the interface model includes a compressive cap, where the compressive inelastic behavior of masonry is lumped.

.. _tension-mode:

Tension mode
~~~~~~~~~~~~

In the tension mode, the exponential softening law is assumed (see :ref:`tensfig`). The yield function has the following form

.. math::

   f_1(\sig, \kappa_1) = \sigma - f_t(\kappa_1)

where the yield value :math:`f_t` is defined as

.. math::
   :label: ft

   f_t = f_{t0} \exp\left(-\frac{f_{t0}}{G^I_f} \kappa_1\right)

.. _tensfig:

.. figure:: /figures/tension.svg
   :width: 70%
   :alt: Tensile behavior of proposed model (:math:`f_t=0.2` MPa, :math:`G_f^I=0.018` N/mm)

   Tensile behavior of proposed model (:math:`f_t=0.2` MPa, :math:`G_f^I=0.018` N/mm)

The :math:`f_{t0}` represents the tensile strength of the joint or interface; and :math:`G^I_f` is the mode-I fracture energy. For the tension mode, the associated flow hypothesis is assumed.

.. _shear-mode:

Shear mode
~~~~~~~~~~

For the shear mode a Coulomb friction envelope is used. The yield function has the form

.. math::

   f_2(\sig,\kappa_2) = \vert\tau\vert+\sigma\tan\phi(\kappa_2)-c(\kappa_2)

According to [Rots]_ the variations of friction angle :math:`\phi` and cohesion :math:`c` are assumed as

.. math::
   :label: c

   c &= c_0\exp\left(-\frac{c_0}{G^{II}_f}\kappa_2\right)\\
   \tan\phi &= \tan\phi_0+(\tan\phi_r-\tan\phi_0)\left(\frac{c_0-c}{c_0}\right)

where :math:`c_0` is initial cohesion of joint, :math:`\phi_0` initial friction angle, :math:`\phi_r` residual friction angle, and :math:`G^{II}_f` fracture energy in mode II failure. A non-associated plastic potential :math:`g_2` is considered as

.. math::

   g_2=\vert\tau\vert+\sigma\tan\Phi-c

.. figure:: /figures/shearconf.svg
   :width: 70%
   :alt: Shear behavior of proposed model for different confinement levels in MPa (:math:`c_0=0.8\ \rm{MPa},\ \tan\phi_0=1.0,\ \tan\phi_r=0.75,{\rm and}\ G_f^{II}=0.05\ {N/mm}`)

   Shear behavior of proposed model for different confinement levels in MPa (:math:`c_0=0.8\ \rm{MPa},\ \tan\phi_0=1.0,\ \tan\phi_r=0.75,{\rm and}\ G_f^{II}=0.05\ {N/mm}`)

.. _coupling-of-tension-shear-modes:

Coupling of tension/shear modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tension and Coulomb friction modes are coupled with isotropic softening. This means that the percentage of softening in the cohesion is assumed to be the same as on the tensile strength

.. math::

   \dot\kappa_1=\lambda_1+\frac{G^I_f}{G^{II}_f}\frac{c_0}{f_{t0}}\lambda_2;\ \ \dot\kappa_2=\frac{G^{II}_f}{G^I_f}\frac{f_{t0}}{c_0}\lambda_1+\lambda_2

This follows from (:eq:`ft`) and (:eq:`c`). However, in the corner region, when both yield surfaces are activated, such approach will lead to a non-acceptable penalty. For this reason a quadratic combination is assumed

.. math::

   \dot\kappa_1=\sqrt{(\lambda_1)^2+\left(\frac{G^I_f}{G^{II}_f}\frac{c_0}{f_{t0}}\lambda_2\right)^2};\ \ \dot\kappa_2=\sqrt{\left(\frac{G^{II}_f}{G^I_f}\frac{f_{t0}}{c_0}\lambda_1\right)^2+(\lambda_2)^2}

.. _cap-mode:

Cap mode
~~~~~~~~

For the cap mode, an ellipsoid interface model is used. The yield condition is assumed as

.. math::

   f_3(\sig, \kappa_3) = C_{nn}\sigma^2+C_{ss}\tau^2 + C_n\sigma-\bar{\sigma}^2(\kappa_3)

where :math:`C_{nn},\ C_{ss}` and :math:`C_n` are material model parameters and :math:`\bar{\sigma}` is yield value, originally assumed in the following form of hardening/softening law [Rots]_

.. math::
   :label: hs3

   \bar{\sigma}_1(\kappa_3) &= \bar{\sigma}_i+(\bar{\sigma}_p-\bar{\sigma}_i)\sqrt{\frac{2\kappa_3}{\kappa_p}-\frac{\kappa_3^2}{\kappa_p^2}};\;\;\kappa_3\in(0,\kappa_p)\\
   \bar{\sigma}_2(\kappa_3) &= \bar{\sigma}_p+(\bar{\sigma}_m-\bar{\sigma}_p)\left(\frac{\kappa_3-\kappa_p}{\kappa_m-\kappa_p}\right)^2;\;\;\kappa_3\in(\kappa_p, \kappa_m)\\
   \bar{\sigma}_3(\kappa_3) &= \bar{\sigma}_r+(\bar{\sigma}_m-\bar{\sigma}_r)\exp\left(m\frac{\kappa_3-\kappa_m}{\bar{\sigma}_m-\bar{\sigma}_r}\right);\;\;\kappa_3\in(\kappa_m, \infty)

with :math:`m=2(\bar{\sigma}_m-\bar{\sigma}_p)/(\kappa_m-\kappa_p)`. The hardening/softening law (:eq:`hs3`) is shown in fig.(:numref:`hs3fig`). Note that the curved diagram is a :math:`C^1` continuous :math:`\sigma-\kappa_3` relation. The energy under the load-displacement diagram can be related to a “compressive fracture energy”.

The original hardening law (:eq:`hs3`.1) exhibits indefinite slope for :math:`\kappa_3=0`, which can cause the problems with numerical implementation. This has been overcomed by replacing this hardening law with parabolic equation given by

.. math::

   \bar{\sigma}_1(\kappa_3) = \bar{\sigma}_i-2*(\bar{\sigma}_i-\bar{\sigma}_p)*\frac{\kappa_3}{\kappa_p}+(\bar{\sigma}_i-\bar{\sigma}_p)\frac{\kappa_3}{\kappa_p}

An associated flow and strain hardening hypothesis are being considered. This yields

.. math::

   \dot\kappa_3=\lambda_3\sqrt{(2C_{nn}\sigma+C_n)*(2C_{nn}\sigma+C_n) + (2C_{ss}\tau)*(2C_{ss}\tau)}

.. figure:: /figures/capmode.png
   :width: 70%
   :alt: Hardening/softening law for cap mode
   :name: hs3fig

   Hardening/softening law for cap mode

The model parameters are summarized in :numref:`compomasonry1_table`.
There is one algorithmic issue, that follows from the model
formulation. Since the cap mode hardening/softening is not coupled to
hardening/softening of shear and tension modes the it may happen that
when the cap and shear modes are activated, the return directions
become parallel for both surfaces. This should be avoided by adjusting
the input parameters accordingly (one can modify dilatancy angle, for example).

.. table:: Composite model for masonry - summary.
   :name: compomasonry1_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Composite plasticity model for masonry                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Masonry02` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`E{rn}`             |
   |                    | :elemparam:`n{rn}` :elemparam:`ft0{rn}` :elemparam:`gfi{rn}` :elemparam:`gfii{rn}`           |
   |                    | :elemparam:`kn{rn}` :elemparam:`ks{rn}` :elemparam:`c0{rn}` :elemparam:`tanfi0{rn}`          |
   |                    | :elemparam:`tanfir{rn}` :elemparam:`tanpsi{rn}` :elemparam:`si{rn}` :elemparam:`sp{rn}`      |
   |                    | :elemparam:`sm{rn}` :elemparam:`sr{rn}` :elemparam:`kp{rn}` :elemparam:`km{rn}`              |
   |                    | :elemparam:`kr{rn}` :elemparam:`cnn{rn}` :elemparam:`css{rn}` :elemparam:`cn{rn}`            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`ft0` tensile strength                                                              |
   |                    | - :param:`gfi` fracture energy for mode I                                                    |
   |                    | - :param:`gfii` fracture energy for mode II                                                  |
   |                    | - :param:`kn` joint elastic property                                                         |
   |                    | - :param:`ks` joint elastic property                                                         |
   |                    | - :param:`c0` initial cohesion                                                               |
   |                    | - :param:`tanfi0` initial friction angle                                                     |
   |                    | - :param:`tanfir` residual friction angle                                                    |
   |                    | - :param:`tanpsi` dilatancy                                                                  |
   |                    | - :param:`{si, sp, sm, sr}` cap parameters :math:`\{\bar{\sigma}_i, \bar{\sigma}_p,          |
   |                    |   \bar{\sigma}_m, \bar{\sigma}_r\}`                                                          |
   |                    | - :param:`{kp, km,kr}` cap parameters :math:`\{\kappa_p, \kappa_m, \kappa_r\}`               |
   |                    | - :param:`cnn`,:param:`css`,:param:`cn` cap mode parametrs                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dInterface                                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+


.. _Rer:

Nonlinear elasto-plastic material model for concrete plates and shells - Concrete2
----------------------------------------------------------------------------------

Nonlinear elasto-plastic material model with hardening.
Takes into account uniaxial stress + transverse shear in concrete
layers with transverse stirrups.
Can be used only for 2d plates and shells with layered cross section
and together with explicit integration method (stiffness matrix is not
provided).
The model description and parameters are summarized
in :numref:`Rer_table`.

.. table:: Nonlinear elasto-plastic material model for concrete - summary.
   :name: Rer_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlinear elasto-plastic material model for concrete plates and shells                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Concrete2` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`E{rn}`             |
   |                    | :elemparam:`n{rn}` :elemparam:`SCCC{rn}` :elemparam:`SCCT{rn}` :elemparam:`EPP{rn}`          |
   |                    | :elemparam:`EPU{rn}` :elemparam:`EOPU{rn}` :elemparam:`EOPP{rn}` :elemparam:`SHEARTOL{rn}`   |
   |                    | :elemparam:`IS_PLASTIC_FLOW{in}` :elemparam:`IFAD{in}` :elemparam:`STIRR_E{rn}`              |
   |                    | :elemparam:`STIRR_Ft{rn}` :elemparam:`STIRR_A{rn}` :elemparam:`STIRR_TOL{rn}`                |
   |                    | :elemparam:`STIRR_EREF{rn}` :elemparam:`STIRR_LAMBDA{rn}`                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`SCCC` pressure strength                                                            |
   |                    | - :param:`SCCT` tension strength                                                             |
   |                    | - :param:`EPP` threshold effective plastic strain for softening in compression               |
   |                    | - :param:`EPU` ultimate eff. plastic strain                                                  |
   |                    | - :param:`EOPP` threshold volumetric plastic strain for softening in tension                 |
   |                    | - :param:`EOPU` ultimate volumetric plastic strain                                           |
   |                    | - :param:`SHEARTOL` threshold value of the relative shear deformation (psi**2/eef) at which  |
   |                    |   shear is considered in layers. For lower relative shear deformations the transverse shear  |
   |                    |   remains elastic decoupled from bending. default value SHEARTOL = 0.01                      |
   |                    | - :param:`IS_PLASTIC_FLOW` indicates that plastic flow (not deformation theory) is used in   |
   |                    |   pressure                                                                                   |
   |                    | - :param:`IFAD` State variables will not be updated, otherwise update state variables        |
   |                    | - :param:`STIRR_E` Young modulus of stirrups                                                 |
   |                    | - :param:`STIRR_R` stirrups uniaxial strength = elastic limit                                |
   |                    | - :param:`STIRR_A` stirrups area/unit length (beam) or /unit area (shell)                    |
   |                    | - :param:`STIRR_TOL` stirrups tolerance of equilibrium in the z direction (=0 no iteration)  |
   |                    | - :param:`STIRR_EREF` stirrups reference strain rate for Peryzna's material                  |
   |                    | - :param:`STIRR_LAMBDA` coefficient for that material (stirrups)                             |
   |                    | - :param:`SHTIRR_H` isotropic hardening factor for stirrups                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dShellLayer, 2dPlateLayer                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+

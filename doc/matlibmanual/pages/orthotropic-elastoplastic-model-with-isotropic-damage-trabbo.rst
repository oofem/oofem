Orthotropic elastoplastic model with isotropic damage - TrabBone3d
==================================================================


.. _local_formulation:

Local formulation
-----------------

The basic equations include an additive decomposition of total strain into elastic (reversible) part and plastic (irreversible) part

.. math::
   :label: additiveDecomposition

   \veps = \veps_{\rm{e}} + \veps_{\rm{p}}

the stress strain law

.. math::
   :label: constitutiveLaw

   \vsig = \left(1 - \omega\right)\bar{\vsig} = \left(1 - \omega\right)\mbf{D} : \veps_{\rm{e}}

the yield function

.. math::

   f(\bar{\vsig}, \kappa) = \sqrt{\bar{\vsig} : \mbf{F} : \bar{\vsig}} - \sigma_Y(\kappa)

loading-unloading conditions

.. math::

   f(\bar{\vsig}, \kappa) \leq 0 \qquad \dot{\lambda} \geq 0 \qquad \dot{\lambda} f(\bar{\vsig}, \kappa) = 0

evolution law for plastic strain

.. math::

   \dot{\veps}_{\rm{p}} = \dot{\lambda} \frac{\partial f}{\partial \bar{\vsig}}

the incremental definition of cumulated plastic strain

.. math::

   \dot{\kappa} = \| \epspd \|

the law governing the evolution of the damage variable

.. math::

   \omega(\kappa) = \omega_c (1 - \mbox{e}^{-a\kappa})

and the hardening law

.. math::

   \sigma_Y(\kappa) = 1 + \sigma_H (1 - \mbox{e}^{-s\kappa})

In the equations above, :math:`\bar{\vsig}` is the effective stress tensor, :math:`\mbf{D}` is the elastic stiffness tensor, :math:`f` is the yield function, :math:`\lambda` is the consistency parameter (plastic multiplier), :math:`\omega` is the damage variable, :math:`\sigma_Y` is the yield stress and :math:`s`, :math:`a`, :math:`\sigma_H` and :math:`\omega_c` are positive material parameters.

Material anisotropy is characterized by the second-order positive definite fabric tensor

.. math::
   :label: fabric

   \mbf{M} = \sum_{i=1}^3 m_i (\mbf{m}_i \otimes \mbf{m}_i)

normalized such that Tr\ :math:`(\mbf{M}) = 3`, :math:`m_i` are the eigenvalues and :math:`\mbf{m}_i` the eigenvectors.

The eigenvectors of the fabric tensor determine the directions of material orthotropy and the components of the elastic stiffness tensor :math:`\mbf{D}` are linked to eigenvalues of the fabric tensor. In the coordinate system aligned with :math:`m_i`, :math:`i=1,2,3`,

the stiffness can be presented in Voigt (engineering) notation as

.. math::

   \mbf{D} = \left[\begin{array}{cccccc}
   \frac{1}{E_1} & -\frac{\nu_{12}}{E_1} & -\frac{\nu_{13}}{E_1} & 0 & 0 & 0 \\
   -\frac{\nu_{21}}{E_2} & \frac{1}{E_2} & -\frac{\nu_{23}}{E_2} & 0 & 0 & 0 \\
   -\frac{\nu_{31}}{E_3} & -\frac{\nu_{32}}{E_3} & \frac{1}{E_3} & 0 & 0 & 0 \\
   0 & 0 & 0 & \frac{1}{G_{23}} & 0 & 0 \\
   0 & 0 & 0 & 0 & \frac{1}{G_{13}} & 0 \\
   0 & 0 & 0 & 0 & 0 & \frac{1}{G_{12}}
   \end{array}\right]^{-1}

where :math:`E_i = E_0 \rho^k m_i^{2l}`, :math:`G_{ij} = G_0 \rho^k m_i^l m_j^l` and :math:`\nu_{ij} = \nu_0 \frac{m_i^l}{m_j^l}`. Here, :math:`E_0`, :math:`G_0` and :math:`\nu_0` are elastic constants characterizing the compact (poreless) material, :math:`\rho` is the volume fraction of solid phase and :math:`k` and :math:`l` are dimensionless exponents.

Similar relations as for the stiffness tensor are also postulated for the components of a fourth-order tensor :math:`\mbf{F}` that is used in the yield condition. The yield condition is divided into tensile and compressive parts. Tensor :math:`\mbf{F}` is different in each part of the effective stress space. This tensor is denoted :math:`\mbf{F}^{+}` in tensile part, characterized by :math:`\hat{\mbf{N}}:\bar{\vsig} \leq 0`, and :math:`\mbf{F}^{-}` in compressive part, characterized by :math:`\hat{\mbf{N}}:\bar{\vsig} \leq 0`, where

.. math::

   \hat{\mbf{N}} = \frac{\sum_{i=1}^{3} m_i^{-2q}}{\sqrt{\sum_{i=1}^3 m_i^{-4q}}}(\mbf{m}_i \otimes \mbf{m}_i)

.. math::
   :label: eq:Fpm

   \mbf{F^{\pm}}=\left[\begin{array}{cccccc}
   \frac{1}{\left({\sigma_{1}^{\pm}}\right)^2} & -\frac{\chi_{12}^{\pm}}{\left({\sigma_{1}^{\pm}}\right)^2}&-\frac{\chi_{13}^{\pm}}{\left({\sigma_{1}^{\pm}}\right)^2}& 0& 0&0\\
   -\frac{\chi_{21}^{\pm}}{\left({\sigma_{2}^{\pm}}\right)^2} & \frac{1}{\left({\sigma_{2}^{\pm}}\right)^2}&-\frac{\chi_{23}^{\pm}}{\left({\sigma_{2}^{\pm}}\right)^2}& 0& 0&0\\
   -\frac{\chi_{31}^{\pm}}{\left({\sigma_{3}^{\pm}}\right)^2} & -\frac{\chi_{32}^{\pm}}{\left({\sigma_{3}^{\pm}}\right)^2}&\frac{1}{\left({\sigma_{3}^{\pm}}\right)^2}& 0& 0&0\\
   0 & 0 & 0 & \frac{1}{\tau_{23}} & 0 & 0\\
   0 & 0 & 0 & 0 & \frac{1}{\tau_{13}} & 0\\
   0 & 0 & 0 & 0 & 0 & \frac{1}{\tau_{12}}
    \end{array}\right].

In the equation above :math:`\sigma_i^{\pm} = \sigma_0^{\pm}\rho^p m_i^{2q}` is uniaxial yield stress along the :math:`i`-th principal axis of orthotropy, :math:`\tau_{ij} = \tau_0 \rho^p m_i^q m_j^q` is the shear yield stress in the plane of orthotropy and :math:`\chi_{ij}^{\pm} = \chi_0^{\pm}\frac{m_i^{2q}}{m_j^{2q}}` is the so-called interaction coefficient, :math:`p` and :math:`q` are dimensionless exponents and parameters with subscript 0 are related to a fictitious material with zero porosity. The yield surface is continuously differentiable if the parameters values are constrained by the condition

.. math::
   :label: eq:yield_condition

   \frac{\chi_0^- +1}{(\sigma_0^-)^2} = \frac{\chi_0^+ +1}{(\sigma_0^+)^2}.

The model description and parameters are summarized in :numref:`trabbone_table`.

.. table:: Anisotropic elastoplastic model with isotropic damage - summary.
   :name: trabbone_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Anisotropic elastoplastic model with isotropic damage                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`TrabBone3d` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`eps0{rn}`              |
   |                    | :elemparam:`nu0{rn}` :elemparam:`mu0{rn}` :elemparam:`expk{rn}` :elemparam:`expl{rn}`        |
   |                    | :elemparam:`m1{rn}` :elemparam:`m2{rn}` :elemparam:`rho{rn}` :elemparam:`sig0pos{rn}`        |
   |                    | :elemparam:`sig0neg{rn}` :elemparam:`chi0pos{rn}` :elemparam:`chi0neg{rn}`                   |
   |                    | :elemparam:`tau0{rn}` :elemparam:`plashardfactor{rn}` :elemparam:`expplashard{rn}`           |
   |                    | :elemparam:`expdam{rn}` :elemparam:`critdam{rn}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`eps0` Young modulus (at zero porosity)                                             |
   |                    | - :param:`nu0` Poisson ratio (at zero porosity)                                              |
   |                    | - :param:`mu0` shear modulus of elasticity (at zero porosity)                                |
   |                    | - :param:`m1` first eigenvalue of the fabric tensor                                          |
   |                    | - :param:`m2` second eigenvalue of the fabric tensor                                         |
   |                    | - :param:`rho` volume fraction of solid phase                                                |
   |                    | - :param:`sig0pos` yield stress in tension                                                   |
   |                    | - :param:`sig0neg` yield stress in compression                                               |
   |                    | - :param:`tau0` yield stress in shear                                                        |
   |                    | - :param:`chi0pos` interaction coefficient in tension                                        |
   |                    | - :param:`plashardfactor` hardening parameter                                                |
   |                    | - :param:`expplashard` exponent in hardening law                                             |
   |                    | - :param:`expdam` exponent in damage law                                                     |
   |                    | - :param:`critdam` critical damage                                                           |
   |                    | - :param:`expk` exponent :math:`k` in the expression for elastic stiffness                   |
   |                    | - :param:`expl` exponent :math:`l` in the expression for elastic stiffness                   |
   |                    | - :param:`expq` exponent :math:`q` in the expression for tensor \mbfF                        |
   |                    | - :param:`expp` exponent :math:`p` in the expression for tensor \mbfF                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+


Nonlocal formulation - TrabBoneNL3d
-----------------------------------

The model is regularized by the over-nonlocal formulation with damage driven by a combination of local and nonlocal cumulated plastic strain

.. math::
   :label: eq:overKappa

   \hat{\kappa} = (1-m)\kappa + m\bar{\kappa},

where :math:`\bar{\kappa}` is the nonlocal contribution to the plastic strain. The nonlocal contribution is given by

.. math::
   :label: eq:nonlocal_contribution

   \bar{\kappa} = \int_{\Omega} \frac{1}{|x-y|} \kappa(y) \, dV(y),

where :math:`\Omega` is the domain of the material and :math:`|x-y|` is the distance between points :math:`x` and :math:`y`. The over-nonlocal parameter :math:`m` controls the degree of nonlocality.

The model description and parameters are summarized in :numref:`trabboneNl_table`.

.. table:: Nonlocal formulation of anisotropic elastoplastic model with isotropic damage -- summary.
   :name: trabboneNl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal anisotropic elastoplastic model with isotropic damage                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`TrabBoneNL3d` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`eps0{rn}`            |
   |                    | :elemparam:`nu0{rn}` :elemparam:`mu0{rn}` :elemparam:`expk{rn}` :elemparam:`expl{rn}`        |
   |                    | :elemparam:`m1{rn}` :elemparam:`m2{rn}` :elemparam:`rho{rn}` :elemparam:`sig0pos{rn}`        |
   |                    | :elemparam:`sig0neg{rn}` :elemparam:`chi0pos{rn}` :elemparam:`chi0neg{rn}`                   |
   |                    | :elemparam:`tau0{rn}` :elemparam:`plashardfactor{rn}` :elemparam:`expplashard{rn}`           |
   |                    | :elemparam:`expdam{rn}` :elemparam:`critdam{rn}` :elemparam:`m{rn}` :elemparam:`R{rn}`       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`eps0` Young modulus (at zero porosity)                                             |
   |                    | - :param:`nu0` Poisson ratio (at zero porosity)                                              |
   |                    | - :param:`mu0` shear modulus (at zero porosity)                                              |
   |                    | - :param:`m1` first eigenvalue of the fabric tensor                                          |
   |                    | - :param:`m2` second eigenvalue of the fabric tensor                                         |
   |                    | - :param:`rho` volume fraction of the solid phase                                            |
   |                    | - :param:`sig0pos` yield stress in tension                                                   |
   |                    | - :param:`tau0` yield stress in shear                                                        |
   |                    | - :param:`chi0pos` interaction coefficient in tension                                        |
   |                    | - :param:`chi0neg` interaction coefficient in compression                                    |
   |                    | - :param:`plashardfactor` hardening parameter                                                |
   |                    | - :param:`expplashard` exponent in the hardening law                                         |
   |                    | - :param:`expdam` exponent in the damage law                                                 |
   |                    | - :param:`critdam` critical damage                                                           |
   |                    | - :param:`expk` exponent :math:`k` in the expression for elastic stiffness                   |
   |                    | - :param:`expl` exponent :math:`l` in the expression for elastic stiffness                   |
   |                    | - :param:`expq` exponent :math:`q` in the expression for tensor \mbfF                        |
   |                    | - :param:`expp` exponent :math:`p` in the expression for tensor \mbfF                        |
   |                    | - :param:`m` over-nonlocal parameter                                                         |
   |                    | - :param:`R` nonlocal interaction radius                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat                                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+

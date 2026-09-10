Material models for tensile failure
===================================


.. _material-models-for-tensile-failure:

Material models for tensile failure
-----------------------------------

.. _nonlinear-elasto-plastic-material-model-for-concrete-plates-and-shells-concrete2:

Nonlinear elasto-plastic material model for concrete plates and shells - Concrete2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The description can be found in section :ref:`Rer`.


.. _smeared-rotating-crack-model-concrete3:

Smeared rotating crack model - Concrete3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of smeared rotating crack model.
Virgin material is modeled as isotropic linear elastic material
(described by Young modulus and Poisson ratio). The onset of cracking begins, when principal stress reaches
tensile strength.
Further behavior is then determined by softening law,
governed by principle of preserving of fracture
energy :math:`G_f`. For large elements, the tension strength can be
artificially reduced
to preserve fracture energy. Multiple cracks are allowed.
The elastic unloading and reloading is assumed.
In compression regime, this model correspond to isotropic linear elastic material.
The model description and parameters are summarized
in :numref:`rcm_table`.

.. table:: Rotating crack model for concrete - summary.
   :name: rcm_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Rotating crack model for concrete                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Concrete3` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`               |
   |                    | :elemparam:`Gf{rn}` :elemparam:`Ft{rn}` :elemparam:`exp_soft{in}` :elemparam:`tAlpha{rn}`    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`Gf` fracture energy                                                                |
   |                    | - :param:`Ft` tension strength                                                               |
   |                    | - :param:`exp_soft` determines the type of softening (0 = linear, 1 = exponential, 2 =       |
   |                    |   Hordijk)                                                                                   |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/sm/concrete_3point.in <https://github.com/oofem/oofem/blob/devel/tests/regression |
   |                    | /benchmark/sm/concrete_3point.in>`_, `benchmark/sm/concrete_3point_direct.in <https://github |
   |                    | .com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/concrete_3point_direct.in>`_       |
   +--------------------+----------------------------------------------------------------------------------------------+


.. _smeared-rotating-crack-model-with-transition-to-scalar-damage-linear-softening-rcsd:

Smeared rotating crack model with transition to scalar damage - linear softening - RCSD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of smeared rotating crack model with transition to
scalar damage with linear softening law.
Improves the classical rotating model (see
section ``rcm``) by introducing the transition to scalar damage
model in  later stages of tension softening.

Traditional smeared-crack models for concrete fracture are known to suffer by stress locking (meaning here spurious stress transfer across widely opening cracks), mesh-induced directional bias, and possible instability at late stages of the loading process. The combined model keeps the anisotropic character of the rotating crack but it does not transfer spurious stresses across widely open cracks. The new model with transition to scalar damage (RC-SD) keeps the anisotropic character of the RCM but it does not transfer spurious stresses across widely open cracks.

Virgin material is modeled as isotropic linear elastic material (described by Young modulus and Poisson ratio). The onset of cracking begins when principal stress reaches tensile strength. Further behavior is then determined by :math:`\bf{linear}` softening law, governed by the principle of preserving fracture energy :math:`G_f`. For large elements, the tension strength can be artificially reduced to preserve fracture energy. The transition to scalar damage model takes place when the softening stress reaches the specified limit. Multiple cracks are allowed. The elastic unloading and reloading is assumed. In compression regime, this model corresponds to isotropic linear elastic material. The model description and parameters are summarized in :numref:`rcsd_table`.

.. table:: RC-SD model for concrete - summary.
   :name: rcsd_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Smeared rotating crack model with transition to scalar damage - linear softening             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RCSD` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`                    |
   |                    | :elemparam:`Gf{rn}` :elemparam:`Ft{rn}` :elemparam:`sdtransitioncoeff{rn}`                   |
   |                    | :elemparam:`tAlpha{rn}`                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`Gf` fracture energy                                                                |
   |                    | - :param:`Ft` tension strength                                                               |
   |                    | - :param:`sdtransitioncoeff` determines the transition from RC to SD model. Transition takes |
   |                    |   plase when ratio of current softening stress to tension strength is less than              |
   |                    |   :param:`sdtransitioncoeff` value                                                           |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer              |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _rcsde:

Smeared rotating crack model with transition to scalar damage - exponential softening - RCSDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of smeared rotating crack model with transition to scalar damage with exponential softening law.
The description and model summary (:numref:`rcsde_table`) are the same as for the RC-SD model with linear softening law (see section ``rcsd``).

.. table:: RC-SD model for concrete - summary.
   :name: rcsde_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Smeared rotating crack model with transition to scalar damage - exponential softening        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RCSDE` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`                   |
   |                    | :elemparam:`Gf{rn}` :elemparam:`Ft{rn}` :elemparam:`sdtransitioncoeff{rn}`                   |
   |                    | :elemparam:`tAlpha{rn}`                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _rcsdnl:

Nonlocal smeared rotating crack model with transition to scalar damage - RCSDNL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of nonlocal version of smeared rotating crack model with transition to scalar damage.
Improves the classical rotating model (see section ``rcm``) by introducing the transition to scalar damage model in later stages of tension softening.
The improved RC-SD (see section ``rcsd``) is further extended to a nonlocal formulation, which not only acts as a powerful localization limiter but also alleviates mesh-induced directional bias. A special type of material instability arising due to negative shear stiffness terms in the rotating crack model is resolved by switching to SD mode. A bell-shaped nonlocal averaging function is used.

Virgin material is modeled as isotropic linear elastic material (described by Young modulus and Poisson ratio). The onset of cracking begins when principal stress reaches tensile strength.
Further behavior is then determined by :bf:`exponential` softening law.

The transition to scalar damage model takes place when the softening stress reaches the specified limit or when the loss of material stability due to negative shear stiffness terms that may arise in the standard RCM formulation, which takes place when the ratio of minimal shear coefficient in stiffness to bulk material shear modulus reaches the limit.

Multiple cracks are allowed.
The elastic unloading and reloading is assumed.
In compression regime, this model corresponds to isotropic linear elastic material.
The model description and parameters are summarized in :numref:`rcsdnl_table`.

.. table:: RCSDNL model for concrete - summary.
   :name: rcsdnl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal smeared rotating crack model with transition to scalar damage for concrete          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`RCSDNL` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`                  |
   |                    | :elemparam:`Ft{rn}` :elemparam:`sdtransitioncoeff{rn}` :elemparam:`sdtransitioncoeff2{rn}`   |
   |                    | :elemparam:`r{rn}` :elemparam:`tAlpha{rn}`                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young modulus                                                                   |
   |                    | - :param:`n` Poisson ratio                                                                   |
   |                    | - :param:`ef` deformation corresponding to fully open crack                                  |
   |                    | - :param:`Ft` tension strength                                                               |
   |                    | - :param:`sdtransitioncoeff` determines the transition from RC to SD model. Transition takes |
   |                    |   place when ratio of current softening stress to tension strength is less than              |
   |                    |   :param:`sdtransitioncoeff` value                                                           |
   |                    | - :param:`sdtransitioncoeff2` determines the transition from RC to SD model. Transition      |
   |                    |   takes place when ratio of current minimal shear stiffness term to virgin shear modulus is  |
   |                    |   less than :param:`sdtransitioncoeff2` value                                                |
   |                    | - :param:`r` parameter specifying the width of nonlocal averaging zone                       |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`regionMap` map indicating the regions (currently region is characterized by cross  |
   |                    |   section number) to skip for nonlocal avaraging. The elements and corresponding IP are not  |
   |                    |   taken into account in nonlocal averaging process if corresponding regionMap value is       |
   |                    |   nonzero.                                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat, 2dPlateLayer, 2dBeamLayer, 3dShellLayer              |
   +--------------------+----------------------------------------------------------------------------------------------+


.. _sec:idmtf:

Isotropic damage model for tensile failure - Idm1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This isotropic damage model assumes that the stiffness degradation is isotropic, i.e., stiffness moduli corresponding to different directions decrease proportionally and independently of the loading direction. The damaged stiffness tensor is expressed as :math:`\mbf{D}=(1-\omega)\mbf{D}_e` where :math:`\omega` is a scalar damage variable and :math:`\mbf{D}_e` is the elastic stiffness tensor. The damage evolution law is postulated in an explicit form, relating the damage variable :math:`\omega` to the largest previously reached equivalent strain level, :math:`\kappa`.

The equivalent strain, :math:`\tilde\varepsilon`, is a scalar measure derived from the strain tensor. The choice of the specific expression for the equivalent strain affects the shape of the elastic domain in the strain space and plays a similar role to the choice of a yield condition in plasticity.

The following definitions of **equivalent strain** are currently supported:

    - **Mazars** (1984) definition based on norm of positive part of strain:

      .. math::
         :label: mazars

         \tilde\varepsilon = \sqrt{\sum_{I=1}^{3} \langle\varepsilon_I\rangle^2}

      where :math:`\langle\varepsilon_I\rangle` are positive parts of principal values of the strain tensor :math:`\mbf{\varepsilon}`.

    - Definitions derived from the **Rankine** criterion of maximum principal stress:

      .. math::
         :label: rankinesmooth

         \tilde\varepsilon &= \frac{1}{E} \sqrt{\sum_{I=1}^{3} \langle\bar{\sigma}_I\rangle^2}
         \\
         \tilde\varepsilon &= \frac{1}{E} \max_{I=1}^{3} \bar{\sigma}_I

      where :math:`\bar{\sigma}_I`, :math:`I=1,2,3`, are the principal values of the effective stress tensor :math:`\bar{\mbf{\sigma}} = \mbf{D}_e:\mbf{\varepsilon}` and :math:`\langle\bar{\sigma}_I\rangle` are their positive parts.

    - **Energy** norm scaled by Young's modulus to obtain a strain-like quantity:

      .. math::
         :label: energynorm

         \tilde\varepsilon = \frac{1}{E} \sqrt{\mbf{\varepsilon}:\mbf{D}_e:\mbf{\varepsilon}}

    - **Modified Mises** definition, proposed by de Vree [Vree:95]_:

      .. math::
         :label: modifiedmises

         \tilde\varepsilon = \frac{(k-1)I_{1\varepsilon}}{2k(1-2\nu)} + \frac{1}{2k} \sqrt{\frac{(k-1)^2}{(1-2\nu)^2} I_{1\varepsilon}^2 + \frac{12 k J_{2\varepsilon}}{(1+\nu)^2}}

      where

      .. math::

         I_{1\varepsilon} = \sum_{I=1}^3 \varepsilon_I

      is the first strain invariant (trace of the strain tensor),

      .. math::

         J_{2\varepsilon} = \frac{1}{2} \sum_{I=1}^3 \varepsilon_I^2 - \frac{1}{6} I_{1\varepsilon}^2

      is the second deviatoric strain invariant,
      and :math:`k` is a model parameter that corresponds to the ratio between the uniaxial compressive strength :math:`f_c` and uniaxial tensile strength :math:`f_t`.

    - **Griffith** definition with a solution on inclined elipsoidal inclusion. This definition handles materials in pure tension and also in compression, where tensile stresses usually appear on specifically oriented elipsoidal inclusion. The derivation of Griffith's criterion is summarized in [Hoek]_. In implementation, first check if Rankine criterion applies

      .. math::
         :label: griffith1

         \tilde\varepsilon = \frac{1}{E} \max_{I=1}^{3} \bar{\sigma}_I

      and if not, use Griffith's solution with ordered principal stresses :math:`\sigma_1 > \sigma_3`. The optional parameter :param:`griff_n` is by default 8 and represents the uniaxial compression/tensile strength ratio.

.. math::

   \tilde\varepsilon = \frac{\partial}{\partial E} \cdot \frac{\partial}{\partial -(\sigma_1 - \sigma_3)^2}{\text{griff\_n}(\sigma_1 + \sigma_3)} 

Note that all these definitions are based on the three-dimensional description of strain (and stress). If they are used in a reduced problem, the strain components that are not explicitly provided by the finite element approximation are computed from the underlying assumptions and used in the evaluation of equivalent strain. For instance, in a plane-stress analysis, the out-of-plane component of normal strain is calculated from the assumption of zero out-of-plane normal stress (using standard Hooke's law).

Since the growth of damage usually leads to softening and may induce localization of the dissipative process, attention should be paid to proper regularization. The most efficient approach is based on a nonlocal formulation; see Section :ref:`sec:nidm`. If the model is kept local, the damage law should be adjusted according to the element size, in the spirit of the crack-band approach. When done properly, this ensures a correct dissipation of energy in a localized band of cracking elements, corresponding to the fracture energy of the material. For various numerical studies, it may be useful to specify the parameters of the damage law directly, independently of the element size. One should be aware that in this case the model would exhibit pathological sensitivity to the size of finite elements if the mesh is changed.

The following **damage laws** are currently implemented:

    :bullet: **Cohesive crack with exponential softening** postulates a relation between the normal stress :math:`\sigma` transmitted by the crack and the crack opening :math:`w` in the form

        .. math::

           \sigma = f_t \exp\left(-\frac{w}{w_f}\right)

    Here, :math:`f_t` is the tensile strength and :math:`w_f` is a parameter with the dimension of length (crack opening), which controls the ductility of the material. In fact, :math:`w_f = G_f / f_t` where :math:`G_f` is the mode-I fracture energy. In the context of the crack-band approach, the crack opening :math:`w` corresponds to the inelastic (cracking) strain :math:`\varepsilon_c` multiplied by the effective thickness :math:`h` of the crack band. The effective thickness :math:`h` is estimated by projecting the finite element onto the direction of the maximum principal strain (and stress) at the onset of damage. The inelastic strain :math:`\varepsilon_c` is the difference between the total strain :math:`\varepsilon` and the elastic strain :math:`\sigma / E`.

    For the damage model, we obtain

    .. math::

       \varepsilon_c = \varepsilon - \frac{\sigma}{E} = \varepsilon - (1 - \omega)\varepsilon = \omega\varepsilon

    and thus :math:`w = h\varepsilon_c = h\omega\varepsilon`. Substituting this into the cohesive law and combining with the stress-strain law for the damage model, we get a nonlinear equation

    .. math::

       (1 - \omega)E\varepsilon = f_t \exp\left(-\frac{h\omega\varepsilon}{w_f}\right)

    For a given strain :math:`\varepsilon`, the corresponding damage variable :math:`\omega` can be solved from this equation by Newton iterations. It can be shown that the solution exists and is unique for every :math:`\varepsilon \ge \varepsilon_0` provided that the element size :math:`h` does not exceed the limit size :math:`h_{max} = w_f / \varepsilon_0`. For larger elements, a local snapback in the stress-strain diagram would occur, which is not admissible. In terms of the material properties, :math:`h_{max}` can be expressed as :math:`EG_f / f_t^2`, which is related to Irwin's characteristic length.

The derivation has been performed for monotonic loading and uniaxial tension. Under general conditions, :math:`\varepsilon` is replaced by the internal variable :math:`\kappa`, which represents the maximum previously reached level of equivalent strain.

In the list of input variables, the tensile strength :math:`f_t` is not specified directly but through the corresponding strain at peak stress, :math:`\varepsilon_0 = f_t / E`, denoted by keyword **e0**. Another input parameter is the characteristic crack opening :math:`w_f`, denoted by keyword **wf**.

The derivative can be expressed explicitly

.. math::
   :label: eq:derivative

   \frac{\partial\omega}{\partial\varepsilon} = - \frac{(\omega \varepsilon_f  - \varepsilon_f) \exp \left ( \frac{\omega \varepsilon}{\varepsilon_f}\right ) - \omega \varepsilon_0}{\varepsilon_f \varepsilon \exp \left ( \frac{\omega \varepsilon}{\varepsilon_f}\right ) - \varepsilon_0 \varepsilon}

.. _cohesive-crack-with-linear-softening:

Cohesive crack with linear softening
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cohesive law is assumed to have a simpler linear form

.. math::

   \sigma = f_t\left(1-\frac{w}{w_f}\right)

The relation between damage and strain can then be derived from the cohesive law and substituting :math:`w=h \omega \varepsilon`

.. math::

   \sigma = (1-\omega)E\varepsilon =  f_t\left(1-\frac{h \omega \varepsilon}{w_f}\right),

which leads to explicit evaluation of the damage variable

.. math::

   \omega=\frac{1-\frac{\varepsilon_0}{\varepsilon}}{1-\frac{h\varepsilon_0}{w_f}}

and no iteration is needed. Parameter :math:`w_f`, denoted again by keyword *wf*, has now the meaning of crack opening at complete failure (zero cohesive stress) and is related to fracture energy by a modified formula :math:`w_f=2G_f/f_t`. The expression for maximum element size, :math:`h_{max}=w_f/\varepsilon_0`, remains the same as for cohesive law with exponential softening, but in terms of the material properties it is now translated as :math:`h_{max}=2EG_f/f_t^2`. The derivative with respect to :math:`\varepsilon` yields

.. math::

   \frac{\partial\omega}{\partial\varepsilon} = \frac{\varepsilon_0}{\varepsilon^2 \left ( 1-\frac{h \varepsilon_0}{w_f}\right )}

.. _cohesive-crack-with-bilinear-softening:

Cohesive crack with bilinear softening
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of properly transforming the crack opening into inelastic strain, the current implementation deals with a stress-strain diagram adjusted such that the areas marked in the right part of :numref:`idm_softening` are equal to the fracture energies :math:`G_f` and :math:`G_{ft}` divided by the element size. The third parameter defining the law is the strain :math:`\eps_k` at which the softening diagram changes slope. Since this strain is considered as fixed, the corresponding stress :math:`\sigma_k` depends on the element size and for small elements gets close to the tensile strength (the diagram then gets close to linear softening with fracture energy :math:`G_{ft}`).

.. _linear-softening-stress-strain-law:

Linear softening stress-strain law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The specified parameters :math:`\varepsilon_0` and :math:`\varepsilon_f`, denoted by keywords *e0* and *ef*, have the meaning of (equivalent) strain at peak stress and at complete failure. The linear relation between stress and strain on the softening branch is obtained with the damage law

.. math::

   \omega =\frac{\varepsilon_f}{\varepsilon_f-\varepsilon_0}\left(1-\frac{\varepsilon_0}{\varepsilon}\right)

Again, to cover general conditions, :math:`\varepsilon` is replaced by :math:`\kappa`.

.. math::

   \frac{\partial\omega}{\partial\varepsilon} = \frac{\varepsilon_0 \varepsilon_f}{\varepsilon^2 ( \varepsilon_f - \varepsilon_0 )}

.. _exponential-softening-stress-strain-law:

Exponential softening stress-strain law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The damage law has a modified dependence of damage on strain:

.. math::

   \omega =1-\frac{\varepsilon_0}{\varepsilon}\exp\left(-\frac{\varepsilon-\varepsilon_0}{\varepsilon_f-\varepsilon_0}\right)

.. math::

   \frac{\partial\omega}{\partial\varepsilon} = \left[ \frac{1}{\varepsilon(\varepsilon_f-\varepsilon_0)} + \frac{1}{\varepsilon^2} \right] \varepsilon_0 \exp \left (\frac{\varepsilon_0-\varepsilon}{\varepsilon_f-\varepsilon_0} \right )

.. _mazars-stress-strain-law:

Mazars stress-strain law
^^^^^^^^^^^^^^^^^^^^^^^^

The dependence of damage on strain is given by

.. math::

   \omega =1-\frac{(1-A_t)\varepsilon_0}{\varepsilon}-A_t\exp\left(B_t(\varepsilon-\varepsilon_0)\right)

.. _smooth-exponential-stress-strain-law:

Smooth exponential stress-strain law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The dependence of damage on strain is given by

.. math::

   \omega = 1-\exp\left(-\left(\frac{\varepsilon}{\varepsilon_0}\right)^{M_d}\right)

This leads to a stress-strain curve that immediately deviates from linearity (has no elastic part) and smoothly changes from hardening to softening, with tensile strength

.. math::

   f_t = E\varepsilon_0\left({\rm e}M_d\right)^{-1/M_d}

.. _extended-smooth-stress-strain-law:

Extended smooth stress-strain law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The damage law has a rather complicated form:

.. math::
   :label: eq:damageEvol

   \omega = \left\{ \begin{array}{ll} 
   1-\exp\left(-\displaystyle\frac{1}{m}\left(\dfrac{\eps}{\varepsilon_{\rm p}}\right)^m\right) & \mbox{if $\eps \leq \varepsilon_1$}\\[2mm]
   1- \dfrac{\varepsilon_3}{\eps} \exp \left(- \dfrac{\eps- \varepsilon_1}{\varepsilon_{\rm f} \left[1+ \left(\frac{\eps-\varepsilon_1}{\varepsilon_2}\right)^n\right]}\right) & \mbox{if $\eps > \varepsilon_1$} \end{array} 
   \right\}

The primary model parameters are the uniaxial tensile strength :math:`f_{\rm t}`, the strain at peak stress (under uniaxial tension) :math:`\varepsilon_{\rm p}`, and additional parameters :math:`\varepsilon_1`, :math:`\varepsilon_2` and :math:`n`, which control the post-peak part of the stress-strain law. In the input record, they are denoted by keywords *ft, ep, e1, e2, nd*. Other parameters that appear in (:eq:`eq:damageEvol`) can be derived from the condition of zero slope of the stress-strain curve at :math:`\kappa=\varepsilon_{\rm p}` and from the conditions of stress and stiffness continuity at :math:`\kappa=\varepsilon_1`:

.. math::
   :label: eqparams1

   m&=\frac{1}{\ln(E\varepsilon_{\rm p}/f_{\rm t})}\\
   \varepsilon_{\rm f}&=\frac{\varepsilon_1}{\left(\varepsilon_1/\varepsilon_{\rm p}\right)^m -1}\\
   \varepsilon_3&=\varepsilon_1\exp\left(-\frac{1}{m}\left(\frac{\varepsilon_1}{\varepsilon_p}\right)^m\right)

.. _trilinear-softening:

Trilinear softening
^^^^^^^^^^^^^^^^^^^

It is formulated for an embedded crack model and it is expressed in terms of crack opening (:math:`w`). This approach uses a trilinear softening diagram that helps reproducing the FRC behaviour; this diagram is defined by four points: :math:`t`, :math:`k`, :math:`r` and :math:`f` that are related to specific properties of concrete, fibres and the fibre proportion. In order to adapt the original formulation to the ``idm1``, which uses damage :math:`\omega` as the driving parameter of fracture, instead of crack opening :math:`w`, the correspondence between both must be taken into account. Figure :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon` shows the trilinear diagram expressed in terms of the crack opening :math:`w` (left) and in terms of the equivalent strain :math:`\varepsilon` (right).

.. figure:: /figures/trilinear_diagram_sigmaVSw_and_sigmaVSepsilon.png
   :align: center
   :alt: Trilinear softening diagram expressed in terms of the crack opening (w) and expressed in terms of strain.
   :name: fig:diagram_sigmaVSw_and_sigmaVSepsilon

   Trilinear softening diagram expressed in terms of the crack opening (:math:`w`) and expressed in terms of strain.

The inelastic strain :math:`\varepsilon_c` is the difference between the total strain :math:`\varepsilon` and the elastic strain :math:`\varepsilon/E`, thus:

.. math::

   \varepsilon_c = \varepsilon - \dfrac{\sigma}{E}

If :math:`\sigma` is obtained by using the damage parameter (:math:`\omega`), it can be expressed as follows:

.. math::
   :label: eq:sigma

   \sigma = (1-\omega) E \varepsilon

Therefore, the inelastic strain can be written as:

.. math::

   \varepsilon_c =  \varepsilon - \dfrac{1}{E}\left( (1-\omega) E \varepsilon\right) \quad \Rightarrow \quad \varepsilon_c = \varepsilon - (1- \omega) \varepsilon = \omega \varepsilon

Thus, crack opening (:math:`w`) is related to damage (:math:`\omega`) through (Note that :math:`\omega` represents damage and :math:`w` represents crack opening.):

.. math::
   :label: eq:w

   w = h \varepsilon_c = h \omega \varepsilon

where :math:`h` stands for the effective thickness of the crack band, which is estimated by projecting the finite element onto the direction of the maximum principal strain at the onset of damage.

Hereafter, the expressions of damage (:math:`\omega`) are obtained for each section of the softening diagram (before damage develops, damage between points :math:`t` and :math:`k`, between points :math:`k` and :math:`r`, between points :math:`r` and :math:`f`, and, finally, after damage has fully developed).

    - **Case 1: :math:`\varepsilon \leq \varepsilon_0`:**
      In this case, damage has not started, thus:

      .. math::

         \omega = 0

      And its derivative:

      .. math::

         \dfrac{\partial \omega}{\partial \varepsilon} = 0

    - **Case 2: :math:`\varepsilon_0 \le \varepsilon \leq \varepsilon_k`:**
      Referred to the :math:`\sigma-w` diagram:

      .. math::

         \sigma = f_t + w \dfrac{f_k - f_t}{w_k}

      Therefore, since :math:`f_t = E \varepsilon_0` and using expressions :eq:`eq:sigma` and :eq:`eq:w`:

      .. math::

         \omega = \dfrac{E}{E + h \left(\dfrac{f_k - f_t}{w_k}\right)} - \dfrac{E \varepsilon_0}{\varepsilon \left[E + h \left(\dfrac{f_k - f_t}{w_k}\right)\right]}

      And its derivative:

      .. math::

         \dfrac{\partial \omega}{\partial \varepsilon} = \dfrac{E \varepsilon_0}{E + h \left(\dfrac{f_k - f_t}{w_k}\right)} \cdot \dfrac{1}{\varepsilon^2}

    - **Case 3: :math:`\varepsilon_k \le \varepsilon \leq \varepsilon_r`:**
      Referred to the :math:`\sigma-w` diagram:

      .. math::

         \sigma = f_k + (w - w_k) \left(\dfrac{f_r - f_k}{w_r - w_k}\right)

      Therefore, using expressions :eq:`eq:sigma` and :eq:`eq:w`:

      .. math::

         \omega = \dfrac{E}{E + h \left(\dfrac{f_r - f_k}{w_r - w_k}\right)} + \dfrac{1}{\varepsilon} \cdot \dfrac{w_k \left(\dfrac{f_r - f_k}{w_r - w_k}\right) - f_k}{E + h \left(\dfrac{f_r - f_k}{w_r - w_k}\right)}

      And its derivative:

      .. math::

         \dfrac{\partial \omega}{\partial \varepsilon} = - \dfrac{w_k \left(\dfrac{f_r - f_k}{w_r - w_k}\right) - f_k}{E + h \left(\dfrac{f_r - f_k}{w_r - w_k}\right)} \cdot \dfrac{1}{\varepsilon^2}

    - **Case 4: :math:`\varepsilon_r \le \varepsilon \leq \varepsilon_f`:**
      Referred to the :math:`\sigma-w` diagram:

      .. math::

         \sigma = f_r + (w - w_r) \left(\dfrac{-f_r}{w_f - w_r}\right)

      Therefore, using expressions :eq:`eq:sigma` and :eq:`eq:w`:

      .. math::

         \omega = \dfrac{E}{E + h \left(\dfrac{-f_r}{w_f - w_r}\right)} + \dfrac{1}{\varepsilon} \cdot \dfrac{w_r \left(\dfrac{-f_r}{w_f - w_r}\right) - f_r}{E + h \left(\dfrac{-f_r}{w_f - w_r}\right)}

      And its derivative:

      .. math::

         \dfrac{\partial \omega}{\partial \varepsilon} = - \dfrac{w_r \left(\dfrac{-f_r}{w_f - w_r}\right) - f_r}{E + h \left(\dfrac{-f_r}{w_f - w_r}\right)} \cdot \dfrac{1}{\varepsilon^2}

    - **Case 5: :math:`\varepsilon \geq \varepsilon_f`:**
      In this case, damage is fully developed:

      .. math::

         \omega = 1

Note that parameter :math:`\it{damlaw}` determines which type of damage law should be used, but the adjustment for element size is done only if parameter :math:`\it{wf}` is specified for :math:`\it{damlaw}=0` or :math:`\it{damlaw}=1`. For other values of :math:`\it{damlaw}`, or if parameter :math:`\it{ef}` is specified instead of :math:`\it{wf}`, the stress-strain curve does not depend on element size and the model would exhibit pathological sensitivity to the mesh size. These cases are intended to be used in combination with a nonlocal formulation. An alternative formulation uses fracture energy to determine fracturing strain.

The model parameters are summarized in :numref:`id_table`. Figure :numref:`idm_softening` shows three modes of a softening law with corresponding variables.

.. figure:: /figures/Damage_material_diag.png
   :width: 99%
   :alt: Implemented stress-strain diagrams for isotropic damage material. Fracturing strain :math:`\varepsilon_f` and crack opening at zero stress :math:`w_f` are interrelated through effective thickness :math:`h` of the crack band. Note that exponential softening approach based on an exponential cohesive law is not exactly equivalent to the approach based on an exponential softening branch of the stress-strain diagram; see the detailed discussion of the damage laws.
   :name: idm_softening

   Implemented stress-strain diagrams for isotropic damage material. Fracturing strain :math:`\varepsilon_f` and crack opening at zero stress :math:`w_f` are interrelated through effective thickness :math:`h` of the crack band. Note that exponential softening approach based on an exponential cohesive law is not exactly equivalent to the approach based on an exponential softening branch of the stress-strain diagram; see the detailed discussion of the damage laws.


.. table:: Isotropic damage model for concrete in tension
   :name: id_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic damage model for concrete in tension                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Idm1` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`    |
   |                    | [:elemparam:`tAlpha{rn}`] [:elemparam:`equivstraintype{in}`] [:elemparam:`k{rn}`]            |
   |                    | [:elemparam:`damlaw{in}`] :elemparam:`e0{rn}` [:elemparam:`wf{rn}`] [:elemparam:`ef{rn}`]    |
   |                    | [:elemparam:`ek{rn}`] [:elemparam:`wk{rn}`] [:elemparam:`sk{rn}`] [:elemparam:`wkwf{rn}`]    |
   |                    | [:elemparam:`skft{rn}`] [:elemparam:`gf{rn}`] [:elemparam:`gft{rn}`] [:elemparam:`At{rn}`]   |
   |                    | [:elemparam:`Bt{rn}`] [:elemparam:`md{rn}`] [:elemparam:`ft{rn}`] [:elemparam:`ep{rn}`]      |
   |                    | [:elemparam:`e1{rn}`] [:elemparam:`e2{rn}`] [:elemparam:`nd{rn}`]                            |
   |                    | [:elemparam:`maxOmega{rn}`] [:elemparam:`checkSnapBack{rn}`]                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`tAlpha` thermal expansion coefficient                                              |
   |                    | - :param:`equivstraintype` allows to choose from different definitions of equivalent strain: |
   |                    |   ; 0 - default = Mazars, eq. (:eq:`mazars`) ; 1 - smooth Rankine, eq.                       |
   |                    |   (:eq:`rankinesmooth`) ; 2 - scaled energy norm, eq. (:eq:`energynorm`) ; 3 -               |
   |                    |   modified Mises, eq. (:eq:`modifiedmises`) ; 4 - standard Rankine, eq.                      |
   |                    |   (``rankinestandard``) ; 5 - elastic energy based on positive stress ; 6 - elastic          |
   |                    |   energy based on positive strain ; 7 - Griffith criterion eq. (``griffith2``)               |
   |                    | - :param:`k` ratio between uniaxial compressive and tensile strength, needed only if         |
   |                    |   equivstraintype=3, default value 1                                                         |
   |                    | - :param:`damlaw` allows to choose from different damage laws: ; 0 - exponential softening   |
   |                    |   (default) with parameters e0 and wf :math:`\vert` ef :math:`\vert` gf ; 1 - linear         |
   |                    |   softening with parameters e0 and wf :math:`\vert` ef :math:`\vert` gf ; 2 - bilinear       |
   |                    |   softening with (e0, gf, gft, ek) :math:`\vert` (e0, wk, sk, wf) :math:`\vert` (e0, wkwf,   |
   |                    |   skft, wf) :math:`\vert` (e0, gf, gft, wk) ; 3 - Hordijk softening (not implemented yet) ;  |
   |                    |   4 - Mazars damage law with parameters At and Bt ; 5 - smooth stress-strain curve with      |
   |                    |   parameters e0 and md ; 6 - disable damage (dummy linear elastic material) ; 7 - extended   |
   |                    |   smooth damage law (:eq:`eq:damageEvol`) with parameters ft, ep, e1, e2, nd ; 11 -          |
   |                    |   trilinear softening diagram with (e0, w_k, w_r, w_f, f_k, f_r)                             |
   |                    | - :param:`e0` strain at peak stress (for damage laws 0,1,2,3), limit elastic strain (for     |
   |                    |   damage law 4), characteristic strain (for damage law 5)                                    |
   |                    | - :param:`wf` parameter controling ductility, has the meaning of crack opening (for damage   |
   |                    |   laws 0 and 1)                                                                              |
   |                    | - :param:`ef` parameter controling ductility, has the meaning of strain (for damage laws 0   |
   |                    |   and 1)                                                                                     |
   |                    | - :param:`ek` strain at knee point in bilinear softening type (for damage law 2)             |
   |                    | - :param:`wk` crack opening at knee point in bilinear softening type (for damage law 2)      |
   |                    | - :param:`sk` stress at knee point in bilinear softening type (for damage law 2)             |
   |                    | - :param:`wkwf` ratio of wk/wf :math:`<0,1>` in bilinear softening type (for damage law 2)   |
   |                    | - :param:`skft` ratio of sk/ft :math:`<0,1>` in bilinear softening type (for damage law 2)   |
   |                    | - :param:`gf` fracture energy (for damage laws 0--2)                                         |
   |                    | - :param:`gft` total fracture energy (for damage law 2)                                      |
   |                    | - :param:`At` parameter of Mazars damage law, used only by law 4                             |
   |                    | - :param:`Bt` parameter of Mazars damage law, used only by law 4                             |
   |                    | - :param:`md` exponent used only by damage law 5, default value 1                            |
   |                    | - :param:`ft` tensile strength, used only by damage law 7                                    |
   |                    | - :param:`ep` strain at peak stress, used only by damage law 7                               |
   |                    | - :param:`e1` parameter used only by damage law 7                                            |
   |                    | - :param:`e2` parameter used only by damage law 7                                            |
   |                    | - :param:`nd` exponent used only by damage law 7                                             |
   |                    | - :param:`griff_n` uniaxial compression/tensile ratio for Griffith's criterion               |
   |                    | - :param:`maxOmega` maximum damage, used for convergence improvement (its value is between 0 |
   |                    |   and 0.999999 (default), and it affects only the secant stiffness but not the stress)       |
   |                    | - :param:`checkSnapBack` parameter for snap back checking, 0 no check, 1 check (default)     |
   |                    | - :param:`w_k` crack opening of point :math:`k` in the trilinear diagram (see Fig.           |
   |                    |   :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)                                         |
   |                    | - :param:`w_r` crack opening of point :math:`r` in the trilinear diagram (see Fig.           |
   |                    |   :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)                                         |
   |                    | - :param:`w_f` crack opening of point :math:`f` in the trilinear diagram (see Fig.           |
   |                    |   :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)                                         |
   |                    | - :param:`f_k` cohesive stress of point :math:`k` in the trilinear diagram (see Fig.         |
   |                    |   :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)                                         |
   |                    | - :param:`f_r` cohesive stress of point :math:`r` in the trilinear diagram (see Fig.         |
   |                    |   :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           | Adaptivity support                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `partests/bar/bar.oofem.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/partests/bar/bar.oofem.in>`_,    |
   |                    | `sm/adapt01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/adapt01.in>`_, |
   |                    | `sm/adapt02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/adapt02.in>`_, |
   |                    | `sm/idm01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm01.in>`_,     |
   |                    | `sm/idm02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm02.in>`_,     |
   |                    | `sm/idm03.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm03.in>`_,     |
   |                    | `sm/idm04.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm04.in>`_,     |
   |                    | `sm/idm05.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm05.in>`_,     |
   |                    | `sm/idm06.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm06.in>`_,     |
   |                    | `sm/idm11.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm11.in>`_ (and |
   |                    | 4 more)                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+

:param k: Ratio between uniaxial compressive and tensile strength, needed only if equivstraintype=3, default value 1
:param damlaw: Allows to choose from different damage laws:

    * 0 - Exponential softening (default) with parameters e0 and wf :math:`|ef| |gf|`
    * 1 - Linear softening with parameters e0 and wf :math:`|ef| |gf|`
    * 2 - Bilinear softening with (e0, gf, gft, ek) :math:`|wk| |sk| |wf| |wkwf| |skft| |wk|`
    * 3 - Hordijk softening (not implemented yet)
    * 4 - Mazars damage law with parameters At and Bt
    * 5 - Smooth stress-strain curve with parameters e0 and md
    * 6 - Disable damage (dummy linear elastic material)
    * 7 - Extended smooth damage law :math:`(``\ eq:damageEvol``)` with parameters ft, ep, e1, e2, nd
    * 11 - Trilinear softening diagram with (e0, w_k, w_r, w_f, f_k, f_r)

:param e0: Strain at peak stress (for damage laws 0,1,2,3), limit elastic strain (for damage law 4), characteristic strain (for damage law 5)
:param wf: Parameter controlling ductility, has the meaning of crack opening (for damage laws 0 and 1)
:param ef: Parameter controlling ductility, has the meaning of strain (for damage laws 0 and 1)
:param ek: Strain at knee point in bilinear softening type (for damage law 2)
:param wk: Crack opening at knee point in bilinear softening type (for damage law 2)
:param sk: Stress at knee point in bilinear softening type (for damage law 2)
:param wkwf: Ratio of wk/wf :math:`<0,1>` in bilinear softening type (for damage law 2)
:param skft: Ratio of sk/ft :math:`<0,1>` in bilinear softening type (for damage law 2)
:param gf: Fracture energy (for damage laws 0--2)
:param gft: Total fracture energy (for damage law 2)
:param At: Parameter of Mazars damage law, used only by law 4
:param Bt: Parameter of Mazars damage law, used only by law 4
:param md: Exponent used only by damage law 5, default value 1
:param ft: Tensile strength, used only by damage law 7
:param ep: Strain at peak stress, used only by damage law 7
:param e1: Parameter used only by damage law 7
:param e2: Parameter used only by damage law 7
:param nd: Exponent used only by damage law 7
:param griff_n: Uniaxial compression/tensile ratio for Griffith's criterion
:param maxOmega: Maximum damage, used for convergence improvement (its value is between 0 and 0.999999 (default), and it affects only the secant stiffness but not the stress)
:param checkSnapBack: Parameter for snap back checking, 0 no check, 1 check (default)
:param w_k: Crack opening of point k in the trilinear diagram (see :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)
:param w_r: Crack opening of point r in the trilinear diagram (see :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)
:param w_f: Crack opening of point f in the trilinear diagram (see :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)
:param f_k: Cohesive stress of point k in the trilinear diagram (see :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)
:param f_r: Cohesive stress of point r in the trilinear diagram (see :numref:`fig:diagram_sigmaVSw_and_sigmaVSepsilon`)

Supported modes: 3dMat, PlaneStress, PlaneStrain, 1dMat

Features: Adaptivity support

Isotropic damage model for tensile failure -- summary.
:name: id_table

.. _sec:nidm:

Nonlocal isotropic damage model for tensile failure - Idmnl1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Nonlocal version of isotropic damage model from Section :ref:`sec:idmtf`.
The nonlocal averaging acts as a powerful localization limiter.

In the standard version of the model, damage is driven by the nonlocal equivalent strain :math:`\bar{\varepsilon}`, defined as a weighted average of the local equivalent strain:

.. math::

   \bar{\varepsilon}(\mbf{x}) = \int_V \alpha(\mbf{x},\mbf{\xi}) \tilde{\varepsilon}(\mbf{\xi}) \; {\rm d}\mbf{\xi}

In the “undernonlocal” formulation, the damage-driving variable is a combination of local and nonlocal equivalent strain, :math:`m\bar{\varepsilon}+(1-m)\tilde{\varepsilon}`, where :math:`m` is a parameter between 0 and 1. (If :math:`m>1`, the formulation is called “overnonlocal”; this case is useful for nonlocal plasticity but not for nonlocal damage.)

Instead of averaging the equivalent strain, one can average the compliance variable :math:`\gamma`, directly related to damage according to the formula :math:`\gamma=\omega/(1-\omega)`.

The weight function :math:`\alpha` contains a certain parameter with the dimension of length, which is in general called the characteristic length. Its specific meaning depends on the type of weight function. The following functions are currently supported:

    - 
      Truncated quartic spline, also called the bell-shaped function,

      .. math::

         \alpha_0(s) = \left\langle 1-\frac{s^2}{R^2}\right\rangle^2

      where :math:`R` is the interaction radius (characteristic length) and :math:`s` is the distance between the interacting points. This function is exactly zero for :math:`s\ge R`, i.e., it has a bounded support.
    - 
      Gaussian function

      .. math::

         \alpha_0(s) = \exp\left(-\frac{s^2}{R^2}\right)

      which is theoretically nonzero for an arbitrary large :math:`s` and thus has an unbounded support. However, in the numerical implementation the value of :math:`\alpha_0` is considered as zero for :math:`s>2.5R`.
    - 
      Exponential function

      .. math::

         \alpha_0(s) = \exp\left(-\frac{s}{R}\right)

      which also has an unbounded support, but is considered as zero for :math:`s>6R`. This function is sometimes called the Green function, because in 1D it corresponds to the Green function of the Helmholtz-like equation used by implicit gradient approaches.
    - 
      Piecewise constant function

      .. math::

         \alpha_0(s) =\left\{\begin{array}{cc} 1 & \mbox{ if } s\le R \\ 0 & \mbox{ if } s> R \end{array}\right.

      which corresponds to uniform averaging over a segment, disc or ball of radius :math:`R`.
    - 
      Function that is constant over the finite element in which point :math:`\mbf{x}` is located, and is zero everywhere else. Of course, this is not a physically objective definition of nonlocal averaging, since it depends on the discretization. However, this kind of averaging was proposed in a boundary layer by Prof.\ Bažant and was implemented into OOFEM for testing purposes.
    - 
      Special function 

      .. math::

         \alpha_0(s) =\int_{-\infty}^\infty \exp\left(-\frac{\sqrt{s^2+t^2}}{R}\right)\mbox{d}t

      obtained by reduction of the exponential function from 2D to 1D. The integral cannot be evaluated in closed form and is computed by OOFEM numerically. This function can be used in one-dimensional simulations of a two-dimensional specimen under uniaxial tension; for more details see [Gra14]_.

The above functions depend only on the distance :math:`s` between the interacting points and are not normalized. If the normalizing condition

.. math::

   \int_{V_\infty} \alpha(\mbf{x},\mbf{\xi})\;{\rm d}\mbf{\xi} = 1

is imposed in an infinite body :math:`V_{\infty}`, it is sufficient to scale :math:`\alpha_0` by a constant and set

.. math::

   \alpha(\mbf{x},\mbf{\xi})=\frac{\alpha_0(\Vert\mbf{x}-\mbf{\xi}\Vert)}{V_{r\infty}}

where

.. math::

   V_{r\infty} = \int_{V_\infty} \alpha_0(\Vert\mbf{\xi}\Vert)\;{\rm d}\mbf{\xi} 

Constant :math:`V_{r\infty}` can be computed analytically depending on the specific type of weight function and the number of spatial dimensions in which the analysis is performed. Since the factor :math:`1/V_{r\infty}` can be incorporated directly in the definition of :math:`\alpha_0`, this case is referred to as “no scaling”.

If the body of interest is finite (or even semi-infinite), the averaging integral can be performed only over the domain filled by the body, and the volume contributing to the nonlocal average at a point :math:`\mbf{x}` near the boundary is reduced as compared to points :math:`\mbf{x}` far from the boundary or in an infinite body. To make sure that the normalizing condition

.. math::

   \int_{V} \alpha(\mbf{x},\mbf{\xi})\;{\rm d}\mbf{\xi} = 1

holds for the specific domain :math:`V`, different approaches can be used. The standard approach defines the nonlocal weight function as

.. math::

   \alpha(\mbf{x},\mbf{\xi})=\frac{\alpha_0(\Vert\mbf{x}-\mbf{\xi}\Vert)}{V_r(\mbf{x})}

where

.. math::

   V_r(\mbf{x}) = \int_{V} \alpha_0(\Vert\mbf{x}-\mbf{\xi}\Vert)\;{\rm d}\mbf{\xi} 

According to the approach suggested by Borino, the weight function is defined as

.. math::

   \alpha(\mbf{x},\mbf{\xi})=\frac{\alpha_0(\Vert\mbf{x}-\mbf{\xi}\Vert)}{V_{r\infty}} + \left(1-\frac{V_r(\mbf{x})}{V_{r\infty}}\right)\delta(\mbf{x}-\mbf{\xi})

where :math:`\delta` is the Dirac distribution. One can also say that the nonlocal variable is evaluated as

.. math::

   \bar{\varepsilon}(\mbf{x}) = \frac{1}{V_{r\infty}}\int_V\alpha_0(\Vert\mbf{x}-\mbf{\xi}\Vert)\tilde\epsilon(\mbf{\xi})\;{\rm d}\mbf{\xi}+\left(1-\frac{V_r(\mbf{x})}{V_{r\infty}}\right)\tilde\varepsilon(\mbf{x})

The term on the right-hand side after the integral is a multiple of the local variable, and so it can be referred to as the local complement.  

In a recent paper [Gra14]_, special techniques that modify the averaging procedure based on the distance from a physical boundary of the domain or on the stress state have been considered. The details are explained in [Gra14]_. These techniques can be invoked by setting the optional parameter :it:`nlVariation` to 1, 2 or 3 and specifying additional parameters :math:`\beta` and :math:`\zeta` for distance-based averaging, or :math:`\beta` for stress-based averaging.  

The model parameters are summarized in :numref:`idnl_table` and :numref:`idnl_table_cont`.

.. table:: Nonlocal isotropic damage model for tensile failure -- summary.
   :name: idnl_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal isotropic damage model for concrete in tension                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`Idmnl1` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`E{rn}` :elemparam:`n{rn}`  |
   |                    | [:elemparam:`tAlpha{rn}`] [:elemparam:`equivstraintype{in}`] [:elemparam:`k{rn}`]            |
   |                    | [:elemparam:`damlaw{in}`] :elemparam:`e0{rn}` [:elemparam:`ef{rn}`] [:elemparam:`At{rn}`]    |
   |                    | [:elemparam:`Bt{rn}`] [:elemparam:`md{rn}`] :elemparam:`r{rn}` [:elemparam:`regionMap{ia}`]  |
   |                    | [:elemparam:`wft{in}`] [:elemparam:`averagingType{in}`] [:elemparam:`m{rn}`]                 |
   |                    | [:elemparam:`scalingType{in}`] [:elemparam:`averagedQuantity{in}`]                           |
   |                    | [:elemparam:`nlVariation{in}`] [:elemparam:`beta{rn}`] [:elemparam:`zeta{rn}`]               |
   |                    | [:elemparam:`maxOmega{rn}`]                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`E` Young's modulus                                                                 |
   |                    | - :param:`n` Poisson's ratio                                                                 |
   |                    | - :param:`tAlpha` thermal expansion coefficient                                              |
   |                    | - :param:`equivstraintype` allows to choose from different definitions of equivalent strain, |
   |                    |   same as for the local model; see Tab. :numref:`id_table`                                   |
   |                    | - :param:`k` ratio between uniaxial compressive and tensile strength, needed only if         |
   |                    |   equivstraintype=3, default value 1                                                         |
   |                    | - :param:`damlaw` allows to choose from different damage laws, same as for the local model;  |
   |                    |   see Tab. :numref:`id_table` (note that parameter *wf* cannot be used for the nonlocal      |
   |                    |   model)                                                                                     |
   |                    | - :param:`e0` strain at peak stress (for damage laws 0,1,2,3), limit elastic strain (for     |
   |                    |   damage law 4), characteristic strain (for damage law 5)                                    |
   |                    | - :param:`ef` strain parameter controling ductility, has the meaning of strain (for damage   |
   |                    |   laws 0 and 1), the tangent modulus just after the peak is                                  |
   |                    |   :math:`E_t=-f_t/(\varepsilon_f-\varepsilon_0)`                                             |
   |                    | - :param:`At` parameter of Mazars damage law, used only by law 4                             |
   |                    | - :param:`Bt` parameter of Mazars damage law, used only by law 4                             |
   |                    | - :param:`md` exponent, used only by damage law 5, default value 1                           |
   |                    | - :param:`r` nonlocal characteristic length :math:`R`; its meaning depends on the type of    |
   |                    |   weight function (e.g., interaction radius for the quartic spline)                          |
   |                    | - :param:`regionMap` map indicating the regions (currently region is characterized by cross  |
   |                    |   section number) to skip for nonlocal avaraging. The elements and corresponding IP are not  |
   |                    |   taken into account in nonlocal averaging process if corresponding regionMap value is       |
   |                    |   nonzero.                                                                                   |
   |                    | - :param:`wft` selects the type of nonlocal weight function: ; 1 - default, quartic spline   |
   |                    |   (bell-shaped function with bounded support) ; 2 - Gaussian function ; 3 - exponential      |
   |                    |   function (Green function in 1D) ; 4 - uniform averaging up to distance :math:`R` ; 5 -     |
   |                    |   uniform averaging over one finite element ; 6 - special function obtained by reducing the  |
   |                    |   2D exponential function to 1D (by numerical integration)                                   |
   |                    | - -- continued in Tab. :numref:`idnl_table_cont` ---                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `partests/barnl/barnl.oofem.in <https://github.com/oofem/oofem/blob/devel/tests/regression/p |
   |                    | artests/barnl/barnl.oofem.in>`_, `partests/barnl/barnl2.oofem.in <https://github.com/oofem/o |
   |                    | ofem/blob/devel/tests/regression/partests/barnl/barnl2.oofem.in>`_,                          |
   |                    | `partests/brazil_2d_nl2/brazil_2d_nl.oofem.in <https://github.com/oofem/oofem/blob/devel/tes |
   |                    | ts/regression/partests/brazil_2d_nl2/brazil_2d_nl.oofem.in>`_, `partests/lb03/lb03.oofem.in  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/partests/lb03/lb03.oofem.in>`_,  |
   |                    | `sm/distancebasedaveraging.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/distancebasedaveraging.in>`_, |
   |                    | `sm/stressbasedaveraging.in                                                                  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/stressbasedaveraging.in>`_    |
   +--------------------+----------------------------------------------------------------------------------------------+

.. table:: Nonlocal isotropic damage model for tensile failure -- continued.
   :name: idnl_table_cont

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlocal isotropic damage model for concrete in tension                                      |
   |                    |                                                                                              |
   |                    | - :param:`averagingType` activates a special averaging procedure, default value 0 does not   |
   |                    |   change anything, value 1 means averaging over one finite element (equivalent to *wft*\ =5, |
   |                    |   but kept here for compatibility with previous version)                                     |
   |                    | - :param:`m` multiplier for overnonlocal or undernonlocal formulation, which use *m*-times   |
   |                    |   the local variable plus :math:`(1-m)`-times the nonlocal variable, default value 1         |
   |                    | - :param:`scalingType` selects the type of scaling of the weight function (e.g. near a       |
   |                    |   boundary): ; 1 - default, standard scaling with integral of weight function in the         |
   |                    |   denominator ; 2 - no scaling (the weight function normalized in an infinite body is used   |
   |                    |   even near a boundary) ; 3 - Borino scaling (local complement)                              |
   |                    | - :param:`averagedQuantity` selects the variable to be averaged, default value 1 corresponds |
   |                    |   to equivalent strain, value 2 activates averaging of compliance variable                   |
   |                    | - :param:`nlVariation` activates a special averaging procedure, default value 0 does not     |
   |                    |   change anything, value 1 means distance-based averaging (the characteristic length is      |
   |                    |   linearly reduced near a physical boundary), value 2 means stress-based averaging (the      |
   |                    |   averaging is anisotropic and the characteristic length is affected by the stress), value 3 |
   |                    |   means distance-based averaging (the characteristic length is exponentially reduced near a  |
   |                    |   physical boundary)                                                                         |
   |                    | - :param:`beta` parameter :math:`\beta`, required only for distance-based and stress-based   |
   |                    |   averaging (i.e., for *nlVariation*\ =1, 2 or 3)                                            |
   |                    | - :param:`zeta` parameter :math:`\zeta`, required only for distance-based averaging (i.e.,   |
   |                    |   for *nlVariation*\ =1 or 3)                                                                |
   |                    | - :param:`maxOmega` maximum damage, used for convergence improvement (its value is between 0 |
   |                    |   and 0.999999 (default), and it affects only the secant stiffness but not the stress)       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dMat, PlaneStress, PlaneStrain, 1dMat                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           | Adaptivity support                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _subsubsection_anisotropic_damage_model_mdm:


Anisotropic damage model - Mdm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



Local formulation
^^^^^^^^^^^^^^^^^


.. _local-formulation:

Local formulation
^^^^^^^^^^^^^^^^^

The concept of isotropic damage is appropriate for materials weakened
by voids, but if the physical source of damage is the initiation and
propagation of microcracks, isotropic stiffness degradation
can be considered only as a first rough
approximation. More refined damage models take into account the highly
oriented nature of cracking, which is reflected by the anisotropic
character of the damaged stiffness or compliance matrices.

A number of anisotropic damage formulations have been proposed
in the literature. Here we use a model outlined by Jirásek
[mdm]_, which is based on the principle of energy
equivalence and on the construction of the inverse integrity
tensor by integration of a scalar over all spatial directions.
Since the model uses certain concepts from the
microplane theory, it is called the microplane-based damage model (MDM).

The general structure of the MDM
model is schematically shown in ``ff4``
and the basic equations are summarized in :numref:`tab2`.
Here, :math:`\veps` and :math:`\vsig` are the (nominal) second-order
strain and stress tensors
with components :math:`\eps_{ij}` and :math:`\sigma_{ij}`; :math:`ě` and :math:`š`
are first-order strain and stress tensors with components :math:`e_i`
and :math:`s_i`, which characterize the strain and stress on “microplanes”
of different orientations given by a unit vector :math:`\mbf{n}`
with components :math:`n_i`;
:math:`\psi` is a dimensionless compliance parameter
that is a scalar but can have different values for different
directions :math:`\mbf{n}`;
the symbol :math:`\delta` denotes a virtual quantity; and a sumperimposed
tilde denotes an effective quantity, which is supposed to characterize the
state of the intact material between defects such as microcracks or voids.

.. math::
   :label: eq1

   \vet = \vepst \cdot \mbf{n}

.. math::
   :label: eq2

   \vst = \psi \vs

.. math::
   :label: eq3

   \vs = \vsig \cdot \mbf{n}

.. math::
   :label: eq4

   \vsigt : \dvepst = \frac{3}{2\pi} \int_{\Omega} \vst \cdot \dvet \; \mbox{d}\Omega

.. math::
   :label: eq5

   \dvs \cdot \ve = d\vst \cdot \vet

.. math::
   :label: eq6

   \dvsig : \veps = \frac{3}{2\pi} \int_{\Omega} \dvs \cdot \ve \; d\Omega

.. math::
   :label: eq7

   \vsigt = \frac{3}{2\pi} \int_\Omega (\vst \otimes \mbf{n}) \sym \; d\Omega

.. math::
   :label: eq8

   \ve = \psi \vet

.. math::
   :label: eq9

   \veps = \frac{3}{2\pi} \int_\Omega (\ve \otimes \mbf{n}) \sym \; d\Omega

.. table:: Basic equations of microplane-based anisotropic damage model
   :name: tab2

   +--------------------------------------------------------------------------------------------+-----------------------------------------+--------------------------------------------------------------------------------+
   |                                                                                            |                                         |                                                                                |
   +--------------------------------------------------------------------------------------------+-----------------------------------------+--------------------------------------------------------------------------------+
   | :math:`\vet = \vepst \cdot \mbf{n}`                                                        | :math:`\vst = \psi š`                   | :math:`š = \vsig \cdot \mbf{n}`                                                |
   +--------------------------------------------------------------------------------------------+-----------------------------------------+--------------------------------------------------------------------------------+
   | :math:`\vsigt : \dvepst = \frac{3}{2\pi} \int_{\Omega} \vst \cdot \dvet \; \mbox{d}\Omega` | :math:`\dvs \cdot ě = d\vst \cdot \vet` | :math:`\dvsig : \veps = \frac{3}{2\pi} \int_{\Omega} \dvs \cdot ě \; d\Omega`  |
   +--------------------------------------------------------------------------------------------+-----------------------------------------+--------------------------------------------------------------------------------+
   | :math:`\vsigt = \frac{3}{2\pi} \int_\Omega (\vst \otimes \mbf{n}) \sym \; d\Omega`         | :math:`ě = \psi \vet`                   | :math:`\veps = \frac{3}{2\pi} \int_\Omega (ě \otimes \mbf{n}) \sym \; d\Omega` |
   +--------------------------------------------------------------------------------------------+-----------------------------------------+--------------------------------------------------------------------------------+

.. figure:: /figures/dm_comp.png
   :width: 80%
   :alt: Structure of microplane-based anisotropic damage model
   :align: center

   Structure of microplane-based anisotropic damage model

Combining the basic equations, it is possible to show that the components of the damaged material compliance tensor are given by

.. math::
   :label: damcom

   C_{ijkl}=M_{pqij}M_{rskl}C^e_{pqrs}

where :math:`C^e_{pqrs}` are the components of the elastic material compliance tensor,

.. math::
   :label: ee27

   M_{ijkl} = \frac{1}{4}\left(
   \psi_{ik}\delta_{jl}+\psi_{il}\delta_{jk}+\psi_{jk}\delta_{il}+\psi_{jl}\delta_{ik}\right)

are the components of the so-called damage effect tensor, and

.. math::
   :label: ee24

   \psi_{ij} = \frac{3}{2\pi}\int_\Omega \psi\, n_i n_j \, d\Omega

are the components of the second-order inverse integrity tensor. The integration domain :math:`\Omega` is the unit hemisphere. In practice, the integral over the unit hemisphere is evaluated by summing the contribution from a finite number of directions, according to one of the numerical integration schemes that are used by microplane models.

The scalar variable :math:`\psi` characterizes the relative compliance
in the direction given by the vector :math:`\mbf{n}`.
If :math:`\psi` is the same in all directions,
the inverse integrity tensor  evaluated from :eq:`ee24`
is equal to the unit second-order tensor (Kronecker delta) multiplied
by :math:`\psi`, the damage effect tensor evaluated from :eq:`ee27`
is equal to the symmetric fourth-order unit tensor multiplied
by :math:`\psi`,
and the damaged
material compliance tensor evaluated from :eq:`damcom` is the
elastic compliance tensor multiplied by :math:`\psi^2`. The factor multiplying
the elastic compliance tensor in the
isotropic damage model is :math:`1/(1-\omega)`, and so :math:`\psi` corresponds
to  :math:`1/\sqrt{1-\omega}`. In the initial undamaged state,
:math:`\psi=1` in all directions.  The evolution of :math:`\psi`
is governed by the history of the projected strain components.
In the simplest case, :math:`\psi` is driven by the normal strain
:math:`e_N=\eps_{ij}n_in_j`. Analogy with the isotropic damage model
leads to the damage law

.. math::

   \psi=f(\kappa)

and loading-unloading conditions

.. math::

   g(e_N,\kappa)\equiv e_N-\kappa\le 0, \hskip 10mm
   \dot{\kappa}\ge 0, \hskip 10mm
   \dot{\kappa}g(e_N,\kappa)=0

in which :math:`\kappa` is a history variable that represents the maximum
level of normal strain in the given direction ever reached in the
previous history of the material. An appropriate modification
of the exponential softening
law leads to the damage law

.. math::
   :label: expsoft2

   f(\kappa)=\left\{
   \begin{array}{ll}
   1 & \mbox{ if } \kappa\le e_0
   \\ 
   \sqrt{\frac{\kappa}{e_0}\exp\left(\frac{\kappa-e_0}{e_f-e_0}\right)}
   & \mbox{ if } \kappa>e_0
   \end{array}
   \right.

where :math:`e_0` is a parameter controlling the elastic limit, and :math:`e_f>e_0`
is another parameter controlling ductility.
Note that softening in a limited number of directions does not necessarily
lead to softening on the macroscopic level, because the response
in the other directions remains elastic. Therefore, :math:`e_0` corresponds
to the elastic limit but not to the state at peak stress.

If the MDM model is used in its basic form described above, the compressive strength turns out to depend on the Poisson ratio and, in applications to concrete, its value is too low compared to the tensile strength. The model is designed primarily for tensile-dominated failure, so the low compressive strength is not considered as a major drawback. Still, it is desirable to introduce a modification that would prevent spurious compressive failure in problems where moderate compressive stresses appear. The desired effect is achieved by redefining the projected strain :math:`e_N` as

.. math::
   :label: ee37

   e_N =  \frac{\eps_{ij}n_in_j}{1-\displaystyle\frac{m}{Ee_0}\sigma_{kk}}

where :math:`m` is a nonnegative parameter that controls the sensitivity to the mean stress, :math:`\sigma_{kk}` is the trace of the stress tensor, and the normalizing factor :math:`Ee_0` is introduced in order to render the parameter :math:`m` dimensionless. Under compressive stress states (characterized by :math:`\sigma_{kk}<0`), the denominator in (:eq:`ee37`) is larger than 1, and the projected strain is reduced, which also leads to a reduction of damage. A typical recommended value of parameter :math:`m` is 0.05.

.. _nonlocal-formulation:

Nonlocal formulation
""""""""""""""""""""

Nonlocal formulation of the MDM model is based on the averaging of the inverse integrity tensor. This roughly corresponds to the nonlocal isotropic damage model with averaging of the compliance variable :math:`\gamma=\omega/(1-\omega)`, which does not cause any spurious locking effects. In equation (:eq:`ee27`) for the evaluation of the damage effect tensor, the inverse integrity tensor is replaced by its weighted average with components

.. math::
   :label: psinl1

   \bar{\psi}_{ij}(\vx)=\int_V \alpha(\vx,\vxi)\psi_{ij}(\vxi)\mbox{d}\vxi

By fitting a wide range of numerical results, it has been found that
the parameters of the nonlocal MDM model can be estimated from the
measurable material properties using the formulas

.. math::

   \lambda_f &= \frac{EG_f}{Rf_t^2}
   \\
   \lambda &= \frac{\lambda_f}{1.47 - 0.0014\lambda_f}
   \\
   e_0 &= \frac{f_t}{(1-m)E(1.56 + 0.006\lambda)}
   \\
   e_f &= e_0[1 + (1-m)\lambda]

where :math:`E` is Young's modulus, :math:`G_f` is the fracture energy, :math:`f_t`
is the uniaxial tensile strength,
:math:`m` is the compressive correction factor, typically chosen
as :math:`m=0.05`, and :math:`R` is the radius of nonlocal interaction reflecting the
internal length of the material.

\paragraph{Input Record}
The model description and parameters are summarized
in :numref:`mdm_table`.

.. table:: MDM model - summary
   :name: mdm_table

   +-----------------------+----------------------------------------------------------------------------------------------+
   | Description           | MDM Anisotropic damage model                                                                 |
   |                       | Common parameters                                                                            |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Record Format         | :descitem:`Mdm` :elemparam:`d{rn}` :elemparam:`nmp{ins}` :elemparam:`talpha{rn}`             |
   |                       | :elemparam:`parmd{rn}` :elemparam:`nonloc{in}` :elemparam:`formulation{in}`                  |
   |                       | :elemparam:`mode{in}`                                                                        |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Parameters            | - :param:`num` material model number                                                         |
   |                       | - :param:`D` material density                                                                |
   |                       | - :param:`nmp` number of microplanes used for hemisphere integration, supported values are   |
   |                       |   21,28, and 61                                                                              |
   |                       | - :param:`talpha` thermal dillatation coeff                                                  |
   |                       | - :param:`parmd`                                                                             |
   |                       | - :param:`nonloc`                                                                            |
   |                       | - :param:`formulation`                                                                       |
   |                       | - :param:`mode`                                                                              |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Nonlocal variant I    |                                                                                              |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Additional params     | :elemparam:`r{rn}` :elemparam:`efp{rn}` :elemparam:`ep{rn}`                                  |
   |                       |                                                                                              |
   |                       | - :param:`r` nonlocal interaction radius                                                     |
   |                       | - :param:`efp` :math:`\varepsilon_fp` is a model parameter that controls the post-peak slope |
   |                       |   :math:`\varepsilon_fp` =:math:`\varepsilon_f-\varepsilon_0`, where :math:`\varepsilon_f`   |
   |                       |   is strain at zero stress level.                                                            |
   |                       | - :param:`ep` max effective strain at peak :math:`\varepsilon_0`                             |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Nonlocal variant II   |                                                                                              |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Additional params     | :elemparam:`r{rn}` :elemparam:`gf{rn}` :elemparam:`ft{rn}`                                   |
   |                       |                                                                                              |
   |                       | - :param:`r` nonlocal intraction radius                                                      |
   |                       | - :param:`gf` fracture energy                                                                |
   |                       | - :param:`ft` tensile strength                                                               |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Local variant I       |                                                                                              |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Additional params     | :elemparam:`efp{rn}` :elemparam:`ep{rn}`                                                     |
   |                       |                                                                                              |
   |                       | - :param:`efp` :math:`\varepsilon_fp` is a model parameter that controls the post-peak slope |
   |                       |   :math:`\varepsilon_fp` =:math:`\varepsilon_f-\varepsilon_0`, where :math:`\varepsilon_f`   |
   |                       |   is strain at zero stress level.                                                            |
   |                       | - :param:`ep` max effective strain at peak :math:`\varepsilon_0`                             |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Local variant II      |                                                                                              |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Additional params     | :elemparam:`gf{rn}` :elemparam:`ep{rn}`                                                      |
   |                       |                                                                                              |
   |                       | - :param:`gf` fracture energy                                                                |
   |                       | - :param:`ep` max effective strain at peak :math:`\varepsilon_0`                             |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Supported modes       | 3dMat, PlaneStress                                                                           |
   +-----------------------+----------------------------------------------------------------------------------------------+
   | Features              | Adaptivity support                                                                           |
   +-----------------------+----------------------------------------------------------------------------------------------+

.. _sec:idmfi:

Isotropic damage model for interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model provides an interface law which can be used to describe a damageable interface between two materials (e.g.\ between steel reinforcement and concrete). The law is formulated in terms of the traction vector and the displacement jump vector. The basic response is elastic, with stiffness :param:`kn` in the normal direction and :param:`ks` in the tangential direction.

Similar to other isotropic damage models, this model assumes that the stiffness degradation is isotropic, i.e., both stiffness moduli decrease proportionally and independently of the loading direction. The damaged stiffnesses are :math:`:param:`\ kn`\times(1-\omega)` and :math:`:param:`\ ks`\times(1-\omega)` where :math:`\omega` is a scalar damage variable.

The damage evolution law is postulated in an explicit form, relating the damage variable :math:`\omega` to the largest previously reached equivalent “strain” level, :math:`\kappa`.

The equivalent “strain”, :math:`\tilde\varepsilon`, is a scalar measure of the displacement jump vector. The choice of the specific expression for the equivalent strain affects the shape of the elastic domain in the strain space and plays a similar role to the choice of a yield condition in plasticity.

Currently, in the present implementation, the equivalent strain is given by

.. math::

   \tilde\varepsilon = \sqrt{\langle w_n\rangle^2 + \beta w_s^2}

where :math:`\langle w_n\rangle` is the positive part of the normal displacement jump (opening of the interface) and :math:`w_s` is the norm of the tangential part of displacement jump (sliding of the interface). Parameter :math:`\beta` is optional and its default value is 0, in which case damage depends on the opening only (not on the sliding).

The dependence of damage :math:`\omega` on maximum equivalent strain :math:`\kappa` is described by the following damage law which corresponds to exponential softening:

.. math::

   \omega = \left\{ \begin{array}{ll}
   0 & \mbox{ for } \kappa\le \varepsilon_0 \\
   1 - \displaystyle\frac{\varepsilon_0}{\kappa}  \exp\left( - \displaystyle\frac{f_t( \kappa - \varepsilon_0 )}{G_f} \right) & \mbox{ for } \kappa> \varepsilon_0
   \end{array}\right.

Here, :math:`\varepsilon_0=f_t/k_n` is the value of equivalent strain at the onset of damage.

Note that if the interface is subjected to shear traction only (with zero or negative normal traction), the propagation of damage starts when the magnitude of the sliding displacement is :math:`\vert w_s\vert=\varepsilon_0/\sqrt{\beta}`, i.e., when the magnitude of the shear traction is equal to

.. math::

   f_s = \frac{k_s \varepsilon_0}{\sqrt{\beta}} = f_t\frac{k_s }{k_n\sqrt{\beta}}

So the ratio between the shear strength and tensile strength of the interface, :math:`f_s/f_t`, is equal to :math:`k_s/k_n\sqrt{\beta}`.

The model parameters are summarized in :numref:`iid_table`. 

.. table:: Isotropic damage model for interface elements -- summary.
   :name: iid_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic damage model for concrete in tension                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`isointrfdm01` :elemparam:`kn{rn}` :elemparam:`ks{rn}` :elemparam:`ft{rn}`         |
   |                    | :elemparam:`gf{rn}` :optelemparam:`maxomega{rn}` :elemparam:`talpha{rn}` :elemparam:`d{rn}`  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`d` material density                                                                |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`kn` elastic stifness in normal direction                                           |
   |                    | - :param:`ks` elastic stifness in tangential direction                                       |
   |                    | - :param:`ft` tensile strength                                                               |
   |                    | - :param:`gf` fracture energy                                                                |
   |                    | - :optparam:`maxomega` maximum damage, used for convergence improvement (its value is        |
   |                    |   between 0 and 0.999999 (default), and it affects only the secant stiffness but not the     |
   |                    |   stress)                                                                                    |
   |                    | - :optparam:`beta` parameter controlling the effect of sliding part of displacement jump on  |
   |                    |   equivalent strain, default value 0                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dInterface, 3dInterface                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           |                                                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/interface3d.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/interface3d.in>`_             |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _sec:idmfiTabulated:

Isotropic damage model for interfaces using tabulated data for damage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The model provides an interface law which can be used to describe a damageable interface between two materials (e.g., between steel reinforcement and concrete). The law is formulated in terms of the traction vector and the displacement jump vector. The basic response is elastic, with stiffness :param:`kn` in the normal direction and :param:`ks` in the tangential direction.

Similar to other isotropic damage models, this model assumes that the stiffness degradation is isotropic, i.e., both stiffness moduli decrease proportionally and independently of the loading direction. The damaged stiffnesses are :param:`kn`\ :math:`\times(1-\omega)` and :param:`ks`\ :math:`\times(1-\omega)` where :math:`\omega` is a scalar damage variable.

The equivalent “strain”, :math:`\tilde\varepsilon`, is a scalar measure derived from the displacement jump vector. The choice of the specific expression for the equivalent strain affects the shape of the elastic domain in the strain space and plays a similar role to the choice of a yield condition in plasticity. Currently, in the present implementation, :math:`\tilde\varepsilon` is equal to the positive part of the normal displacement jump (opening of the interface).

The damage evolution law is postulated in a separate file that should have the following format. Each line should contain one strain, damage pair separated by a whitespace character. The exception to this is the first line which should contain a single integer stating how many strain, damage pairs that the file will contain. The strains given in the file is defined as the equivalent strain minus the limit of elastic deformation.
To find the damage for arbitrary strains linear interpolation between the tabulated values is used. If a strain larger than one in the given table is achieved the respective damage for the largest tabulated strain will be used. Both the strains and damages must be given in a strictly increasing order.

The model parameters are summarized in :numref:`iidTabulated_table`. 

.. table:: Isotropic damage model for interface elements using tabulated data for damage -- summary.
   :name: iidTabulated_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Isotropic damage model for concrete in tension                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`isointrfdm02` :elemparam:`kn{rn}` :elemparam:`ks{rn}` :elemparam:`ft{rn}`         |
   |                    | :elemparam:`tablename{rn}` :optelemparam:`maxomega{rn}` :elemparam:`talpha{rn}`              |
   |                    | :elemparam:`d{rn}`                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`d` material density                                                                |
   |                    | - :param:`tAlpha` thermal dilatation coefficient                                             |
   |                    | - :param:`kn` elastic stifness in normal direction                                           |
   |                    | - :param:`ks` elastic stifness in tangential direction                                       |
   |                    | - :param:`ft` tensile strength                                                               |
   |                    | - :param:`tablename` file name of the table with the strain damage pairs                     |
   |                    | - :param:`maxomega` maximum damage, used for convergence improvement (its value is between 0 |
   |                    |   and 0.999999 (default), and it affects only the secant stiffness but not the stress)       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dInterface, 3dInterface                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Features           |                                                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+

.. math::
   :label: eq:tabular_example

   \begin{array}{|c|c|}
   \hline
   \text{Strain} & \text{Damage} \\
   \hline
   0 & 0 \\
   0.1 & 0.01 \\
   0.2 & 0.04 \\
   0.3 & 0.09 \\
   \hline
   \end{array}

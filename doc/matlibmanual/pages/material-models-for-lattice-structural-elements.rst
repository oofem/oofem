Material models for lattice structural elements
===============================================


.. _lattice_material_models:

Material models for lattice structural elements
-----------------------------------------------

.. _lattice_linear_elastic_model:

Linear elastic lattice model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a linear elastic material used together with lattice elements.
It uses a stress strain law of the form

.. math::

   \boldsymbol{\sigma} = \mathbf{D}_{\rm e} \boldsymbol{\varepsilon}

where :math:`\boldsymbol{\sigma}` is a vector of tractions and rotational components, and :math:`\boldsymbol{\varepsilon}` is a vector of strains obtained from displacement jumps smeared over the element length and rotational components.
Furthermore, :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix which is based on the elastic modulus of the lattice material :math:`E`, and a parameter :math:`a_1` which is the ratio of the modulus of the shear and normal direction. There is the option to include a rotational stiffness. The model parameters are summarised in :numref:`latticelinearelastic_table`.

.. table:: Linear elastic material model for lattice elements -- summary.
   :name: latticelinearelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Linear elastic model for lattice elements                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticelinearelastic` :elemparam:`in` :elemparam:`d{rn}`                          |
   |                    | :optelemparam:`talpha{rn}` :optelemparam:`calpha{rn}` :elemparam:`e{rn}`                     |
   |                    | :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal strain expansion coefficient. Default is 0.                     |
   |                    | - :optparam:`calpha` Thermal displacement expansion coefficient. Default is 0.               |
   |                    | - :param:`e` Young's modulus of the equivalent lattice material                              |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dlattice, 3dlattice                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice3delastic.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3delastic.in>`_,       |
   |                    | `lm/lattice3dshell_combined_n1_lle.in <https://github.com/oofem/oofem/blob/devel/tests/regre |
   |                    | ssion/lm/lattice3dshell_combined_n1_lle.in>`_, `lm/lattice3dshell_combined_n4_lle.in <https: |
   |                    | //github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell_combined_n4_lle.in>`_ |
   |                    | , `lm/lattice3dshell_lle.in                                                                  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell_lle.in>`_,     |
   |                    | `lm/lattice3dshell_lle_twist_n4.in <https://github.com/oofem/oofem/blob/devel/tests/regressi |
   |                    | on/lm/lattice3dshell_lle_twist_n4.in>`_, `lm/latticedyn1.in                                  |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticedyn1.in>`_             |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _scalar-damage-lattice-model:

Scalar damage lattice model
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a scalar damage material model used together with lattice elements.
It uses a scalar damage model, which results in a stress-strain law of the form

.. math::

   \boldsymbol{\sigma} = \left(1-\omega\right) \mathbf{D}_{\rm e} \boldsymbol{\varepsilon}

where :math:`\boldsymbol{\sigma}` is a vector of tractions and rotational components, and :math:`\boldsymbol{\varepsilon}` is a vector of strains obtained from displacement jumps smeared over the element length and rotational components.
Furthermore, :math:`\omega` is the damage variable varying from 0 (undamaged) to 1 (fully damaged). 
Also, :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix which is based on the elastic modulus of the lattice material :math:`E`, and a parameter :math:`a_1` which is the ratio of the modulus of the shear and normal direction.
The strength envelope (onset of damage) is elliptic and determined by three parameters, :math:`f_{\rm t}`, :math:`f_{\rm q}`, and :math:`f_{\rm c}`. The evolution of the damage variable :math:`\omega` is controlled by normal stress-normal crack opening law. The three possible laws are linear, bilinear, and exponential.

This simple damage material model for lattice elements has been used in many recent articles. For instance, the two-dimensional version has been used in *P. Grassl and M. Jirásek. "Meso-scale approach to modelling the fracture process zone of concrete subjected to uniaxial tension". International Journal of Solids and Structures. Volume 47, Issues 7-8, pp. 957-968, 2010*. An example of an application of the three-dimensional version of the model is found in *\ P. Grassl, J. Bolander. "Three-Dimensional Network Model for Coupling of Fracture and Mass Transport in Quasi-Brittle Geomaterials", Materials, 9, 782, 2016*.

The model parameters are summarised in :numref:`latticedamage_table`.

.. table:: Scalar damage model for lattice elements -- summary.
   :name: latticedamage_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Scalar damage model for lattice elements                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticedamage` :elemparam:`in` :elemparam:`d{rn}` :optelemparam:`talpha{rn}`      |
   |                    | :elemparam:`e{rn}` :optelemparam:`a1{rn}` :optelemparam:`a2{rn}` :optelemparam:`e0{rn}`      |
   |                    | :optelemparam:`coh{rn}` :optelemparam:`ec{rn}` :optelemparam:`stype{rn}`                     |
   |                    | :optelemparam:`wf{rn}` :optelemparam:`wf1{rn}`                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal strain expansion coefficient. Default i 0.                      |
   |                    | - :optparam:`calpha` thermal displacement expansion coefficient. Default is 0.               |
   |                    | - :param:`e` normal modulus of lattice material                                              |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   |                    | - :param:`e0` strain at tensile strength: :math:`f_{\rm t}/E`                                |
   |                    | - :param:`coh` ratio of shear and tensile strength: :math:`f_{\rm q}/f_{\rm t}`              |
   |                    | - :param:`ec` ratio of compressive and tensile strength: :math:`f_{\rm c}/f_{\rm t}`         |
   |                    | - :optparam:`stype` softening types: 1-linear, 2-bilinear and 3-exponential. Default is 1.   |
   |                    | - :param:`wf` displacement threshold related to fracture energy used in all three softening  |
   |                    |   types.                                                                                     |
   |                    | - :optparam:`wf1` displacement threshold related to softening type 2. Default is wf1=0.15    |
   |                    |   wf.                                                                                        |
   |                    | - :optparam:`e01` strain threshold related to softening type 2. Default is wf1=0.15 wf       |
   |                    | - :optparam:`bio` Biot's coefficient.                                                        |
   |                    | - :optparam:`btype` Type to consider how Biot's coefficient changes with crack opening.      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dlattice, 3dlattice                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice2dboundary1.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2dboundary1.in>`_,     |
   |                    | `lm/lattice2drandom.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2drandom.in>`_,        |
   |                    | `lm/lattice3d1.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d1.in>`_,             |
   |                    | `lm/lattice3d2.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d2.in>`_,             |
   |                    | `lm/lattice3d3.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d3.in>`_,             |
   |                    | `lm/lattice3d4.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d4.in>`_,             |
   |                    | `lm/lattice3dboundarytruss.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dboundarytruss.in>`_, |
   |                    | `lm/lattice3dnlshell_dmg_bend_n7.in <https://github.com/oofem/oofem/blob/devel/tests/regress |
   |                    | ion/lm/lattice3dnlshell_dmg_bend_n7.in>`_, `lm/lattice3drandom.in                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3drandom.in>`_,        |
   |                    | `lm/lattice3dshell_dmg_bend_n7.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                    | n/lm/lattice3dshell_dmg_bend_n7.in>`_                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _plasticity_damage_lattice_model:

Plasticity damage lattice model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This combined plasticity damage model for lattice elements has been developed to describe the failure process of geomaterials such as concrete and rock. The plasticity part is based on the effective stress and damage is evaluated by the plastic strain. The model was first introduced in [GraDav11]_ without hardening and then extended to hardening in [AthWheGra18]_.

The main stress-strain law has the form

.. math::

   \boldsymbol{\sigma} = \left(1-\omega\right) \mathbf{D}_{\rm e} \left(\boldsymbol{\varepsilon}-\boldsymbol{\varepsilon}_{\rm p} \right)

where :math:`\boldsymbol{\sigma}` is a vector of tractions and rotational components, :math:`\omega` is the damage variable and :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix, :math:`\boldsymbol{\varepsilon}` is a vector of strains obtained from displacement jumps smeared over the element length and rotational components and :math:`\boldsymbol{\varepsilon}_{\rm p}` are the plastic strains.

The plasticity part is based on the effective stresses and uses only a subset of the strains and stresses which are the normal stress :math:`\sigma_{\rm n}` and the two shear stresses :math:`\sigma_{\rm s}` and :math:`\sigma_{\rm t}`. The yield surface is composed of two ellipse which are arranged so that the transition of the two ellipses is smooth (Figure ``plastdamyieldfig``).

.. figure:: /figures/plastdamyieldfig.svg
   :width: 60%
   :alt: Latticeplastdam yield surface, which is composed of ellipse with a smooth transition.

   Latticeplastdam yield surface, which is composed of ellipse with a smooth transition.

The resulting yield function :math:`f` is

.. math::

   f = \left\{ \begin{array}{l} \alpha^2\bar{\sigma}_{\rm n}^2 + 2 \dfrac{\alpha^2\left(f_{\rm c} - \alpha \beta f_{\rm t}\right)}{\left(1+\alpha \beta\right)} q \bar{\sigma}_{\rm n} + \bar{\sigma}_{\rm q}^2 - \dfrac{2 \alpha^2 f_{\rm c}f_{\rm t} + \alpha^2 \left(1-\alpha \beta \right) f_{\rm t}^2}{1+\alpha \beta} q^2 \\ \mbox{if $\bar{\sigma}_{\rm n} \geq -\dfrac{f_{\rm c} - \alpha \beta f_{\rm t}}{1+\alpha \beta} q$} \vspace{0.5cm}\\
     \dfrac{\bar{\sigma}_{\rm n}^2}{\beta^2} + 2 \dfrac{f_{\rm c} - \alpha \beta f_{\rm t}}{\beta^2 \left(1+\alpha \beta \right)} q \bar{\sigma}_{\rm n} + \bar{\sigma}_{\rm q}^2 + \dfrac{\left(1-\alpha \beta\right) f_{\rm c}^2 -2 \alpha \beta f_{\rm c} f_{\rm t} }{\beta^2 \left(1+\alpha \beta \right)} q^2\\ \mbox{if $\bar{\sigma}_{\rm n} < - \dfrac{f_{\rm c} - \alpha \beta f_{\rm t}}{1+\alpha \beta} q$} \end{array} \right.

Here parameters :math:`\alpha` and :math:`\beta` are the slopes, and :math:`f_{\rm c}` and :math:`f_{\rm t}` are the strength parameters shown in Figure ``plastdamyieldfig``. The hardening function :math:`q` controls the shape of the yield function.

The model parameters are summarised in :numref:`latticeplastdam_table`.

.. table:: Combined plasticity damage model for lattice elements
   :name: latticeplastdam_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Combined plasticity damage model for lattice elements                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeplastdam` :elemparam:`in` :elemparam:`d{rn}` :optelemparam:`talpha{rn}`    |
   |                    | :optelemparam:`calpha{rn}` :elemparam:`e{rn}` :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`  |
   |                    | :elemparam:`ft{rn}` :elemparam:`fc{rn}` :optelemparam:`angle1{rn}`                           |
   |                    | :optelemparam:`angle2{rn}` :optelemparam:`flow{rn}` :optelemparam:`stype{rn}`                |
   |                    | :elemparam:`wf{rn}` :optelemparam:`ft1{rn}` :optelemparam:`wf1{rn}`                          |
   |                    | :optelemparam:`ahard{rn}` :optelemparam:`damage{in}` :optelemparam:`sub{in}`                 |
   |                    | :optelemparam:`tol{rn}` :optelemparam:`iter{in}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal expansion coefficient. Default is 0.                            |
   |                    | - :optparam:`calpha` thermal displacement expansion coefficient. Optional parameter. Default |
   |                    |   is 0.                                                                                      |
   |                    | - :param:`e` modulus of the equivalent lattice material                                      |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   |                    | - :param:`ft` tensile strength                                                               |
   |                    | - :param:`fc` compressive strength                                                           |
   |                    | - :optparam:`angle1` ratio of compressive and tensile strength. Default is 0.5               |
   |                    | - :optparam:`angle2` ratio of compressive and tensile strength. Default is 0.5               |
   |                    | - :optparam:`flow` ratio of compressive and tensile strength. Default is 0.25                |
   |                    | - :optparam:`stype` softening types: 0-exponential, 1-bilinear. Default is 0                 |
   |                    | - :param:`wf` displacement threshold related to fracture energy used in both softening       |
   |                    |   types.                                                                                     |
   |                    | - :optparam:`wf1` displacement threshold related to softening type 1. Default is wf1=0.1 wf. |
   |                    | - :optparam:`ft1` stress threshold related to softening type 1. Optional parameter. Default  |
   |                    |   is ft1=0.15 ft.                                                                            |
   |                    | - :optparam:`ahard` displacement threshold related to fracture energy used in all three      |
   |                    |   softening types. Default value is 0.1                                                      |
   |                    | - :optparam:`damage` flag to switch on and off damage. Default value is 1, which considers   |
   |                    |   damage.                                                                                    |
   |                    | - :optparam:`sub` maximum number of subincrementations in the plasticity part. Default is    |
   |                    |   10.                                                                                        |
   |                    | - :optparam:`tol` tolerance for the newton iteration of the plasticity part. Default value   |
   |                    |   is 1.e-6.                                                                                  |
   |                    | - :optparam:`iter` maximum number of iterations for the stress return of the plasticity      |
   |                    |   model. Default is 100.                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice3ddamplast1.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3ddamplast1.in>`_,     |
   |                    | `lm/lattice3ddamplast2.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3ddamplast2.in>`_,     |
   |                    | `lm/lattice3ddamplast3.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3ddamplast3.in>`_,     |
   |                    | `lm/latticeplastdam3drandom.in                                                               |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticeplastdam3drandom.in>`_ |
   +--------------------+----------------------------------------------------------------------------------------------+



The theory of the model is described in the paper *P. Grassl and T. Davies. "Lattice modelling of corrosion induced cracking and bond in reinforced concrete". Cement and Concrete Composites, Volume 33, pp. 918-924, 2011*.

The model parameters are summarised in :numref:`latticebond_table`.

.. table:: Bond plasticity model for lattice elements -- summary.
   :name: latticebond_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Bond model for lattice elements                                                              |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticebondplast` :elemparam:`in` :elemparam:`d{rn}` :optelemparam:`talpha{rn}`   |
   |                    | :optelemparam:`calpha{rn}` :elemparam:`e{rn}` :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`  |
   |                    | :elemparam:`fc{rn}` :elemparam:`angle1{rn}` :elemparam:`ef{rn}` :elemparam:`sub{in}`         |
   |                    | :elemparam:`iter{in}` :elemparam:`tol{rn}`                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal expansion coefficient. Default is 0.                            |
   |                    | - :optparam:`calpha` thermal displacement expansion coefficient. Optional parameter. Default |
   |                    |   is 0.                                                                                      |
   |                    | - :param:`e` modulus of the equivalent lattice material                                      |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   |                    | - :param:`fc` compressive strength                                                           |
   |                    | - :param:`angle1` friction angle                                                             |
   |                    | - :optparam:`ef` strain threshold to control hardening                                       |
   |                    | - :optparam:`sub` maximum number of subincrementations                                       |
   |                    | - :optparam:`iter` maximum number of newton iterations                                       |
   |                    | - :optparam:`tol` tolerance for newton method                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice3dbondplast1.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dbondplast1.in>`_,    |
   |                    | `lm/lattice3dbondplast2.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dbondplast2.in>`_     |
   +--------------------+----------------------------------------------------------------------------------------------+

\subsubsection{Viscoelastic lattice model}

This model combines a viscoelastic model, described in previous sections, with a linear elastic lattice material. Two material entries are used. One is for the viscoelastic lattice model and the other is for the viscoelastic model of choice. 

The model description and parameters are summarized in :numref:`latticeviscoelastic_table`.

.. table:: Visco elastic material model for lattice elements -- summary.
   :name: latticeviscoelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Viscoelastic lattice model                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeviscoelastic` :elemparam:`in` :elemparam:`viscomat{in}` :elemparam:`d{rn}` |
   |                    | :optelemparam:`talpha{rn}` :optelemparam:`calpha{rn}` :elemparam:`e{rn}`                     |
   |                    | :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`                                                |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`viscomat` material number for viscoelastic model                                   |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal strain expansion coefficient. Default is 0.                     |
   |                    | - :optparam:`calpha` Thermal displacement expansion coefficient. Default is 0.               |
   |                    | - :param:`e` Young's modulus of the equivalent lattice material                              |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dlattice, 3dlattice                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/latticeviscoelastic_mps_1.in <https://github.com/oofem/oofem/blob/devel/tests/regression |
   |                    | /lm/latticeviscoelastic_mps_1.in>`_                                                          |
   +--------------------+----------------------------------------------------------------------------------------------+


.. table:: Model description and parameters for viscoelastic damage lattice model.

   +-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Description | Combines a viscoelastic model with a damage lattice material model. Two material entries are used. One is for the viscoelastic extension of the damage lattice model and the other is for the viscoelastic model of choice. |
   +-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table:: Damage visco elastic material model for lattice elements -- summary.
   :name: latticedamageviscoelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Viscoelastic damage lattice model                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticedamageviscoelastic` :elemparam:`in` :elemparam:`viscomat{in}`              |
   |                    | :elemparam:`timefactor{rn}` :elemparam:`d{rn}` :optelemparam:`talpha{rn}`                    |
   |                    | :optelemparam:`calpha{rn}` :elemparam:`e{rn}` :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`  |
   |                    | :optelemparam:`e0{rn}` :optelemparam:`coh{rn}` :optelemparam:`ec{rn}`                        |
   |                    | :optelemparam:`stype{rn}` :optelemparam:`wf{rn}` :optelemparam:`wf1{rn}`                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`viscomat` material number for viscoelastic model                                   |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal strain expansion coefficient. Default is 0.                     |
   |                    | - :optparam:`calpha` Thermal displacement expansion coefficient. Default is 0.               |
   |                    | - :param:`e` Young's modulus of the equivalent lattice material                              |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   |                    | - :param:`e0` strain at tensile strength: :math:`f_{\rm t}/E`                                |
   |                    | - :param:`coh` ratio of shear and tensile strength: :math:`f_{\rm q}/f_{\rm t}`              |
   |                    | - :param:`ec` ratio of compressive and tensile strength: :math:`f_{\rm c}/f_{\rm t}`         |
   |                    | - :optparam:`stype` softening types: 1-linear, 2-bilinear and 3-exponential. Default is 1.   |
   |                    | - :param:`wf` displacement threshold related to fracture energy used in all three softening  |
   |                    |   types.                                                                                     |
   |                    | - :optparam:`wf1` displacement threshold related to softening type 2. Default is wf1=0.15    |
   |                    |   wf.                                                                                        |
   |                    | - :optparam:`e01` strain threshold related to softening type 2. Default is wf1=0.15 wf       |
   |                    | - :optparam:`bio` Biot's coefficient.                                                        |
   |                    | - :optparam:`btype` Type to consider how Biot's coefficient changes with crack opening.      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dlattice, 3dlattice                                                                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/latticedamagevisco_mps_1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/ |
   |                    | lm/latticedamagevisco_mps_1.in>`_                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+


Viscoelastic plasticity damage lattice model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This model combines a viscoelastic model, described in previous sections, with a plasticity damage lattice material model. Two material entries are used. One is for the viscoelastic extension of the plasticity damage lattice model and the other is for the viscoelastic model of choice. 

The model description and parameters are summarized in :numref:`latticeplasticdamageviscoelastic_table`.

.. table:: Plasticity-damage viscoelastic material model for lattice elements -- summary
   :name: latticeplasticdamageviscoelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Plasticity-damage viscoelastic lattice model                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeplasticitydamageviscoelastic` :elemparam:`in` :elemparam:`viscomat{in}`    |
   |                    | :elemparam:`timefactor{rn}` :elemparam:`d{rn}` :optelemparam:`talpha{rn}`                    |
   |                    | :optelemparam:`calpha{rn}` :elemparam:`e{rn}` :optelemparam:`a1{rn}` :optelemparam:`a2{rn}`  |
   |                    | :elemparam:`ft{rn}` :elemparam:`fc{rn}` :optelemparam:`angle1{rn}`                           |
   |                    | :optelemparam:`angle2{rn}` :optelemparam:`flow{rn}` :optelemparam:`stype{rn}`                |
   |                    | :elemparam:`wf{rn}` :optelemparam:`ft1{rn}` :optelemparam:`wf1{rn}`                          |
   |                    | :optelemparam:`ahard{rn}` :optelemparam:`damage{in}` :optelemparam:`sub{in}`                 |
   |                    | :optelemparam:`tol{rn}` :optelemparam:`iter{in}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`viscomat` material number for viscoelastic material                                |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal strain expansion coefficient. Default is 0.                     |
   |                    | - :optparam:`calpha` Thermal displacement expansion coefficient. Default is 0.               |
   |                    | - :param:`e` Young's modulus of the equivalent lattice material                              |
   |                    | - :optparam:`a1` ratio of shear and normal modulus. Default is 1.                            |
   |                    | - :optparam:`a2` ratio of rotational and normal modulus. Default is 1.                       |
   |                    | - :param:`ft` tensile strength                                                               |
   |                    | - :param:`fc` compressive strength                                                           |
   |                    | - :optparam:`angle1` ratio of compressive and tensile strength. Default is 0.5               |
   |                    | - :optparam:`angle2` ratio of compressive and tensile strength. Default is 0.5               |
   |                    | - :optparam:`flow` ratio of compressive and tensile strength. Default is 0.25                |
   |                    | - :optparam:`stype` softening types: 0-exponential, 1-bilinear. Default is 0                 |
   |                    | - :param:`wf` displacement threshold related to fracture energy used in both softening       |
   |                    |   types.                                                                                     |
   |                    | - :optparam:`wf1` displacement threshold related to softening type 1. Default is wf1=0.1 wf. |
   |                    | - :optparam:`ft1` stress threshold related to softening type 1. Optional parameter. Default  |
   |                    |   is ft1=0.15 ft.                                                                            |
   |                    | - :optparam:`ahard` displacement threshold related to fracture energy used in all three      |
   |                    |   softening types. Default value is 0.1                                                      |
   |                    | - :optparam:`damage` flag to switch on and off damage. Default value is 1, which considers   |
   |                    |   damage.                                                                                    |
   |                    | - :optparam:`sub` maximum number of subincrementations in the plasticity part. Default is    |
   |                    |   10.                                                                                        |
   |                    | - :optparam:`tol` tolerance for the newton iteration of the plasticity part. Default value   |
   |                    |   is 1.e-6.                                                                                  |
   |                    | - :optparam:`iter` maximum number of iterations for the stress return of the plasticity      |
   |                    |   model. Default is 100.                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/latticeplastdamvisco3drandom.in <https://github.com/oofem/oofem/blob/devel/tests/regress |
   |                    | ion/lm/latticeplastdamvisco3drandom.in>`_, `lm/latticeplastdamvisco_mps_1.in <https://github |
   |                    | .com/oofem/oofem/blob/devel/tests/regression/lm/latticeplastdamvisco_mps_1.in>`_,            |
   |                    | `lm/latticeplastdamvisco_mps_2.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                    | n/lm/latticeplastdamvisco_mps_2.in>`_, `lm/latticeplastdamvisco_mps_3.in <https://github.com |
   |                    | /oofem/oofem/blob/devel/tests/regression/lm/latticeplastdamvisco_mps_3.in>`_                 |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _latticeframeelastic:

Elastic material model for lattice based frame elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an elastic material model used together with 3D lattice frame elements. It has been specifically developed to describe the failure of steel frame members using the lattice elements lattice3d and lattice3dnl. However, it can also be used with any other 3D lattice element.

The stress-strain law has the form

.. math::

   \mathbf{s} = \mathbf{D}_{\rm e} \mathbf{e}

where :math:`\mathbf{e}` is a vector of internal forces and moments, :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix, and :math:`\mathbf{e}` is a vector of generalized strains obtained from displacement jumps smeared over the element length and rotational components. The model parameters are summarized in Table :numref:`latticeframesteelplastic_table`.

.. table:: Elastic model for 3D lattice based frame elements -- summary.
   :name: latticeframeelastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Plasticity lattice model for steel                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeframeelastic` :elemparam:`in` :elemparam:`d{rn}`                           |
   |                    | :optelemparam:`talpha{rn}` :elemparam:`e{rn}` :elemparam:`n{rn}`                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal expansion coefficient. Default is 0.                            |
   |                    | - :param:`e` Young's modulus of the lattice material.                                        |
   |                    | - :param:`n` Poisson's ratio of the material that the beam element is made of.               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice3dnl_logrot.in                                                                    |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dnl_logrot.in>`_,     |
   |                    | `lm/lattice3dnlshell.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dnlshell.in>`_,       |
   |                    | `lm/lattice3dshell.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell.in>`_,         |
   |                    | `lm/lattice3dshell_combined_n1.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                    | n/lm/lattice3dshell_combined_n1.in>`_, `lm/lattice3dshell_combined_n4.in <https://github.com |
   |                    | /oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell_combined_n4.in>`_,                |
   |                    | `lm/latticeframe3dnl.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticeframe3dnl.in>`_,       |
   |                    | `lm/latticeframeelastictemp1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/ |
   |                    | lm/latticeframeelastictemp1.in>`_, `lm/latticeframesteelplastictemp1.in <https://github.com/ |
   |                    | oofem/oofem/blob/devel/tests/regression/lm/latticeframesteelplastictemp1.in>`_,              |
   |                    | `lm/latticelink3dnl_assembly.in <https://github.com/oofem/oofem/blob/devel/tests/regression/ |
   |                    | lm/latticelink3dnl_assembly.in>`_                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+

The linear-elastic response of the frame elements is given directly as a relationship between generalized stress (internal forces) and generalized strain.

.. _latticeframesteelplastic:

Steel Plasticity material model for lattice based frame elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a plasticity material used together with 3D lattice frame elements. It has been specifically developed to describe the failure process of steel frame members using the lattice elements lattice3d and lattice3dnl.

The generalised stress-strain law has the form

.. math::

   \mathbf{s} = \mathbf{D}_{\rm e} \left(\mathbf{e} - \mathbf{e}_{\rm p} \right)

where :math:`\mathbf{e}` is a vector of internal forces and moments, :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix, :math:`\mathbf{e}` is a vector of generalised strains obtained from displacement jumps smeared over the element length and rotational components, and :math:`\mathbf{e}_{\rm p}` are the plastic strains. The model parameters are summarised in Table :numref:`latticeframesteelplastic_table`.

.. table:: Plasticity steel model for 3D lattice frame elements -- summary
   :name: latticeframesteelplastic_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Plasticity lattice model for steel                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeframesteelplastic` :elemparam:`in` :elemparam:`d{rn}`                      |
   |                    | :optelemparam:`talpha{rn}` :elemparam:`e{rn}` :elemparam:`n{rn}` :elemparam:`nx0{rn}`        |
   |                    | :elemparam:`mx0{rn}` :elemparam:`my0{rn}` :elemparam:`mz0{rn}` :elemparam:`sub{in}`          |
   |                    | :elemparam:`iter{in}` :elemparam:`tol{rn}`                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal expansion coefficient. Default is 0.                            |
   |                    | - :param:`e` Young's modulus of the lattice material.                                        |
   |                    | - :param:`n` Poisson's ratio of the material that the beam element is made of.               |
   |                    | - :param:`nx0` ultimate capacity under pure axialloads.                                      |
   |                    | - :param:`mx0` ultimate capacity under pure moment about :math:`x`.                          |
   |                    | - :param:`my0` ultimate capacity under pure moment about :math:`y`.                          |
   |                    | - :param:`mz0` ultimate capacity under pure moment about :math:`z`.                          |
   |                    | - :optparam:`sub` maximum number of subincrementations.                                      |
   |                    | - :optparam:`iter` maximum number of newton iterations.                                      |
   |                    | - :optparam:`tol` tolerance for newton method.                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/latticeframesteelplastic1.in <https://github.com/oofem/oofem/blob/devel/tests/regression |
   |                    | /lm/latticeframesteelplastic1.in>`_, `lm/latticeframesteelplastic2.in <https://github.com/oo |
   |                    | fem/oofem/blob/devel/tests/regression/lm/latticeframesteelplastic2.in>`_                     |
   +--------------------+----------------------------------------------------------------------------------------------+

The nonlinear response of the beam element is given directly as a relationship between internal force and strain. The basic equations include an additive decomposition of total strain into elastic part and plastic part of the form

.. math::

   \mathbf{e} = \mathbf{e}_{\rm e} + \mathbf{e}_{\rm p}

The trial stress is calculated as

.. math::

   \mathbf{s}_{n+1}^{tr}  = \mathbf{D}^{\mathrm{e}} (\mathbf{e}_{n+1}  - \mathbf{e}_n^p)

The yield function is given as

.. math::

   \left(\dfrac{N_x}{N_0}\right)^2+ \left(\dfrac{M_x}{M_{x0}}\right)^2+ \left(\dfrac{M_y}{M_{y0}}\right)^2+ \left(\dfrac{M_z}{M_{z0}}\right)^2 - 1 = 0

where :math:`N_x`, :math:`M_x`, :math:`M_y` and :math:`M`, are the two components of bending moment, an axial force and a torsional moment, respectively. The corresponding components with the subscript 0 indicate a fully plastic value under the condition that each component of resultant forces acts independently.

The loading-unloading conditions are

.. math::

   f(\bar{\mathbf{s}},\kappa)\le 0 \qquad \dot{\lambda}\geq 0 \qquad \dot{\lambda}f(\bar{\mathbf{s}},\kappa)=0

The flow rule is

.. math::

   \dot{\mathbf{e}}_{\rm{p}} = \dot{\lambda} \dfrac{\partial f}{\partial \bar{\mathbf{s}}}


Reinforced concrete plasticity-damage material model for lattice based frame elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a plasticity-damage material used together with 3d lattice frame elements. It has been specifically developed to describe the failure process of reinforced concrete frame members using the lattice elements lattice3d and lattice3dnl.

The generalised stress-strain law has the form

.. math::

   \mathbf{s} = \mathbf{D}_{\rm e} \left(\mathbf{e}-\mathbf{e}_{\rm p} \right)

where :math:`\mathbf{e}` is a vector of internal forces and moments, :math:`\mathbf{D}_{\rm e}` is the elastic stiffness matrix, :math:`\mathbf{e}` is a vector of translation and rotational jumps smeared over the element length and :math:`\mathbf{e}_{\rm p}` are the plastic smeared jumps. The model parameters are summarised in Table :numref:`latticeframesteelplastic_table`.

.. table:: Damage Plasticity concrete model for 3D lattice frame elements -- summary

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Damage-plasticity lattice model for reinforced concrete                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticeframeconcreteplastic` :elemparam:`in` :elemparam:`d{rn}`                   |
   |                    | :optelemparam:`talpha{rn}` :elemparam:`e{rn}` :elemparam:`n{rn}` :elemparam:`nx0{rn}`        |
   |                    | :elemparam:`nx1{rn}` :elemparam:`mx0{rn}` :elemparam:`mx1{rn}` :elemparam:`my0{rn}`          |
   |                    | :elemparam:`my1{rn}` :elemparam:`mz0{rn}` :elemparam:`mz1{rn}` :elemparam:`wu{rn}`           |
   |                    | :elemparam:`wf{rn}` :elemparam:`sub{in}` :elemparam:`iter{in}` :elemparam:`tol{rn}`          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :optparam:`talpha` Thermal expansion coefficient. Default is 0.                            |
   |                    | - :param:`e` Young's modulus of the lattice material.                                        |
   |                    | - :param:`n` Poisson's ratio of the material that the beam element is made of.               |
   |                    | - :param:`nx0` ultimate positive capacity under pure axialloads.                             |
   |                    | - :param:`mx0` ultimate positive capacity under pure moment about :math:`x`.                 |
   |                    | - :param:`my0` ultimate positive capacity under pure moment about :math:`y`.                 |
   |                    | - :param:`mz0` ultimate positive capacity under pure moment about :math:`z`.                 |
   |                    | - :param:`nx1` ultimate negative capacity under pure axialloads.                             |
   |                    | - :param:`mx1` ultimate negative capacity under pure moment about :math:`x`.                 |
   |                    | - :param:`my1` ultimate negative capacity under pure moment about :math:`y`.                 |
   |                    | - :param:`mz1` ultimate negative capacity under pure moment about :math:`z`.                 |
   |                    | - :param:`wu` plastic generalised displacement jump norm (computed from rotational           |
   |                    |   components only) threshold.                                                                |
   |                    | - :param:`wf` threshold controlling the evolution of the damage.                             |
   |                    | - :optparam:`sub` maximum number of subincrementations.                                      |
   |                    | - :optparam:`iter` maximum number of newton iterations.                                      |
   |                    | - :optparam:`tol` tolerance for newton method.                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 3dlattice                                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/latticeframeconcreteplastic.in <https://github.com/oofem/oofem/blob/devel/tests/regressi |
   |                    | on/lm/latticeframeconcreteplastic.in>`_                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+

.. |latticeframesteelplastic_table| replace:: latticeframesteelplastic_table

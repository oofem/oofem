Material models for interfaces
==============================


.. _material_models_for_interfaces:

Material models for interfaces
------------------------------

Interface elements have to be used with material models describing
the constitutive behavior of interfaces between two materials 
(e.g.\ between steel reinforcement and concrete),
or between two bodies in contact.
Such interface laws are formulated in terms of the traction vector :math:`t` and the displacement jump vector :math:`\delta`.

.. _cohesive_interface_material_cohInt:

Cohesive interface material - cohInt
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A simple interface material with generally different stiffness in compression, tension and shear. It is intended for simulating a contact with low stiffness in tension and high in compression. The parameter :param:`transitionOpening` specifies initial gap embedded in an element which causes transition from tension to compression. For example, :param:`transitionOpening`\ =0.01~m means that there exists the embedded gap which needs displacement jump 0.01~m to close. Traction-separation law for normal direction takes the following form:

.. math::
   :label: eq:cohint1

   t_n &=& \left (\frac{\pi}{2} + \arctan (s_m \bar\delta)\right) \frac{k_t}{\pi} \bar\delta + \left (\frac{\pi}{2} - \arctan (s_m \bar\delta)\right) \frac{k_n}{\pi} \bar\delta\\
   s_m &=& \mathit{smoothMag}\\
   \delta_0 &=& \mathit{transitionOpening}\\
   \bar\delta &=& \delta_n + \delta_0\\
   k_t &=& k_n \cdot \mathit{stiffCoeffKn}

where :param:`smoothMag` controls smoothing magnitude. Tangential stiffness for normal direction is found by differentiating Eq.~(:eq:`eq:cohint1`)

.. math::

   k_{nn} &=& \left (\frac{\pi}{2} + \arctan (s_m \bar\delta)\right) \frac{k_t}{\pi} + \left (\frac{\pi}{2} - \arctan (s_m \bar\delta)\right) \frac{k_n}{\pi} +\\
   &+& \frac{s_m \cdot k_t \bar\delta}{\pi(s_m^2\bar\delta^2+1)} - \frac{s_m \cdot k_n \bar\delta}{\pi(s_m^2\bar\delta^2+1)}

Shear stiffness remains constant during all possible loadings and there is no influence of normal direction.

.. table:: Cohesive interface material -- summary.
   :name: cohint_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Cohesive interface material                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`cohInt` :elemparam:`in` :elemparam:`kn{rn}` :elemparam:`ks{rn}`                   |
   |                    | :optelemparam:`stiffCoeffKn{rn}` :optelemparam:`smoothMag{rn}`                               |
   |                    | :optelemparam:`transitionOpening{rn}`                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`kn` (penalty) stiffness in compression                                             |
   |                    | - :param:`ks` stiffness in shear                                                             |
   |                    | - :param:`stiffCoeffKn` ratio (tensile stiffness / compression stiffness)                    |
   |                    | - :param:`smoothMag` smoothing parameter for transition between tensile/compressive behavior |
   |                    | - :param:`transitionOpening` embedded gap in the element                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _1dInterface,_2dInterface,_3dInterface                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/InterfaceEL_Line1.in                                                                     |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Line1.in>`_,      |
   |                    | `sm/InterfaceEL_Point3D_01.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_01.in>`_, |
   |                    | `sm/InterfaceEL_Point3D_02.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_02.in>`_, |
   |                    | `sm/InterfaceEL_Point3D_03.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_03.in>`_  |
   +--------------------+----------------------------------------------------------------------------------------------+


Isotropic damage model for interfaces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This model is described in Section :ref:`sec:idmfi`.

Simple interface material - obsolete

This model provides a simple interface law with penalty-type contact and friction. 
In the normal direction, the response is linear elastic, but with different stiffnesses in tension and in compression (stiffness :param:`kn` in compression, :math:`:param:`\ kn`*:param:`stiffcoeff`\ ` in tension). 
By setting :param:`kn` to a high value, the penetration (overlap) can be reduced,
in the sense of the penalty approach. By setting  :param:`stiffcoeff` to 0, free opening
of the gap can be allowed. 

The shear response is elastoplastic, with the yield limit dependent on the normal traction. Crosssection's width, height or area have no influence on the results.
The magnitude of the shear traction :math:`\vsig_T` must not exceed the yield limit computed 
according to a Coulomb-like friction law as the product of the (negative part of)
normal traction :math:`\sigma_N` and a dimensionless coefficient of friction :math:`fc`:

.. math::

   ||\vsig_T|| \leq \left\langle-\sigma_N\cdot fc\right\rangle 

:math:`\sigma_N` is computed multiplying a known normal strain and a known stiffness in 
tension or in compression.

If :param:`regularized` is set to true, the regularized version of the model is used. The regularization parameter :param:`m` allows to control how close to the original formulation the regularized model is. The default value is 15; the higher the value, the closer the responses are. The regularized version has the following form:

.. math::

   \sigma_N &=& 0.5 k_n(\delta_0+\delta_n)-
               \del{0.5 k_n}{m}\log(|(\cosh(m(\delta_0+\delta_n))|)+\\
            &+&c\left(0.5 k_n(\delta_0+\delta_n)+
               \del{0.5 k_n}{m}\log(|(\cosh(m(\delta_0+\delta_n)))|)\right)

where :math:`c` is the ratio of tensile/compressive stiffness (:param:`stiffcoeff`), :math:`\delta_0` is the normal clearance, :math:`\delta_n` is the normal strain (jump), and :math:`m` denotes the regularization coefficient. Consistent tangent stiffness is computed by differentiating the above relation.

The shear elastic stiffness is assumed to be :param:`kn` (i.e., equal to the normal stiffness). If the normal traction is tensile or zero, no shear traction can be transmitted by the interface. In the future, this law will be enriched by a cohesion-like term.

If used in connection with two nodes with :tt:`interface1delement`, the element-level parameter :param:`normal` determines what is compression and what is tension (from the 1st node to the 2nd node in :param:`normal` direction, it is tension). If :param:`normalClearance` is defined, the element has tensile stiffness until the gap within the element is closed. This feature allows simulating two faces with a gap. When the gap is smaller than :param:`normalClearance`, the element stiffness changes to compression, see Figure :numref:`SimpleInterfaceMat`.

The model description and parameters are summarized in :numref:`simpleinterfacemat_table`.

.. table:: Simple interface material
   :name: simpleinterfacemat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Simple interface material                                                                    |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`simpleintermat` :elemparam:`in` :elemparam:`kn{rn}` :optelemparam:`fc{rn}`        |
   |                    | :optelemparam:`stiffcoeff{rn}` :optelemparam:`normalClearance{rn}`                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`kn` (penalty) stiffness in compression                                             |
   |                    | - :param:`fc` friction coefficient                                                           |
   |                    | - :param:`stiffcoeff` ratio (tensile stiffness / compression stiffness)                      |
   |                    | - :param:`normalClearance` free distance within element                                      |
   |                    | - :param:`regularized` when true use regularized formulation (Default is False)              |
   |                    | - :param:`m` regularization coefficient :math:`m` (default value is 15)                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _1dInterface                                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/interface01.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/interface01.in>`_             |
   +--------------------+----------------------------------------------------------------------------------------------+

:name: `simpleinterfacemat_table`

.. figure:: /figures/Simple_interface_material_diag.png
   :width: 70%
   :alt: Working diagram of SimpleInterfaceMat.
   :name: SimpleInterfaceMat

   Working diagram of SimpleInterfaceMat.


Bond-slip model for reinforced concrete - BondCEB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


The complex steel/concrete interaction can be modelled with a bond-slip relation, which describes the bond stress (tangential interface traction), :math:`\tau`, in terms of the relative reinforcement slip (tangential interface jump), :math:`s`. This interface material model is based on the local bond--slip relationship for reinforced concrete under good bond conditions outlined in the the *fib* Model Code for Concrete Structures 2010 [fib:2010]_.

In general, the model is formulated in terms of the interface traction vector :math:`\mbf{t}` and the interface jump vector :math:`\mbf{\delta}`. In the model, the interface is assumed to have only elastic stiffness :math:`k_{\text{n}}` in the normal direction. The normal interface traction is then evaluated elastically from the normal interface jump, i.e. :math:`t_{\text{n}} =  k_{\text{n}} \delta_{\text{n}}`. The tangential traction :math:`t_{\text{t}} = \tau` is evaluated from the user-specified function :math:`\tau (s)`, see Figure :numref:`BondCEB`:

.. math::

   \tau (s) = \begin{cases}
       \tau_\text{max} \left( \dfrac{s}{s_1} \right) ^ \alpha & \text{for } 0 \leq s \leq s_1, \\
       \tau_\text{max} & \text{for } s_1 < s \leq s_2, \\
       \tau_\text{max} - \dfrac{\left( \tau_\text{max} - \tau_\text{f} \right) \left(s-s_2 \right)}{\left( s_3 - s_2 \right)} & \text{for } s_2 < s \leq s_3, \\
       \tau_\text{f} & \text{for } s > s_3.
   \end{cases}

For an exponent :math:`\alpha` smaller than 1 (default value 0.4), the initial tangential stiffness of the interface, :math:`k_{\text{s}}`, is undefined and cannot be used in numerical computations. Therefore, this stiffness needs to be specified manually. Should the provided stiffness be smaller than :math:`k_n = \dfrac{\tau_\text{max}}{s_1}`, it will be automatically adjusted to this value. Note that only elastic tangent stiffness is supported. Hence, in 2D the tangent stiffness matrix takes the form

.. math::

   \dfrac{\partial \mbf{t}}{\partial \mbf{\delta}}= \left[ \begin{array}{cc}
   \dfrac{\partial t_\text{n}}{\partial \delta_\text{n}} &
   \dfrac{\partial t_\text{n}}{\partial s} \\
   \dfrac{\partial \tau}{\partial \delta_\text{n}} & \dfrac{\partial \tau}{\partial s}
   \end{array} \right] =
   \left[ \begin{array}{cc}
   k_\text{n} & 0 \\ 0 & k_\text{s}
   \end{array} \right]

Model description and the input parameters are summarized in :numref:`bondceb_table`.

.. figure:: /figures/bondceb_diag.png
   :width: 50%
   :alt: Bond stress (:math:`\tau`) - reinforcement slip (:math:`s`) diagram for BondCEB model.
   :name: BondCEB

   Bond stress (:math:`\tau`) - reinforcement slip (:math:`s`) diagram for BondCEB model.

.. table:: Bond-slip model for reinforced concrete -- summary.
   :name: bondceb_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Bond-slip model for reinforced concrete                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`bondceb` :elemparam:`in` :elemparam:`kn{rn}` :elemparam:`ks{rn}`                  |
   |                    | :elemparam:`s1{rn}` :elemparam:`s2{rn}` :elemparam:`s3{rn}` :elemparam:`taumax{rn}`          |
   |                    | :optelemparam:`tauf{rn}` :optelemparam:`alpha{rn}`                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`kn` interface elastic normal stiffness                                             |
   |                    | - :param:`ks` interface elastic tangential (shear) stiffness                                 |
   |                    | - :param:`s1` characteristic slip value                                                      |
   |                    | - :param:`s2` characteristic slip value                                                      |
   |                    | - :param:`s3` characteristic slip value                                                      |
   |                    | - :param:`taumax` maximum bond stress                                                        |
   |                    | - :param:`tauf` residual bond stress at reinforcement pull-out                               |
   |                    | - :param:`alpha` bond-slip curve parameter (exponent)                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dInterface, _3dInterface                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/bondceb01.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bondceb01.in>`_,              |
   |                    | `sm/bondceb02.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bondceb02.in>`_               |
   +--------------------+----------------------------------------------------------------------------------------------+

.. _linkslip:

LinkSlip - Slip model for trusses and beams
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is an elasto-plastic model for slip between reinforcement modelled as three-dimensional beam or truss elements and matrix elements. The theory of the model is described in the paper [SciGraLarRun20]_ for the case of beams. The model parameters are summarised in :numref:`linkslip_table`.

There are three types of constitutive relations supported which are shown in :ref:`linkslip`, which are chosen by means of the input parameter :param:`type`.

.. figure:: /figures/linkslip.png
   :width: 80%
   :alt: Bond stress (:math:`\tau`) - slip (:math:`s`) diagram for link slip model.

   Bond stress (:math:`\tau`) - slip (:math:`s`) diagram for link slip model.


.. table:: Slip model for truss and beam elements -- summary
   :name: linkslip_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Slip model                                                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`linkslip` :elemparam:`in` :elemparam:`kn{rn}` :elemparam:`kl{rn}`                 |
   |                    | :optelemparam:`type{in}` :elemparam:`s1{rn}` :elemparam:`s2{rn}` :elemparam:`s3{rn}`         |
   |                    | :elemparam:`t0{rn}` :optelemparam:`tf{rn}` :optelemparam:`alpha{rn}`                         |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`kn` interface elastic axial stiffness (in the direction of slip)                   |
   |                    | - :param:`kl` interface elastic lateral stiffness                                            |
   |                    | - :param:`type` type of bond law (default is zero)                                           |
   |                    | - :param:`s1` characteristic slip value                                                      |
   |                    | - :param:`s2` characteristic slip value                                                      |
   |                    | - :param:`s3` characteristic slip value                                                      |
   |                    | - :param:`t0` maximum bond stress                                                            |
   |                    | - :param:`tf` residual bond stress at reinforcement pull-out                                 |
   |                    | - :param:`alpha` bond-slip curve parameter (exponent)                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _3dInterface                                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/bond_link_1.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_1.in>`_,            |
   |                    | `sm/bond_link_2.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_2.in>`_,            |
   |                    | `sm/bond_link_3.in                                                                           |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_3.in>`_,            |
   |                    | `sm/linkslip01.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip01.in>`_,             |
   |                    | `sm/linkslip02.in                                                                            |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip02.in>`_              |
   +--------------------+----------------------------------------------------------------------------------------------+

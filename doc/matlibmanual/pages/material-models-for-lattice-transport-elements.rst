Material models for lattice transport elements
==============================================


.. _latticetransmat:

Material for unsaturated flow in lattice models -- LatticeTransMat
------------------------------------------------------------------

Lattice transport elements have to be used with this material model. A positive sign is assumed for liquid tension, unlike the convention of soil mechanics which assumes compression positive.

These transport elements are idealised as one-dimensional conductive pipes. The gradient of hydraulic head, which governs flow rate along each transport element, is determined from the capillary pressures :math:`P_{c}` at the two nodes.

The mass balance equation describes the change in moisture inside a porous element as a consequence of liquid flow and solid-liquid retention. It leads to the following partial differential equation

.. math::
   :label: latticetransmat1

   c\frac{\partial P_{c}}{\partial t} + k\mbox{div}(\nabla(\frac{P_{c}}{g}-\rho z)) =0

where :math:`P_{c}` is the capillary pressure, :math:`c` is the mass capacity function(:math:`s^{2}m^{-2}`), :math:`k` is the Darcy hydraulic conductivity (:math:`ms^{-1}`), :math:`\rho` is the fluid mass density, :math:`g` is the acceleration of gravity, :math:`z` is the capillary height and :math:`t` is the time.

The hydraulic conductivity :math:`k` consists of

.. math::
   :label: latticetransmat11

   k=k_{0} + k_{c}

where :math:`k_{0}` is the hydraulic conductivity of the intact material and :math:`k_{c}` is the additional conductivity due to cracking.

Darcy hydraulic conductivity :math:`k_{0}` is defined as

.. math::
   :label: latticetransmat12

   k_{0}=\frac{\rho g}{\mu}\kappa

where :math:`\mu` is the dynamic viscosity (:math:`Pa.s`), :math:`\kappa` is the permeability also called intrinsic conductivity(:math:`m^{2}`), and :math:`\kappa_{r}` is the relative conductivity. :math:`\kappa_{r}` is a function of the effective degree of saturation.

The cracking part is

.. math::
   :label: latticetransmat13

   k_{c}=\xi\frac{\rho g}{\mu}\frac{w^{3}_{c}}{12h} \kappa_{r}

where :math:`\xi` is a tortuosity factor taking into account the roughness of the crack surface, :math:`w_{c}` is the equivalent crack opening of the dual mechanical lattice and :math:`h` is the length of the dual mechanical element.

.. _one-dimensional-transport-element:

One-dimensional transport element
---------------------------------

The discrete form of the differential equation for mass transport for a one-dimensional transport element is

.. math::
   :label: latticetransmat3

   \alpha_{e}P_{c} + C_{e}\frac{\partial P_{c}}{\partial t}  = f_{e}

where :math:`P_{c}` is a vector containing the nodal values of the capillary pressure, :math:`\alpha_{e}` is the conductivity matrix, :math:`C_{e}` is the capacity matrix and :math:`f_{e}` is the nodal flow rate vector (:math:`kg.m^{3}`). 

The capacity matrix is

.. math::
   :label: latticetransmat4

   C_{e} = \frac{Al}{12}c \left(  \begin{array}{ c c }  2 & 1 \\  1 & 2 \end{array}  \right)

where :math:`c` is the capacity of the material, :math:`l` is the length of the transport element and :math:`A` is the cross-sectional area of the transport element. 

The conductivity matrix is defined as

.. math::
   :label: latticetransmat5

   \alpha_{e} = \frac{A}{l}\frac{k}{g} \left(  \begin{array}{ c c }  1 & -1 \\  -1 & 1 \end{array}  \right)

The mass transport equation is based on the constitutive laws for the capacity :math:`c` and the hydraulic conductivity :math:`k`.

The capacity :math:`c` is defined as 

.. math::
   :label: latticetransmat6

   c=-\rho\frac{\partial \theta}{\partial P_{c}}

where :math:`\theta` is the volumetric water content (:math:`\theta=\frac{V_{w}}{V_{T}}` with :math:`V_{w}` the volume of water, and :math:`V_{T}` the total volume) which is calculated by a modified version of van Genuchten’s retention model. Note that the presence of a crack in an element does not influence the capacity in the present model.

The volumetric water content is 

.. math::
   :label: latticetransmat7

   \theta = S_{e} (\theta_{s} - \theta_{r})+\theta_{r}

where :math:`\theta_{r}` and :math:`\theta_{s}` are the residual and saturated water contents corresponding to effective saturation values of :math:`S_{e} = 0` and :math:`S_{e} = 1`, respectively.

The effective degree of saturation :math:`S_{e}` is defined as

.. math::
   :label: latticetransmat8

   S_{e}  = \left\{
       \begin{array}{ll}
            \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}\left(1+\left(\frac{ P_{c}}{a}\right)^{\frac{1}{1-m}}\right)^{-m}  & \mbox{if } P_{c} \ge P_{c(aev)} \\
            1 & \mbox{if } P_{c}<P_{c(aev)}
       \end{array}
   \right.

where :math:`\theta_{m}` is an additional model parameter and :math:`P_{c(aev)}` is the air-entry value of capillary pressure which separates saturated (:math:`P_{c}<P_{c(aev)}`) from unsaturated states (:math:`P_{c} \ge P_{c(aev)}`). It is intuitive that the smaller the pore size of the material, the larger the value of :math:`P_{c(aev)}` will be.

The relative conductivity :math:`\kappa_{r}` is a function of the effective degree of saturation and is defined as

.. math::
   :label: latticetransmat9

   \kappa_{r}  = \sqrt{S_{e}}\left(\frac{1-\left[1-\left(\frac{S_{e}}{  \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}}\right)^{\frac{1}{m}}\right]^{m}}{1-\left[1-\left(\frac{1}{  \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}}\right)^{\frac{1}{m}}\right]^{m}}\right)^{2}

If :math:`\theta_{m} = \theta_{s}`, we have :math:` \frac{\theta_{m}-\theta_{r}}{\theta_{s}-\theta_{r}}=1` : the equation reduces to the expression of the relative conductivity of the original van Genuchten model.

The model parameters are summarized

.. table:: Model Parameters
   :widths: 15 85

   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Description   | Material for fluid transport in lattice models                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Record Format |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
   | Parameters    | - :param:`num` material model number - :param:`d` fluid mass density - :param:`k` permeability (:math:`m^{2}`) - :param:`vis` dynamic viscosity (:math:`Pa.s`) - :param:`contype` unsaturated flow allowed when :param:`contype`\ =1 - :param:`thetas` saturated water content - :param:`thetar` residual water content - :param:`paev` air-entry value of capillary pressure - :param:`m` van Genuchten parameter - :param:`a` van Genuchten parameter - :param:`thetam` additional model parameter for the modified version of van Genuchten’s retention model - :param:`ctor` coefficient of tortuosity (:math:`ctor=\frac{1}{\tau}\le1`) |
   +---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. table:: Material for unsaturated flow in lattice models - summary.
   :name: Iatticetransmat_table
   :widths: 20 94

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Material for fluid transport in lattice models                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`latticetransmat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`k{rn}`       |
   |                    | :elemparam:`vis{rn}` :elemparam:`contype{in}` :elemparam:`thetas{rn}`                        |
   |                    | :elemparam:`thetar{rn}` :elemparam:`paev{rn}` :elemparam:`m{rn}` :elemparam:`a{rn}`          |
   |                    | :elemparam:`thetam{rn}` :optelemparam:`ctor{rn}`                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` fluid mass density                                                              |
   |                    | - :param:`k` permeability (:math:`m^{2}`)                                                    |
   |                    | - :param:`vis` dynamic viscosity (:math:`Pa.s`)                                              |
   |                    | - :param:`contype` unsaturated flow allowed when contype=1                                   |
   |                    | - :param:`thetas` saturated water content                                                    |
   |                    | - :param:`thetar` residual water content                                                     |
   |                    | - :param:`paev` air-entry value of capillary pressure                                        |
   |                    | - :param:`m` van Genuchten parameter                                                         |
   |                    | - :param:`a` van Genuchten parameter                                                         |
   |                    | - :param:`thetam` additional model parameter for the modified version of van Genuchten’s     |
   |                    |   retention model                                                                            |
   |                    | - :param:`ctor` coefficient of tortuosity (:math:`ctor=\frac{1}{\tau}\le1`)                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dMassLatticeTransport                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `lm/lattice2dcrackinput.in                                                                   |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2dcrackinput.in>`_,    |
   |                    | `lm/lattice3d_mt1.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt1.in>`_,          |
   |                    | `lm/lattice3d_mt2.in                                                                         |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt2.in>`_,          |
   |                    | `lm/latticetransmat.in                                                                       |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticetransmat.in>`_         |
   +--------------------+----------------------------------------------------------------------------------------------+

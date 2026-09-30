Isotropic material for moisture transport based on Ba zant and Najjar - BazantNajjarMoisture
============================================================================================


.. _sec:BazantNajjarMoistureMat:

Isotropic material for moisture transport based on Bažant and Najjar -- BazantNajjarMoisture
--------------------------------------------------------------------------------------------

This is a specific model for nonlinear moisture transport in isotropic cementitious materials, based on [Bazant:72]_.

The governing equation

.. math::
   :label: BNmodel:governing

   \frac{\partial h}{\partial t} = \nabla \cdot \left( C(h) \nabla h \right)

is a special case of :eq:`nlisomoisture:governing`, valid under the assumption that the slope of the sorption isotherm is linear, i.e.\ the moisture capacity is constant. In :eq:`BNmodel:governing`, :math:`h` is the relative humidity and :math:`C(h)` is the humidity-dependent diffusivity approximated by

.. math::
   :label: BNmodel:diffusivity

   C (h) = C_1 \left( \alpha_0
   + \frac{1-\alpha_0}{1+\left(\frac{1-h}{1-h_c}\right)^n} \right)

where :math:`C_1` is the diffusivity at saturation (typical value for concrete :math:`\approx 30` mm\ :math:`^2`/day), :math:`\alpha_0` is the dimensionless ratio of diffusivity at low humidity to diffusivity at saturation (typical value :math:`\approx 0.05`), :math:`h_c` is the humidity “in the middle” of the transition between low and high diffusivity (typical value :math:`\approx 0.8`), and :math:`n` is dimensionless exponent (high values of :math:`n`, e.g.\ 12, lead to a rapid transition between low and high diffusivity). Optionally, it is possible to specify the moisture capacity. This property is not needed for solution of the diffusion equation  :eq:`BNmodel:governing`, but it is needed if the computed change of relative humidity is transformed into water content loss (mass of lost water per unit volume).

The model parameters are summarized in :numref:`BazantNajjarMoistureMat`.

.. table:: Model parameters
   :name: BazantNajjarMoistureMat

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Nonlinear isotropic material for moisture transport                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`BazantNajjarMoistureMat` :elemparam:`num{in}` :elemparam:`d{rn}`                  |
   |                    | :elemparam:`c1{rn}` :elemparam:`n{rn}` :elemparam:`alpha0{rn}` :elemparam:`hc{rn}`           |
   |                    | :optelemparam:`capa{rn}`                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`c1` moisture diffusivity at full saturation [m\ :math:`^2` s\ :math:`^{-1}`]       |
   |                    | - :param:`n` exponent [-]                                                                    |
   |                    | - :param:`alpha0` ratio between minimum and maximum diffusivity [-]                          |
   |                    | - :param:`hc` relative humidity at which the diffusivity is exactly between its minimum and  |
   |                    |   maximum value [-]                                                                          |
   |                    | - :param:`capa` moisture capacity (default value is 1.0)                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/bazantnajjar.in                                                                          |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/bazantnajjar.in>`_            |
   +--------------------+----------------------------------------------------------------------------------------------+

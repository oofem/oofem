
.. _BinghamFluidMaterial:

Bingham fluid - BinghamFluid
============================

Constitutive model of Bingham fluid. This is a constitutive model of
non-Newtonian type. The model parameters are summarized
in :numref:`BinghamFluidMaterial_table`.

In the Bingham model the flow is characterized by following
constitutive equation

.. math::
   :label: eq:bigham-model

   \mbf{\tau} &=& \mbf{\tau}_0 + \mu \mbf{\dot\gamma} \;\;\;\text{if}\ \dot\tau \ge \tau_0 \\
   \dot\gamma &=& \mbf{0} \;\;\;\;\;\;\;\;\;\;\;\;\;\;\text{if}\ \dot\tau \le \tau_0

where :math:`\mbf{\tau}` is the shear stress applied to material, :math:`\dot\tau = \sqrt{\mbf{\tau} : \mbf{\tau}}`
is the shear stress measure, :math:`\mbf{\dot\gamma}` is
the shear rate, :math:`\mbf{\tau}_0` is the yield stress, and :math:`\mu` is the plastic
viscosity.
The parameters for the model can be in general determined using two
possibilities: (i) stress controlled rheometer, when the stress is applied
to material and shear rate is measured, and (ii) shear rate controlled
rheometer, where concrete is sheared and stress is measured. However,
most of the widely used tests are unsatisfactory in the sense, that
they measure only one parameter. These one-factor tests include slump
test, penetrating rod test, and Ve-Be test. Recently, some tests
providing two parameters on output have been designed (BTRHEOM, IBB,
and BML rheometers). Also a refined version of the
standard slump test has been developed for estimating yield stress and
plastic viscosity. The test is based on measuring the time necessary
for the upper surface of the concrete cone in the slump to fall a
distance 100 mm. Semi-empirical models are then proposed for estimating
yield stress and viscosity based on measured results. The advantage
is, that this test does not require any special equipment, provided that
the one for the standard version is available.

In order to avoid numerical difficulties caused by the existence of
the sharp angle in material model
response at :math:`\tau = \tau_0`, the numerical implementation uses
following smoothed relation for viscosity

.. math::
   :label: eq:smooth-bigham-model

   \mu = \mu_0 + \frac{\tau_0}{\dot\gamma}(1 - e^{-m\dot\gamma})

where :math:`m` is so called stress growth parameter. The higher value of
parameter :math:`m`, the closer approximation of the original
constitutive equation :eq:`eq:bigham-model` is obtained.

.. table:: Bingham Fluid material - summary.
   :name: BinghamFluidMaterial_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Bingham fluid material                                                                       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`BinghamFluid` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`mu0{rn}`        |
   |                    | :elemparam:`tau0{rn}`                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`mu0` viscosity                                                                     |
   |                    | - :param:`tau0` Yield stress                                                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2d, 3d flow                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+

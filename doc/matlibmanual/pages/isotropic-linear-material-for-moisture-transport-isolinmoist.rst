
.. _IsoLinMoistureMat:

Isotropic linear material for moisture transport -- IsoLinMoisture
==================================================================

Linear isotropic material model for moisture transport problems described
by the linear diffusion equation. Note that the symbols :math:`k` and :math:`c`
in :eq:`lindiff` have a different meaning than in :eq:`lindiffheat`.
The reason is that the nonlinear model for moisture transport described
in Section :ref:`sec:NlIsoMoistureMat` traditionally
uses :math:`c` for permeability and :eq:`lindiff` should be obtained
as a special case of :eq:`nlisomoisture:governing`. 
On the other hand, the heat conduction model
from Section ``IsoLET`` was implemented earlier and the input parameters
are directly called :math:`k` and :math:`c`, so changing this notation now could lead
to confusion for some older input files.

.. math::
   :label: lindiff

   k\frac{\partial h}{\partial t} = \nabla \cdot \left( c \nabla h \right)

where :math:`h` is the pore relative humidity (dimensionless, between 0 and 1), 
:math:`k` is the moisture capacity [kg/m\ :sup:`3`],
and :math:`c` is the moisture permeability [kg/m\ :sup:`\cdot`\ s].

The model parameters are summarized in :numref:`IsoLinmoistureMat_table`.

.. table:: Linear isotropic material for moisture transport - summary.
   :name: IsoLinmoistureMat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Linear isotropic material for moisture transport                                             |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`IsoLinMoistureMat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`perm{rn}`  |
   |                    | :elemparam:`capa{rn}`                                                                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`perm` moisture permeability                                                        |
   |                    | - :param:`capa` moisture capacity                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/isolinmoisture.in                                                                        |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/isolinmoisture.in>`_          |
   +--------------------+----------------------------------------------------------------------------------------------+

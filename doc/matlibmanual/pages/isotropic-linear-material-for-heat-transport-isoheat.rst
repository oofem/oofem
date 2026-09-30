
Isotropic linear material for heat transport -- IsoHeat
=======================================================

Linear isotropic material model for heat transport problems described
by the linear diffusion equation

.. math::
   :label: lindiffheat

   c\frac{\partial T}{\partial t} = \nabla \cdot \left( k \nabla T \right)

where :math:`T` is the temperature, :math:`c` is the specific heat capacity,
and :math:`k` is the conductivity.

The model parameters are summarized
in :numref:`Isoheat_table`.

.. table:: Linear isotropic material for heat transport - summary.
   :name: Isoheat_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Linear isotropic elastic material                                                            |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`IsoHeat` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`k{rn}`               |
   |                    | :elemparam:`c{rn}`                                                                           |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`k` Conductivity                                                                    |
   |                    | - :param:`c` Specific heat capacity                                                          |
   |                    | - :param:`maturityT0` Baseline value for material maturity, default 0                        |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | _2dHeat                                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `tm/hydratingConcreteMat07.in                                                                |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat07.in>`_, |
   |                    | `tm/line01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/line01.in>`_,   |
   |                    | `tm/line02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/line02.in>`_,   |
   |                    | `tm/qbrick_01.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_01.in>`_,              |
   |                    | `tm/qbrick_02.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_02.in>`_,              |
   |                    | `tm/qbrick_03.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_03.in>`_,              |
   |                    | `tm/qquad01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qquad01.in>`_, |
   |                    | `tm/tmpatch31.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch31.in>`_,              |
   |                    | `tm/tmpatch40.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch40.in>`_,              |
   |                    | `tm/tmquad12.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmquad12.in>`_ (and 37 more)  |
   +--------------------+----------------------------------------------------------------------------------------------+


Newtonian fluid - NewtonianFluid
================================

Constitutive model of Newtonian fluid. The model parameters are summarized
in :numref:`NewtonianFluidMaterial_table`.

.. table:: Newtonian Fluid material - summary.
   :name: NewtonianFluidMaterial_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Newtonian Fluid material                                                                     |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`NewtonianFluid` :elemparam:`num{in}` :elemparam:`d{rn}` :elemparam:`mu{rn}`       |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - :param:`num` material model number                                                         |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`mu` viscosity                                                                      |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2d, 3d flow                                                                                  |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `benchmark/fm/axi01.oofem.in                                                                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/fm/axi01.oofem.in>`_,  |
   |                    | `benchmark/fm/axi02.oofem.in                                                                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/fm/axi02.oofem.in>`_,  |
   |                    | `benchmark/fm/axi04.oofem.in                                                                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/fm/axi04.oofem.in>`_,  |
   |                    | `benchmark/fm/bdam7.oofem.in                                                                 |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/fm/bdam7.oofem.in>`_,  |
   |                    | `fm/cbs1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs1.in>`_,       |
   |                    | `fm/cbs2.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs2.in>`_,       |
   |                    | `fm/cbs3.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs3.in>`_,       |
   |                    | `fm/scctest01.in                                                                             |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/scctest01.in>`_,              |
   |                    | `fmpfem/pfemFreeFall.in                                                                      |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/fmpfem/pfemFreeFall.in>`_,       |
   |                    | `fmpfem/pfemPrescribedVelocity.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                    | n/fmpfem/pfemPrescribedVelocity.in>`_ (and 2 more)                                           |
   +--------------------+----------------------------------------------------------------------------------------------+

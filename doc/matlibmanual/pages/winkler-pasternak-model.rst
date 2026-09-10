
.. _WinklerPasternak_model:

Winkler-Pasternak model
=======================

Implementation of 2D Winkler-Pasternak model for plate (and potentially beam) subsoil model.

.. table:: Winkler Pasternak material - summary.
   :name: WinklerPasternak_table

   +--------------------+----------------------------------------------------------------------------------------------+
   | Description        | Winkler Pasternak isotropic material for subsoil interaction                                 |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Record Format      | :descitem:`winklerpasternak` :elemparam:`in` :elemparam:`d{rn}` :elemparam:`c1{rn}`          |
   |                    | :elemparam:`c2{rn}`                                                                          |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Parameters         | - material number                                                                            |
   |                    | - :param:`d` material density                                                                |
   |                    | - :param:`c1` C1 Winkler parameter (foundation modulus)                                      |
   |                    | - :param:`c2` C2 shear interaction modulus                                                   |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Supported modes    | 2dPlateSubSoil                                                                               |
   +--------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples     | `sm/test_wp1.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp1.in>`_,               |
   |                    | `sm/test_wp2.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp2.in>`_,               |
   |                    | `sm/test_wp3.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp3.in>`_,               |
   |                    | `sm/test_wp4.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp4.in>`_,               |
   |                    | `sm/test_wp5.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp5.in>`_,               |
   |                    | `sm/test_wp6.in                                                                              |
   |                    | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp6.in>`_                |
   +--------------------+----------------------------------------------------------------------------------------------+

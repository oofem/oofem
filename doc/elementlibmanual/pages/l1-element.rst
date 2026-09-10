
.. _L1:

L1 element
==========

Represents 1D isoparametric line element in 3D space. The element geometry is defined
using  2 nodes. The element supports linear and quadratic approximation.

.. table:: l1 element summary
   :name: l1summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | l1                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | l1 element                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Depending on configured variables                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear and quadratic unknown approximation; Linear approximation of geometry.                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | As required by individual terms                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | None, need to be defined using terms                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `mpm/mpms05.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms05.in>`_, |
   |                          | `mpm/mpms06.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms06.in>`_, |
   |                          | `mpm/mpms08.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms08.in>`_, |
   |                          | `mpm/mpms01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms01.in>`_, |
   |                          | `mpm/mpms02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms02.in>`_, |
   |                          | `mpm/mpms07.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms07.in>`_, |
   |                          | `mpm/cook2_u1p0.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u1p0.in>`_,            |
   |                          | `mpm/cook2_u2p1.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u2p1.in>`_,            |
   |                          | `mpm/cook2_u1p0_2.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/cook2_u1p0_2.in>`_,          |
   |                          | `mpm/mpms03.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms03.in>`_  |
   |                          | (and 2 more)                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+

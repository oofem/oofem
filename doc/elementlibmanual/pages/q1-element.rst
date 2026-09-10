
.. _Q1:

Q1 element
==========

Represents 2D isoparametric quad element. The element geometry in (x,y) plane is defined
using  4 corner nodes. The element supports linear and quadratic approximation.

.. table:: q1 element summary
   :name: q1summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | q1                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | q1 element                                                                                   |
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
   | Tests/Examples           | `mpm/t02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/t02.in>`_,       |
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
   |                          | `mpm/mpms03.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms03.in>`_, |
   |                          | `mpm/mpms_cook2_u2p1.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms_cook2_u2p1.in>`_ (and 1 |
   |                          | more)                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+

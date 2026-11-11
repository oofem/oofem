
.. _L1:

L1 element
==========

Represents 1D line element in the (x,y) plane. The element supports linear and quadratic
approximation. It is given either by its 2 end nodes, or by 3 nodes: the 2 end nodes followed by
the inner node of the quadratic approximation (``l1 n nodes 3 a b c``).

* Given by its end nodes, the element has a linear geometry; a quadratic approximation generates
  its inner node at the midpoint.
* Given with its inner node, the element has a quadratic (isoparametric) geometry through all 3
  nodes, and the inner node may lie anywhere between the end nodes. The gradients of all
  unknowns, linear ones included, and the integration follow this geometry. Being an input node,
  the inner node can be referenced by sets, boundary conditions and loads, which a generated node
  cannot.

See the Cell records section of the Input Manual for the node specification.

.. table:: l1 element summary
   :name: l1summary
   :widths: 26 94

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | l1                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | l1 element                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Depending on configured variables                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear and quadratic unknown approximation; linear approximation of geometry (2 nodes),      |
   |                          | quadratic approximation of geometry (3 nodes, inner node given).                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact; approximate for a quadratic geometry with the inner node off the midpoint, where the  |
   |                          | integrands are rational functions                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | \-                                                                                           |
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
   |                          | `mpm/mpms03.in <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms03.in>`_, |
   |                          | `mpm/mpms_l1inner.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/mpm/mpms_l1inner.in>`_           |
   |                          | (and 2 more)                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+

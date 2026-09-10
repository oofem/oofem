Free warping analysis elements
==============================


Free warping analysis elements
------------------------------

TrWarp
~~~~~~

Implements 2D linear triangular three-node finite element for :param:`FreeWarping`
analysis. Each node has 1 degree of freedom. The node numbering is anti-clockwise. The
element features are summarized in :numref:`TrWarpsummary`.

.. table:: TrWarp element summary
   :name: TrWarpsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | TrWarp                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular warping element                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | One dof (the value of deplanation function) is required in each node.                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using one point gauss integration formula.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | This type of element is supported in :param:`FreeWarping` analysis only.                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Only :param:`warpingCS` is supported.                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Edge loads corresponding to free warping problem are generated automatically. Additional     |
   |                          | edge loads are not supported.                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/freewarpingtest2.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/freewarpingtest2.in>`_        |
   +--------------------------+----------------------------------------------------------------------------------------------+

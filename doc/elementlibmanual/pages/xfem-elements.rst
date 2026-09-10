XFEM elements
=============


XFEM elements
-------------

XFEM elements allow simulations where the unknown fields are enriched through the
partition of unity concept. Two elements are currently available for 2D XFEM
simulations: *TrPlaneStress2dXFEM* (subclass of *TrPlaneStress2d*) and
*PlaneStress2dXfem* (subclass of PlaneStress2d).

TrPlaneStress2dXFEM element
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. table:: TrPlaneStress2dXFEM element summary
   :name: TrPlaneStress2dXFEMsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | TrPlaneStress2dXFEM                                                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Two dimensional 3-node triangular XFEM element                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`czmaterial{in}` :optelemparam:`nipcz{in}`:optelemparam:`useplanestrain{in}`   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`czmaterial`: Interface material for the cohesive zone. (If no material is specified, |
   |                          | traction free crack surfaces are assumed.) :param:`nipcz`: Number of integration points used |
   |                          | on each segment of the cohesive zone. :param:`useplanestrain`: If plane strain or plane      |
   |                          | stress should be assumed. 0 implies plane stress and 1 implies plane strain.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two continuous (standard) DOFs (u-displacement, v-displacement) and a variable number of     |
   |                          | enriched DOFs (can be continuous or discontinuous).                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Elements cut by an XFEM interface are divided into subtriangles.                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/xFemCrackVal.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xFemCrackVal.in>`_            |
   +--------------------------+----------------------------------------------------------------------------------------------+

PlaneStress2dXfem element
~~~~~~~~~~~~~~~~~~~~~~~~~

.. table:: PlaneStress2dXfem element summary
   :name: PlaneStress2dXfemsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | PlaneStress2dXfem                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Two dimensional 4-node quad XFEM element                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`czmaterial{in}` :optelemparam:`nipcz{in}`:optelemparam:`useplanestrain{in}`   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`czmaterial`: Interface material for the cohesive zone. (If no material is specified, |
   |                          | traction free crack surfaces are assumed.) :param:`nipcz`: Number of integration points used |
   |                          | on each segment of the cohesive zone. :param:`useplanestrain`: If plane strain or plane      |
   |                          | stress should be assumed. 0 implies plane stress and 1 implies plane strain.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two continuous (standard) DOFs (u-displacement, v-displacement) and a variable number of     |
   |                          | enriched DOFs (can be continuous or discontinuous).                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Elements cut by an XFEM interface are divided into subtriangles.                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/plasticRemap1.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/plasticRemap1.in>`_,          |
   |                          | `sm/xfemCohesiveZone1.in                                                                     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xfemCohesiveZone1.in>`_,      |
   |                          | `sm/xfemMultipleCracks1.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xfemMultipleCracks1.in>`_,    |
   |                          | `sm/xFemCrackValBranchCZ.in                                                                  |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xFemCrackValBranchCZ.in>`_,   |
   |                          | `sm/xFemIntersectingCracks.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xFemIntersectingCracks.in>`_, |
   |                          | `sm/xfemCrackPropMatForce.in                                                                 |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xfemCrackPropMatForce.in>`_,  |
   |                          | `sm/xFemCrackValBranch.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/xFemCrackValBranch.in>`_      |
   +--------------------------+----------------------------------------------------------------------------------------------+

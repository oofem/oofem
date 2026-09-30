Lattice transport elements
==========================


Lattice transport elements
~~~~~~~~~~~~~~~~~~~~~~~~~~

These are two-node lattice transport elements intended for coupled and stand-alone mass-
or heat-transport analyses on Delaunay/Voronoi lattice networks. Each node carries one
pressure (or temperature) degree of freedom. Flow along the element follows a linear
pressure interpolation; the cross-section polygon is the dual Voronoi facet and enters
strength/stiffness-like quantities through :param:`polycoords`.

The element assembles a capacity matrix that can optionally be switched from the default
consistent form to a diagonal (row-sum-preserving) lumped form. Lumping makes the scheme
equivalent to a two-point-flux-approximation (TPFA) finite-volume discretisation on the
Voronoi dual and is recommended when the capacity :math:`c(p)` is strongly nonlinear
(e.g. drying, Richards-type problems), where the consistent form is known to produce
pressure oscillations and an apparent loss of conductivity near sharp fronts. Lumping is
opt-in via :param:`lumpedcapacity`; the default preserves existing behaviour.

latticemt2d element
^^^^^^^^^^^^^^^^^^^

Two-node 2D lattice mass-transport element. The element features are summarised in Table
:numref:`latticemt2dsummary`.

.. table:: latticemt2d element summary
   :name: latticemt2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticemt2d                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2d lattice mass-transport element                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`thick{rn}` :elemparam:`width{rn}` :elemparam:`gpcoords{ra}`                      |
   |                          | :optelemparam:`dim{rn}` :optelemparam:`crackwidth{rn}` :optelemparam:`couplingflag{in}`      |
   |                          | :optelemparam:`couplingnumber{in}` :optelemparam:`lumpedcapacity{in}`                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`thick`: out-of-plane thickness. :param:`width`: Voronoi-facet width (cross-section   |
   |                          | length). :param:`gpcoords`: integration-point coordinates in the global system.              |
   |                          | :param:`dim`: dimension factor (optional, default 2). :param:`crackwidth`: crack width used  |
   |                          | in conductivity scaling (optional). :param:`couplingflag`: flag (optional, default 0)        |
   |                          | activating coupling with a mechanical lattice element. :param:`couplingnumber`: number of    |
   |                          | the coupled mechanical element. :param:`lumpedcapacity`: flag (optional, default 0). If set  |
   |                          | to 1, the capacity matrix is replaced by its diagonal (row-sum-preserving) lumped form,      |
   |                          | yielding a TPFA-equivalent, monotone scheme. Recommended for strongly nonlinear              |
   |                          | :math:`c(p)`.                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (P_f - pressure / moisture) per node.                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear pressure along the element.                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | One Gauss point at :param:`gpcoords`.                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | simplecs.                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice2dcrackinput.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2dcrackinput.in>`_,    |
   |                          | `lm/latticetransmat.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticetransmat.in>`_         |
   +--------------------------+----------------------------------------------------------------------------------------------+

latticemt3d element
^^^^^^^^^^^^^^^^^^^

Two-node 3D lattice mass-transport element. The cross-section is a polygonal Voronoi
facet supplied via :param:`polycoords`. The element features are summarised in Table
:numref:`latticemt3dsummary`.

.. table:: latticemt3d element summary
   :name: latticemt3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticemt3d                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice mass-transport element                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`polycoords{ra}` :optelemparam:`dim{rn}` :optelemparam:`area{rn}`                 |
   |                          | :optelemparam:`mlength{rn}` :optelemparam:`crackwidths{ra}` :optelemparam:`couplingflag{in}` |
   |                          | :optelemparam:`couplingnumber{ia}` :optelemparam:`lumpedcapacity{in}`                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`polycoords`: coordinates of the mid-cross-section (Voronoi facet) vertices in the    |
   |                          | global system. :param:`dim`: dimension factor (optional, default 3). :param:`area`:          |
   |                          | prescribed cross-section area (optional; otherwise computed from :param:`polycoords`).       |
   |                          | :param:`mlength`: minimum element length threshold (optional). :param:`crackwidths`: per-    |
   |                          | vertex crack widths used in conductivity scaling (optional). :param:`couplingflag`: flag     |
   |                          | (optional, default 0) activating coupling with a mechanical lattice element.                 |
   |                          | :param:`couplingnumber`: array of coupled mechanical element numbers.                        |
   |                          | :param:`lumpedcapacity`: flag (optional, default 0). If set to 1, the capacity matrix is     |
   |                          | replaced by its diagonal (row-sum-preserving) lumped form, yielding a TPFA-equivalent,       |
   |                          | monotone scheme. Recommended for strongly nonlinear :math:`c(p)`.                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (P_f - pressure / moisture) per node.                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear pressure along the element.                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | One Gauss point at the element midpoint.                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | simplecs.                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice3d_mt2.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt2.in>`_,          |
   |                          | `lm/lattice3d_mt1.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt1.in>`_           |
   +--------------------------+----------------------------------------------------------------------------------------------+

latticemt3dboundary element
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three-noded 3D lattice mass-transport element for boundaries of 3D periodic cells. The
first two nodes carry the single pressure / moisture dof as in latticemt3d; the third is
a control node that supplies the macro gradients driving the periodic cell. The dof of
the node that lies outside the periodic cell is reconstructed from its periodic image
inside the cell plus the control-node dofs, using the same translation pattern as
lattice3dboundary (:numref:`lattice3dboundaryfig`). The element features are
summarised in :numref:`latticemt3dboundarysummary`.

.. table:: latticemt3dboundary element summary
   :name: latticemt3dboundarysummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticemt3dboundary                                                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice mass-transport boundary element                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`polycoords{ra}` :elemparam:`location{ia}` :optelemparam:`dim{rn}`                |
   |                          | :optelemparam:`area{rn}` :optelemparam:`mlength{rn}` :optelemparam:`crackwidths{ra}`         |
   |                          | :optelemparam:`couplingflag{in}` :optelemparam:`couplingnumber{ia}`                          |
   |                          | :optelemparam:`lumpedcapacity{in}`                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`polycoords`: coordinates of the mid-cross-section (Voronoi facet) vertices in the    |
   |                          | global system. :param:`location`: array of two integers between 1 and 26 specifying the      |
   |                          | location of the two transport nodes with respect to the 3D periodic cell — same encoding as  |
   |                          | lattice3dboundary. :param:`dim`: dimension factor (optional, default 3). :param:`area`:      |
   |                          | prescribed cross-section area (optional; otherwise computed from :param:`polycoords`).       |
   |                          | :param:`mlength`: minimum element length threshold (optional). :param:`crackwidths`: per-    |
   |                          | vertex crack widths used in conductivity scaling (optional). :param:`couplingflag`: flag     |
   |                          | (optional, default 0) activating coupling with a mechanical lattice element.                 |
   |                          | :param:`couplingnumber`: array of coupled mechanical element numbers.                        |
   |                          | :param:`lumpedcapacity`: flag (optional, default 0). Replaces the capacity matrix with its   |
   |                          | diagonal lumped form for monotone TPFA-equivalent behaviour.                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (P_f) per transport node. The third (control) node carries the periodic-cell dofs |
   |                          | as for lattice3dboundary but reduced to the transport-relevant components.                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear pressure along the element.                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | One Gauss point at the element midpoint.                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | simplecs.                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice3d_mt2.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d_mt2.in>`_           |
   +--------------------------+----------------------------------------------------------------------------------------------+

References
----------

.. [RobertCook1989] R. D. Cook and D. S. Malkus and M. E. Plesha, “Concepts and Applications of Finite Element Analysis”, Third Edition, isbn: 0-471-84788-7, 1989.

.. [BittnarSejnoha1996] Z. Bittnar and J. Sejnoha, “Numerical Methods in Structural Mechanics”,Thomas Telford,isbn:978-0784401705, 1996.

.. [RagnarLarsson2011] R. Larsson and J. Mediavilla and M. Fagerström, “Dynamic fracture modeling in shell structures based on XFEM”, International Journal for Numerical Methods in Engineering, vol. 86, no. 4-5, 499--527, 2011.

.. [GraJir10] P. Grassl and M. Jirásek, “Meso-scale approach to modelling the fracture process zone of concrete subjected to uniaxial tension”, International Journal of Solids and Structures, vol. 47, iss. 7-8, pp. 957-968, 2010..

.. [GraBol16] P. Grassl, J. Bolander, “Three-Dimensional Network Model for Coupling of Fracture and Mass Transport in Quasi-Brittle Geomaterials”, Materials, 9, 782, 2016

.. [AthWheGra18] I. Athanasiadis, S. Wheeler and P. Grassl. “Hydro-mechanical network modelling of particulate composites”, International Journal of Solids and Structures, vol. 130-131, pp. 49-60, 2018.

.. [GraAnt19] P. Grassl and A. Antonelli. “3D network modelling of fracture processes in fibre-reinforced geomaterials”, International Journal of Solids and Structures, vol. 156-157, Pages 234-242, 2019.

.. [AbdGra24] G. Abdelrhim and P. Grassl. “A simple frame element for large rotations”, Available at SSRN 4850146, 2024.

.. [SciGraLarRun20] A. Sciegaj and P. Grassl and F. Larsson and K. Runesson and K. Lundgren. “Upscaling of three-dimensional reinforced concrete representative volume elements to effective beam and plate models”, International Journal of Solids and Structures, vol. 202, pp. 835-853, 2020.

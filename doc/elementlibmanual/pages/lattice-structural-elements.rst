Lattice structural elements
===========================


Lattice structural elements
---------------------------

Lattice2d element
~~~~~~~~~~~~~~~~~

Represents two-node lattice element. Each node has 3 degrees of freedom. The element is
based on the Rigid Body Spring Model originally developed by Kawai and later delveloped
by Bolander for modelling fracture in concrete. The main idea is to model the elastic
and inelastic response of a connection of two nodes by a set of springs located at the
contact facet of two rigid bodies, which is the mid-cross-section of the element.
Displacement jumps are computed at the mid-cross-section, which are smeared out over the
element length in the form of strains. The element is defined in x,y plane (see Figure
:ref:`lattice2dfig`). The element features are summarized in Table
:numref:`lattice2dsummary`.

.. _lattice2dfig:

.. tikz:: Lattice2d element. Node numbering, DOF numbering and definition of integration point :math:`C`.

   \tikzstyle{elemnode} = [draw,thin,circle,inner sep=1,fill=white]
   \tikzstyle{dofstyle} = [red]
   \tikzstyle{nodestyle} = [black]

   \begin{tikzpicture}[scale=7,>=stealth]
    \coordinate (b) at (0.15,0.1);
    \newcommand{\beamlength}{0.6};

    \draw[->] (-0.05,0) -- (0.6,0) node[below,at end] {$x_g$};
    \draw[->] (0,-0.05) -- (0,0.5) node[right,at end] {$y_g$};

    \draw[very thick] (b) -- +(30:\beamlength)
       coordinate[midway] (bmid)
       coordinate[at end] (bend);

    \draw[dotted,<-] (bmid)++(-60:0.05) node[below left] {$y_{\mathrm{l}}$} -- +(-60:0.1) coordinate (bmid2);
    \draw[dotted,->] (bmid2) -- +(30:0.1) node[below right] {$x_{\mathrm{l}}$};
    \node[elemnode] at (bmid2) {};
    \node[below] at (bmid2) {$C$};

    \draw[thin,->] (b) -- +(0:0.1) node[right, dofstyle] {1};
    \draw[thin,->] (b) -- +(90:0.1) node[above, dofstyle] {2};
    \draw[thin,->] (bend) -- +(0:0.1) node[right, dofstyle] {4};
    \draw[thin,->] (bend) -- +(90:0.1) node[above, dofstyle] {5};
    \draw[thin,->] (b)++(-60:0.05) arc (-60:150:0.05) node[left, dofstyle] {3};
    \draw[thin,->] (bend)++(-60:0.05) arc (-60:150:0.05) node[left, dofstyle] {6};

    \node[elemnode] at (b) {}; \node[below, nodestyle] at (b) {1};
    \node[elemnode] at (bend) {}; \node[below, nodestyle] at (bend) {2};
   \end{tikzpicture}

.. table:: lattice2d element summary
   :name: lattice2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lattice2d                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2d lattice element                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`thick{rn}` :elemparam:`width{rn}` :elemparam:`gpCoords{ra}`                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`thick`: defines the out of plane (:math:`z`-direction) thickness :param:`width`:     |
   |                          | defines the width of the midpoint cross-section in the :math:`x`-:math:`y` plane with the    |
   |                          | point :math:`C` at its centre :param:`gpCooords`: array of the coordinates of the            |
   |                          | integration point :math:`C` in the global coordinate system                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-rotation) are required |
   |                          | in each node.                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice2drandom.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2drandom.in>`_,        |
   |                          | `lm/latticeplastdamvisco_mps_1.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                          | n/lm/latticeplastdamvisco_mps_1.in>`_, `lm/latticeplastdamvisco_mps_2.in <https://github.com |
   |                          | /oofem/oofem/blob/devel/tests/regression/lm/latticeplastdamvisco_mps_2.in>`_,                |
   |                          | `lm/latticeplastdamvisco_mps_3.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                          | n/lm/latticeplastdamvisco_mps_3.in>`_, `lm/latticedamagevisco_mps_1.in <https://github.com/o |
   |                          | ofem/oofem/blob/devel/tests/regression/lm/latticedamagevisco_mps_1.in>`_,                    |
   |                          | `lm/latticeplastdamvisco3drandom.in <https://github.com/oofem/oofem/blob/devel/tests/regress |
   |                          | ion/lm/latticeplastdamvisco3drandom.in>`_, `lm/latticeviscoelastic_mps_1.in <https://github. |
   |                          | com/oofem/oofem/blob/devel/tests/regression/lm/latticeviscoelastic_mps_1.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

The theory of lattice2d is described in the paper “P. Grassl and M. Jirásek. Meso-scale
approach to modelling the fracture process zone of concrete subjected to uniaxial
tension. International Journal of Solids and Structures. Volume 47, Issues 7-8, pp.
957-968, 2010.”

latticeboundary2D element
~~~~~~~~~~~~~~~~~~~~~~~~~

Represents three-node lattice element for boundary of 2d periodic cells. The first two
nodes have 3 degrees of freedom as for the element lattice2d. The third node is used to
control the loading of the periodic cell. It has three components which are
displacements, which are produces of the macroscopic (average) strain components and
length of the peridodic cell as :math:`aE_{xx}`, :math:`bE_{xx}` and :math:`bG_{xy}` and
the length of the periodic cell. The DOFs of the node that lies outside the periodic are
computed from those of the periodic image inside the cell and the DOFs at the third node
(:numref:`lattice2dboundaryfig`). The coordinates :math:`x` and :math:`y` of the
third node are the lengths :math:`a` and :math:`b` of the periodic cell, respectively.
The element is defined in x,y plane. The strain components at the additional node have
the meaning of average strains in the periodic cell. The specific input parameters for
this element in addition to those used for lattice2d are shown in Table
``lattice2dboundarysummary``.

.. figure:: /figures/lattice2dboundaryfig.svg
   :width: 40%
   :name: lattice2dboundaryfig

   latticeboundary2D element. Elements with cross boundary of periodic cell use DOFs of node inside cell and average strain values of periodic cell.

.. table:: latticeboundary2D element summary
   :name: latticeboundary2D summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticeboundary2D                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2d lattice boundary element                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`location{in}`                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`location`: number between 1 and 8 which specifies the location of the node with      |
   |                          | respect to the periodic cell.                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Nine DOFs are required, which are :math:`u`-displacement, :math:`v`-displacement and         |
   |                          | :math:`w`-rotation at nodes 1 and 2, and :math:`aE_{xx}`, :math:`bE_{yy}` and                |
   |                          | :math:`bG_{xy}` at node 3.                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [GraJir10]_                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice2dboundary1.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice2dboundary1.in>`_      |
   +--------------------------+----------------------------------------------------------------------------------------------+

Lattice3d element
~~~~~~~~~~~~~~~~~

Lattice3d represents a two-node 3d lattice element. Each node has six degrees of freedom
as shown in :numref:`lattice3dfig`. The element is based on the Rigid Body Spring
Model originally developed by Kawai and later delveloped by Bolander for modelling
fracture in concrete. The main idea is to model the elastic and inelastic response of a
connection of two nodes by a set of springs located at the contact facet of two rigid
bodies, which is the mid-cross-section of the element. The properties of the
mid-cross-section are internally computed from its vertices which are given as input in
the global coordinate sytem. Displacement jumps are computed at the mid-cross-section,
which are smeared out over the element length in the form of strains. The input
parameters for this element are shown in :numref:`lattice3dsummary`.

.. figure:: /figures/lattice3d.svg
   :width: 30%
   :name: lattice3dfig

   Lattice3d element. Node numbering, DOF numbering, cross-section vertices and local coordinate system at integration point :math:`C`.

.. table:: lattice3d element summary
   :name: lattice3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lattice3d                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice element                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`polycoords{ra}` :elemparam:`couplingflag{in}` :elemparam:`couplingnumbers{ra}`   |
   |                          | :elemparam:`pressures{ra}` :elemparam:`mlength{rn}`                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`polycooords`: array of the coordinates of the vertices of the mid-cross-section of   |
   |                          | the lattice element in the global coordinate system. :param:`couplingflag`: flag (optional   |
   |                          | parameter. Default is 0) which activates coupling with a transport lattice element.          |
   |                          | :param:`couplingnumbers`: array of numbers of transport lattice elements (optional           |
   |                          | parameter), which are coupled with the 3d lattice element. :param:`pressures`: array of      |
   |                          | pressure values (optional parameter), which are used to consider influence of fluid pressure |
   |                          | on mechanical response. :param:`mlength`: minimum length (optional parameter) is used to     |
   |                          | check if the cross- section of the element is not too small. Default value is 1.e-20.        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in each node.    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [GraBol16]_                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice3dshell.in                                                                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell.in>`_,         |
   |                          | `lm/lattice3dshell_lle.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dshell_lle.in>`_,     |
   |                          | `lm/lattice3dshell_lle_twist_n4.in <https://github.com/oofem/oofem/blob/devel/tests/regressi |
   |                          | on/lm/lattice3dshell_lle_twist_n4.in>`_, `lm/lattice3d3.in                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d3.in>`_,             |
   |                          | `lm/lattice3d4.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d4.in>`_,             |
   |                          | `lm/lattice3drandom.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3drandom.in>`_,        |
   |                          | `lm/lattice3ddamplast3.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3ddamplast3.in>`_,     |
   |                          | `lm/lattice3dshell_dmg_bend_n7.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                          | n/lm/lattice3dshell_dmg_bend_n7.in>`_, `lm/lattice3dbondplast1.in                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dbondplast1.in>`_,    |
   |                          | `lm/lattice3delastic.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3delastic.in>`_ (and   |
   |                          | 17 more)                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Lattice3dBoundary element
~~~~~~~~~~~~~~~~~~~~~~~~~

This element represents a three-noded 3d lattice element for boundaries of 3d periodic
cells. The first two nodes have 6 degrees of freedom as for the element lattice3d. The
third node is used to control the loading of the periodic cell using three normal
(:math:`aE_{xx}`, :math:`bE_{yy}` and :math:`zE_{zz}`) and three shear strain
(:math:`cG_{yz}`, :math:`cG_{xz}`, :math:`bG_{xy}`) components. Here, :math:`a`,
:math:`b` and :math:`c` are the three dimensions of the periodic cell. The DOFs of the
node that lies outside the periodic are computed from those of the periodic image inside
the cell and the DOFs at the third node (:numref:`lattice3dboundaryfig`). The
connection between periodic nodes is defined as

.. math::

   \mathbf{x}^{'} = \mathbf{M} \mathbf{x}

Here, :math:`\mathbf{x}^{'}` and :math:`\mathbf{x}` are the nodes inside and outside,
respectively, and :math:`\mathbf{M}` is the translation matrix, for which the input is
provided in the form of a location parameter as shown in Table
:numref:`lattice3dboundarysummary`.

.. figure:: /figures/lattice3dboundaryfig.svg
   :width: 60%
   :name: lattice3dboundaryfig

   Lattice3dboundary element. Elements with cross boundary of periodic cell use DOFs of node inside cell and average strain values of periodic cell.

.. table:: lattice3dboundary element summary
   :name: lattice3dboundarysummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lattice3dboundary                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice boundary element                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`location{in}`                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`location`: array of two numbers between number between 1 and 26, which specifies the |
   |                          | location of the two nodes with respect to the 3d periodic cell.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in each of the   |
   |                          | first two node. Node 3 requires the 6 quantities to control the periodic cell                |
   |                          | :math:`aE_{xx}`, :math:`bE_{yy}`, :math:`zE_{zz}`, :math:`cG_{yz}`, :math:`cG_{xz}` and      |
   |                          | :math:`bG_{xy}`                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [AthWheGra18]_                                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice3d1.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3d1.in>`_              |
   +--------------------------+----------------------------------------------------------------------------------------------+

Latticelink3d
~~~~~~~~~~~~~

This element represents a two-node 3d link element connecting 3d beam and 3d lattice
elements. Each node has six degrees of freedom. The node ordering is important: the
first node carries the rigid arm and must therefore have rotational degrees of freedom,
while the slip (displacement jump) is evaluated at the second node. Which physical
entity is placed first is the user's modelling choice; it only has to be the node from
which the rigid arm is offset. For a lattice model both matrix and reinforcement nodes
usually have rotational degrees of freedom, so the matrix node is typically given first
and the reinforcement node second, evaluating the slip at the reinforcement.  The input
parameters for this element are shown in :numref:`latticelink3dsummary`.

.. table:: latticelink3d element summary
   :name: latticelink3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticelink3d                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice link element                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`length{rn}` :elemparam:`diameter{rn}` :elemparam:`dirvector{ra}`                 |
   |                          | :elemparam:`l_end{rn}`                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`length`: bond length :param:`diameter`: diameter :param:`dirvector`: direction       |
   |                          | vector in which bond-slip occurs. :param:`l_end`: array of pressure values (optional         |
   |                          | parameter), which are used to consider influence of fluid pressure on mechanical response.   |
   |                          | :param:`mlength`: minimum length (optional parameter) is used to check if the cross- section |
   |                          | of the element is not too small. Default value is 1.e-20.                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in each node.    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [GraAnt19]_                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/latticelink1.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticelink1.in>`_            |
   +--------------------------+----------------------------------------------------------------------------------------------+

Latticelink3dboundary
~~~~~~~~~~~~~~~~~~~~~

Represents three-node 3d boundary link element connecting 3d beam and 3d lattice
elements. The first two nodes have the same meaning as for latticelink3d. The third node
is used to control the loading of the periodic cell using three normal (xx, yy and zz)
and three shear strain (yz, xz, xy) components. The specific input parameters for this
element in addition of those for latticelink3d are shown in Table
:numref:`latticelink3dboundarysummary`.

.. table:: latticelink3dboundary element summary
   :name: latticelink3dboundarysummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | latticelink3dboundary                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d lattice link boundary element                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`location{in}` rn                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`location`:                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in each node.    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [GraAnt19]_                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+

Lattice3dnl element
~~~~~~~~~~~~~~~~~~~

Lattice3dnl is the geometrically nonlinear (large-rotation) counterpart of the lattice3d
element. It accepts the same input and the same cross-section options (facet/shell via
:param:`polycoords`, or beam via the cross-section properties together with
:param:`zaxis` and :param:`s`), but the displacement jumps and nodal forces are
evaluated in the deformed configuration, so the element is suitable for large rotations
and catenary action. The rotational part of the generalised strain is the exact relative
rotation of the two nodes, i.e. the axial vector of :math:`\mathbf{R}_2
\mathbf{R}_1^{\mathrm{T}}`, which reduces to the difference of the nodal rotation
vectors for rotation about a single axis.

.. table:: lattice3dnl element summary
   :name: lattice3dnlsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lattice3dnl                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d geometrically nonlinear lattice element                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`polycoords{ra}` :elemparam:`shellnormal{ra}` :optelemparam:`zaxis{ra}`           |
   |                          | :optelemparam:`s{rn}`                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`polycoords`, :param:`shellnormal`, :param:`zaxis`, :param:`s`: as for the lattice3d  |
   |                          | element.                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Cross-section            | As for the lattice3d element.                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in each node.    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [AbdGra24]_                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `lm/lattice3dnlshell.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dnlshell.in>`_,       |
   |                          | `lm/lattice3dnl_logrot.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/lattice3dnl_logrot.in>`_,     |
   |                          | `lm/lattice3dnlshell_dmg_bend_n7.in <https://github.com/oofem/oofem/blob/devel/tests/regress |
   |                          | ion/lm/lattice3dnlshell_dmg_bend_n7.in>`_, `lm/latticelink3dnl_assembly.in <https://github.c |
   |                          | om/oofem/oofem/blob/devel/tests/regression/lm/latticelink3dnl_assembly.in>`_,                |
   |                          | `lm/latticeframe3dnl.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/lm/latticeframe3dnl.in>`_        |
   +--------------------------+----------------------------------------------------------------------------------------------+

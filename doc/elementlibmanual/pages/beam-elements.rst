Beam Elements
=============


Beam Elements
-------------

Beam2d element
~~~~~~~~~~~~~~

Beam element for 2D analysis, based on Timoshenko hypothesis. Structure should be
defined in x,z plane. The internal condensation of arbitrary DOF is supported and is
performed in local coordinate system. On output, the local end displacement and local
end forces are printed. The element features are summarized in Table
:numref:`beam2dsummary`.

.. _beam2dfig:

.. tikz:: Beam2d element. Definition of local c.s.(a) and definition of
 local end forces and local element dofs (b).

   \begin{tikzpicture}[scale=7,>=stealth]
    \tikzstyle{elemnode} = [draw,thin,circle,inner sep=1,fill=white]
    \coordinate (a) at (0.1,-0.2);
    \coordinate (b) at (0.6,-0.1);
    \newcommand{\beamlength}{0.4};

    \draw[->] (-0.05,0) -- (1.1,0) node[below,at end] {$x_g$};
    \draw[->] (0,0.05) -- (0,-0.5) node[right,at end] {$z_g$};
    \draw[very thick] (a) --  +(-30:0.4)
       node[at start,elemnode] {} node[at start,above right] {1}
       node[at end,elemnode] {} node[at end,above right] {2}
       node[midway,below left,inner sep=2] {(a)};
    \draw[dotted,->] (a)++(-30:0.42) -- +(-30:0.1) node[below] {$X_1$};
    \draw[dotted,->] (a)++(-30-90:0.02) -- +(-120:0.1) node[right] {$Z_1$};

    \draw[very thick] (b) -- +(-30:0.4)
       node[midway,below left, inner sep=2] {(b)}
       coordinate[at end] (bend);
    \draw[dotted,->] (b)++(-30-90:0.1) -- +(-120:0.1);
    \draw[dotted,->] (bend)++(-30:0.1) -- +(-30:0.1);

    \draw[thin,<-] (b) -- +(-30:-0.1) node[below left,midway] {1};
    \draw[thin,->] (b) -- +(-120:0.1) node[below right,at end] {2};
    \draw[thin,->] (bend) -- +(-30:0.1) node[below left,at end] {4};
    \draw[thin,->] (bend) -- +(-120:0.1) node[below right,at end] {5};
    \draw[thin,->] (b)++(-20:0.05) arc (-30:120:0.05); \node at (0.68,-0.05) {3};
    \draw[thin,->] (bend)++(-20:0.05) arc (-30:120:0.05); \node at (1.03,-0.27) {6};
   \end{tikzpicture}

.. table:: beam2d element summary
   :name: beam2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | beam2d                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D beam element                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`dofstocondense{ia}`                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`dofstocondense`: allows to specify local element dofs that will be condensed. The    |
   |                          | numbering of local element dofs is shown in fig. :ref:`beam2dfig`. The size of this array    |
   |                          | should be equal to number of local element dofs (6) and nonzero value indicates the          |
   |                          | corresponding dof will be condensed.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, w-displacement, y-rotation) are required in each node.           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Cubic approximations of lateral displacement and rotation are used. For longitudinal         |
   |                          | displacement the linear one is assumed.                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Full dynamic analysis support. Linear stability analysis support.                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Area,inertia moment along y-axis (:param:`iy` parameter) and equivalent shear area           |
   |                          | (:param:`shearareaz` parameter) should be specified.                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant and linear edge loads are supported, shear influence is taken into account. Edge    |
   |                          | number should be equal to 1. Temperature load is supported, the first coefficient of         |
   |                          | temperature load represent mid-plane temperature change, the second one represent difference |
   |                          | between temperature change of local z+ and local z- surfaces of beam (in local coordinate    |
   |                          | system). Temperature load require that the “thick” property of cross section model is        |
   |                          | defined.                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/pdelta01.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/pdelta01.in>`_,               |
   |                          | `sm/pdelta02.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/pdelta02.in>`_,               |
   |                          | `sm/rigarm04.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rigarm04.in>`_,               |
   |                          | `sm/beam2d_5.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam2d_5.in>`_,               |
   |                          | `sm/nlstatic01.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/nlstatic01.in>`_,             |
   |                          | `sm/spring05.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring05.in>`_,               |
   |                          | `sm/rigarm03.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rigarm03.in>`_,               |
   |                          | `sm/layeredcs03.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layeredcs03.in>`_,            |
   |                          | `sm/rigarm01.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rigarm01.in>`_,               |
   |                          | `sm/rigarm02.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rigarm02.in>`_ (and 7 more)   |
   +--------------------------+----------------------------------------------------------------------------------------------+

Beam3d element
~~~~~~~~~~~~~~

Beam element for 3D **linear** analysis, based on Timoshenko hypothesis. The internal
condensation of arbitrary DOF is supported and is performed in local coordinate system.
On output, the local end-displacement and local end-forces are printed. Requires the
local coordinate system to be chosen according to main central axes of inertia. Local
element  coordinate system is determined by the following rules:

1. let first element node has following coordinates :math:`(x_i, y_i, z_i)` and
   the second one :math:`(x_j, y_j, z_j)`,

2. direction vector of local x-axis is then :math:`\mathbf{a}_1 = (x_j-x_i,
   y_j-y_i, z_j-z_i)`,

3. local y-axis direction vector lies in plane defined by local x-axis
   direction vector (:math:`\mathbf{a}_1`) and given point (k-node with
   coordinates :math:`(x_k, y_k, z_k)`) - so called reference node,

4. local z-axis is then determined as vector product of local x-axis direction
   vector (:math:`\mathbf{a}_1`) by vector :math:`(x_k-x_i, y_k-y_i, z_k-z_i)`,

5. local y-axis is then determined as vector product of local z-axis direction
   vector by local x-axis direction vector.

The element features are summarized in :numref:`beam3dsummary`.

.. _beam3dfig:

.. tikz:: Beam3d element. Definition of local c.s., local end forces
 and local element dofs numbering.

   \begin{tikzpicture}[scale=6,>=stealth,
     x={(0.8cm,0cm)}, y={(-0.5cm,-0.5cm)}, z={(0cm,-1cm)}]
    \tikzstyle{elemnode} = [draw,thin,circle,inner sep=1,fill=white]
    %\coordinate (a) at (0.1,-0.2);
    %\coordinate (b) at (0.6,-0.1);
    \newcommand{\beamlength}{0.4};

    \begin{scope}[yshift=-3]
    \draw[->] (-0.05,0,0) -- (0.5,0,0) node[at end, below] {$x_g$};
    \draw[->] (0,-0.05,0) -- (0,0.5,0) node[at end, below] {$y_g$};
    \draw[->] (0,0,-0.05) -- (0,0,0.25) node[at end, right] {$z_g$};
    \end{scope}

    \coordinate (i) (0,0,0.5);
    \coordinate (j) (0.7,0,0.25);

    \draw[very thick,-] (0,0,0.5) coordinate (i) -- (0.7,0,0.25) coordinate(j)
      node[elemnode,at start] {} node[at start,above left] {1,i}
      node[elemnode,at end] {} node[at end,above left] {2,j};
    \draw[dashed,->] (0,0,0.5) -- (0.4,0,0.5)
      node[elemnode,at end] {} node[at end,right] {k};
    \draw[->] (j) -- +(0.7*0.4,0,-0.25*0.4) node[right] {$x_1$};
    \draw[->] (i) -- +(0.3,0.0,0.1) node[above] {$y_1$};
    \draw[->] (i) -- +(0,0,0.3) node[right] {$z_1$};

    \begin{scope}[xshift=20,yshift=0]
     \draw[very thick,-] (0,0,0.5) coordinate (i) -- (0.7,0,0.25) coordinate(j)
      node[elemnode,at start] {}
      node[elemnode,at end] {};
     \draw[<-] (i)++(-0.7*0.02,0,0.25*0.02) -- +(-0.7*0.18,0,0.25*0.18) node[above,midway] {1};
     \draw[<<-] (i)++(-0.7*0.21,0,0.25*0.21) -- +(-0.7*0.2,0,0.25*0.2) node[above,midway] {4};
     \draw[->] (i)+(0.015,0.0,0.005) -- +(0.15,0.0,0.05) node[above,near end] {2};
     \draw[->>] (i)++(0.16,0.0,0.0533) -- +(0.15,0.0,0.05) node[above,near end] {5};
     \draw[->] (i)+(0,0,0.02) -- ++(0,0,0.15) node[left,midway] {3};
     \draw[->>] (i)++(0,0,0.16) -- +(0,0,0.16) node[left,midway] {6};

     \draw[->] (j)++(0.7*0.02,0,-0.25*0.02) -- +(0.7*0.2,0,-0.25*0.2) node[above,midway] {7};
     \draw[->>] (j)++(0.7*0.23,0,-0.25*0.23) -- +(0.7*0.18,0,-0.25*0.18) node[above,midway] {4};
     \draw[->] (j)+(0.015,0.0,0.005) -- +(0.15,0.0,0.05) node[above,near end] {8};
     \draw[->>] (j)++(0.16,0.0,0.0533) -- +(0.15,0.0,0.05) node[above,near end] {11};
     \draw[->] (j)+(0,0,0.02) -- ++(0,0,0.15) node[left,midway] {9};
     \draw[->>] (j)++(0,0,0.16) -- +(0,0,0.16) node[left,midway] {12};

    \end{scope}

   \end{tikzpicture}

.. table:: beam3d element summary
   :name: beam3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | beam3d                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D beam element                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`refnode{in}` :optelemparam:`dofstocondense{ia}`                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`refnode`: sets reference node to determine the local coordinate system of element.   |
   |                          | :param:`dofstocondense`: allows to specify local element dofs that will be condensed. The    |
   |                          | numbering of local element dofs is shown in fig. :ref:`beam3dfig`. The size of this array    |
   |                          | should be equal to number of local element dofs (12) and nonzero value indicates the         |
   |                          | corresponding dof will be condensed.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and x,y,z-rotations) are required in each node.                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Cubic approximations of lateral displacement and rotation (along local y,z axes) are used.   |
   |                          | For longitudinal displacement and the rotation along local x-axis (torsion) the linear       |
   |                          | approximations are assumed.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Full dynamic analysis support. Linear stability analysis support.                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Area, inertia moment along y and z axis (:param:`iy` and :param:`iz` parameters), torsion    |
   |                          | inertia moment (:param:`ik` parameter) and either cross section area shear correction factor |
   |                          | (:param:`beamshearcoeff` parameter) or equivalent shear areas (:param:`shearareay` and       |
   |                          | :param:`shearareaz` parameters) are required. These cross section properties are assumed to  |
   |                          | be defined in local coordinate system of element.                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant and linear edge loads are supported. Edge number should be equal to 1. Temperature  |
   |                          | load is supported, the first coefficient of temperature load represent mid- plane            |
   |                          | temperature change, the second one represent difference between temperature change of local  |
   |                          | z+ surface and local z- surface surface of beam and the third one represent difference       |
   |                          | between temperature change of local y+ surface and local y- surface of beam. Requires the    |
   |                          | “thick” (measured in direction of local z axis) and “width” (measured in direction of local  |
   |                          | y axis) cross section model properties to be defined.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/eigen02_beam3d.in                                                                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/eigen02_beam3d.in>`_,         |
   |                          | `sm/beam3d_3.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam3d_3.in>`_,               |
   |                          | `sm/beam3dsubsoil01.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam3dsubsoil01.in>`_,        |
   |                          | `sm/fiberedcs01.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/fiberedcs01.in>`_,            |
   |                          | `sm/pdelta03.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/pdelta03.in>`_,               |
   |                          | `sm/beam3d_2.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam3d_2.in>`_,               |
   |                          | `sm/spring06.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring06.in>`_,               |
   |                          | `sm/beam3d_1.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam3d_1.in>`_,               |
   |                          | `sm/Buckling02.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/Buckling02.in>`_,             |
   |                          | `sm/eigen_beam3d.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/eigen_beam3d.in>`_            |
   +--------------------------+----------------------------------------------------------------------------------------------+

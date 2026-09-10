Plane Strain Elements
=====================


Plane Strain Elements
---------------------

Quad1PlaneStrain
~~~~~~~~~~~~~~~~

Represents isoparametric four-node quadrilateral plane-strain finite element. Each node
has 2 degrees of freedom. Structure should be defined in x,y plane.  The nodes should be
numbered anti-clockwise (positive rotation around z-axis). The element features are
summarized in :numref:`quad1planestrainsummary`.

.. _Quad1PlaneStrainfig:

.. tikz:: Quad1PlaneStrain element. Node numbering, Side numbering and
 definition of local edge c.s.(a).

   \begin{tikzpicture}[scale=6,>=stealth]
    \tikzstyle{elemnode} = [draw,circle,inner sep=1,fill=white]
    \newcommand{\lcoordsys}[2]{
      \begin{scope}[transform canvas={shift={#2},scale=0.5,rotate=#1}]
       \draw[->] (0,0.05) ++(-0.02,0) -- ++(0.2,0) node[above] {$x$};
       \draw[->] (0,0.05) ++(0,-0.02) -- ++(0,0.1) node[right] {$y$};
      \end{scope}
    }

    \draw[->] (-0.05,0) -- (0.8,0) node[above] {$x_g$};
    \draw[->] (0,-0.05) -- (0,0.5) node[right] {$y_g$};

    \draw[thick,xshift=-2]
        (0.2,0.1)
     -- (0.7,0.15) coordinate[midway] (e1) node[below,midway,blue] {1} node[elemnode] {} node[below right] {2}
     -- (0.8,0.45) coordinate[midway] (e2) node[right,midway,blue] {2} node[elemnode] {} node[above] {3}
     -- (0.3,0.5)  coordinate[midway] (e3) node[above,midway,blue] {3} node[elemnode] {} node[above] {4}
     -- (0.2,0.1)  coordinate[midway] (e4) node[left,midway,blue] {4} node[elemnode] {} node[below] {1};

    \lcoordsys{  7}{(e1)};
    \lcoordsys{ 70}{(e2)};
    \lcoordsys{175}{(e3)};
    \lcoordsys{255}{(e4)};
   \end{tikzpicture}

.. table:: quad1planestrain element summary
   :name: quad1planestrainsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | quad1planestrain                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear quadrilateral plane-strain element                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using gauss integration formula in 4 (the default), 9   |
   |                          | or 16 integration points. The default number of integration points used can be overloaded    |
   |                          | using :param:`NIP` parameter. Reduced integration for shear terms is employed. Shear terms   |
   |                          | are always integrated using 1 point integration rule.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Nonlocal constitutive support, Adaptivity support.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 4). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`Quad1PlaneStrainfig`)).                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/DruckerPrager_01.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/DruckerPrager_01.in>`_,       |
   |                          | `sm/ogden1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ogden1.in>`_,   |
   |                          | `sm/ogden2.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/ogden2.in>`_,   |
   |                          | `sm/mooneyrivlin1.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/mooneyrivlin1.in>`_,          |
   |                          | `sm/timestepreduction.in                                                                     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/timestepreduction.in>`_,      |
   |                          | `sm/mooneyrivlin2.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/mooneyrivlin2.in>`_,          |
   |                          | `sm/contact2d_sliding.in                                                                     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d_sliding.in>`_,      |
   |                          | `sm/contact2d_generalized_vertex_direct.in <https://github.com/oofem/oofem/blob/devel/tests/ |
   |                          | regression/sm/contact2d_generalized_vertex_direct.in>`_, `sm/contact2d_projection_outside.in |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/contact2d_projection_outside. |
   |                          | in>`_, `sm/contact2d_projection_boundary.in <https://github.com/oofem/oofem/blob/devel/tests |
   |                          | /regression/sm/contact2d_projection_boundary.in>`_ (and 9 more)                              |
   +--------------------------+----------------------------------------------------------------------------------------------+

TrplaneStrain
~~~~~~~~~~~~~

Implements an triangular three-node constant strain plane-strain  finite element. Each
node has 2 degrees of freedom. The node numbering is anti-clockwise. The element
features are summarized in :numref:`trplanestrainsummary`.

.. _TrplaneStrain:

.. tikz:: TrplaneStrain element - node and side numbering.

   \begin{tikzpicture}[scale=6,>=stealth]
    \tikzstyle{elemnode} = [draw,circle,inner sep=1,fill=white]
    \newcommand{\lcoordsys}[2]{
      \begin{scope}[transform canvas={shift={#2},scale=0.5,rotate=#1}]
       \draw[->] (0,0.05) ++(-0.02,0) -- ++(0.2,0) node[above] {$x$};
       \draw[->] (0,0.05) ++(0,-0.02) -- ++(0,0.1) node[right] {$y$};
      \end{scope}
    }
    \draw[->] (-0.05,0) -- (0.8,0) node[above] {$x_g$};
    \draw[->] (0,-0.05) -- (0,0.5) node[right] {$y_g$};

    \draw[thick]
        (0.2,0.1) node[elemnode] {} node[below] {1}
     -- (0.7,0.2) node[elemnode] {} node[below] {2} node[blue,midway,below] {1} coordinate[midway] (e1)
     -- (0.4,0.5) node[elemnode] {} node[above] {3} node[blue,midway,above right] {2} coordinate[midway] (e2)
     -- (0.2,0.1) node[blue,midway,above left] {3} coordinate[midway] (e3);

    \lcoordsys{ 12}{(e1)};
    \lcoordsys{135}{(e2)};
    \lcoordsys{243}{(e3)};
   \end{tikzpicture}

.. table:: trplanestrain element summary
   :name: trplanestrainsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestrain                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular plane-strain element                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using one point gauss integration formula.              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Nonlocal constitutive support. Edge load support, Adaptivity support.                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and are computed using numerical      |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 3). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`TrplaneStrain`)).                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch107.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch107.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

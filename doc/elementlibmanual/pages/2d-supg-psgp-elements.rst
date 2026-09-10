2D SUPG/PSGP Elements
=====================


2D SUPG/PSGP Elements
---------------------

.. _Tr1SUPG:

Tr1SUPG element
~~~~~~~~~~~~~~~

Represents the linear triangular finite element for transient incompressible flow
analysis using SUPG/PSPG stabilization with equal order approximation of velocity and
pressure fields. Each node has 3 degrees of freedoms (two components of velocity and
pressure). The node numbering is anti-clockwise. The element features are summarized in
:numref:`Tr1SUPGsummary`.

.. _Tr1SUPG2fig:

.. tikz:: Tr1SUPG element. Node numbering, Side numbering and
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

    \draw[thick]
        (0.2,0.1) node[elemnode] {} node[below] {1}
     -- (0.7,0.2) node[elemnode] {} node[below] {2} node[blue,midway,below] {1} coordinate[midway] (e1)
     -- (0.4,0.5) node[elemnode] {} node[above] {3} node[blue,midway,above right] {2} coordinate[midway] (e2)
     -- (0.2,0.1) node[blue,midway,above left] {3} coordinate[midway] (e3);

    \lcoordsys{ 12}{(e1)};
    \lcoordsys{135}{(e2)};
    \lcoordsys{243}{(e3)};
   \end{tikzpicture}

.. table:: Tr1SUPG element summary
   :name: Tr1SUPGsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr1SUPG                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | linear triangular finite element for transient incompressible flow analysis using SUPG/PSPG  |
   |                          | algorithm                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`vof{rn}` :optelemparam:`pvof{rn}`                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two velocity components (V_u and V_v) and pressure (P_f) are required in each node.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of velocity and pressure fields.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | exact                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant boundary tractions are supported. Body loads representing the self-weight load are  |
   |                          | supported.                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Multi-fluid analysis     | The element has support for solving problems with two immiscible fluids in a fixed spatial   |
   |                          | domain. In the present implementation, a VOF and LevelSet tracking algorithms are used to    |
   |                          | track the position of interface. In case of VOF tracking, an initial VOF fraction (volume    |
   |                          | fraction of reference fluid) can be specified using :param:`vof` (default is zero). Element  |
   |                          | can also be marked as allways filled with reference fluid (some form of source) using        |
   |                          | parameter :param:`pvof` which specifies the permanent VOF value. In case of LevelSet         |
   |                          | tracking, the initial levelset is specified using reference polygon (see corresponding       |
   |                          | levelset record in oofem input manual). The material model should be of type **Keyword**:    |
   |                          | :param:`twofluidmat`, that supports modelling of two immiscible fluids.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _Tr21SUPG:

Tr21SUPG element
~~~~~~~~~~~~~~~~

Implementation of P2P1 Taylor Hood element for transient incompressible flow analysis
using SUPG and LSIC stabilization. It consists of globally continuous, piecewise
quadratic functions for  approximation in velocity space and globally continuous,
piecewise linear functions for approximation in pressure space. LBB condition is
satisfied. There are 3 degrees of freedom in vertices (two components of velocity and
pressure), and 2 degrees of freedom in edge nodes (two components of velocity only). The
node numbering is anti-clockwise, vertices are numbered first. The element features are
summarized in :numref:`Tr21SUPGsummary`.

.. _Tr21SUPGfig:

.. tikz:: Tr21SUPG element - node and side numbering.

   \tikzstyle{elemnode} = [draw,circle,inner sep=1,fill=white]
   \begin{tikzpicture}[scale=6,>=stealth]
    \draw[->] (-0.05,0) -- (0.8,0) node[above] {$x_g$};
    \draw[->] (0,-0.05) -- (0,0.5) node[right] {$y_g$};

    \draw[thick,xshift=-2]
        (0.2,0.1) node[elemnode] {} node[below] {1}
     to[out=-10,in=210] node[elemnode,midway] {} node[below,midway] {4} node[blue,below,near end] {1} (0.7,0.2)  node[elemnode] {} node[below right] {2}
     to[out=110,in=-20] node[elemnode,midway] {} node[above right,midway] {5} node[blue,above right,near end] {2} (0.4,0.5) node[elemnode] {} node[above] {3}
     to[out=-90,in=40] node[elemnode,midway] {} node[above left,midway] {6} node[blue,above left,near end] {3} (0.2,0.1);
   \end{tikzpicture}

.. table:: Tr21SUPG element summary
   :name: Tr21SUPGsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr21SUPG                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | P2P1 Taylor Hood element                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two velocity components (V_u and V_v) and pressure (P_f) in vertices and two velocity        |
   |                          | components (V_u and V_v) in edge nodes are required.                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of velocity and linear approximation of pressure fields.             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration is exact, each submatrix of element stiffness matrix is evaluated in proper      |
   |                          | number of Gauss points. Submatrices connected with velocity are evaluated in 7 or 13 points, |
   |                          | mixed velocity-pressure submatrices in 3 or 7 points, submatrices connected with pressure in |
   |                          | 3 points.                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant boundary tractions are supported. Body loads representing the self-weight load are  |
   |                          | supported.                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Multi-fluid analysis     | The element has no support for solving problems with two immiscible fluids in a fixed        |
   |                          | spatial domain.                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _Tr1SUPGAxi:

Tr1SUPGAxi element
~~~~~~~~~~~~~~~~~~

Represents the linear triangular finite element for transient incompressible flow
analysis using SUPG/PSPG stabilization with equal order approximation of velocity and
pressure fields in 2d-axisymmetric setting. Each node has 3 degrees of freedoms (two
components of velocity and pressure). The y-axis is axis of ratational symmetry. The
node numbering is anti-clockwise. The element features are summarized in Table
:numref:`Tr1SUPGAxisummary`.

.. _Tr1SUPGAxifig:

.. tikz:: Tr1SUPGAxi element. Node numbering, Side numbering and
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

    \draw[thick]
        (0.2,0.1) node[elemnode] {} node[below] {1}
     -- (0.7,0.2) node[elemnode] {} node[below] {2} node[blue,midway,below] {1} coordinate[midway] (e1)
     -- (0.4,0.5) node[elemnode] {} node[above] {3} node[blue,midway,above right] {2} coordinate[midway] (e2)
     -- (0.2,0.1) node[blue,midway,above left] {3} coordinate[midway] (e3);

    \lcoordsys{ 12}{(e1)};
    \lcoordsys{135}{(e2)};
    \lcoordsys{243}{(e3)};
   \end{tikzpicture}

.. table:: Tr1SUPGAxi element summary
   :name: Tr1SUPGAxisummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr1SUPGAxi                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | linear equal order approximation axisymmetric element                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`vof{rn}`:optelemparam:`pvof{rn}`                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two velocity components (V_u and V_v) and pressure (P_f) are required in each node.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of velocity and pressure fields.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Gauss integration in seven point employed.                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant boundary tractions are supported. Body loads representing the self-weight load are  |
   |                          | supported.                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Multi-fluid analysis     | The element has support for solving problems with two immiscible fluids in a fixed spatial   |
   |                          | domain. In the present implementation, a VOF tracking algorithm is used to track the         |
   |                          | position of interface. An initial VOF fraction (volume fraction of reference fluid) can be   |
   |                          | specified using :param:`vof` (default is zero). Element can also be marked as always filled  |
   |                          | with reference fluid (some form of source) using parameter :param:`pvof` which specifies the |
   |                          | permanent VOF value. In this case, the material model should be of type **Keyword**:         |
   |                          | :param:`twofluidmat`, that supports modelling of two immiscible fluids.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+

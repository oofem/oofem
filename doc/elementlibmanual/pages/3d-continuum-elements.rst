3D Continuum Elements
=====================


3D Continuum Elements
---------------------

This section contains description of continuum elements.

.. _lspacesect:

LSpace element
~~~~~~~~~~~~~~

Implementation of Linear 3d eight - node  finite element. Each node has 3 degrees of
freedom. The element features are summarized in :numref:`lspacesummary`.

.. tikz:: LSpace element (Node numbers in black, side numbers in blue,
 and surface numbers in red).

   \begin{tikzpicture}[scale=2,>=stealth,
     x={(-0.4cm,-0.3cm)}, y={ (1cm,0cm) }, z={(0cm,1cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]
    \newcommand{\fs}{0.25}

   % Coord.sys. (shifted for readability)
    \draw[->,xshift=-5] (-0.1,0,0) -- (1.5,0,0) node[at end,below] {$\xi$};
    \draw[->,xshift=-5] (0,-0.1,0) -- (0,1,0) node[at end,below] {$\eta$};
    \draw[->,xshift=-5] (0,0,-0.1) -- (0,0,0.8) node[at end,right] {$\zeta$};

   % Can't use rectangle in 3d
    \draw[thick] (-1,-1,1) -- (-1,1,1) -- (1,1,1) -- (1,-1,1) -- cycle;
    \draw[thick,background] (-1,-1,-1) -- (-1,1,-1);
    \draw[thick] (-1,1,-1)-- (1,1,-1);
    \draw[thick] (1,1,-1) -- (1,-1,-1);
    \draw[thick,background] (1,-1,-1) -- (-1,-1,-1);
    \draw[thick,background] (-1,-1,1) -- (-1,-1,-1);
    \draw[thick] (-1,1,1) -- (-1,1,-1);
    \draw[thick] (1,-1,1) -- (1,-1,-1);
    \draw[thick] (1,1,1) -- (1,1,-1);

   % Faces
    \draw[red] (1,-\fs,-\fs) -- (1,-\fs,\fs) -- (1,\fs,\fs) -- (1,\fs,-\fs) -- cycle;
    \node[red] at (1,0,0) {5};
    \draw[red] (-\fs,1,-\fs) -- (-\fs,1,\fs) -- (\fs,1,\fs) -- (\fs,1,-\fs) -- cycle;
    \node[red] at (0,1,0) {4};
    \draw[red] (-\fs,-\fs,1) -- (-\fs,\fs,1) -- (\fs,\fs,1) -- (\fs,-\fs,1) -- cycle;
    \node[red] at (0,0,1) {1};
    \draw[red!50!black,background] (-1,-\fs,-\fs) -- (-1,-\fs,\fs) -- (-1,\fs,\fs) -- (-1,\fs,-\fs) -- cycle;
    \node[red!50!black] at (-1,0,0) {3};
    \draw[red!50!black,background] (-\fs,-1,-\fs) -- (-\fs,-1,\fs) -- (\fs,-1,\fs) -- (\fs,-1,-\fs) -- cycle;
    \node[red!50!black] at (0,-1,0) {6};
    \draw[red!50!black,background] (-\fs,-\fs,-1) -- (-\fs,\fs,-1) -- (\fs,\fs,-1) -- (\fs,-\fs,-1) -- cycle;
    \node[red!50!black] at (0,0,-1) {2};

   % Nodes
    \node[elemnode] (n1) at (-1,-1, 1) {}; \node[above left] at (n1) {1};
    \node[elemnode] (n2) at (-1, 1, 1) {}; \node[above right] at (n2) {2};
    \node[elemnode] (n3) at ( 1, 1, 1) {}; \node[above left] at (n3) {3};
    \node[elemnode] (n4) at ( 1,-1, 1) {}; \node[above left] at (n4) {4};
    \node[elemnode] (n5) at (-1,-1,-1) {}; \node[below right] at (n5) {5};
    \node[elemnode] (n6) at (-1, 1,-1) {}; \node[below right] at (n6) {6};
    \node[elemnode] (n7) at ( 1, 1,-1) {}; \node[below right] at (n7) {7};
    \node[elemnode] (n8) at ( 1,-1,-1) {}; \node[below right] at (n8) {8};

   % Edges
    \node[blue,above] at (-1,0,1) {1};
    \node[blue,above left] at (0,1,1) {2};
    \node[blue,above left] at (1,0,1) {3};
    \node[blue,above left] at (0,-1,1) {4};

    \node[blue!50!black,left] at (-1,-1,0) {5};
    \node[blue,left] at (-1,1,0) {6};
    \node[blue,left] at (1,1,0) {7};
    \node[blue,left] at (1,-1,0) {8};

    \node[blue!50!black,above] at (-1,0,-1) {9};
    \node[blue,above left] at (0,1,-1) {10};
    \node[blue,above left] at (1,0,-1) {11};
    \node[blue!50!black,above left] at (0,-1,-1) {12};
   \end{tikzpicture}

.. table:: lspace element summary
   :name: lspacesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | lspace                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Linear isoparametric brick element                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points (possible completions are 1, 8  |
   |                          | (default), or 27).                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacement and geometry.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components.                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Supports adaptivity, geometric nonlinearity, and layered cross section support               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/python/usrdefboundaryload01.in <https://github.com/oofem/oofem/blob/devel/tests/regressi |
   |                          | on/sm/python/usrdefboundaryload01.in>`_, `smmfront/mfront01.in                               |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/smmfront/mfront01.in>`_,         |
   |                          | `sm/deadweight02.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/deadweight02.in>`_,           |
   |                          | `sm/brick_nlgeo_1.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_1.in>`_,          |
   |                          | `sm/blatzko1.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/blatzko1.in>`_,               |
   |                          | `sm/brick_nlgeo_2.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_2.in>`_,          |
   |                          | `sm/brick_nlgeo_3.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_3.in>`_,          |
   |                          | `sm/brick_nlgeo_4.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_4.in>`_,          |
   |                          | `sm/brick_nlgeo_5.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_5.in>`_,          |
   |                          | `sm/brick_nlgeo_6.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/brick_nlgeo_6.in>`_ (and 30   |
   |                          | more)                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+

LSpaceBB element
~~~~~~~~~~~~~~~~

Implementation of 3d brick eight - node  linear approximation element with selective
integration of deviatoric and volumetric  strain contributions (B-bar formulation) for
incompressible problems.  Features and description identical to conventional lspace
element, see section :ref:`lspacesect`.

.. _QSpace_element:

QSpace element
~~~~~~~~~~~~~~

Implementation of quadratic 3d 20-node  finite element. Each node has 3 degrees of
freedom. The element features are summarized in :numref:`qspacesummary`.

.. tikz:: QSpace element.

   \begin{tikzpicture}[scale=1.75,>=stealth,
     x={(-0.4cm,-0.3cm)}, y={ (1cm,0cm) }, z={(0cm,1cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]

   % Coord.sys. (shifted for readability)
    \draw[->,xshift=-5] (-0.1,0,0) -- (1.5,0,0) node[at end,below] {$\xi$};
    \draw[->,xshift=-5] (0,-0.1,0) -- (0,1,0) node[at end,below] {$\eta$};
    \draw[->,xshift=-5] (0,0,-0.1) -- (0,0,0.8) node[at end,right] {$\zeta$};

   % Can't use rectangle in 3d
    \draw[thick] (-1,-1,1) -- (-1,1,1) -- (1,1,1) -- (1,-1,1) -- cycle;
    \draw[thick,background] (-1,-1,-1) -- (-1,1,-1);
    \draw[thick] (-1,1,-1)-- (1,1,-1);
    \draw[thick] (1,1,-1) -- (1,-1,-1);
    \draw[thick,background] (1,-1,-1) -- (-1,-1,-1);
    \draw[thick,background] (-1,-1,1) -- (-1,-1,-1);
    \draw[thick] (-1,1,1) -- (-1,1,-1);
    \draw[thick] (1,-1,1) -- (1,-1,-1);
    \draw[thick] (1,1,1) -- (1,1,-1);

   % Nodes
    \node[elemnode] (n1) at (-1,-1, 1) {}; \node[above left] at (n1) {1};
    \node[elemnode] (n2) at (-1, 1, 1) {}; \node[above right] at (n2) {2};
    \node[elemnode] (n3) at ( 1, 1, 1) {}; \node[above left] at (n3) {3};
    \node[elemnode] (n4) at ( 1,-1, 1) {}; \node[above left] at (n4) {4};
    \node[elemnode] (n5) at (-1,-1,-1) {}; \node[below right] at (n5) {5};
    \node[elemnode] (n6) at (-1, 1,-1) {}; \node[below right] at (n6) {6};
    \node[elemnode] (n7) at ( 1, 1,-1) {}; \node[below right] at (n7) {7};
    \node[elemnode] (n8) at ( 1,-1,-1) {}; \node[below right] at (n8) {8};

    \node[elemnode] (n9)  at (-1, 0, 1) {}; \node[above] at (n9) {9};
    \node[elemnode] (n10) at ( 0, 1, 1) {}; \node[above left] at (n10) {10};
    \node[elemnode] (n11) at ( 1, 0, 1) {}; \node[above] at (n11) {11};
    \node[elemnode] (n12) at ( 0,-1, 1) {}; \node[above left] at (n12) {12};

    \node[elemnode] (n13) at (-1, 0,-1) {}; \node[below] at (n13) {13};
    \node[elemnode] (n14) at ( 0, 1,-1) {}; \node[below right] at (n14) {14};
    \node[elemnode] (n15) at ( 1, 0,-1) {}; \node[below] at (n15) {15};
    \node[elemnode] (n16) at ( 0,-1,-1) {}; \node[below right] at (n16) {16};

    \node[elemnode] (n17) at (-1,-1, 0) {}; \node[right] at (n17) {17};
    \node[elemnode] (n18) at (-1, 1, 0) {}; \node[right] at (n18) {18};
    \node[elemnode] (n19) at ( 1, 1, 0) {}; \node[right] at (n19) {19};
    \node[elemnode] (n20) at ( 1,-1, 0) {}; \node[right] at (n20) {20};

   \end{tikzpicture}

.. table:: qspace element summary
   :name: qspacesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | qspace                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Quadratic isoparametric brick element                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points (possible completions are 1, 8  |
   |                          | (default), or 27).                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacement and geometry.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components.                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/layered_cube_lcs.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube_lcs.in>`_,       |
   |                          | `sm/layered_cube.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layered_cube.in>`_,           |
   |                          | `sm/vitrification01.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/vitrification01.in>`_,        |
   |                          | `sm/cantilever_Qspace.in                                                                     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/cantilever_Qspace.in>`_       |
   +--------------------------+----------------------------------------------------------------------------------------------+

LTRSpace element
~~~~~~~~~~~~~~~~

Implementation of tetrahedra four-node finite element.  Each node has 3 degrees of
freedom. The element features are summarized in :numref:`LTRSpacesummary`.
Following node numbering convention is adopted (see also Fig.
:ref:`lintetrahedron_fig`):

- Select a face that will contain the first three corners. The excluded corner
  will be the last one.

- Number these three corners in a counterclockwise sense when looking at the
  face from the  excluded corner.

.. _lintetrahedron_fig:

.. tikz:: LTRSpace element. Definition and node numbering convention.

   \begin{tikzpicture}[scale=4,>=stealth,x={(1cm,0cm)}, y={ (0.4cm,-0.3cm) }, z={(0.4cm,0.8cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]
    \newcommand{\fs}{0.23}

   % Can't use rectangle in 3d
    \draw[thick,background] (0,0,0) -- (1,0,0) node[midway,above,blue!50!black] {3};
    \draw[thick] (0,0,1) -- (0,0,0) node[midway,above left,blue] {4}
                         -- (0,1,0) node[midway,below left,blue] {1}
                         -- (1,0,0) node[midway,below right,blue] {2}
                         -- (0,0,1) node[midway,above right,blue] {6}
                         -- (0,1,0) node[near start,left,blue] {5};

   % Nodes
    \node[elemnode] (n1) at (0,0,0) {}; \node[below left] at (n1) {1};
    \node[elemnode] (n2) at (0,1,0) {}; \node[below left] at (n2) {2};
    \node[elemnode] (n3) at (1,0,0) {}; \node[below right] at (n3) {3};
    \node[elemnode] (n4) at (0,0,1) {}; \node[above left] at (n4) {4};

   % Faces
    \draw[red!50!black,background] (\fs,\fs,0) -- (\fs,1-2*\fs,0) -- (1-2*\fs,\fs,0) -- cycle;
    \node[red!50!black] at (1/3,1/3,0) {1};
    \draw[red!50!black,background] (\fs,0,\fs) -- (\fs,0,1-2*\fs) -- (1-2*\fs,0,\fs) -- cycle;
    \node[red!50!black] at (1/3,0,1/3) {4};

    \draw[red] (0,\fs,\fs) -- (0,\fs,1-2*\fs) -- (0,1-2*\fs,\fs) -- cycle;
    \node[red] at (0,1/3,1/3) {2};
    \draw[red] (1-2*\fs,\fs,\fs) -- (\fs,1-2*\fs,\fs) -- (\fs,\fs,1-2*\fs) -- cycle;
    \node[red] at (1/3,1/3,1/3) {3};

   \end{tikzpicture}

.. table:: LTRSpace element summary
   :name: LTRSpacesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | LTRSpace                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Linear tetrahedra element                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points (possible completions are 1, 8  |
   |                          | (default), or 27).                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry using linear volume coordinates.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using four point Gauss integration formula.        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Adaptivity support, Geometric nonlinearity support.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Surface and Edge loadings supported.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/con2dpm9.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm9.in>`_,               |
   |                          | `sm/con2dpm10.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm10.in>`_,              |
   |                          | `sm/con2dpm11.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm11.in>`_,              |
   |                          | `sm/con2dpm12.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm12.in>`_,              |
   |                          | `sm/con2dpm6.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm6.in>`_,               |
   |                          | `sm/con2dpm7.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm7.in>`_,               |
   |                          | `sm/con1dpm1.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con1dpm1.in>`_,               |
   |                          | `sm/con1dpm3.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con1dpm3.in>`_,               |
   |                          | `sm/con2dpm1.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm1.in>`_,               |
   |                          | `sm/con2dpm2.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm2.in>`_ (and 15 more)  |
   +--------------------------+----------------------------------------------------------------------------------------------+

QTRSpace element
~~~~~~~~~~~~~~~~

Implementation of tetrahedra ten-node finite element.  Each node has 3 degrees of
freedom. The element features are summarized in :numref:`QTRSpacesummary`.
Following node numbering convention is adopted (see also Fig.
:ref:`qtetrahedron_fig`):

.. _qtetrahedron_fig:

.. tikz:: QTRSpace element. Definition and node numbering convention.

   \begin{tikzpicture}[scale=4,>=stealth,x={(1cm,0cm)}, y={ (0.4cm,-0.3cm) }, z={(0.4cm,0.8cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]
    \newcommand{\fs}{0.23}

   % Can't use rectangle in 3d
    \draw[thick,background] (0,0,0) -- (1,0,0) node[midway,above,blue!50!black]{};
    \draw[thick] (0,0,1) -- (0,0,0) node[midway,above left,blue] {}
                         -- (0,1,0) node[midway,below left,blue] {}
                         -- (1,0,0) node[midway,below right,blue]{}
                         -- (0,0,1) node[midway,above right,blue]{}
                         -- (0,1,0) node[near start,left,blue] {};
   % Nodes
    \node[elemnode] (n1) at (0,0,0) {}; \node[below left] at (n1) {1};
    \node[elemnode] (n2) at (0,1,0) {}; \node[below left] at (n2) {2};
    \node[elemnode] (n3) at (1,0,0) {}; \node[below right] at (n3) {3};
    \node[elemnode] (n4) at (0,0,1) {}; \node[above left] at (n4) {4};
    \node[elemnode] (n5) at (0,0.5,0) {}; \node[below left] at (n5) {5};
    \node[elemnode] (n6) at (0.5,0.5,0) {}; \node[below right] at (n6) {6};
    \node[elemnode] (n7) at (0.5,0.,0) {}; \node[above right] at (n7) {7};
    \node[elemnode] (n8) at (0,0,0.5) {}; \node[above left] at (n8) {8};
    \node[elemnode] (n9) at (0,0.5,0.5) {}; \node[above right] at (n9) {9};
    \node[elemnode] (n10) at (0.5,0,0.5) {}; \node[above right] at (n10) {10};

   \end{tikzpicture}

.. table:: QTRSpace element summary
   :name: QTRSpacesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | QTRSpace                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D tetrahedra element with quadratic interpolation                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to alter the default integration formula (possible completions are 1, 4 |
   |                          | (default), 5, 11, 15, 24, and 45 point intergartion formulas).                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry using linear volume coordinates.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using four point Gauss integration formula.        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

LWedge element
~~~~~~~~~~~~~~

Implementation of wedge six-node finite element.  Each node has 3 degrees of freedom.
The element features are summarized in :numref:`LWedgesummary`. Following node
numbering convention is adopted (see also Fig. :ref:`linwedge_fig`):

.. _linwedge_fig:

.. tikz:: LWedge element. Node numbering convention in black, edge numbering in blue and face numbering in red.

   \begin{tikzpicture}[scale=4,>=stealth,x={(1cm,0cm)}, y={ (0.4cm,-0.3cm) }, z={(0.0cm,0.8cm)}]
    %\tikzstyle{elemnode} = [fill,circle,inner sep=2]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]
    \newcommand{\fs}{0.20}
    \newcommand{\fsb}{0.35}

     \coordinate (n1) at (0,0,0);
     \coordinate (n2) at (0,1,0);
     \coordinate (n3) at (1,0,0);
     \coordinate (n4) at (0,0,1);
     \coordinate (n5) at (0,1,1);
     \coordinate (n6) at (1,0,1);

     % Faces behind
     \draw[red!50!black,background] (\fs,\fs,0) -- (\fs,1-2*\fs,0) -- (1-2*\fs,\fs,0) -- cycle;
     \node[red!50!black] at (1/3,1/3,0) {1};
     \draw[red!50!black,background] (\fsb,0,\fsb) -- (\fsb,0,1-\fsb) -- (1-\fsb,0,1-\fsb) -- (1-\fsb,0,\fsb) -- cycle;
     \node[red!50!black] at (1/2,0,1/2) {5};

     % Can't use rectangle in 3d
     \draw[thick,background] (n1) -- (n3) node[midway,above,blue!50!black] {3};
     \draw[thick] (n2) -- (n1) node[midway,below,blue] {1};
     \draw[thick] (n2) -- (n3) node[midway,below,blue] {2};
     \draw[thick] (n5) -- (n4) node[midway,below,blue] {4};
     \draw[thick] (n5) -- (n6) node[midway,below,blue] {5};
     \draw[thick] (n4) -- (n6) node[midway,above,blue] {6};
     \draw[thick] (n1) -- (n4) node[midway,left,blue] {7};
     \draw[thick] (n2) -- (n5) node[midway,above left,blue] {8};
     \draw[thick] (n3) -- (n6) node[midway,right,blue] {9};

     % Nodes
     \node[elemnode] at (n1) {}; \node[below left] at (n1) {1};
     \node[elemnode] at (n2) {}; \node[below] at (n2) {2};
     \node[elemnode] at (n3) {}; \node[below right] at (n3) {3};
     \node[elemnode] at (n4) {}; \node[left] at (n4) {4};
     \node[elemnode] at (n5) {}; \node[above] at (n5) {5};
     \node[elemnode] at (n6) {}; \node[right] at (n6) {6};

     % Faces
     \draw[red,background] (\fs,\fs,1) -- (\fs,1-2*\fs,1) -- (1-2*\fs,\fs,1) -- cycle;
     \node[red] at (1/3,1/3,1) {2};

     \draw[red] (0,\fsb,\fsb) -- (0,\fsb,1-\fsb) -- (0,1-\fsb,1-\fsb) -- (0,1-\fsb,\fsb) -- cycle;
     \node[red] at (0,1/2,1/2) {3};

     \draw[red] (\fsb,1-\fsb,\fsb) -- (1-\fsb,\fsb,\fsb) -- (1-\fsb,\fsb,1-\fsb) -- (\fsb,1-\fsb,1-\fsb) -- cycle;
     \node[red] at (1/2,1/2,1/2) {4};

   \end{tikzpicture}

.. table:: LWedge element summary
   :name: LWedgesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | LWedge                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D wedge six-node finite element with linear interpolation                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to alter the default integration formula (possible completions are 2    |
   |                          | (default) and 9 point integration formulas).                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using four point Gauss integration formula.        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

QWedge element
~~~~~~~~~~~~~~

Implementation of wedge fifteen-node finite element.  Each node has 3 degrees of
freedom. The element features are summarized in :numref:`QWedgesummary`. Following
node numbering convention is adopted (see also Fig. :ref:`qwedge_fig`):

.. _qwedge_fig:

.. tikz:: QWedge element. Node numbering convention in black, edge numbering in blue and face numbering in red.

   \begin{tikzpicture}[scale=4,>=stealth,x={(1cm,0cm)}, y={ (0.4cm,-0.3cm) }, z={(0.0cm,0.8cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [densely dashed]

     \coordinate (n1) at (0,0,0);
     \coordinate (n2) at (0,1,0);
     \coordinate (n3) at (1,0,0);
     \coordinate (n4) at (0,0,1);
     \coordinate (n5) at (0,1,1);
     \coordinate (n6) at (1,0,1);

     % Can't use rectangle in 3d
     \draw[thick,background] (n1) -- (n3) coordinate[midway] (e3);
     \draw[thick] (n2) -- (n1) coordinate[midway] (e1);
     \draw[thick] (n2) -- (n3) coordinate[midway] (e2);
     \draw[thick] (n5) -- (n4) coordinate[midway] (e4);
     \draw[thick] (n5) -- (n6) coordinate[midway] (e5);
     \draw[thick] (n4) -- (n6) coordinate[midway] (e6);
     \draw[thick] (n1) -- (n4) coordinate[midway] (e7);
     \draw[thick] (n2) -- (n5) coordinate[midway] (e8);
     \draw[thick] (n3) -- (n6) coordinate[midway] (e9);

     % Nodes
     \node[elemnode] at (n1) {}; \node[below left ] at (n1) {1};
     \node[elemnode] at (n2) {}; \node[below      ] at (n2) {2};
     \node[elemnode] at (n3) {}; \node[below right] at (n3) {3};
     \node[elemnode] at (n4) {}; \node[      left ] at (n4) {4};
     \node[elemnode] at (n5) {}; \node[below left ] at (n5) {5};
     \node[elemnode] at (n6) {}; \node[      right] at (n6) {6};

     \node[elemnode] at (e1) {}; \node[below left ] at (e1) {7};
     \node[elemnode] at (e2) {}; \node[below right] at (e2) {8};
     \node[elemnode] at (e3) {}; \node[below      ] at (e3) {9};
     \node[elemnode] at (e4) {}; \node[below left ] at (e4) {10};
     \node[elemnode] at (e5) {}; \node[below right] at (e5) {11};
     \node[elemnode] at (e6) {}; \node[above      ] at (e6) {12};
     \node[elemnode] at (e7) {}; \node[      left ] at (e7) {13};
     \node[elemnode] at (e8) {}; \node[above left ] at (e8) {14};
     \node[elemnode] at (e9) {}; \node[      right] at (e9) {15};

   \end{tikzpicture}

.. table:: QWedge element summary
   :name: QWedgesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | QWedge                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D wedge six-node finite element with quadratic interpolation                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to alter the default integration formula (possible completions are 2    |
   |                          | (default) and 9 point integration formulas).                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using four point Gauss integration formula.        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0,1,2.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/qwedge_01.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/qwedge_01.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

Layer stacking sequence definition for 3D elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Selected 3D elements (bricks and wedge geometries) support using the LayeredCrossSection
model to define layer stack as a sequence of individual layers. Individual layers are
assumed to lie in element parametric :math:`\xi-\eta` plane and are stacked along
parametric :math:`\zeta` coordinate of the element. The direction of parametric
coordinates is determined by element node numbering convention, see figures with element
geometries above. Note, that the stacking direction is in general the function of
element geometry.

It is important to understand concept of element and material coordinate systems.

The element coordinate system (elemCS) coincides, by default, with the global coordinate
system.  The user-defined element coordinate system can be defined using lcs parameter.
The lcs parameter defines an array of size 6, where the first 3 components define
direction of  local element-axis and remaining 3 components define direction of local
element y-axis.  The local element z axis is computed using vector product
:math:`e_z=e_x\times e_y`.

The material properties of each layer are defined in material coordinate system (matCS).
Also the solver output for individual layers is done in matCS. By default, the material
coordinate system coincides with global coordinate system. Additionally, the material
coordinate system for individual layer can be rotated around material CS z-axis by
angle, defined by layered cross section :param:`rotations` keyword. This array parameter
allows to define rotation angle for individual layers and should be defined in degrees
not radians. If :param:`matcs` element keyword is present, but no lcs element record is
defined, then the following definition of elemCS is assumed:
:math:`e_x=\{\frac{dx(\xi,\eta,\zeta)}{d\xi},\frac{dy(\xi,\eta,\zeta)}{d\xi},\frac{dz(\xi,\eta,\zeta)}{d\xi}\},
h=\{\frac{dx(\xi,\eta,\zeta)}{d\eta},\frac{dy(\xi,\eta,\zeta)}{d\eta},\frac{dz(\xi,\eta,\zeta)}{d\eta}\},
e_z=e_x\times h, e_y=e_z\times e_x`, where :math:`\xi,\eta,\zeta` are parametric element
coordinates.  The strains and stresses in individual layers are always reported in
material coordinate system. LayeredCS integration  The layered cross section integration
can be set up using number of integrations points in layer plane
(:param:`nintegrationpoints` parameter ) and using number of integration points per
layer  thickness (:param:`layerintegrationpoints` parameter).

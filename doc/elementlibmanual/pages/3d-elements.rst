3D Elements
===========


3D Elements
-----------

Tetrah1ht - tetrahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Represents isoparametric four-node tetrahedral element. Each node has 1 degree of
freedom. The same numbering convection is adopted as in mechanics, see Fig.
:ref:`lintetrahedron_fig`. The element features are summarized in Table
:numref:`Tetrah1htsummary`.

.. table:: Tetrah1ht element summary
   :name: Tetrah1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tetrah1ht                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Isoparametric, four-node tetrahedral element with linear approximation for heat transfer     |
   |                          | problems                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (T_f - temperature) is required in each node.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of temperature.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using gauss integration formula in 1 (the default), or 4 integration points. The |
   |                          | default number of integration point used can be overloaded using :param:`NIP` parameter.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side and surface numbering is shown in Fig. :ref:`lintetrahedron_fig`.      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tm/tmpatch40.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch40.in>`_,              |
   |                          | `tm/tmpatch41.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch41.in>`_,              |
   |                          | `tm/tmpatch42.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch42.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _Brick1ht:

Brick1ht - hexahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Represents isoparametric eight-node brick/hexahedron finite element for heat transfer
problems. Each node has 1 degree of freedom. The element features are summarized in
:numref:`Brick1htsummary`.

.. _Brick1htfig:

.. tikz:: Brick1ht element. Node numbers are in black, side numbers are in blue,
 and surface numbers are in red.

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

.. table:: Brick1ht element summary
   :name: Brick1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Brick1ht                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Isoparametric, hexahedral 3D element with linear approximation for heat transfer problems    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (T_f - temperature) is required in each node.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of temperature.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using gauss integration formula in 8 (the default), or 27 integration points.    |
   |                          | The default number of integration point used can be overloaded using :param:`NIP` parameter. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side and surface numbering is shown in fig. (:ref:`Brick1htfig`)).          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tm/tmpatch46.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch46.in>`_,              |
   |                          | `tmcemhyd/cemhyd02.in                                                                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tmcemhyd/cemhyd02.in>`_,         |
   |                          | `tm/tmpatch31.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch31.in>`_,              |
   |                          | `tm/tmpatch34.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch34.in>`_,              |
   |                          | `tm/tmpatch32.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch32.in>`_,              |
   |                          | `tm/tmpatch33.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch33.in>`_,              |
   |                          | `tm/tmpatch37.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch37.in>`_,              |
   |                          | `tm/tmpatch36.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch36.in>`_,              |
   |                          | `tm/tmpatch39.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch39.in>`_,              |
   |                          | `tm/tmpatch35.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch35.in>`_ (and 2 more)  |
   +--------------------------+----------------------------------------------------------------------------------------------+

Brick1hmt - hexahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Represents isoparametric eight-node quadrilateral finite element for heat and mass (one
constituent) transfer problems.  Two dofs (T_f - temperature and C_1 - concentration)
are required in each node. Linear approximation of temperature and mass concentration.
Other features are similar to Brick1 element, see section :ref:`Brick1ht`.

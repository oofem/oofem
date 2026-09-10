Stokes' Flow Elements
=====================


Stokes' Flow Elements
---------------------

Stokes' flow elements neglect acceleration, and thus requires no additional
stabilization.

Tr21Stokes element
~~~~~~~~~~~~~~~~~~

Standard 6 node triangular element for stokes flow, with quadratic geometry, velocity
and linear pressure. Both compressible and incompressible material behavior is supported
(and also the seamless transition between the two). The element features are summarized
in :numref:`Tr121Stokessummary`.

.. table:: Tr121Stokes element summary
   :name: Tr121Stokessummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr121Stokes                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Standard 6 node triangular element for stokes flow, with quadratic geometry, velocity and    |
   |                          | linear pressure                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Unknown pressure in nodes 1--3 with unknown velocity (V_u and V_v) in all 6 nodes.           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of geometry and velocity, linear pressure approximation.             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Tet21Stokes element
~~~~~~~~~~~~~~~~~~~

Standard 10 node tetrahedral element for stokes flow, with quadratic geometry, velocity
and linear pressure. The element features are summarized in Table
:numref:`Tet21Stokessummary`.

.. table:: Tet21Stokes element summary
   :name: Tet21Stokessummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tet21Stokes                                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Standard 10 node tetrahedral element for stokes flow, with quadratic geometry, velocity and  |
   |                          | linear pressure                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Unknown pressure (P_f) in nodes 1--4, and unknown velocity (V_u, V_v, V_w) in all nodes.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of geometry and velocity, linear pressure approximation.             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Untested                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Hexa21Stokes element
~~~~~~~~~~~~~~~~~~~~

Standard 27 node hexahedral element for stokes flow, with quadratic geometry, velocity
and linear pressure. The element features are summarized in Table
:numref:`Hexa21Stokessummary`.

.. _Hexa21Stokes:

.. tikz:: Hexa21Stokes element. Node numbering and face numbering.

   \begin{tikzpicture}[scale=1.75,>=stealth,
     x={(-0.4cm,-0.3cm)}, y={ (1cm,0cm) }, z={(0cm,1cm)}]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \tikzstyle{background} = [draw=black!50]
    \tikzstyle{background2} = [black!50]

   % Coord.sys. (shifted for readability)
    \draw[->,thin,color=blue!50] (-0.1,0,0) -- (1.25,0,0) node[at end,below] {$\xi$};
    \draw[->,thin,color=blue!50] (0,-0.1,0) -- (0,0.8,0) node[at end,above] {$\eta$};
    \draw[->,thin,color=blue!50] (0,0,-0.1) -- (0,0,0.82) node[at end,right] {$\zeta$};

   % Can't use rectangle in 3d
    \draw[thick,background] (-1,-1,-1) -- (-1,1,-1);
    \draw[thick,background] (1,-1,-1) -- (-1,-1,-1);
    \draw[thick,background] (-1,-1,1) -- (-1,-1,-1);
    \draw[thick] (-1,-1,1) -- (-1,1,1) -- (1,1,1) -- (1,-1,1) -- cycle;
    \draw[thick] (-1,1,-1)-- (1,1,-1);
    \draw[thick] (1,1,-1) -- (1,-1,-1);
    \draw[thick] (-1,1,1) -- (-1,1,-1);
    \draw[thick] (1,-1,1) -- (1,-1,-1);
    \draw[thick] (1,1,1) -- (1,1,-1);

   % Nodes
    \node[elemnode] (n1) at (-1,-1, 1) {}; \node[above left] at (n1) {1};
    \node[elemnode] (n2) at (-1, 1, 1) {}; \node[above right] at (n2) {2};
    \node[elemnode] (n3) at ( 1, 1, 1) {}; \node[above left] at (n3) {3};
    \node[elemnode] (n4) at ( 1,-1, 1) {}; \node[above left] at (n4) {4};
    \node[elemnode,background] (n5) at (-1,-1,-1) {}; \node[below right,background2] at (n5) {5};
    \node[elemnode] (n6) at (-1, 1,-1) {}; \node[below right] at (n6) {6};
    \node[elemnode] (n7) at ( 1, 1,-1) {}; \node[below right] at (n7) {7};
    \node[elemnode] (n8) at ( 1,-1,-1) {}; \node[below left] at (n8) {8};

    \node[elemnode] (n9)  at (-1, 0, 1) {}; \node[above] at (n9) {9};
    \node[elemnode] (n10) at ( 0, 1, 1) {}; \node[above left] at (n10) {10};
    \node[elemnode] (n11) at ( 1, 0, 1) {}; \node[above] at (n11) {11};
    \node[elemnode] (n12) at ( 0,-1, 1) {}; \node[above left] at (n12) {12};

    \node[elemnode,background] (n13) at (-1, 0,-1) {}; \node[below,background2] at (n13) {13};
    \node[elemnode] (n14) at ( 0, 1,-1) {}; \node[below right] at (n14) {14};
    \node[elemnode] (n15) at ( 1, 0,-1) {}; \node[below] at (n15) {15};
    \node[elemnode,background] (n16) at ( 0,-1,-1) {}; \node[below right,background2] at (n16) {16};

    \node[elemnode,background] (n17) at (-1,-1, 0) {}; \node[right,background2] at (n17) {17};
    \node[elemnode] (n18) at (-1, 1, 0) {}; \node[right] at (n18) {18};
    \node[elemnode] (n19) at ( 1, 1, 0) {}; \node[right] at (n19) {19};
    \node[elemnode] (n20) at ( 1,-1, 0) {}; \node[left] at (n20) {20};

    \node[elemnode] (n21) at ( 0, 0, 1) {}; \node[right] at (n21) {21};
    \node[elemnode,background] (n22) at ( 0, 0,-1) {}; \node[right,background2] at (n22) {22};
    \node[elemnode,background] (n23) at (-1, 0, 0) {}; \node[above,background2] at (n23) {23};
    \node[elemnode] (n24) at ( 0, 1, 0) {}; \node[right] at (n24) {24};
    \node[elemnode] (n25) at ( 1, 0, 0) {}; \node[right] at (n25) {25};
    \node[elemnode,background] (n26) at ( 0,-1, 0) {}; \node[right,background2] at (n26) {26};

    \node[elemnode,background] (n27) at ( 0,0, 0) {}; \node[above right,background2] at (n27) {27};

   \end{tikzpicture}

.. table:: Hexa21Stokes element summary
   :name: Hexa21Stokessummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Hexa21Stokes                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Standard 10 node tetrahedral element for stokes flow, with quadratic geometry, velocity and  |
   |                          | linear pressure                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Unknown pressure (P_f) in nodes 1--8, and unknown velocity (V_u, V_v, V_w) in all nodes.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of geometry and velocity, linear pressure approximation.             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Untested                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Tr1BubbleStokes element
~~~~~~~~~~~~~~~~~~~~~~~

So called “Mini” element in 2D. A 3 node triangular element for stokes flow, with linear
geometry and pressure. Velocity is enriched by a bubble function. Should not be used
with materials that have memory (which is uncommon for flow problems). The element
features are summarized in :numref:`Tr1BubbleStokessummary`.

.. table:: Tr1BubbleStokes element summary
   :name: Tr1BubbleStokessummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr1BubbleStokes                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | So called “Mini” 2D element                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Unknown pressure (P_f) in all nodes, and unknown velocity (V_u, V_v) in all nodes and one    |
   |                          | internal dof manager.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear geometry and pressure. Velocity is enriched by a bubble function.                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Untested                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Tet1BubbleStokes element
~~~~~~~~~~~~~~~~~~~~~~~~

So called “Mini” element in 3D. A 4 node tetrahedral element for stokes flow, with
linear geometry and pressure. Velocity is enriched by a bubble function. Should not be
used with materials that have memory (which is uncommon for flow problems). The element
features are summarized in :numref:`Tet1BubbleStokessummary`.

.. table:: Tet1BubbleStokes element summary
   :name: Tet1BubbleStokessummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tet1BubbleStokes                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | So called “Mini” 3D element                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Unknown pressure (P_f) in all nodes, and unknown velocity (V_u, V_v) in all nodes and one    |
   |                          | internal dof manager.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear geometry and pressure. Velocity is enriched by a bubble function.                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Untested                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

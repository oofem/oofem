3D SUPG/PSGP Elements
=====================


3D SUPG/PSGP Elements
---------------------

.. _PY1_3D_SUPG:

Tet1_3D_SUPG element
~~~~~~~~~~~~~~~~~~~~

Represents 3D linear pyramid element for transient incompressible flow analysis using
SUPG/PSPG stabilization with equal order approximation of velocity and pressure fields.
Each node has 3 degrees of freedoms (two components of velocity and pressure). The
element features are summarized in :numref:`TET1SUPGsummary`.

.. _PY1_3D_SUPGfig:

.. tikz:: Tet1_3D_SUPG element.

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

.. table:: TET1SUPG element summary
   :name: TET1SUPGsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | TET1SUPG                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | linear equal order approximation axisymmetric element                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`vof{rn}`:optelemparam:`pvof{rn}`                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three velocity components (V_u, V_v, and V_w) and pressure (P_f) are required in each node.  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of velocity and pressure fields.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | exact                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Constant boundary tractions are supported. Body loads representing the self-weight load are  |
   |                          | supported.                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Multi-fluid analysis     | The element has support for solving problems with two immiscible fluids in a fixed spatial   |
   |                          | domain. In the present implementation, a LevelSet tracking algorithm is used to track the    |
   |                          | position of interface. The material model should be of type **Keyword**:                     |
   |                          | :param:`twofluidmat`, that supports modelling of two immiscible fluids.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+

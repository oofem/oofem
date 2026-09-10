Interface elements
==================


Interface elements
------------------

Interface elements represents an interaction between points, edges or surfaces. They are
used to model debonding between surfaces or more general fracture processes through the
use of cohesive zone (cz) models. They can also be used to model contact between
elements. Specific interface material models needs to be used - consult the
*matlibmanual* for supported models.

Ordering convention:
~~~~~~~~~~~~~~~~~~~~

All inerface elements have *plus*-side and a *minus*-side and all nodes should first be
specified for the minus-side and then the plus-side. The normal to the interface is
defined to point from the minus-side to the-plus side. Direction of normal vector on an
element specifies normal stress/cohesion across the element. It is assumed that normal
jump, normal traction are at the first position of corresponding vectors. Stiffness
matrix in local coordinates has always on position 1,1 normal stiffness and then shear
stiffness.

IntElPoint, Interface1d elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of one dimensional (slip) interface element.  This element connects two
separate nodes and the interaction is governed by a one-dimensional slip law.  This law
determines the force acting between the nodes as a function on their relative
displacement in the slip direction. The element can be used in 1D, 2D, and 3D (default)
and its features are summarized in :numref:`Interface1dsummary`.

.. table:: IntElPoint element summary
   :name: Interface1dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | IntElPoint, Interface1d (deprecated)                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | One dimensional (slip) interface element                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`refnode{in}` :optelemparam:`normal{ra}`                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`refnode`: determines the reference node, which is used to specify a reference        |
   |                          | direction (the direction vector is obtained by subtracting the coordinates of the first node |
   |                          | from the reference node). :param:`normal`: The reference direction can be directly specified |
   |                          | by the optional parameter :param:`normal`. Although both :param:`refnode` and                |
   |                          | :param:`normal` are optional, at least one of them must be specified.                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | One, two, or three DOFs (u-displacement, v-displacement, w-displacement) are required in     |
   |                          | each node, according to element mode (determined from domain type).                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element requires material model with _1dInterface support.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/InterfaceEL_Point2D_01.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point2D_01.in>`_, |
   |                          | `sm/interface01.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/interface01.in>`_,            |
   |                          | `sm/InterfaceEL_Point3D_01.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_01.in>`_, |
   |                          | `sm/InterfaceEL_Point3D_02.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_02.in>`_, |
   |                          | `sm/InterfaceEL_Point3D_03.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Point3D_03.in>`_, |
   |                          | `sm/bondceb02.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bondceb02.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

Bondlink3dtruss element
~~~~~~~~~~~~~~~~~~~~~~~

This class implements a bond link for connecting trusses and continuum elements for
meshes in which the nodes of the trusses and continuum elements coincide. The bond area
for this element is calculated from the length and diameter of the truss element that is
bonded to the continuum as :math:`A_{\mathrm{bond}} = \pi d l`, where
:math:`A_{\mathrm{bond}}` is the bond area, :math:`d` is the reinforcement diameter and
:math:`l` is the length of the element. The input parameters for this element are shown
in :numref:`bondlink3dtrusssummary`.

.. table:: bondlink3dtruss element summary
   :name: bondlink3dtrusssummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | bondlink3dtruss                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d bond link element for trusses                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`length{rn}` :elemparam:`diameter{rn}` :elemparam:`dirvector{ra}`                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`length`: bond length :param:`diameter`: diameter of rebar modelled by truss element  |
   |                          | :param:`dirvector`: direction vector in which bond-slip occurs.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement) are      |
   |                          | required in each node.                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element is used together with the :param:`linkslip` bond-slip material.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/linkslip02.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip02.in>`_              |
   +--------------------------+----------------------------------------------------------------------------------------------+

IntElLine1, Interface2dlin elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of a two dimensional line element with a linear approximation of the
displacement jump.  The element can be used to tie together two element edges and is
defined by four nodes - two on each edge. Note that, the nodes along the interface are
doubled, each couple with identical coordinates. Nodes on the negative side are numbered
first, followed by nodes on the positive part. Requires material model with _2dInterface
support. The element features are summarized in :numref:`IntElLine1`.

.. _interf2d_lin_fig:

.. tikz:: Interface2dlin element with linear interpolation. Definition and node numbering convention

   \begin{tikzpicture}[scale=3,>=stealth]
    \tikzstyle{elemnode} = [draw,thin,circle,inner sep=1,fill=white]
    \newcommand{\eoff}{0.03}
    \draw[dashed, xshift=2.5, yshift=1]
       (0,0) -- (1,0) node[at start,elemnode] {}
       to node[at start,elemnode] {}
       (0.5,1) -- (0,0) node[at start,elemnode] {} ;

    \draw[dashed,xshift=12,yshift=4]
       (1,0) to node[at start,elemnode] {}
       (0.5,1) -- (1.5,1) node[at start,elemnode] {}
       to node[at start,elemnode] {} (1,0);

    \begin{scope}[xshift=7.5,yshift=2.5]
     \draw[<->] (0.8,1.2) -- (0.5,1.0) node[at start,left] {$\eta$}
       to (0.99,0) -- (1.1,-0.2) node[at end,right] {$\xi$};

     \draw[thick] (0.5,1)+(0:-\eoff) to
       node[at start,elemnode] (n1) {}
       node[at end,elemnode] (n2) {}
       (1-\eoff,-\eoff)
       (0.5,1)+(60:\eoff) to
       node[at start,elemnode] (n3) {}
       node[at end,elemnode] (n4) {}
       (1+\eoff,0);
      \node[yshift=3,left] at (n1) {1};
      \node[yshift=-4,left] at (n2) {2};
      \node[xshift=2,yshift=6] at (n3) {3}; % Fine adjustments
      \node[yshift=-4,right] at (n4) {4};

    \end{scope}
   \end{tikzpicture}

.. table:: IntElLine1 element summary
   :name: IntElLine1

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | IntElLine1, Interface2dlin(deprecated)                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D interface element with linear approximation                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`axisymmode`                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`axisymmode`: Flag controlling axisymmetric mode (integration over unit               |
   |                          | circumferential angle).                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using two point integration formula.               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element requires material model with _2dInterface support.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/InterfaceEL_Line1.in                                                                     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_Line1.in>`_,      |
   |                          | `sm/bondceb01.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bondceb01.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

IntElLine2, Interface2dquad elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of a two dimensional interface element with quadratic approximation of
displacement field. Can be used to glue together two elements with quadratic
displacement approximation along the shared edge. Note, that the nodes along the
interface are doubled, each couple with identical coordinates. Nodes on the negative
side are numbered first, followed by nodes on the positive part. Requires material model
with _2dInterface support. The element features are summarized in Table
:numref:`IntElLine2`.

.. _interf2d_quad_fig:

.. tikz:: Interface2dquad element with quadratic interpolation. Definition and node numbering convention

   \begin{tikzpicture}[scale=3,>=stealth]
    \tikzstyle{elemnode} = [draw,thin,circle,inner sep=1,fill=white]
    \newcommand{\eoff}{0.03}
    \draw[dashed]
       (0,0) -- (1,0) node[at start,elemnode] {} node[midway,elemnode] {}
       to[out=90,in=-30] node[at start,elemnode] {} node[midway,elemnode] {}
       (0.5,1) -- (0,0) node[at start,elemnode] {} node[midway,elemnode] {};

    \draw[dashed,xshift=15,yshift=5]
       (1,0) to[out=90,in=-30] node[at start,elemnode] {} node[midway,elemnode] {}
       (0.5,1) -- (1.5,1) node[at start,elemnode] {} node[midway,elemnode] {}
       to[out=-90,in=45] node[at start,elemnode] {} node[midway,elemnode] {} (1,0);

    \begin{scope}[xshift=7.5,yshift=2.5]
     \draw[<->] (0.5,1)+(60:0.3) -- (0.5,1) node[at start,left] {$\eta$}
       to[out=-30,in=90] (1,0) -- (1,-0.3) node[at end,right] {$\xi$};

     \draw[thick] (0.5,1)+(60:-\eoff) to[out=-30,in=90]
       node[at start,elemnode] (n1) {}
       node[at end,elemnode] (n2) {}
       node[midway,elemnode] (n3) {}
       (1-\eoff,0)
       (0.5,1)+(60:\eoff) to[out=-30,in=90]
       node[at start,elemnode] (n4) {}
       node[at end,elemnode] (n5) {}
       node[midway,elemnode] (n6) {}
       (1+\eoff,0);
      \node[left] at (n1) {1};
      \node[left] at (n2) {2};
      \node[left] at (n3) {3};
      \node[yshift=2,right] at (n4) {4}; % Fine adjustments
      \node[right] at (n5) {5};
      \node[right] at (n6) {6};
    \end{scope}
   \end{tikzpicture}

.. table:: IntElLine2 element summary
   :name: IntElLine2

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | IntElLine2, Interface2dquad (deprecated)                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D interface element with quadratic approximation                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`axisymmode`                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`axisymmode`: Flag controlling axisymmetric mode (integration over unit               |
   |                          | circumferential angle).                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all strain components using four point integration formula.              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element requires material model with _2dInterface support.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+

IntElSurfTr1, Interface3dtrlin elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of a three dimensional interface element with linear approximation of
displacement field. Can be used to glue together two elements with linear displacement
approximation along the shared triangular surface. Note, that the nodes along the
interface are doubled, each couple with identical coordinates. Nodes on the negative
surface are numbered first, followed by nodes on the positive part. The numbering of
surface nodes on positive surface (+) determines the positive normal (right hand rule).
Requires material model with _3dInterface support. The element features are summarized
in :numref:`Interface3dtrlinsummary`.

.. _interf3d_lin_fig:

.. tikz:: Interface3dtrlin element with linear interpolation. Definition and node numbering convention

   \begin{tikzpicture}[scale=5,>=stealth,
     x={(1cm,0cm)}, y={(-0.5cm,-0.5cm)}, z={(0cm,1cm)}]
    \tikzstyle{elemnode} = [solid,draw,thin,circle,inner sep=1,fill=white]

    \begin{scope}
    %\draw[->] (-0.05,0,0) -- (0.5,0,0) node[at end, below] {$x_g$};
    %\draw[->] (0,-0.05,0) -- (0,0.5,0) node[at end, below] {$y_g$};
    %\draw[->] (0,0,-0.05) -- (0,0,0.25) node[at end, right] {$z_g$};
    \end{scope}

    \begin{scope}[xshift=-2.5]
     \draw[dotted] (0,0,0) -- (-1,0,0) (0,0.5,1) -- (0,0,0) -- (0,1,0);
     \draw[dashed] (0,1,0) -- (-1,0,0) -- (0,0.5,1) -- cycle;
    \end{scope}

    \begin{scope}[xshift=5]
     \draw[dotted] (0,0,0) -- (0.5,0,0) (0,0.5,1) -- (0,0,0) -- (0,1,0);
     \draw[dashed] (0,1,0) -- (0.5,0,0) -- (0,0.5,1) -- cycle;
    \end{scope}

    \draw[fill=red,opacity=0.05] (0,0,0) -- (0,1,0) -- (0,0.5,1) -- cycle;
    \draw[thick,red] (0,0,0) node[elemnode] {}  node[left] {1}
     -- (0,1,0) node[elemnode] {} node[below] {3}
     -- (0,0.5,1) node[elemnode] {} node[above] {2}
     -- cycle node[draw,circle,inner sep=0] at (0,0.5,0.33) {\scriptsize $-$};

    \draw[fill=black,opacity=0.1,xshift=2.5] (0,0,0) -- (0,1,0) -- (0,0.5,1) -- cycle;
    \draw[thick,xshift=2.5] (0,0,0) node[elemnode] {} node[right] {4}
     -- (0,1,0) node[elemnode] {} node[below] {6}
     -- (0,0.5,1) node[elemnode] {} node[above] {5}
     -- cycle node[draw,circle,inner sep=0] at (0,0.5,0.33) {\scriptsize $+$};

   \end{tikzpicture}

.. table:: IntElSurfTr1 element summary
   :name: Interface3dtrlinsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | IntElSurfTr1, Interface3dtrlin (deprecated)                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D interface element with linear approximation                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`refnode`: determines the reference node, which is used to specify a reference        |
   |                          | direction (the direction vector is obtained by subtracting the coordinates of the first node |
   |                          | from the reference node).                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all components using one point integration formula.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element requires material model with _3dInterface support.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+

IntElSurfQuad1 element
~~~~~~~~~~~~~~~~~~~~~~

Implementation of a three dimensional interface element with linear approximation of
displacement field. Can be used to glue together two elements with linear displacement
approximation along the shared quad surface. Note, that the nodes along the interface
are doubled, each couple with identical coordinates. Nodes on the negative surface are
numbered first, followed by nodes on the positive part. The numbering of surface nodes
on positive surface (+) determines the positive normal (right hand rule). Requires
material model with _3dInterface support. The element features are summarized in Table
:numref:`Interface3dquadlinsummary`.

.. _interf3dquad_lin_fig:

.. tikz:: IntElSurfQuad1 element with linear interpolation. Definition and node numbering convention

   \begin{tikzpicture}[scale=5,>=stealth,
     x={(1cm,0cm)}, y={(-0.5cm,-0.5cm)}, z={(0cm,1cm)}]
    \tikzstyle{elemnode} = [solid,draw,thin,circle,inner sep=1,fill=white]

    \begin{scope}
    %\draw[->] (-0.05,0,0) -- (0.5,0,0) node[at end, below] {$x_g$};
    %\draw[->] (0,-0.05,0) -- (0,0.5,0) node[at end, below] {$y_g$};
    %\draw[->] (0,0,-0.05) -- (0,0,0.25) node[at end, right] {$z_g$};
    \end{scope}

    \begin{scope}[xshift=-2.5]
     \draw[dotted] (0,0,0) -- (-1,0,0) (0,0,1) -- (0,0,0) -- (0,1,0) ;
     \draw[dotted] (-1,0,1) -- (-1,0,0) -- (-1,1,0) ;
     \draw[dashed] (-1,1,0) -- (0,1,0) -- (0,1,1) -- (-1,1,1) -- cycle;
     \draw[dashed] (-1,1,1) -- (-1,0,1) -- (0,0,1) -- (0,1,1) -- cycle;
    \end{scope}

    \begin{scope}[xshift=5]
     \draw[dotted] (0,0,0) -- (1,0,0) (0,0,1) -- (0,0,0) -- (0,1,0);
     \draw[dotted] (1,0,1) -- (1,0,0) -- (1,1,0) ;
     \draw[dashed] (1,1,0) -- (0,1,0) -- (0,1,1) -- (1,1,1) -- cycle;
     \draw[dashed] (1,1,1) -- (1,0,1) -- (0,0,1) -- (0,1,1) -- cycle;
    \end{scope}

    \draw[fill=red,opacity=0.05] (0,0,0) -- (0,1,0) -- (0,1,1) -- (0,1,1) -- (0,0,1) --cycle;
    \draw[thick,red] (0,0,0) node[elemnode] {}  node[left] {1}
     -- (0,1,0) node[elemnode] {} node[below] {4}
     -- (0,1,1) node[elemnode] {} node[above] {3}
     -- (0,0,1) node[elemnode] {} node[above] {2}
     -- cycle node[draw,circle,inner sep=0] at (0,0.5,0.5) {\scriptsize $-$};

    \draw[fill=black,opacity=0.1,xshift=2.5] (0,0,0) -- (0,1,0) -- (0,1,1) -- (0,0,1) -- cycle;
    \draw[thick,xshift=2.5] (0,0,0) node[elemnode] {} node[right] {5}
     -- (0,1,0) node[elemnode] {} node[below] {8}
     -- (0,1,1) node[elemnode] {} node[above] {7}
     -- (0,0,1) node[elemnode] {} node[above] {6}
     -- cycle node[draw,circle,inner sep=0] at (0,0.5,0.5) {\scriptsize $+$};

   \end{tikzpicture}

.. table:: IntElSurfQuad1 element summary
   :name: Interface3dquadlinsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | IntElSurfQuad1                                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D interface element with linear approximation                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`refnode`: determines the reference node, which is used to specify a reference        |
   |                          | direction (the direction vector is obtained by subtracting the coordinates of the first node |
   |                          | from the reference node).                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration of all components using one point integration formula.                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Element requires material model with _3dInterface support.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/InterfaceEL_SurfQuad2.in                                                                 |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/InterfaceEL_SurfQuad2.in>`_   |
   +--------------------------+----------------------------------------------------------------------------------------------+

Bondlink3d element
~~~~~~~~~~~~~~~~~~

This class implements a bond link for connecting beam (frame) and continuum elements in
unstructured meshes. The node ordering is important: the first node carries the rigid
arm and must therefore have rotational degrees of freedom, while the slip (displacement
jump) is evaluated at the second node. For continuum (matrix) element nodes that have no
rotational degrees of freedom, the reinforcement (beam/rebar) node must be given first
and the continuum node second. The slip is then evaluated at the matrix. The main idea
is to use the rotation of the beam element and the rigid arm from the beam node to the
continuum element node to compute the displacement jump at the matrix node (and two
components, which are perpendicular to each other and lie in a plane for which the
direction along the rebar is normal to. This element differs from the other interface
elements since the first stiffness component is in the direction of the slip. The bond
area for this element is calculated from the length and diameter of the truss element
that is bonded to the continuum as :math:`A_{\mathrm{bond}} = \pi d l`, where
:math:`A_{\mathrm{bond}}` is the bond area, :math:`d` is the reinforcement diameter and
:math:`l` is the length of the element. Each node has six degrees of freedom. The input
parameters for this element are shown in :numref:`bondlink3dsummary`.

.. table:: bondlink3d element summary
   :name: bondlink3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | bondlink3d                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3d bond link element for beams                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`length{rn}` :elemparam:`diameter{rn}` :elemparam:`dirvector{ra}`                 |
   |                          | :elemparam:`l_end{rn}`                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`length`: bond length :param:`diameter`: diameter :param:`dirvector`: direction       |
   |                          | vector in which bond-slip occurs.                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement,            |
   |                          | :math:`u`-rotation, :math:`v`-rotation and :math:`w`-rotation) are required in the beam node |
   |                          | and three dofs (:math:`u`-displacement, :math:`v`-displacement, :math:`w`-displacement) in   |
   |                          | the continnum node. The order of input of nodes must be beam node first and continuum node   |
   |                          | second.                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [SciGraLarRun20]_                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/bond_link_3.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_3.in>`_,            |
   |                          | `sm/bond_link_2.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_2.in>`_,            |
   |                          | `sm/linkslip01.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip01.in>`_,             |
   |                          | `sm/bond_link_1.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/bond_link_1.in>`_             |
   +--------------------------+----------------------------------------------------------------------------------------------+

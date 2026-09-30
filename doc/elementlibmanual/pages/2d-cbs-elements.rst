2D CBS Elements
===============


2D CBS Elements
---------------

.. _Tr1CBS:

Tr1CBS element
~~~~~~~~~~~~~~

Represents the linear triangular finite element for transient incompressible flow
analysis using cbs algorithm with equal order approximation of velocity and pressure
fields. Each node has 3 degrees of freedoms (two components of velocity and pressure).
The node numbering is anti-clockwise. The element features are summarized in Table
:numref:`Tr1CBSsummary`.

.. _Tr1CBSfig:

.. tikz:: Tr1CBS element. Node numbering, Side numbering and
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

.. table:: Tr1CBS element summary
   :name: Tr1CBSsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr1CBS                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | linear triangular finite element for transient incompressible flow analysis using cbs        |
   |                          | algorithm                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`bsides{ia}` :optelemparam:`bcodes{ia}`                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | Since the problem formulation requires to evaluate some boundary terms, the element boundary |
   |                          | edges should be specified as well as the types of boundary conditions applied at these       |
   |                          | boundary edges. The boundary edges (their numbers) are specified using :param:`bsides`       |
   |                          | array. The type of boundary condition(s) applied to corresponding boundary side is           |
   |                          | determined by :param:`bcodes` array. The available/supported boundary codes are following: 1 |
   |                          | for prescribed traction, 2 for prescribed normal velocity, 4 for prescribed tangential       |
   |                          | velocity, and 8 for prescribed pressure. If the element side is subjected to a combination   |
   |                          | of these fundamental types boundary conditions, the corresponding code is obtained by        |
   |                          | summing up the corresponding codes.                                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two velocity components (V_u and V_v) and pressure (P_f) are required in each node.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Equal order approximation of velocity and pressure fields.                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | exact                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Constant boundary tractions are supported\footnoteIn CBS algorithm formulation the           |
   |                          | prescribed traction boundary condition leads indirectly to pressure boundary condition in    |
   |                          | nodes associated to loaded edge. Such boundary condition is represented by                   |
   |                          | PrescribedTractionPressureBC. See section on boundary conditions in OOFEM input manual..     |
   |                          | Body loads representing the self-weight load are supported.                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Untested                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `fm/cbs3.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs3.in>`_,       |
   |                          | `fm/cbs1.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs1.in>`_,       |
   |                          | `fm/cbs2.in <https://github.com/oofem/oofem/blob/devel/tests/regression/fm/cbs2.in>`_        |
   +--------------------------+----------------------------------------------------------------------------------------------+

Truss Elements
==============


Truss Elements
--------------

.. _Truss1d:

Truss 1D element
~~~~~~~~~~~~~~~~

Represents linear isoparametric truss element in 1D. The elements are assumed to be
located along the x-axis. Requires cross section area to be specified. The element
features are summarized in :numref:`truss1dsummary`.

.. table:: truss1d element summary
   :name: truss1dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | truss1d                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 1D truss element                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (u-displacement) is required in each node                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacement and geometry                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Full dynamic analysis support, Full nonlocal constitutive support, Adaptivity support        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Area is required                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are not supported in current implementation         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/con2dpm8.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm8.in>`_,               |
   |                          | `sm/con2dpm5.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/con2dpm5.in>`_,               |
   |                          | `sm/isoasymm01.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/isoasymm01.in>`_,             |
   |                          | `partests/dyn_bar03/dyn_bar03.oofem.in <https://github.com/oofem/oofem/blob/devel/tests/regr |
   |                          | ession/partests/dyn_bar03/dyn_bar03.oofem.in>`_, `partests/dyn_bar03/bar03.statics.oofem.in  |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/partests/dyn_bar03/bar03.statics |
   |                          | .oofem.in>`_, `sm/Mises02.in                                                                 |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/Mises02.in>`_, `sm/idm02.in   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm02.in>`_, `sm/idm03.in     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm03.in>`_, `sm/idm04.in     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm04.in>`_, `sm/idm05.in     |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm05.in>`_ (and 9 more)      |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _Truss2d:

Truss 2D element
~~~~~~~~~~~~~~~~

Two node linear isoparametric truss element for 2D analysis. The element geometry can be
specified in (x,z), (x,y), or (y,z) plane. The element features are summarized in Table
:numref:`truss2dsummary`.

.. tikz:: Truss2d element in (x,z) plane.

   \begin{tikzpicture}[scale=7,>=stealth]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \newcommand{\trusslength}{0.5};

    \draw[->] (-0.05,0) -- (0.7,0) node[above left,at end] {$x_g$};
    \draw[->] (0,-0.05) -- (0,0.5) node[below right,at end] {$z_g$};
    \draw[very thick] (0.1,0.1) --  +(30:\trusslength)
       node[elemnode,at start] {} node[at start,yshift=2,above left] {1}
       node[elemnode,at end] {} node[at end,yshift=2,above left] {2};
    %\draw[dotted,->] (0.1,0.1)++(-30:0.42) -- +(-30:0.1) node[below] {$X_1$};
    %\draw[dotted,->] (a)++(-30-90:0.02) -- +(-120:0.1) node[right] {$Z_1$};
   \end{tikzpicture}

.. table:: truss2d element summary
   :name: truss2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | truss2d                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D truss element                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`cs{in}`                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`cs`: this parameter can be used to change default definition plane. The supported    |
   |                          | values of :param:`cs` are following: 0 for (x,z) plane (default), 1 for (x,y) plane, and 3   |
   |                          | for (y,z) plane.                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs representing displacements in definition plane are required in each node. The       |
   |                          | element can be used in different planes, default definition plane is (x,z). The parameter    |
   |                          | :param:`cs` can be used to change default definition plane.                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Full dynamic analysis support. Full nonlocal constitutive support.                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | cross section area should be provided.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Edge loads are supported, Edge number should be equal to 1                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/linear_constraint_1.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linear_constraint_1.in>`_,    |
   |                          | `sm/linear_constraint_2.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linear_constraint_2.in>`_,    |
   |                          | `sm/linear_constraint_3.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linear_constraint_3.in>`_,    |
   |                          | `sm/patch160.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch160.in>`_,               |
   |                          | `partests/dyn_bar01/dyn_bar01.oofem.in <https://github.com/oofem/oofem/blob/devel/tests/regr |
   |                          | ession/partests/dyn_bar01/dyn_bar01.oofem.in>`_, `sm/trussb3_creep.in                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/trussb3_creep.in>`_,          |
   |                          | `sm/trussb3_relax.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/trussb3_relax.in>`_,          |
   |                          | `sm/steelRelaxMat2.in                                                                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/steelRelaxMat2.in>`_,         |
   |                          | `sm/rotated_1.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rotated_1.in>`_,              |
   |                          | `sm/rotated_2.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/rotated_2.in>`_ (and 7 more)  |
   +--------------------------+----------------------------------------------------------------------------------------------+

truss2d element summarytruss2dsummary

.. _Truss3d:

Truss 3D element
~~~~~~~~~~~~~~~~

Two node linear isoparametric truss element for 3D analysis. The element geometry is
specified in (x,y,z) space. The element features are summarized in Table
:numref:`truss3dsummary`.

.. tikz:: Truss3d element in (x,y,z) space.

   \begin{tikzpicture}[scale=7,>=stealth]
    \tikzstyle{elemnode} = [draw=black,thick,fill=white,circle,inner sep=1]
    \newcommand{\trusslength}{0.5};

    \begin{scope}
    \draw[->] (-0.05,0,0) -- (0.5,0,0) node[at end, below] {$x_g$};
    \draw[->] (0,-0.05,0) -- (0,0.5,0) node[at end, below right] {$y_g$};
    \draw[->] (0,0,-0.05) -- (0,0,0.5) node[at end, right] {$z_g$};
    \end{scope}

    \draw[very thick] (0.1,0.1) --  +(30:\trusslength)
       node[elemnode,at start] {} node[at start,yshift=2,above left] {1}
       node[elemnode,at end] {} node[at end,yshift=2,above left] {2};
    %\draw[dotted,->] (0.1,0.1)++(-30:0.42) -- +(-30:0.1) node[below] {$X_1$};
    %\draw[dotted,->] (a)++(-30-90:0.02) -- +(-120:0.1) node[right] {$Z_1$};
   \end{tikzpicture}

.. table:: truss3d element summary
   :name: truss3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | truss3d                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D truss element                                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three displacement DOFs (in x, y, and z directions) are required in each node.               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Exact.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Full dynamic analysis support. Full nonlocal constitutive support.                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | cross section area should be provided.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/nldeidynamic1.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/nldeidynamic1.in>`_,          |
   |                          | `sm/compoDamMat.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/compoDamMat.in>`_,            |
   |                          | `sm/linkslip02.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/linkslip02.in>`_              |
   +--------------------------+----------------------------------------------------------------------------------------------+

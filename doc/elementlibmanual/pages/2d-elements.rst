2D Elements
===========


2D Elements
-----------

.. _Tr1ht:

Tr1ht element
~~~~~~~~~~~~~

Implements the linear triangular finite element for heat transfer problems. Each node
has 1 degree of freedom. The cross section thickness property is requested form cross
section model. The node numbering is anti-clockwise. The element features are summarized
in :numref:`Tr1htsummary`.

.. _Tr1htfig:

.. tikz:: Tr1ht element - node and side numbering.

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

.. table:: Tr1ht element summary
   :name: Tr1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Tr1ht                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | triangular finite element with linear approximation for heat transfer problems               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (T_f - temperature) is required in each node.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of temperature.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using one point gauss integration formula.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and are computed using numerical      |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 3). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`Tr1htfig`)).                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tm/hydratingConcreteMat04.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat04.in>`_, |
   |                          | `tm/hydratingConcreteMat03.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat03.in>`_, |
   |                          | `tm/tmpatch45-1.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch45-1.in>`_,            |
   |                          | `tm/hydratingConcreteMat01.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat01.in>`_, |
   |                          | `tm/hydratingConcreteMat06.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat06.in>`_, |
   |                          | `tm/hydratingConcreteMat02.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat02.in>`_, |
   |                          | `tm/hydratingConcreteMat07.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat07.in>`_, |
   |                          | `tm/hydratingConcreteMat05.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat05.in>`_, |
   |                          | `tm/tmpatch43-1.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch43-1.in>`_,            |
   |                          | `tm/tmpatch43-3.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch43-3.in>`_ (and 1      |
   |                          | more)                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+

Tr1mt element
~~~~~~~~~~~~~

Isoparametric triangular finite element with linear approximation of moisture. Other
features are the same as for Tr1ht in Section :ref:`Tr1ht`.

Tr1hmt element
~~~~~~~~~~~~~~

Isoparametric triangular finite element with linear approximations of  temperature and
moisture.  Other features are the same as for Tr1ht in Section :ref:`Tr1ht`.

Tests/Examples: `tests/regression/tm/HeMoKunzel_1.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/tm/HeMoKunzel_1.in>`_

.. _Quad1ht:

Quad1ht element
~~~~~~~~~~~~~~~

Represents isoparametric four-node quadrilateral finite element for heat transfer
problems. Each node has 1 degree of freedom. Problem should be defined in x,y plane. The
cross section thickness property is requested form cross section model. The nodes should
be numbered anti-clockwise (positive rotation around z-axis). The element features are
summarized in :numref:`Quad1htsummary`.

.. _Quad1htfig:

.. tikz:: Quad1ht element. Node numbering, Side numbering and
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

.. table:: Quad1ht element summary
   :name: Quad1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Quad1ht                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Isoparametric four-node quadrilateral linear interpolation element for heat transfer         |
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
   | Integration              | Integration using gauss integration formula in 4 (the default), 9, or 16 integration points. |
   |                          | The default number of integration point used can be overloaded using :param:`NIP` parameter. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 4). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`Quad1htfig`)).                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tmcemhyd/cemhyd01.in                                                                        |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tmcemhyd/cemhyd01.in>`_,         |
   |                          | `tm/tmpatch02.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch02.in>`_,              |
   |                          | `tm/tmpatch01.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch01.in>`_,              |
   |                          | `tm/tmpatch06.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch06.in>`_,              |
   |                          | `tm/tmpatch11.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch11.in>`_,              |
   |                          | `tm/tmpatch11dtf.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch11dtf.in>`_,           |
   |                          | `tm/tmpatch04.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch04.in>`_,              |
   |                          | `tm/tmpatch05.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch05.in>`_,              |
   |                          | `tm/tmpatch17.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmpatch17.in>`_,              |
   |                          | `tm/tmquad12.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/tmquad12.in>`_ (and 6 more)   |
   +--------------------------+----------------------------------------------------------------------------------------------+

Quad1mt element
~~~~~~~~~~~~~~~

Isoparametric four-node quadrilateral finite element. Other features are the same as for
Quad1ht in Section :ref:`Quad1ht`.

Tests/Examples: `tests/regression/tm/nlisomoisture01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/tm/nlisomoisture01.in>`_,
`tests/regression/tm/nlisomoisture02.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/tm/nlisomoisture02.in>`_

Quad1hmt element
~~~~~~~~~~~~~~~~

Represents isoparametric four-node quadrilateral finite element for heat and mass (one
constituent) transport problems.  Two dofs (T_f - temperature and C_1 - concentration)
are required in each node. Linear approximation of temperature and mass concentration.
Other features are similar to Quad1 element, see section :ref:`Quad1ht`.

.. _QQuad1ht:

QQuad1ht element
~~~~~~~~~~~~~~~~

Represents isoparametric quadratic eight-node quadrilateral finite element for heat
transfer problems. Each node has 1 degree of freedom. Problem should be defined in x,y
plane. The cross section thickness property is requested form the cross section model.
The nodes should be numbered anti-clockwise (positive rotation around z-axis), see fig.
:ref:`qplanstrssfig`. The element has the same features as in Table
:numref:`Quad1htsummary`.

Tests/Examples: `tests/regression/tm/qquad01.in
<https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qquad01.in>`_

QQuad1mt element
~~~~~~~~~~~~~~~~

Element for mass transport problems, see the parent element in sec. :ref:`QQuad1ht`.

QQuad1hmt element
~~~~~~~~~~~~~~~~~

Element for heat and mass transport problems, see the parent element in sec.
:ref:`QQuad1ht`.

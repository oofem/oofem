1D Elements
===========


1D Elements
-----------

.. _Line1ht:

Line1ht element
~~~~~~~~~~~~~~~

Two node linear isoparametric element for for heat transfer problems. Each node has 1
degree of freedom. The cross section property “area” is requested from the cross
section. The element geometry is specified in (x,y,z) plane. The element features are
summarized in :numref:`Line1htsummary`. Stabilization through lumped capacity
matrix is suggested in transient transport problems, using **lumped**.

.. tikz:: Line1ht element in (x,y,z) space.

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

.. table:: Line1ht element summary
   :name: Line1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Line1ht                                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Line (truss-like) finite element with linear approximation for heat transfer problems        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used.                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (T_f - temperature) is required in each node.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of temperature.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using one point gauss integration formula.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are unsupported.                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Area                                                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tm/line01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/line01.in>`_,   |
   |                          | `tm/line02.in <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/line02.in>`_,   |
   |                          | `tm/hydratingConcreteMat07.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/hydratingConcreteMat07.in>`_  |
   +--------------------------+----------------------------------------------------------------------------------------------+

Line1mt element
~~~~~~~~~~~~~~~

Line (truss-like) finite element with linear approximation of moisture. Other features
are the same as for Line1ht in Section :ref:`Line1ht`.

Line1hmt element
~~~~~~~~~~~~~~~~

Line (truss-like) finite element with linear approximations of temperature and moisture.
Other features are the same as for Tr1ht in Section :ref:`Line1ht`.

Plane Stress Elements
=====================


Plane Stress Elements
---------------------

PlaneStress2d
~~~~~~~~~~~~~

Represents isoparametric four-node quadrilateral plane-stress finite element. Each node
has 2 degrees of freedom. Structure should be defined in x,y plane.  The nodes should be
numbered anti-clockwise (positive rotation around z-axis). The element features are
summarized in :numref:`planestress2dsummary`.

The generalization of this element, that can be positioned arbitrarily in space is
:param:`linquad3dplanestress` element. This element requires 3 displacement degrees of
freedon in each node and assumes, that the element geometry is flat, i.e. all nodes are
in the same plane. The element features are summarized in Table
:numref:`linquad3dplanestresssummary`.

.. _Planestress2dfig:

.. tikz:: PlaneStress2d element. Node numbering, edge numbering and definition of local edge c.s.(a).

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

.. table:: planestress2d element summary
   :name: planestress2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | planestress2d                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D quadrilateral element for plane stress analysis                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using Gauss integration formula in 1, 4 (default), 9 or |
   |                          | 16 integration points. The default number of integration points used can be overloaded using |
   |                          | :param:`NIP` parameter. Reduced integration for shear terms is employed. Shear terms are     |
   |                          | always integrated using the 1-point integration rule.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Nonlocal constitutive support, Geometric nonlinearity support.                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 4). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`Planestress2dfig`)).                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0, 1.                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/idm08.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm08.in>`_,     |
   |                          | `sm/idm11.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/idm11.in>`_,     |
   |                          | `sm/materOrient01.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/materOrient01.in>`_,          |
   |                          | `sm/concrete_fcm_visco.in                                                                    |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/concrete_fcm_visco.in>`_,     |
   |                          | `sm/EC2creep_casting.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/EC2creep_casting.in>`_,       |
   |                          | `sm/control_switch_1.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/control_switch_1.in>`_,       |
   |                          | `sm/control_switch_2.in                                                                      |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/control_switch_2.in>`_,       |
   |                          | `sm/hangingnode01.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/hangingnode01.in>`_,          |
   |                          | `sm/hangingnode02.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/hangingnode02.in>`_,          |
   |                          | `sm/hangingnode03.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/hangingnode03.in>`_ (and 23   |
   |                          | more)                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: linquad3dplanestress element summary
   :name: linquad3dplanestresssummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | linquad3dplanestress                                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D quadrilateral element for plane stress analysis                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, w-displacement) are required in each node.       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using Gauss integration formula in 1, 4 (default), 9 or |
   |                          | 16 integration points. The default number of integration points used can be overloaded using |
   |                          | :param:`NIP` parameter. Reduced integration for shear terms is employed. Shear terms are     |
   |                          | always integrated using the 1-point integration rule.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Nonlocal constitutive support, Geometric nonlinearity support.                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and computed using numerical          |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 4). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`Planestress2dfig`)).                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0, 1.                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Basic functionality tested, element loads need further testing.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch303.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch303.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

QPlaneStress2d
~~~~~~~~~~~~~~

Implementation of quadratic isoparametric eight-node quadrilateral plane-stress finite
element. Each node has 2 degrees of freedom. The node numbering is anti-clockwise and is
explained in fig. (:ref:`qplanstrssfig`). The element features are summarized in
:numref:`qplanestress2dsummary`.

.. _qplanstrssfig:

.. tikz:: QPlaneStress2d element - node numbering.

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
        (0.2,0.1) node[elemnode] {} node[below] {1}
     to[out=-10,in=210] coordinate[midway] (e1) node[elemnode,midway] {} node[below,midway] {5} (0.7,0.15) node[elemnode] {} node[below right] {2}
     to[out=110,in=-100]  coordinate[midway] (e2)  node[elemnode,midway] {} node[right,midway] {6} (0.7,0.4) node[elemnode] {} node[above] {3}
     to[out=190,in=-10] coordinate[midway] (e3) node[elemnode,midway] {} node[above,midway] {7} (0.2,0.4) node[elemnode] {} node[above] {4}
     to[out=-100,in=100] coordinate[midway] (e4) node[elemnode,midway] {} node[left,midway] {8} (0.2,0.1);

    \lcoordsys{  7}{(e1)};
    \lcoordsys{ 90}{(e2)};
    \lcoordsys{180}{(e3)};
    \lcoordsys{275}{(e4)};
   \end{tikzpicture}

.. table:: qplanestress2d element summary
   :name: qplanestress2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | qplanestress2d                                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D quadratic isoparametric plane stress element                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration using Gauss integration formula in 4 (the default), 9 or 16 integration     |
   |                          | points. The default number of integration points used can be overloaded using :param:`NIP`   |
   |                          | parameter.                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Adaptivity support.                                                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body and boundary loads are supported.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0, 1.                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Stable                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/deadweight01.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/deadweight01.in>`_            |
   +--------------------------+----------------------------------------------------------------------------------------------+

TrPlaneStress2d
~~~~~~~~~~~~~~~

Implements an triangular three-node constant strain plane-stress  finite element. Each
node has 2 degrees of freedom. The node numbering is anti-clockwise. The element
features are summarized in :numref:`trplanestress2dsummary`.

.. _TrPlanestressfig:

.. tikz:: TrPlaneStress2d element - node and side numbering.

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

.. table:: trplanestress2d element summary
   :name: trplanestress2dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestress2d                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular isoparametric plane stress element                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using one point gauss integration formula.              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Nonlocal constitutive support, Edge load support, Geometric nonlinearity support, Adaptivity |
   |                          | support.                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported and are computed using numerical      |
   |                          | integration. The side numbering is following. Each i-th element side begins in i-th element  |
   |                          | node and ends on next element node (i+1-th node or 1-st node, in the case of side number 3). |
   |                          | The local positive edge x-axis coincides with side direction, the positive local edge y-axis |
   |                          | is rotated 90 degrees anti-clockwise (see fig. (:ref:`TrPlanestressfig`)).                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0, 1.                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/deadweight03.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/deadweight03.in>`_,           |
   |                          | `sm/adapt01.in <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/adapt01.in>`_, |
   |                          | `sm/patch104.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch104.in>`_,               |
   |                          | `sm/setprops01.in                                                                            |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/setprops01.in>`_,             |
   |                          | `sm/patch130.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch130.in>`_,               |
   |                          | `sm/patch102.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch102.in>`_,               |
   |                          | `sm/distancebasedaveraging.in                                                                |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/distancebasedaveraging.in>`_, |
   |                          | `sm/stressbasedaveraging.in                                                                  |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/stressbasedaveraging.in>`_,   |
   |                          | `sm/layeredcs01.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/layeredcs01.in>`_             |
   +--------------------------+----------------------------------------------------------------------------------------------+

QTrPlStr
~~~~~~~~

Implementation of quadratic six-node plane-stress finite element. Each node has 2
degrees of freedom. Node numbering is anti-clockwise and is shown in fig.
(:ref:`qtrplanstressfig`). The element features are summarized in Table
:numref:`qtrplstrsummary`.

.. _qtrplanstressfig:

.. tikz:: QTrPlStr element - node and side numbering.

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

.. table:: qtrplstr element summary
   :name: qtrplstrsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | qtrplstr                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D quadratic triangular plane stress element                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacements and geometry.                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Full integration using gauss integration formula in 4 points (the default) or in 7 points    |
   |                          | (using :param:`NIP` parameter).                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Adaptivity support (error indicator).                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Boundary loads are supported.                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0, 1.                                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch140.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch140.in>`_ (and 1 more)   |
   +--------------------------+----------------------------------------------------------------------------------------------+

TrPlaneStrRot
~~~~~~~~~~~~~

Implementation of triangular three-node plane-stress  finite element with independent
rotation field. Each node has 3 degrees of freedom. The element features are summarized
in :numref:`trplanestrrotsummary`.

The generalization of this element, that can be positioned arbitrarily in space is
:param:`trplanestrrot3d` element. This element requires 6 degrees of freedon in each
node. The element features are summarized in :numref:`trplanestrrot3dsummary`.

The implementation is based on the following paper: Ibrahimbegovic, A., Taylor, R.L.,
Wilson, E. L.: A robust membrane qudritelar element with rotational degrees of freedom,
Int. J. Num. Meth. Engng., 30, 445-457, 1990. The rotation field is defined as
:math:`\omega = \del{1}{2}(\der{v}{x}-\der{u}{y}) = \nabla_u\mbf{u}`. The following form
of potential energy functial is assumed:

.. math::

   \Pi = \del{1}{2}\int_{\Omega}\mbf{\sigma}^T\mbf{\varepsilon}\ d\Omega + \int_{\Omega}\mbf{\tau}^T(\nabla_u\mbf{u}-\omega)\ d\Omega - \int_\Omega \mbf{X}^T\mbf{u}\ d\Omega

where :math:`\mbf{\tau}` is pseudo-stress (component of anti-symmetric stress tensor)
working on dislocation :math:`(\nabla_u\mbf{u}-\omega)`; the following constitutive
relation foris assumed: :math:`\mbf{\tau} = G(\nabla_u\mbf{u}-\omega)`, where :math:`G`
is elasticity modulus in shear.

.. table:: trplanestrrot element summary
   :name: trplanestrrotsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestrrot                                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular plane stress element with rotational DOFs                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`NIPRot{in}`                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms. :param:`NIPRot`: allows to set the number of integration points for integration of    |
   |                          | terms associated to rotational field.                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, z-rotation) are required in each node.           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using gauss integration formula in 4 points (default)   |
   |                          | or using 1 or 7 points (using :param:`NIP` parameter). Integration of strains associated     |
   |                          | with rotational field integration using 1 point is default (4 and 7 points rules can be      |
   |                          | specified using :param:`NIPRot` parameter).                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch150.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch150.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: trplanestrrot3d element summary
   :name: trplanestrrot3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestrrot3d                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 3D linear triangular plane stress element with rotational DOFs                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`NIPRot{in}`                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms. :param:`NIPRot`: allows to set the number of integration points for integration of    |
   |                          | terms associated to rotational field.                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u-displacement, v-displacement, w-displacement, x-rotation, y-rotation,            |
   |                          | z-rotation) are required in each node.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and geometry.                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using gauss integration formula in 4 points (default)   |
   |                          | or using 1 or 7 points (using :param:`NIP` parameter). Integration of strains associated     |
   |                          | with rotational field integration using 1 point is default (4 and 7 points rules can be      |
   |                          | specified using :param:`NIPRot` parameter).                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+

TrPlaneStressRotAllman
~~~~~~~~~~~~~~~~~~~~~~

Implementation of triangular three-node plane-stress  with nodal rotations. Each node
has 3 degrees of freedom. The element features are summarized in Table
:numref:`TrPlaneStressRotAllmansummary`.

The generalization of this element, that can be positioned arbitrarily in space is
:param:`trplanestressrotallman3d` element. This element requires 6 degrees of freedon in
each node. The element features are summarized in Table
:numref:`trplanestressrotallman3dsummary`.

The implementation is based on the following paper: Allman, D.J.: A compatible
triangular element including vertex rotations for plane elasticity analysis, Computers &
Structures, vol. 19, no. 1-2, pp. 1-8, 1984. The element is based on plane stress
element with quadratic interpolation. The displacements in midside nodes are expressed
using vertex displacements and vertex rotations (for edge normal displacement
component);  the tangential component is interpolated from vertex values. For particular
element side starting at i-th vertex and ending in j-th vertex the normal and tangential
displacements at edge midpoint can be expressed as

.. math::

   u_n\vert_{l/2} &=& \del{u_{ni}+u_{nj}}{2}+\del{l}{8}(\omega_i-\omega_j)\\
   u_t\vert_{l/2} &=& \del{u_{ti}+u_{tj}}{2}

where :math:`l` is edge length. This allows to express global displacements in element
midside nodes using vertex displacements and rotations. For a single edge, one obtains:

.. math::

   u\vert_{l/2} &=&-\del{u_{ni}+u_{nj}}{2}+\del{l}{8}(\omega_i-\omega_j)\del{\Delta y_{ji}}{l}+(\del{u_{t1}+u_{t2}}{2})\del{\Delta x_{ji}}{l}\\
   v\vert_{l/2} &=& \del{u_{ni}+u_{nj}}{2}+\del{l}{8}(\omega_i-\omega_j)\del{\Delta x_{ji}}{l}+(\del{u_{t1}+u_{t2}}{2})\del{\Delta y_{ji}}{l}\\

.. table:: trplanestressrotallman element summary
   :name: TrPlaneStressRotAllmansummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestressrotallman                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular plane stress element with rotational DOFs                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (u-displacement, v-displacement, z-rotation) are required in each node.           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of geometry, quadratic interpolation of displacements.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using gauss integration formula in 4 points.            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Zero energy mode         | The zero energy mode (equal rotations) is handled by adding additional energy term           |
   |                          | preventing spurious modes.                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch108.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch108.in>`_, `sm/beam44.in |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/beam44.in>`_                  |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: trplanestressrotallman3d element summary
   :name: trplanestressrotallman3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | trplanestressrotallman3d                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D linear triangular plane stress element with rotational DOFs                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (D_u, D_v, D_w, R_x, R_y, R_z) are required in each node.                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of geometry, quadratic interpolation of displacements.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of membrane strain terms using gauss integration formula in 4 points.            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Zero energy mode         | The zero energy mode (equal rotations) is handled by adding additional energy term           |
   |                          | preventing spurious modes.                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+

Plate & Shell Elements
======================


Plate & Shell Elements
----------------------

\subsubsection DKT Element  Implementation of Discrete Kirchhoff Triangle (DKT) plate
element. This element is suitable for thin plates, as the traswerse shear strain energy
is neglected. The structure should be defined in x,y plane, nodes should be numbered
anti-clockwise (positive rotation around z-axis). The element features are summarized in
:numref:`dktplatesummary`.

.. table:: DKTplate element summary
   :name: dktplatesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | dktplate                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D Discrete Kirchhoff Triangular plate element                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (w-displacement, u and v - rotations) are required in each node.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of rotations, cubic approximation of displacement along the edges.   |
   |                          | Note: there is no need to define interpolation for displacement on the element.              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Default integration of all terms using three point formula.                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary load support is beta.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the generalized shell strain/force momentum vectors in global coordinate system   |
   |                          | are printed, with the following meaning: :math:`s_{\varepsilon}  = \{\varepsilon_x,          |
   |                          | \varepsilon_y, \varepsilon_{xz}, \kappa_x, \kappa_y, \kappa_{xy}, \gamma_{xz},               |
   |                          | \gamma_{yz}\},\; s_{\sigma}  = \{n_x, n_y, n_{xy}, m_x, m_y, m_z, m_{xy}, q_{xz}, q_{yz}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_{xy}` are membrane in plane normal    |
   |                          | deformations, :math:`\gamma_{zx}, \gamma_{xz}` are (out of plane and in plane) shear         |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_{xy}` are curvatures, :math:`n_x, n_y, n_{xy},  |
   |                          | q_{xz}, q_{yz}` are integral forces (normal and shear forces), and :math:`m_x, m_y, m_{xy}`  |
   |                          | are bending moments. Please note, for example, that bending moment :math:`m_x` is defined as |
   |                          | :math:`m_x=\int \sigma_x z\ dz`, so it acts along the y-axis and positive value causes       |
   |                          | tension in bottom layer.                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | J.L.Batoz, K.J.Bathe, L.W.Ho: A study of three-node triangular plate bending elements,       |
   |                          | IJNME, 15(12):1771-1812, 1980                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/dkt_twist01.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/dkt_twist01.in>`_,            |
   |                          | `sm/dkt_bending01.in                                                                         |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/dkt_bending01.in>`_,          |
   |                          | `sm/dkt_twist02.in                                                                           |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/dkt_twist02.in>`_,            |
   |                          | `sm/dkt_rhombic_cantilever_4x4.in <https://github.com/oofem/oofem/blob/devel/tests/regressio |
   |                          | n/sm/dkt_rhombic_cantilever_4x4.in>`_                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection QDKT Element  Implementation of Discrete Kirchhoff Theory plate quad
element (QDKT). This element is suitable for thin plates, as the traswerse shear strain
energy is neglected. The structure should be defined in x,y plane, nodes should be
numbered anti-clockwise (positive rotation around z-axis). The element features are
summarized in :numref:`qdktplatesummary`.

.. table:: QDKTplate element summary
   :name: qdktplatesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | qdktplate                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D Discrete Kirchhoff Quad plate element                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (w-displacement, u and v - rotations) are required in each node.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of rotations, cubic approximation of displacement along the edges.   |
   |                          | Note: there is no need to define interpolation for displacement on the element.              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Default integration of all bending terms using four point formula.                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported.                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the generalized shell strain/force momentum vectors in global coordinate system   |
   |                          | are printed, with the following meaning: :math:`s_{\varepsilon}  = \{\varepsilon_x,          |
   |                          | \varepsilon_y, \varepsilon_{xz}, \kappa_x, \kappa_y, \kappa_{xy}, \gamma_{xz},               |
   |                          | \gamma_{yz}\},\; s_{\sigma}  = \{n_x, n_y, n_{xy}, m_x, m_y, m_z, m_{xy}, q_{xz}, q_{yz}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_{xy}` are membrane in plane normal    |
   |                          | deformations, :math:`\gamma_{zx}, \gamma_{xz}` are (out of plane and in plane) shear         |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_{xy}` are curvatures, :math:`n_x, n_y, n_{xy},  |
   |                          | q_{xz}, q_{yz}` are integral forces (normal and shear forces), and :math:`m_x, m_y, m_{xy}`  |
   |                          | are bending moments. Please note, for example, that bending moment :math:`m_x` is defined as |
   |                          | :math:`m_x=\int \sigma_x z\ dz`, so it acts along the y-axis and positive value causes       |
   |                          | tension in bottom layer.                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | J.L.Batoz, K.J.Bathe, L.W.Ho: A study of three-node triangular plate bending elements,       |
   |                          | IJNME, 15(12):1771-1812, 1980                                                                |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection CCT Element  Implementation of constant curvature triangular element for
plate analysis. Formulation based on Mindlin hypothesis. The structure should be defined
in x,y plane.  The nodes should be numbered anti-clockwise (positive rotation around
z-axis). The element features are summarized in :numref:`cctplatesummary`.

.. table:: cctplate element summary
   :name: cctplatesummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | cctplate                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | 2D constant curvature triangular plate element                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (w-displacement, u and v - rotations) are required in each node.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of rotations, quadratic approximation of displacement.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of all terms using one point formula.                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are not supported now.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the generalized shell strain/force momentum vectors in global coordinate system   |
   |                          | are printed, with the following meaning: :math:`s_{\varepsilon}  =\{\varepsilon_x,           |
   |                          | \varepsilon_y, \varepsilon_{xz}, \kappa_x, \kappa_y, \kappa_{xy}, \gamma_{xz},               |
   |                          | \gamma_{yz}\},\; s_{\sigma}  =\{n_x, n_y, n_{xy}, m_x, m_y, m_z, m_{xy}, q_{xz}, q_{yz}\}`   |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_{xy}` are membrane in plane normal    |
   |                          | deformations, :math:`\gamma_{zx}, \gamma_{xz}` are (out of plane and in plane) shear         |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_{xy}` are curvatures, :math:`n_x, n_y, n_{xy},  |
   |                          | q_{xz}, q_{yz}` are integral forces (normal and shear forces), and :math:`m_x, m_y, m_{xy}`  |
   |                          | are bending moments. Please note, for example, that bending moment :math:`m_x` is defined as |
   |                          | :math:`m_x=\int \sigma_x z\ dz`, so it acts along the y-axis and positive value causes       |
   |                          | tension in bottom layer.                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/patch_cct_02.in                                                                          |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch_cct_02.in>`_,           |
   |                          | `sm/patch200.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/patch200.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection CCT3D Element Implementation of constant curvature triangular element for
plate analysis. Formulation based on Mindlin hypothesis. The element could be
arbitrarily oriented in space.  The nodes should be numbered anti-clockwise (positive
rotation around element normal).  The element features are summarized in Table
:numref:`cctplate3dsummary`.

.. table:: cctplate3d element summary
   :name: cctplate3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | cctplate3d                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Constant curvature triangular plate element in arbitray position                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and u,v,w rotations) are in general required in each node.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of ratations, quadratic approximation of displacement.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of all terms using one point formula.                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are not supported now.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the shell force (:math:`s_f`), shell strain (:math:`s_s`), shell momentum         |
   |                          | (:math:`s_m`), and shell curvature (:math:`s_c`) tensors in global coordinate system are     |
   |                          | printed as vector form with 6 components, with the following meaning: :math:`s_f  =\{n_x,    |
   |                          | n_y, n_z, v_{yz}, v_{xz}, v_{xy}\},\; s_s  =\{\varepsilon_x, \varepsilon_y, \varepsilon_z,   |
   |                          | \gamma_{yz}, \gamma_{xz}, \gamma_{xy}\},\; s_m  =\{m_x, m_y, m_z, m_{yz}, m_{xz},            |
   |                          | m_{xy}\},\; s_c  =\{\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}\}`   |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_z` are membrane normal deformations,  |
   |                          | :math:`\gamma_{zy}, \gamma_{zx}, \gamma_{xy}` are (out of plane and in plane) shear          |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}` are   |
   |                          | curvatures, :math:`n_x, n_y, n_z, v_{yz}, v_{xz}, v_{xy}` are integral forces (normal and    |
   |                          | shear forces), and :math:`m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}` are bending moments. Please |
   |                          | note, for example, that bending moment :math:`m_x` is defined as :math:`m_x=\int \sigma_x z\ |
   |                          | dz`, so it acts along the y-axis and positive value causes tension in bottom layer.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection RerShell Element Combination of CCT plate element (Mindlin hypothesis)
with triangular plane stress element for membrane behavior. The element curvature can be
specified.  Although element requires generally six DOFs per node, no stiffness to local
rotation along z-axis (rotation around element normal) is supplied.  The element
features are summarized in :numref:`rershellsummary`.

.. table:: rershell element summary
   :name: rershellsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | rershell                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Simple shell based on combination of CCT plate element (Mindlin hypothesis) with triangular  |
   |                          | plane stress element. element can be arbitrary positioned in space.                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and u,v,w rotations) are in general required in each node.     |
   |                          | Note, that although element it requires generally six DOFs per node, no stiffness to local   |
   |                          | rotation along z-axis (rotation around element normal) is supplied.                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of ratations, quadratic approximation of displacement.                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of all terms using one point formula.                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are not supported now.                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the shell force (:math:`s_f`), shell strain (:math:`s_s`), shell momentum         |
   |                          | (:math:`s_m`), and shell curvature (:math:`s_c`) tensors in **global coordinate system** are |
   |                          | printed as vector form with 6 components, with the following meaning: :math:`s_f  = \{n_x,   |
   |                          | n_y, n_z, v_{yz}, v_{xz}, v_{xy}\},\; s_s  = \{\varepsilon_x, \varepsilon_y, \varepsilon_z,  |
   |                          | \gamma_{yz}, \gamma_{xz}, \gamma_{xy}\},\; s_m  = \{m_x, m_y, m_z, m_{yz}, m_{xz},           |
   |                          | m_{xy}\},\; s_c  = \{\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_z` are membrane normal deformations,  |
   |                          | :math:`\gamma_{zy}, \gamma_{zx}, \gamma_{xy}` are (out of plane and in plane) shear          |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}` are   |
   |                          | curvatures, :math:`n_x, n_y, n_z, v_{yz}, v_{xz}, v_{xy}` are integral forces (normal and    |
   |                          | shear forces), and :math:`m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}` are bending moments. Please |
   |                          | note, for example, that bending moment :math:`m_x` is defined as :math:`m_x=\int \sigma_x z\ |
   |                          | dz`, so it acts along the y-axis and positive value causes tension in bottom layer.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection tr_shell11 element Combination of CCT3D plate element (Mindlin
hypothesis) with triangular plane stress element for membrane behavior. It comes with
complete set of 6 DOFs per node.  The element features are summarized in Table
:numref:`trshell01summary`.

.. tikz:: Geometry of tr_shell11 element.

   \tikzstyle{elemnode} = [solid,draw,thin,circle,inner sep=1,fill=white]

   \begin{tikzpicture}[scale=6,>=stealth,
     x={(1cm,0cm)}, y={(0.5cm,0.5cm)}, z={(0cm,1cm)}]

    \begin{scope}
    \draw[->] (-0.05,0,0) -- (0.5,0,0) node[at end, below] {$x_g$};
    \draw[->] (0,-0.05,0) -- (0,0.5,0) node[at end, below right] {$y_g$};
    \draw[->] (0,0,-0.05) -- (0,0,0.5) node[at end, right] {$z_g$};
    \end{scope}

    \draw[very thick,-]
       (0.50, 0.10, 0.25) node[elemnode] {} node[below left] {1} -- node[midway, blue, below right] {1}
       (1.00, 0.25, 0.25) node[elemnode] {} node[right] {2}      -- node[near end, blue, above] {2}
       (0.10, 0.30, 0.50) node[elemnode] {} node[above left] {3} -- node[midway, blue, below left] {3}
       (0.50, 0.10, 0.25);

    % Shadow helps to visualize the depth
    \fill[fill=black!10]
       (0.50, 0.10, 0.) --
       (1.00, 0.25, 0.) --
       (0.10, 0.30, 0.) --
       (0.50, 0.10, 0.);

    % Draw normal (using latex arrow head, looks nicest)
    \draw[-latex] (0.5333,0.22, 0.33) -- +(0.1,-0.125,0.25) node[below right] {$\mathbf{n}$};

   \end{tikzpicture}

.. table:: tr_shell01 element summary
   :name: trshell01summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | tr_shell11                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Triangular shell element combining CCT3D plate element (Mindlin hypothesis) with triangular  |
   |                          | plane stress element with rotational DOFs                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and u,v,w rotations) are in general required in each node.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | See description of cct and trplanstrrot elements                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Four point integration, reduced integration of strain associated to normal rotation and      |
   |                          | shear terms.                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported (only surface loads).                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the shell force (:math:`s_f`), shell strain (:math:`s_s`), shell momentum         |
   |                          | (:math:`s_m`), and shell curvature (:math:`s_c`) tensors in **global coordinate system** are |
   |                          | printed as vector form with 6 components, with the following meaning: :math:`s_f  = \{n_x,   |
   |                          | n_y, n_z, v_{yz}, v_{xz}, v_{xy}\},\; s_s  = \{\varepsilon_x, \varepsilon_y, \varepsilon_z,  |
   |                          | \gamma_{yz}, \gamma_{xz}, \gamma_{xy}\},\; s_m  = \{m_x, m_y, m_z, m_{yz}, m_{xz},           |
   |                          | m_{xy}\},\; s_c  = \{\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_z` are membrane normal deformations,  |
   |                          | :math:`\gamma_{zy}, \gamma_{zx}, \gamma_{xy}` are (out of plane and in plane) shear          |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}` are   |
   |                          | curvatures, :math:`n_x, n_y, n_z, v_{yz}, v_{xz}, v_{xy}` are integral forces (normal and    |
   |                          | shear forces), and :math:`m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}` are bending moments. Please |
   |                          | note, for example, that bending moment :math:`m_x` is defined as :math:`m_x=\int \sigma_x z\ |
   |                          | dz`, so it acts along the y-axis and positive value causes tension in bottom layer.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection tr_shell02 element Combination of thin-plate DKT plate element with plane
stress element (TrPlanestressRotAllman). This element comes with complete set of 6 DOFs
per node.  The element features are summarized in :numref:`trshell02summary`.

.. table:: tr_shell02 element summary
   :name: trshell02summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | tr_shell02                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Triangular shell element combining DKT plate element with triangular plane stress element    |
   |                          | with rotational DOFs                                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and u,v,w rotations) are in general required in each node.     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | See description of cct and trplanstrrot elements                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | 4 integration points necessary, use "NIP 4" in element record.                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body loads are supported. Boundary loads are supported (only surface loads).                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the shell force (:math:`s_f`), shell strain (:math:`s_s`), shell momentum         |
   |                          | (:math:`s_m`), and shell curvature (:math:`s_c`) tensors in **global coordinate system** are |
   |                          | printed as vector form with 6 components, with the following meaning: :math:`s_f  = \{n_x,   |
   |                          | n_y, n_z, v_{yz}, v_{xz}, v_{xy}\},\; s_s  = \{\varepsilon_x, \varepsilon_y, \varepsilon_z,  |
   |                          | \gamma_{yz}, \gamma_{xz}, \gamma_{xy}\},\; s_m  = \{m_x, m_y, m_z, m_{yz}, m_{xz},           |
   |                          | m_{xy}\},\; s_c  = \{\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_z` are membrane normal deformations,  |
   |                          | :math:`\gamma_{zy}, \gamma_{zx}, \gamma_{xy}` are (out of plane and in plane) shear          |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}` are   |
   |                          | curvatures, :math:`n_x, n_y, n_z, v_{yz}, v_{xz}, v_{xy}` are integral forces (normal and    |
   |                          | shear forces), and :math:`m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}` are bending moments. Please |
   |                          | note, for example, that bending moment :math:`m_x` is defined as :math:`m_x=\int \sigma_x z\ |
   |                          | dz`, so it acts along the y-axis and positive value causes tension in bottom layer.          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Works only with linear material models, as bending and membrane actions are uncoupled        |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _quad1mindlin:

Quad1Mindlin Element
~~~~~~~~~~~~~~~~~~~~

This class implements an quadrilateral, bilinear, four-node Mindlin plate. This type of
element exhibit strong shear locking (thin plates exhibit almost no bending). Implements
the lumped mass matrix. The element features are summarized in Table
:numref:`quad1mindlinsummary`.

.. table:: quad1mindlin element summary
   :name: quad1mindlinsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | quad1mindlin                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Quadrilateral, bilinear, four-node Mindlin plate                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Three dofs (w-displacement, u and v - rotation) are required in each node.                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear for all unknowns.                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Default uses 4 integration points. No reduced integration is used, as it causes numerical    |
   |                          | problems.                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Dead weight loads, and edge loads are supported.                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the generalized shell strain/force momentum vectors in global coordinate system   |
   |                          | are printed, with the following meaning: :math:`s_{\varepsilon}  = \{\varepsilon_x,          |
   |                          | \varepsilon_y, \varepsilon_{xz}, \kappa_x, \kappa_y, \kappa_{xy}, \gamma_{xz},               |
   |                          | \gamma_{yz}\},\; s_{\sigma}  = \{n_x, n_y, n_{xy}, m_x, m_y, m_z, m_{xy}, q_{xz}, q_{yz}\}`  |
   |                          | where :math:`\varepsilon_x, \varepsilon_y, \varepsilon_{xy}` are membrane in plane normal    |
   |                          | deformations, :math:`\gamma_{zx}, \gamma_{xz}` are (out of plane and in plane) shear         |
   |                          | componets, :math:`\kappa_x, \kappa_y, \kappa_{xy}` are curvatures, :math:`n_x, n_y, n_{xy},  |
   |                          | q_{xz}, q_{yz}` are integral forces (normal and shear forces), and :math:`m_x, m_y, m_{xy}`  |
   |                          | are bending moments. Please note, for example, that bending moment :math:`m_x` is defined as |
   |                          | :math:`m_x=\int \sigma_x z\ dz`, so it acts along the y-axis and positive value causes       |
   |                          | tension in bottom layer.                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [RobertCook1989]_                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Experimental                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _tr2shell7:

Tr2Shell7 Element
~~~~~~~~~~~~~~~~~

This class implements a triangular, quadratic, six-node shell element. The element is a
so-called seven parameter shell with seven dofs per node -- a displacement field (3
dofs), an extensible director field (3 dofs) and a seventh dof representing inhomogenous
thickness strain. This last parameter is included in the model in order to deal with
volumetric/Poisson lock effects.

The element features are summarized in :numref:`quad1mindlinsummary`.

.. table:: tr2shell7 element summary
   :name: tr2shell7summary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | tr2shell7                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Triangular, quadratic, six-node shell with 7 dofs/node                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Seven dofs (displacement in u, v and w-direction; change in director field in u, v and       |
   |                          | w-direction; and inhomgenous thickness stretch) are required in each node.                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic for all unknowns.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Default uses 6 integration points in the midsurface plane. Number of integration points in   |
   |                          | the thickness direction is determined by the Layered cross section.                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Layered cross section support.                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | This element must be used with a Layered cross section.                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Edge loads, constant pressure loads and surface loads are supported.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | Not applicable. The implementation is for large defomrations and hence geometrical           |
   |                          | nonlinearities will always be present, regardless the value of Nlgeo.                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [RagnarLarsson2011]_                                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | Experimental                                                                                 |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `benchmark/sm/tr2shell7.in                                                                   |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/tr2shell7.in>`_,    |
   |                          | `sm/staggeredsolver.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/staggeredsolver.in>`_         |
   +--------------------------+----------------------------------------------------------------------------------------------+

\subsubsection MITC4Shell Element A four-node quadrilateral shell element formulated
using three-dimensional continuum mechanics theory degenerated to shell behaviour. The
element is applicable to thick and thin shells as the “mixed interpolation of tensorial
components” (MITC) approach is used to remove shear locking. The implementation is based
on the following paper: Dvorkin, E.N., Bathe, K.J.: A continuum mechanics based
four-node shell element for general non-linear analysis, Eng.Comput., Vol.1, 77-88,
1984.

Although element requires generally six DOFs per node, no stiffness to local rotation
along z-axis (rotation around director vector) is supplied. The element features are
summarized in :numref:`mitc4shellsummary`.

.. table:: mitc4shell element summary
   :name: mitc4shellsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | mitc4shell                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Quadrilateral, bilinear, four-node shell element using the MITC technique.                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`NIPZ{in}` :optelemparam:`directorType{in}`            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points in local x-y plane (default 4). |
   |                          | :param:`NIPZ`: allows to set the number of integration points in local z-direction (default  |
   |                          | 2). :param:`directorType`: allows to set director vectors. Director vectors can be set as    |
   |                          | normal to the plane (:param:`directorType` = 0, default), or calculated for each node as an  |
   |                          | average of neighbouring elements of same crosssection (:param:`directorType` = 1), or can be |
   |                          | specified at crosssection (:param:`directorType` =2).                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Six dofs (u,v,w-displacements and u,v,w rotations) are in general required in each node.     |
   |                          | Note, that although element requires generally six DOFs per node, no stiffness to local      |
   |                          | rotation along z-axis (rotation around director vector) is supplied.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacements and rotations.                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration of all terms using Gauss integration formula in 8 points (default) or specified  |
   |                          | using :param:`NIP` and :param:`NIPZ` parameters.                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | Variable cross section support.                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | Cross section thickness is required (measured along director vector). Director vectors       |
   |                          | components may be specified                                                                  |
   |                          | :optelemparam:`directorx{in}`:optelemparam:`directory{in}`:optelemparam:`directorz{in}` in   |
   |                          | case of :param:`directorType 2`.                                                             |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Body and boundary loads are supported.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Output                   | On output, the shell force (:math:`s_f`), shell momentum (:math:`s_m`), shell strain         |
   |                          | (:math:`s_s`), shell curvature (:math:`s_c`), strain (:math:`\varepsilon`), and stress       |
   |                          | (:math:`\sigma`) tensors in **global coordinate system** are printed as vector form with 6   |
   |                          | components, with the following meaning: :math:`s_f  = \{n_x, n_y, n_z, v_{yz}, v_{xz},       |
   |                          | v_{xy}\},\; s_m  = \{m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}\},\; s_s  = \{\varepsilon_x,      |
   |                          | \varepsilon_y, \varepsilon_z, \gamma_{yz}, \gamma_{xz}, \gamma_{xy}\},\; s_c  = \{\kappa_x,  |
   |                          | \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}\}\; \varepsilon = \{\varepsilon_x, |
   |                          | \varepsilon_y, \varepsilon_z, \gamma_{yz}, \gamma_{zx}, \gamma_{xy}\},\; \sigma  =           |
   |                          | \{\sigma_x, \sigma_y, \sigma_z, \sigma_{yz}, \sigma_{xz}, \sigma_{xy}\}.` where :math:`n_x,  |
   |                          | n_y, n_z, v_{yz}, v_{xz}, v_{xy}` are integral forces (normal and shear forces), and         |
   |                          | :math:`m_x, m_y, m_z, m_{yz}, m_{xz}, m_{xy}` are bending moments, :math:`\varepsilon_x,     |
   |                          | \varepsilon_y, \varepsilon_z` are membrane normal deformations, :math:`\gamma_{zy},          |
   |                          | \gamma_{zx}, \gamma_{xy}` are (out of plane and in plane) shear componets, :math:`\kappa_x,  |
   |                          | \kappa_y, \kappa_z, \kappa_{yz}, \kappa_{xz}, \kappa_{xy}` are curvatures. Please note, for  |
   |                          | example, the bending moment :math:`m_x` is defined as :math:`m_x=\int \sigma_x z\ dz`, so it |
   |                          | acts along the y-axis and positive value causes tension in bottom layer (positive            |
   |                          | z-coordinate). The shell force (:math:`s_f`), shell momentum (:math:`s_m`), shell strain     |
   |                          | (:math:`s_s`), and shell curvature (:math:`s_c`) tensors are evaluated at the midplane of    |
   |                          | the element (thus are constant along the thickness) while the strain (:math:`\varepsilon`),  |
   |                          | and stress (:math:`\sigma`) tensors are evaluated at each Gausspoint.                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/scordelis_mitc4.in                                                                       |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/scordelis_mitc4.in>`_         |
   +--------------------------+----------------------------------------------------------------------------------------------+

Sub-soil Elements
~~~~~~~~~~~~~~~~~

.. _quad1platesubsoil:

quad1plateSubsoil Element
~~~~~~~~~~~~~~~~~~~~~~~~~

This class implements an quadrilateral, bilinear, four-node plate subsoil element.
Typically this element is combined with suitable plate element with quadrilateral
geometry to model plate element on  (elastic) subsoill foundation, but it can be used
alone. The element geometry should be define in xy plane. The element features are
summarized in :numref:`quad1platesubsoilsummary`.

.. table:: quad1platesubsoil element summary
   :name: quad1platesubsoilsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | quad1plateSubsoil                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Quadrilateral, bilinear, four-node sub-soil plate element                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of membrane     |
   |                          | terms.                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | One dof (w-displacement) is required in each node.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear for transwersal displacement.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | 4 integration points.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Surface load support.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Requires material model with 2dPlateSubSoil mode support.                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [BittnarSejnoha1996]_                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/test_wp5.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp5.in>`_,               |
   |                          | `sm/test_wp1.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp1.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. _tria1platesubsoil:

Tria1PlateSubSoil Element
~~~~~~~~~~~~~~~~~~~~~~~~~

This class implements an quadrilateral, bilinear, four-node plate subsoil element.
Typically this element is combined with suitable plate element with quadrilateral
geometry to model plate element on  (elastic) subsoill foundation, but it can be used
alone. The element geometry should be define in xy plane. The element features are
summarized in :numref:`quad1platesubsoilsummary`.

.. table:: tria1platesubsoil element summary
   :name: tria1platesubsoilsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | tria1platesubsoil                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Tringular, three-node sub-soil plate element with linear interpolation                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | One dof (w-displacement) is required in each node.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear for transwersal displacement.                                                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | 1 integration points.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Surface load support.                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | Requires material model with 2dPlateSubSoil mode support.                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Reference                | [BittnarSejnoha1996]_                                                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/test_wp5.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp5.in>`_,               |
   |                          | `sm/test_wp2.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/test_wp2.in>`_                |
   +--------------------------+----------------------------------------------------------------------------------------------+

Axisymmetric Elements
=====================


Axisymmetric Elements
---------------------

Implementation relies on elements located exclusively in :math:`x,y` plane. The
coordinate :math:`x` corresponds to radius, :math:`y` is the axis of rotation.
Approximation of displacement functions :math:`u,v` is carried out on a particular
finite element. Nonzero strains read

.. math::

   \varepsilon_{x}=\varepsilon_{r}&=&\frac{\partial u}{\partial x}\\
   \varepsilon_{y}=\varepsilon_{z}&=&\frac{\partial v}{\partial y}\\
   \varepsilon_{\theta}&=&\frac{u}{r}\\
   \gamma_{xy}=\gamma_{rz}&=&\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}

Stress components can be computed from elasticity matrix. Note that this matrix
corresponds to a submatrix of the full 3D elasticity matrix.

.. math::

   \left\{\begin{array}{c} \sigma_x \\ \sigma_y \\ \sigma_\theta \\ \sigma_{xy} \end{array} \right\} =
   \frac{E}{(1+\nu)(1-2\nu)}\left[ \begin{array}{cccc} 1-\nu & \nu & \nu & 0 \\ \nu & 1-\nu & \nu & 0 \\ \nu & \nu & 1-\nu & 0 \\ 0 & 0 & 0 & (1-2\nu)/2 \end{array} \right] =
   \left\{\begin{array}{c} \varepsilon_x \\ \varepsilon_y \\ \varepsilon_\theta \\ \gamma_{xy} \end{array} \right\}

In OOFEM, the strain vector is arranged as :math:`\{\varepsilon_{x}, \varepsilon_{y},
\varepsilon_{\theta}, 0, 0, \gamma_{xy}\}^T` and the stress vector :math:`\{\sigma_{x},
\sigma_{y}, \sigma_{\theta}, 0, 0, \tau_{xy}\}^T`. Implementation assumes a segment of 1
rad.

Axisymm3d element
~~~~~~~~~~~~~~~~~

Implementation of triangular three-node finite element  for axisymmetric continuum. Each
node has 2 degrees of freedom.  Node numbering and edge position is the same as in Fig.
:ref:`TrPlanestressfig`.  The element features are summarized in Table
:numref:`axisymm3dsummary`.

.. table:: Axisymm3d element summary
   :name: axisymm3dsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Axisymm3d                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Triangular axisymmetric linear element                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`NIPfish{in}`                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points (possible completions are 1     |
   |                          | (default), 4 and 7 point integration rule).                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacement and geometry.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | The integration can be altered using :param:`NIP` paramter (default is 1 point integration). |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Boundary and body loads are supported.                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/axisymm06.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm06.in>`_,              |
   |                          | `sm/axisymm04.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm04.in>`_,              |
   |                          | `sm/axisymm05.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm05.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

Q4axisymm element
~~~~~~~~~~~~~~~~~

Implementation of quadratic isoparametric eight-node quadrilateral - finite element for
axisymmetric 3d continuum.  Each node has 2 degrees of freedom. The element features are
summarized in :numref:`q4axisymmsummary`.

.. table:: Q4axisymm element summary
   :name: q4axisymmsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Q4axisymm                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Quadratic isoparametric eight-node quadrilateral for axisymmetric analysis                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}` :optelemparam:`NIPfish{in}`                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of terms        |
   |                          | corresponding to :math:`\varepsilon_x` and :math:`\varepsilon_y` strains (possible           |
   |                          | completions are 1, 4 (default), 9, and 16). :param:`NIPfish`: allows to set the number of    |
   |                          | integration points for integration of remain terms (corresponding to                         |
   |                          | :math:`\varepsilon_\theta` and :math:`\gamma_{rz}`) (Supported values include 1 (default),   |
   |                          | 4, 9, and 16 integration point formula).                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of displacement and geometry.                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | The integration of terms corresponding to :math:`\varepsilon_x` and :math:`\varepsilon_y`    |
   |                          | strains can be altered using :param:`NIP` parameter (default is 4 point formula). The        |
   |                          | remaining terms (creesponding to :math:`\varepsilon_\theta` and :math:`\gamma_{rz}`) are     |
   |                          | integrated by default using 1 point formula (see :param:`NIPfish` parameter).                |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | No boundary and body loads are supported.                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+

L4axisymm element
~~~~~~~~~~~~~~~~~

Implementation of isoparametric four-node quadrilateral axisymmetric finite element with
linear interpolations of displacements :math:`u, v`. Node numbering and edge position is
the same as in Fig. :ref:`Planestress2dfig`. The element features are summarized in
:numref:`l4axisymmsummary`.

.. table:: L4axisymm element summary
   :name: l4axisymmsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | L4axisymm                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Isoparametric four-node quadrilateral element for axisymmetric analysis                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to set the number of integration points for integration of terms        |
   |                          | corresponding to :math:`\varepsilon_x` and :math:`\varepsilon_y` strains (possible           |
   |                          | completions are 1, 4 (default), 9, and 16).                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Two dofs (u-displacement, v-displacement) are required in each node.                         |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Linear approximation of displacement and geometry.                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | The integration of :math:`\varepsilon_x` and :math:`\varepsilon_y` strains can be altered    |
   |                          | using :param:`NIP` parameter (possible completions are 1, 4 (default), 9 or 16 point         |
   |                          | integration rule). The remaining strain components (:math:`\varepsilon_\theta` and           |
   |                          | :math:`\gamma_{rz}`) are integrated using one point integration formula.                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | Boundary and body loads supported.                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Nlgeo                    | 0.                                                                                           |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/axisymm03.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm03.in>`_,              |
   |                          | `sm/axisymm01.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm01.in>`_,              |
   |                          | `sm/axisymm02.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/axisymm02.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

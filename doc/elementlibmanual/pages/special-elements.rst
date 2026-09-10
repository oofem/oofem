Special elements
================


Special elements
----------------

LumpedMass element
~~~~~~~~~~~~~~~~~~

This element, defined by a single node, allows to introduce additional concentrated mass
and/or rotational inertias in a node. A different mass and rotary inertia may be
assigned to each coordinate direction. At present, individual mass/inertia components
can be specified for every degree of freedom of element node. Only displacement and
rotational degrees od freedom are considered. The element features are summarized in
:numref:`LumpedMasssummary`.

.. table:: LumpedMass element summary
   :name: LumpedMasssummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | LumpedMass                                                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Lumped mass element                                                                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`components{ra}`                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`components`: allows to specify additional concentrated mass components (Force*Time\  |
   |                          | :math:`^2`/Length) and rotary inertias (Force*Length*Time\ :math:`^2`) about the nodal       |
   |                          | coordinate axes. :param:`dofs`: dofs to which the components apply.                          |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | As specified by :param:`dofs`.                                                               |
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
   | Status                   | Reliable                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+

Spring element
~~~~~~~~~~~~~~

This element represent longitudial or torsional spring element. It is defined by two
nodes, orientation and a spring constant. The spring element has no mass associated, the
mass can be added using LumpedMass element. The spring is linear and works the same way
in tension or in compression. The element features are summarized in Table
:numref:`Springsummary`.

.. table:: Spring element summary
   :name: Springsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | Spring                                                                                       |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Spring element                                                                               |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :elemparam:`mode{in}` :elemparam:`k{rn}` :optelemparam:`m{rn}` :elemparam:`orientation{ra}`  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`mode`: defines the type of spring element (see :numref:`spring_mode_table`).         |
   |                          | :param:`k`: determines the spring constant, corresponding units are [Force/Length] for       |
   |                          | longitudinal spring and [Force*Length/Radian] for torsional spring.                          |
   |                          | :param:`orientation`:defines orientation vector of spring element (of size 3) - for          |
   |                          | longitudinal spring it defines the direction of spring, for torsional spring it defines the  |
   |                          | axis of rotation. :param:`m`: determines optional mass of the element, zero value assumed by |
   |                          | default.                                                                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Note                     | the spring element nodes doesn't need to be coincident, but the spring orientation is always |
   |                          | determined by :param:`orientation` vector.                                                   |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `sm/spring01.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring01.in>`_,               |
   |                          | `sm/spring02.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring02.in>`_,               |
   |                          | `sm/spring03.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring03.in>`_,               |
   |                          | `sm/spring04.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring04.in>`_,               |
   |                          | `sm/spring05.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring05.in>`_,               |
   |                          | `sm/spring06.in                                                                              |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/sm/spring06.in>`_,               |
   |                          | `benchmark/sm/contact/friction_wip/visual/contact2d_inclined_plane_stick.in <https://github. |
   |                          | com/oofem/oofem/blob/devel/tests/regression/benchmark/sm/contact/friction_wip/visual/contact |
   |                          | 2d_inclined_plane_stick.in>`_,                                                               |
   |                          | `benchmark/sm/contact/friction_wip/visual/contact2d_inclined_plane_slip.in <https://github.c |
   |                          | om/oofem/oofem/blob/devel/tests/regression/benchmark/sm/contact/friction_wip/visual/contact2 |
   |                          | d_inclined_plane_slip.in>`_                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+

.. table:: Supported spring element modes
   :name: spring_mode_table

   +------+--------------------------------------------------------------+
   | mode | description                                                  |
   +======+==============================================================+
   | 0    | 1D spring element along x-axis,                              |
   +------+--------------------------------------------------------------+
   |      | requires D_u DOF in each node, orientation vector is {1,0,0} |
   +------+--------------------------------------------------------------+
   | 1    | 2D spring element in xy plane,                               |
   +------+--------------------------------------------------------------+
   |      | requires D_u and D_v DOFs in each node                       |
   +------+--------------------------------------------------------------+
   |      | (orientation vector should be in xy plane)                   |
   +------+--------------------------------------------------------------+
   | 2    | 2D spring element in xz plane,                               |
   +------+--------------------------------------------------------------+
   |      | requires D_u and D_w DOFs in each node                       |
   +------+--------------------------------------------------------------+
   |      | (orientation vector should be in xz plane)                   |
   +------+--------------------------------------------------------------+
   | 3    | 2D torsional spring element in xz plane,                     |
   +------+--------------------------------------------------------------+
   |      | requires R_v DOFs in each node                               |
   +------+--------------------------------------------------------------+
   | 4    | 3D spring element in space,                                  |
   +------+--------------------------------------------------------------+
   |      | requires D_u, D_v, and D_w DOFs in each node                 |
   +------+--------------------------------------------------------------+
   | 5    | 3D torsional spring in space,                                |
   +------+--------------------------------------------------------------+
   |      | requires R_u, R_v, and R_w DOFs in each node                 |
   +------+--------------------------------------------------------------+

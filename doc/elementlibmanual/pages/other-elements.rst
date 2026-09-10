Other Elements
==============


Other Elements
--------------

.. _QBrick1ht_element:

QBrick1ht - quadratic hexahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of quadratic 3d 20-node  finite element. Each node has 1 degree of
freedom. See section :ref:`QSpace_element`  for node numbering order and order of
faces. The element features are summarized in :numref:`QBrick1htsummary`.

.. table:: QBrick1ht element summary
   :name: QBrick1htsummary

   +--------------------------+----------------------------------------------------------------------------------------------+
   | Keyword                  | QBrick1ht                                                                                    |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Description              | Isoparametric, hexahedral 3D element with quadratic approximation for heat transfer problems |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Specific parameters      | :optelemparam:`NIP{in}`                                                                      |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Parameters               | :param:`NIP`: allows to change the default number of integration point used, possible values |
   |                          | are 8, 27 (default) and 64.                                                                  |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Unknowns                 | Single dof (T_f - temperature) is required in each node.                                     |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Approximation            | Quadratic approximation of temperature and geometry..                                        |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Integration              | Integration using gauss integration formula in 8, 27 (default), or 64 integration points.    |
   |                          | The default number of integration point used can be overloaded using :param:`NIP` parameter. |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Loads                    | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Features                 | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | CS properties            | -                                                                                            |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Status                   |                                                                                              |
   +--------------------------+----------------------------------------------------------------------------------------------+
   | Tests/Examples           | `tm/qbrick_02.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_02.in>`_,              |
   |                          | `tm/qbrick_03.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_03.in>`_,              |
   |                          | `tm/qbrick_01.in                                                                             |
   |                          | <https://github.com/oofem/oofem/blob/devel/tests/regression/tm/qbrick_01.in>`_               |
   +--------------------------+----------------------------------------------------------------------------------------------+

QBrick1mt - quadratic hexahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same element as QBrick1ht for mass transfer problems, see
:ref:`QBrick1ht_element`.  Linear approximation of mass concentration.

QBrick1hmt - quadratic hexahedral 3D element
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same element as QBrick1ht for heat and mass (one constituent) transfer problems.
Two dofs (T_f - temperature and C_1 - concentration) are required in each node. Linear
approximation of temperature and mass concentration. Other features are similar to
QBrick1ht element, see section :ref:`QBrick1ht_element`.

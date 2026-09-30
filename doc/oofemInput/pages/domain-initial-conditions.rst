.. _InitialConditions:

Initial condition records
=========================

These records describe the initial conditions.  They are counted by
:param:`nic` in the :ref:`ComponentsSizeRecord`, and the numbering must start
at one and end at *n*, where *n* is the number of records.

An initial condition works in one of three mutually exclusive modes, selected
by which of :param:`conditions`, :param:`f` and :param:`field` is present —
tested in that order.  Giving none of them is an error.

Expression mode
---------------

The preferred form.  The initial value, velocity and acceleration are given as
expressions, which may depend on position through the variables ``x``, ``y``
and ``z``, the coordinates of the node.

.. record::

   :descitem:`InitialCondition` :elemparam:`num{in}` :elemparam:`f{expr}`
   :elemparam:`set{in}` :elemparam:`dofs{ia}`
   :optelemparam:`dfdt{expr}` :optelemparam:`d2fdt2{expr}`
   :optelemparam:`valType{in}`

**Parameters**

:param:`f{expr}`
    Expression for the initial value of the unknown.  Defaults to ``0.0``.

:optparam:`dfdt{expr}`
    Expression for the initial velocity, the first time derivative.  Defaults
    to ``0.0``.

:optparam:`d2fdt2{expr}`
    Expression for the initial acceleration, the second time derivative.
    Defaults to ``0.0``.

:param:`dofs{ia}`
    The individual DOFs the condition applies to, using the ``DofIDItem``
    values listed under :ref:`NodeElementSideRecords`.

:param:`set{in}`
    The set of nodes the initial condition is applied to; see
    :ref:`SetRecords`.

:optparam:`valType{in}`
    Physical meaning of the value; see :param:`valType` under
    :ref:`LoadBoundaryInitialConditions`.

Field mode
----------

Initialises the unknowns from an externally registered field, for example one
exported by a preceding sub-problem of a staggered analysis; see
:ref:`exportfields`.

.. record::

   :descitem:`InitialCondition` :elemparam:`num{in}`
   :elemparam:`field{in}` :elemparam:`set{in}` :elemparam:`dofs{ia}`
   :optelemparam:`valType{in}`

**Parameters**

:param:`field{in}`
    Identifier of the field to read the initial values from, as a
    ``FieldType`` enumeration value; see ``src/core/field.h``.

Dictionary mode (deprecated)
----------------------------

.. warning::

   Retained only for compatibility with old input files.  Use the expression
   mode instead.

.. record::

   :descitem:`InitialCondition` :elemparam:`num{in}`
   :elemparam:`conditions{dc}` :optelemparam:`valType{in}`

**Parameters**

:param:`conditions{dc}`
    Dictionary of initial values, listing keys followed by their values.  The
    key ``v`` is the velocity and ``a`` the acceleration.

In this mode the condition is attached to dof managers through their ``ic``
array rather than through a set; see :ref:`NodeElementSideRecords`.

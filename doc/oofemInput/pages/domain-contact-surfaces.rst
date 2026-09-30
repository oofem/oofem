.. _ContactSurfaceRecords:

Contact surface records
=======================

A contact surface is a geometric entity taking part in contact interactions.
Each surface is formed by a set of contact elements, and can act as the *master*
or the *slave* within a contact boundary condition.  Contact surfaces are the
link between the geometric description of the contacting interfaces and the
boundary conditions that enforce their interaction.

The number of contact surface records is declared by :param:`ncontactsurf` in
the :ref:`ComponentsSizeRecord`.

.. record::

   :descitem:`FEContactSurfaceType` :elemparam:`num{in}`
   :elemparam:`ce_set{in}`

**Parameters**

:param:`num{in}`
    Unique surface identifier.

:param:`ce_set{in}`
    Number of the element set containing the contact elements that form this
    surface; see :ref:`SetRecords`.

**Supported surface types**

``StructuralFEContactSurface``
    Contact surface for structural contact, made of
    ``StructuralContactElement_*`` elements.  Used by
    :ref:`structural-penalty-contact`.

``ThermalFEContactSurface``
    Contact surface for thermal contact.  Used by
    :ref:`thermal-s2s-contact`.

Description
-----------

Each contact surface groups one or more contact elements into a single logical
entity.  These surfaces are then referenced by the contact boundary conditions
(see :ref:`ContactBoundaryConditions`) as their :param:`mastersurface` or
:param:`slavesurface`.  Both surfaces must exist before the contact boundary
condition that pairs them is defined.

Contact surfaces introduce no new DOFs; they provide only the geometric
organisation used when the contact is evaluated.  Each can serve as a master or
as a slave surface, depending on how the contact boundary condition is written.

**Example**

.. code-block:: none

   StructuralFEContactSurface 1 ce_set 3
   StructuralFEContactSurface 2 ce_set 4

Structural finite element surface 1 is built from the contact elements of set 3
and surface 2 from those of set 4.  The two can then be paired by a
``structuralpenaltycontactbc`` record to define a contact interaction.  See
``tests/regression/sm/contact2d_projection_outside.in`` for a complete input
file.

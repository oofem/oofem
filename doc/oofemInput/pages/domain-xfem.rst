.. _XFEMManagerRecords:

XFEM manager and associated records
===================================

The XFEM manager record specifies the number of enrichment items and the
simulation options common to all of them.  It is created when
:param:`nxfemman` is ``1``; see :ref:`ComponentsSizeRecord`.

The functions used for enrichment — Heaviside, absolute value, branch functions
and so on — are *not* specified here, nor is the geometrical representation of
an enrichment item, such as a polygon line or a circle: both are specified for
each enrichment item separately.

.. note::

   OOFEM currently supports XFEM simulations of cracks and material interfaces
   in 2D.

XFEM manager record
-------------------

.. record::

   :descitem:`XfemManager` :elemparam:`numberofenrichmentitems{in}`
   :optelemparam:`numberofnucleationcriteria{in}`
   :optelemparam:`numberofgppertri{in}`
   :optelemparam:`numberoftrirefs{in}`
   :optelemparam:`enrdofscalefac{rn}` :optelemparam:`debugvtk{in}`
   :optelemparam:`vtkexport{in}` :optelemparam:`exportfields{ia}`

**Parameters**

:param:`numberofenrichmentitems{in}`
    Number of enrichment items, whose records follow.

:optparam:`numberofnucleationcriteria{in}`
    Number of nucleation criteria, whose records follow the enrichment items.

:optparam:`numberofgppertri{in}`
    Number of Gauss points in each subtriangle of a cut element.  Default
    ``12``.

:optparam:`numberoftrirefs{in}`
    How many times a subtriangle is refined.  Default ``0``.

:optparam:`enrdofscalefac{rn}`
    Scale factor applied to the enrichment DOFs, which can improve the
    conditioning of the system matrix.  Default ``1.0``.

:optparam:`debugvtk{in}`
    ``1`` writes additional debug VTK files; ``0`` is the default.

:optparam:`vtkexport{in}`
    A nonzero value turns on the XFEM VTK export.  When it is on,
    :param:`exportfields` becomes **mandatory**.

:optparam:`exportfields{ia}`
    Which XFEM fields to export.  Mandatory when :param:`vtkexport` is
    nonzero, and used by the ``vtkxmlxfem`` export module; see
    :ref:`ExportModulesSec`.

Enrichment item records
-----------------------

The specification of an enrichment item consists of several consecutive lines.
See the test ``tests/regression/sm/xFemCrackValBranch.in`` for a complete
example.

First line — the enrichment item
    .. record:: none

       :descitem:`EntType` :elemparam:`num{in}`
       :optelemparam:`enrichmentfront{in}`
       :optelemparam:`propagationlaw{in}`

    :optparam:`enrichmentfront{in}`
        Number of an enrichment front.  The front of the enrichment is treated
        separately: we may, for example, use branch functions at a crack tip
        and Heaviside enrichment along the rest of the crack.  Optional.

    :optparam:`propagationlaw{in}`
        Number of a rule for crack propagation.  Optional.

        .. warning::

           Crack propagation is still highly experimental.

Second line — the enrichment function
    .. record:: none

       :descitem:`EntType` :elemparam:`num{in}`

Third line — the geometric description
    :descitem:`EntType` :elemparam:`num{in}` plus the attributes the chosen
    geometry needs.  The number and type of those attributes vary with the
    geometry: a centre and a radius for a circle, a number of points for a
    polygon line, and so on.

Fourth line — the enrichment front
    :descitem:`EntType` :elemparam:`num{in}` plus its attributes.  Present
    only if an :param:`enrichmentfront` was specified on the first line.

Fifth line — the propagation law
    :descitem:`EntType` :elemparam:`num{in}` plus its attributes.  Present
    only if a :param:`propagationlaw` was specified on the first line.

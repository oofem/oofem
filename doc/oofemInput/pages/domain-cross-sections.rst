.. _CrossSectionRecords:

Cross section records
=====================

These records describe the cross section models.

.. record::

   :descitem:`CrossSectType` :elemparam:`num{in}` :optelemparam:`set{in}`

The order of the records is arbitrary and the cross section number is given by
:param:`num`.  Unlike nodes and elements, the numbering must start at one and
end at *n*, where *n* is the number of records.

**Parameter common to all cross sections**

:optparam:`set{in}`
    Number of an element set (see :ref:`SetRecords`) to which this cross
    section assigns itself.  This is the usual way of binding a cross section —
    and through it a material — to elements, and it is why an element record
    often carries neither ``crossSect`` nor ``mat``.

All cross section properties are understood to be defined in the local
coordinate system of the element.

.. _SimpleCS:

Integral cross section with constant properties
-----------------------------------------------

The integral type of cross section model.  For some problems, 3D ones in
particular, the volume and the cross section dimensions follow from the element
geometry, in which case some or all of these parameters may be omitted; consult
the Element Library Manual for the properties an individual element requires.

.. record::

   :descitem:`SimpleCS` :elemparam:`num{in}` :optelemparam:`material{in}`
   :optelemparam:`set{in}` :optelemparam:`thick{rn}`
   :optelemparam:`width{rn}` :optelemparam:`area{rn}`
   :optelemparam:`iy{rn}` :optelemparam:`iz{rn}` :optelemparam:`ik{rn}`
   :optelemparam:`beamshearcoeff{rn}` :optelemparam:`shearareay{rn}`
   :optelemparam:`shearareaz{rn}` :optelemparam:`drillstiffness{rn}`
   :optelemparam:`reldrillstiffness{rn}` :optelemparam:`drilltype{in}`
   :optelemparam:`directorx{rn}` :optelemparam:`directory{rn}`
   :optelemparam:`directorz{rn}`

**Parameters**

:optparam:`material{in}`
    Number of the material record for the bulk material.

:optparam:`thick{rn}`, :optparam:`width{rn}`
    Cross section thickness and width.

:optparam:`area{rn}`
    Cross section area.

:optparam:`iy{rn}`, :optparam:`iz{rn}`
    Inertia moments about the local *y* and *z* axes.

:optparam:`ik{rn}`
    Saint-Venant torsional constant.

:optparam:`beamshearcoeff{rn}`
    Shear correction factor.  Alternatively give the equivalent shear areas
    :param:`shearareay` and :param:`shearareaz`.

:optparam:`shearareay{rn}`, :optparam:`shearareaz{rn}`
    Equivalent shear areas in the *y* and *z* directions.

:optparam:`drillstiffness{rn}`, :optparam:`reldrillstiffness{rn}`, :optparam:`drilltype{in}`
    Artificial stiffness added for drilling DOFs: an absolute penalty term, a
    relative penalty term, and the type of the added stiffness.  Relevant for
    shell elements with drilling degrees of freedom.

:optparam:`directorx{rn}`, :optparam:`directory{rn}`, :optparam:`directorz{rn}`
    Components of the director vector, for shell formulations that need one.

.. _VariableCS:

Integral cross section with variable properties
-----------------------------------------------

The same integral cross section model, but with the individual properties given
as arbitrary functions of the global coordinates *x*, *y*, *z*.  As with
``SimpleCS``, for some problems the volume and cross section dimensions follow
from the element geometry and many parameters may be omitted.

.. record::

   :descitem:`VariableCS` :elemparam:`num{in}` :optelemparam:`material{in}`
   :optelemparam:`set{in}` :optelemparam:`thick{expr}`
   :optelemparam:`width{expr}` :optelemparam:`area{expr}`
   :optelemparam:`iy{expr}` :optelemparam:`iz{expr}`
   :optelemparam:`ik{expr}` :optelemparam:`shearareay{expr}`
   :optelemparam:`shearareaz{expr}`

The parameters have the same meaning as in :ref:`SimpleCS`, but each takes an
expression (see :ref:`attribute-types`) rather than a number.
:param:`shearareay` and :param:`shearareaz` determine the shear area required
by beam and plate elements.

.. _LayeredCS:

Layered cross section
---------------------

The layered cross section model, based on the geometrical hypothesis that cross
sections remain planar after deformation.

Elements using this cross section model must implement the layered cross
section extension; see the Element Library Manual.

.. record::

   :descitem:`LayeredCS` :elemparam:`num{in}` :elemparam:`nLayers{in}`
   :elemparam:`LayerMaterials{ia}` :elemparam:`Thicks{ra}`
   :elemparam:`Widths{ra}` :optelemparam:`set{in}`
   :optelemparam:`midSurf{rn}` :optelemparam:`nintegrationpoints{in}`
   :optelemparam:`layerintegrationpoints{ia}`
   :optelemparam:`beamshearcoeffxz{rn}`
   :optelemparam:`interfacematerials{ia}`
   :optelemparam:`rotations{ra}` :optelemparam:`initiationlimits{ra}`

**Parameters**

:param:`nLayers{in}`
    Number of layers.

:param:`LayerMaterials{ia}`
    Material number for each layer.

:param:`Thicks{ra}`, :param:`Widths{ra}`
    Thickness and width of each layer.

:optparam:`midSurf{rn}`
    Position of the mid surface, given as its distance from the bottom of the
    cross section; normal and moment forces are computed with respect to it.
    By default it lies at the average thickness position.

:optparam:`nintegrationpoints{in}`
    Number of integration points per layer.  Default is one.  The Gauss
    integration rule is used within each layer.

:optparam:`layerintegrationpoints{ia}`
    Number of integration points for each individual layer; its size should
    equal the number of layers.  This **overrides**
    :param:`nintegrationpoints`.

:optparam:`beamshearcoeffxz{rn}`
    Shear correction factor for 2D beam sections, controlling the effective
    shear area used to evaluate the shear force.  Default ``1.0``.

:optparam:`interfacematerials{ia}`
    Material numbers of the interfaces between the layers, for delamination
    modelling.

:optparam:`rotations{ra}`
    Rotation of each layer about the cross section normal, for laminates whose
    plies have different fibre orientations.

:optparam:`initiationlimits{ra}`
    Limits controlling the initiation of interface failure.

.. _FiberedCS:

Fibered cross section
---------------------

A cross section represented as a set of rectangular fibers.  It rests on the
geometrical hypothesis that cross sections remain planar after deformation, and
is the 3D generalization of the layered approach for beams.

.. important::

   The geometry of the cross section is determined from the fiber dimensions
   and fiber positions, all given in the local coordinate system of the beam,
   that is in the *yz* plane.

.. record::

   :descitem:`FiberedCS` :elemparam:`num{in}` :elemparam:`nfibers{in}`
   :elemparam:`fibermaterials{ia}` :elemparam:`thicks{ra}`
   :elemparam:`widths{ra}` :elemparam:`thick{rn}` :elemparam:`width{rn}`
   :elemparam:`fiberycentrecoords{ra}`
   :elemparam:`fiberzcentrecoords{ra}` :optelemparam:`set{in}`

**Parameters**

:param:`nfibers{in}`
    Number of fibers forming the overall cross section.

:param:`fibermaterials{ia}`
    Material model number for each fiber.  The material model specified at
    element level has no meaning in this case.

:param:`thicks{ra}`, :param:`widths{ra}`
    Thickness and width of each fiber.

:param:`thick{rn}`, :param:`width{rn}`
    Overall thickness and width of the cross section.

:param:`fiberycentrecoords{ra}`, :param:`fiberzcentrecoords{ra}`
    *y* and *z* coordinates of the centre of each fiber.

.. _WarpingCS:

Warping cross section
---------------------

Cross section for free warping analysis; see :ref:`FreeWarping`.

.. record::

   :descitem:`WarpingCS` :elemparam:`num{in}`
   :elemparam:`WarpingNode{in}` :optelemparam:`set{in}`

**Parameters**

:param:`WarpingNode{in}`
    Number of the external node carrying the prescribed boundary condition
    that corresponds to the relative twist of the warping cross section.

.. _LatticeCS:

Lattice cross section
---------------------

.. note::

   Part of the lattice module; compile with ``USE_LM``.

Cross section for the 3D lattice elements ``lattice3d`` and ``lattice3dnl``.

.. record::

   :descitem:`LatticeCS` :elemparam:`num{in}` :elemparam:`material{in}`
   :optelemparam:`set{in}` :optelemparam:`shape{in}`
   :optelemparam:`radius{rn}` :optelemparam:`area{rn}`
   :optelemparam:`thickness{rn}` :optelemparam:`iy{rn}`
   :optelemparam:`iz{rn}` :optelemparam:`ik{rn}`
   :optelemparam:`shearcoeff{rn}` :optelemparam:`shearareay{rn}`
   :optelemparam:`shearareaz{rn}` :optelemparam:`nlayers{in}`

**Parameters**

:param:`material{in}`
    Number of the material record.

:optparam:`shape{in}`
    How the section properties are obtained.  Default ``0``.

    ``0``
        A general polygon or facet built from the element vertex coordinates
        (``polycoords``).
    ``1``
        A circle, from :param:`radius`.
    ``2``
        A rectangle, from four vertices.
    ``3``
        A property-defined section, given directly with no geometry.

    For the geometric shapes ``0``, ``1`` and ``2`` the properties are computed
    from the geometry, and any individual property may be overridden by
    supplying its value here.

:optparam:`radius{rn}`
    Radius, used by ``shape 1``.

:optparam:`area{rn}`
    Cross section area.  **Mandatory** for the property-defined section
    (``shape 3``).

:optparam:`thickness{rn}`
    Cross section thickness.

:optparam:`iy{rn}`, :optparam:`iz{rn}`
    Second moments of area about the local *y* and *z* axes.

:optparam:`ik{rn}`
    Saint-Venant torsional constant.

:optparam:`shearcoeff{rn}`
    Shear correction factor.  Alternatively give :param:`shearareay` and
    :param:`shearareaz`.

:optparam:`shearareay{rn}`, :optparam:`shearareaz{rn}`
    Equivalent shear areas.

:optparam:`nlayers{in}`
    Number of through-thickness layers of a shell.  Shell mode is activated by
    the element's ``shellnormal``; the default is a single integration point at
    the centroid.

For the property-defined section (``shape 3``) the element represents a two-node
frame member — rod, beam and so on.  :param:`area` is then mandatory and the
remaining properties optional; the transverse axes are oriented by the element's
``zaxis``, and the integration point position along the element by its ``s``
parameter.  See the Element Library Manual.

Other cross section types
-------------------------

``SimpleTransportCS``
    Cross section for transport problems.  Takes
    :elemparam:`mat{in}`, :optelemparam:`thickness{rn}` and
    :optelemparam:`area{rn}`.

``FluidCS``
    Cross section for fluid dynamic problems.  Takes :elemparam:`mat{in}`.

``InterfaceCS``
    Cross section for structural interface elements.  Takes
    :elemparam:`material{in}` and :optelemparam:`thickness{rn}`.

``EmptyCS``, ``DummyCS``
    Placeholder cross sections carrying no properties, used where an element
    requires a cross section but none of its properties are needed.

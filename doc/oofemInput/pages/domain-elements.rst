.. _ElementsRecords:

Element records
===============

These records describe the individual elements.

.. record::

   :descitem:`ElementType` :elemparam:`num{in}` :elemparam:`nodes{ia}`
   :elemparam:`crossSect{in}` :optelemparam:`mat{in}`
   :optelemparam:`bodyLoads{ia}` :optelemparam:`boundaryLoads{ia}`
   :optelemparam:`activityltf{in}` :optelemparam:`lcs{ra}`
   :optelemparam:`nip{in}`
   ``<``\ :optelemparam:`remote{}`\ ``>``
   ``<``\ :optelemparam:`partitions{ia}`\ ``>``

The order of the records is arbitrary and the element number is given by
:param:`num`.  Numbering may be non-contiguous: think of the number as a label
assigned to an individual element, by which it is referenced.

.. note::

   The element types themselves, and the parameters specific to each, are
   documented in the separate **Element Library Manual**.  This section covers
   only the attributes common to all elements.

**Parameters**

:param:`nodes{ia}`
    The dof managers — nodes, sides and so on — defining the element geometry.

:param:`crossSect{in}`
    Number of the cross section record for the element.  Required — but it need
    not be given in the element record itself: a cross section record that
    names a ``set`` assigns itself to every element of that set, which is the
    usual style in current input files.  If the element ends up with no cross
    section from either source, the analysis stops.

:optparam:`mat{in}`
    Number of the material record for the element.  Usually supplied by the
    cross section instead (see :ref:`CrossSectionRecords`), so it is optional
    here.

:optparam:`nip{in}`
    Number of integration points.  The permitted values, and the default,
    depend on the element type.  Note that adaptive error estimation often
    needs a higher integration order than the default; see
    :ref:`errorestimators`.

:optparam:`lcs{ra}`
    Orientation of the local coordinate system.  The array holds six numbers:
    the first three are a direction vector of the local *x*-axis, the next
    three a direction vector of the local *y*-axis, and the local *z*-axis
    follows from the vector product.  This is particularly useful for
    orthotropic materials that follow the element orientation.  On a beam or
    truss element :param:`lcs` has no effect, and the 1D element orientation is
    aligned with the global :math:`xx` component.

Loads
-----

:optparam:`bodyLoads{ia}`
    Body loads acting on the element volume, given as load record numbers.  The
    records must be of a body load type, or an error is raised.

:optparam:`boundaryLoads{ia}`
    Boundary loads acting on the element boundary.  The array has the layout

    .. math:: 2\cdot size \; lnum(1)~id(1)~\dots~lnum(size)~id(size),

    where :math:`size` is the total number of loads applied to the element,
    :math:`lnum(i)` is the applied load number, and :math:`id(i)` the
    corresponding entity number to which the load is applied, for example a
    side or surface number.  The entity numbering is element dependent and is
    described in the Element Library Manual.  The records must be of a boundary
    load type, or an error is raised.

Element activity
----------------

:optparam:`activityltf{in}`
    Number of a time function controlling whether the element is active.  A
    nonzero value of the function means the element is active — the default — and
    zero means it is inactive at that solution step.

This supports the insertion and removal of elements during the analysis, and is
tested for structural and transport elements.  It makes it possible to model,
for instance, the temperature evolution of layered casting of concrete, where
certain layers must be inactive before they are cast.  See the OOFEM tests for
an example of a hydrating material model, boundary conditions and element
activity acting together.

.. _parallel-elements:

Parallel mode
-------------

.. note::

   This section applies only to the parallel build.  See
   :ref:`parallel-model` for the overall picture.

:optparam:`remote{}`
    Forces the element to be a remote element.  A remote element does not
    contribute to the governing equations of the local partition.  Remote
    elements exist to implement the band of elements involved in the
    computation of nonlocal variables — see :numref:`fig-nodecut-nlm`, which
    illustrates the approach for node-cut partitioning.  Their role is to
    provide a local mirror of the integration point values of the
    corresponding element on the remote partition, which then undergo nonlocal
    averaging on the local partition.  When this flag is not used, the element
    is assumed to be a local partition element.

:optparam:`partitions{ia}`
    When :param:`remote` is used, this should contain the number of the remote
    partition on which the corresponding element is local; the array then has
    size one.

The parallel version also requires a unique global numbering of the elements
across all partitions; see :ref:`parallel-model`.

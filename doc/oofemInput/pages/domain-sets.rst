.. _SetRecords:

Set records
===========

A set specifies a region of the geometry as a combination of volumes, surfaces,
edges and nodes.  Sets are used above all to attach a region of elements to a
cross section, or to apply a boundary condition, but they have many other uses.

.. record::

   :descitem:`Set` :elemparam:`num{in}` :optelemparam:`elements{ia}`
   :optelemparam:`elementranges{rl}` :optelemparam:`allElements{}`
   :optelemparam:`nodes{ia}` :optelemparam:`noderanges{rl}`
   :optelemparam:`allNodes{}` :optelemparam:`elementboundaries{ia}`
   :optelemparam:`elementedges{ia}` :optelemparam:`elementsurfaces{ia}`
   :optelemparam:`internalelementnodes{ia}`
   ``<ver. 1.6>`` :optelemparam:`dofmanprops{s}`
   :optelemparam:`elemprops{s}`

**Volumes and nodes**

:optparam:`elements{ia}`, :optparam:`nodes{ia}`
    Explicit lists of element, respectively node, numbers.

:optparam:`elementranges{rl}`, :optparam:`noderanges{rl}`
    The same, given as range lists; see :ref:`attribute-types`.

:optparam:`allElements{}`, :optparam:`allNodes{}`
    Flags selecting every element, respectively every node, in the domain.

**Edges and surfaces**

:optparam:`elementedges{ia}`, :optparam:`elementboundaries{ia}`, :optparam:`elementsurfaces{ia}`
    Interleaved lists in which every other number is the element number and the
    one after it the edge, boundary or surface number.  The total length of the
    list is therefore twice the number of entities.  The internal numbering of
    edges and surfaces is given in the **Element Library Manual**.

:optparam:`internalelementnodes{ia}`
    Interleaved list of element number and internal node number, for elements
    that manage internal dof managers.

Which entity kind to use
------------------------

The entity kind must match the dimensionality of what is applied to it:

* bulk loads, in newtons per volume, belong on :param:`elements`;
* surface loads, in newtons per area, belong on :param:`elementboundaries`;
* edge loads — singular loads in newtons per length, or the equivalent — belong
  on :param:`elementedges`.

**Example 1.**  A deadweight (gravity) load is applied to the
:param:`elements` of a set, while a distributed line load is applied to the
midline "edge" of a beam element and therefore belongs to an
:param:`elementedges` set.  In the latter case the midline of the beam is
defined as the first — and only — "edge" of the beam.

**Example 2.**  In an axisymmetric structural analysis a deadweight load is
applied to the :param:`elements` of a set, while an external pressure is
defined as a surface load and applied to the :param:`elementboundaries` of a
set.  The element integrates the load analytically around the axis, so it still
counts as a surface load.

Properties set through a set
----------------------------

Arbitrary dof manager and element properties can be assigned through a set
(supported from version 1.6).

:optparam:`dofmanprops{s}`
    A string of dof manager parameters, applied to every dof manager in the
    set.

:optparam:`elemprops{s}`
    A string of element parameters, applied to every element in the set.

The syntax inside these strings is the ordinary input record syntax.

.. note::

   Parameters defined this way have the **lowest** priority: they are
   overridden by the parameters given in the individual dof manager or element
   records.  See ``tests/regression/sm/setprops01.in`` for an example.

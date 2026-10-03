Material and nonlocal barrier records
=====================================

.. _MaterialTypeRecords:

Material type records
---------------------

These records describe the material models.

.. record::

   :descitem:`MaterialType` :elemparam:`num{in}` :elemparam:`d{rn}`

The order of the records is arbitrary and the material number is given by
:param:`num`.  The numbering must start at one and end at *n*, where *n* is the
number of records.

**Parameters**

:param:`d{rn}`
    Material density.  Compulsory for every material model.

.. note::

   The available material models, their formulation and their own parameters
   are documented in the separate **Material Library Manual**.  This section
   covers only what is common to every material record.

.. _NonlocalBarrierRecords:

Nonlocal barrier records
------------------------

Nonlocal material models of the integral type replace a suitable local quantity
in the local constitutive law by its nonlocal counterpart, obtained as a
weighted average over some characteristic volume.  The weighted average is a sum
of the remote value multiplied by the value of a weight function, which
typically depends on the distance between the remote and the receiver point and
decreases as that distance grows.

Sometimes the mutual interaction between two points has to be disregarded — for
example when they lie on opposite sides of a thin notch, which prevents the
nonlocal interaction from taking place.  Barriers are the way to introduce such
a constraint.  A barrier is a curve in 2D or a surface in 3D; when the line
connecting the receiver and the remote point intersects a barrier, the barrier
is activated and the corresponding interaction is not taken into account.

The number of barrier records is declared by :param:`nbarrier` in the
:ref:`ComponentsSizeRecord`.

Polyline barrier
~~~~~~~~~~~~~~~~

A polyline barrier for 2D problems, defined as a sequence of nodes representing
its vertices.

.. record::

   :descitem:`polylinebarrier` :elemparam:`num{in}`
   :elemparam:`vertexnodes{ia}` :optelemparam:`xcoordindx{in}`
   :optelemparam:`ycoordindx{in}`

**Parameters**

:param:`vertexnodes{ia}`
    Numbers of the nodes representing the polyline vertices.

:optparam:`xcoordindx{in}`, :optparam:`ycoordindx{in}`
    Select the plane — *xy*, *yz* or *xz* — in which the barrier is defined.
    :param:`xcoordindx` is the first coordinate index and
    :param:`ycoordindx` the second.  The defaults are ``1`` and ``2``,
    giving a barrier in the *xy* plane.

Symmetry barrier
~~~~~~~~~~~~~~~~

An implementation of a symmetry barrier, allowing up to three mutually
orthogonal planes of symmetry.  It models the symmetry of the averaged field on
the boundary without the other part of the structure across the plane of
symmetry having to be modelled, by modifying the integration weights of the
source points to account for the symmetry.

The potential symmetry planes are given by an orthogonal, right-handed
coordinate system whose axes are the normals of the corresponding symmetry
planes.

.. record::

   :descitem:`symmetrybarrier` :elemparam:`num{in}`
   :elemparam:`origin{ra}` :elemparam:`normals{ra}`
   :elemparam:`activemask{ia}`

**Parameters**

:param:`origin{ra}`
    Origin of the coordinate system.

:param:`normals{ra}`
    Three components of the *x*-axis direction vector followed by three
    components of the *y*-axis direction vector, both in the global coordinate
    system.  The *z*-axis follows from the orthogonality conditions.

:param:`activemask{ia}`
    Which symmetry planes are active.  The *i*-th nonzero value activates the
    symmetry barrier for the plane whose normal is the corresponding coordinate
    axis: *x* = 1, *y* = 2, *z* = 3.

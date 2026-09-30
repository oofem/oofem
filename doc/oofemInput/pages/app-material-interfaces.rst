.. _materialinterfaces:

Material interfaces
===================

Material interfaces represent and track the position of various interfaces on
fixed grids — a free surface, or an evolving interface between two materials,
for example.  The representation is selected by the :param:`miflag` parameter
of the analysis record; see :ref:`supgIncomp`.

.. table:: Available material interface representations.
   :name: materialinterfacetable

   +----------+--------------+---------------+
   | MI       | ``miflag``   | Compatibility |
   +==========+==============+===============+
   | LEPlic   | 0            | 2D triangular |
   +----------+--------------+---------------+
   | LevelSet | 1            | 2D triangular |
   +----------+--------------+---------------+

LEPlic
------

A representation based on the Volume-Of-Fluid approach.  The initial
distribution of VOF fractions should be specified for each element; see the
Element Library Manual.

.. record:: Syntax

   :optelemparam:`refvol{rn}`

**Parameters**

:optparam:`refvol{rn}`
    Initial volume of the reference fluid.  When it is given, the reference
    volume is computed and printed in every step, so that accuracy and mass
    conservation can be monitored.

LevelSet
--------

A level-set based representation.

.. record:: Syntax

   (:elemparam:`levelset{ra}` |
   (:elemparam:`refmatpolyx{ra}` :elemparam:`refmatpolyy{ra}`))
   :optelemparam:`lsra{in}` :optelemparam:`rdt{rn}`
   :optelemparam:`rerr{rn}`

**Parameters**

:param:`levelset{ra}`
    Initial level set values for all nodes, given directly.  Its size should
    equal the total number of nodes in the domain.

:param:`refmatpolyx{ra}`, :param:`refmatpolyy{ra}`
    Alternatively, initialize the level set from the interface geometry given
    as a 2D polygon: the polygon describes the initial zero level set, and the
    level set values are then the signed distance from it, positive on the left
    side when walking along the polygon.  :param:`refmatpolyx` holds the *x*
    coordinates of the polygon vertices and :param:`refmatpolyy` the *y*
    coordinates.

    .. important::

       The level set must be initialized, either with :param:`levelset` or
       with :param:`refmatpolyx` and :param:`refmatpolyy`.

:optparam:`lsra{in}`
    Level set reinitialization algorithm: ``0`` no reinitialization; ``1``
    (default) reinitializes the level set representation by solving
    :math:`d_{\tau} = S(\phi)(1-\vert\boldsymbol{\nabla}d\vert)` to steady
    state; ``2`` uses the fast marching method to build a signed-distance level
    set representation.

:optparam:`rdt{rn}`
    Time step of the reinitialization integration algorithm used by
    ``lsra 1``.  Default ``0.0``, for which the algorithm chooses the step
    itself; any positive value overrides that choice.

:optparam:`rerr{rn}`
    Error limit used to detect the steady state of the ``lsra 1``
    reinitialization algorithm.  Default ``1.e-6``.

:optparam:`nsd{in}`
    Number of spatial dimensions of the level set.

:optparam:`ci1{rn}`, :optparam:`ci2{rn}`
    Coefficients of the interface reconstruction.

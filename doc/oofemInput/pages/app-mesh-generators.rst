.. _meshpackages:

Mesh generator interfaces
=========================

The mesh generator interface provides the link to a specific mesh generator.  It
is selected by the :param:`meshpackage` parameter of an adaptive analysis; see
:ref:`AdaptiveLinearStatic` and :ref:`AdaptiveNonLinearStatic`.

``MPT_T3D`` — :param:`meshpackage` = 0
    T3d mesh interface.  The default.  Supports 1D, 2D (triangular) and 3D
    (tetrahedral) meshes.  Reliable.

``MPT_TARGE2`` — :param:`meshpackage` = 1
    Interface to the Targe2 2D mesh generator.

``MPT_SUBDIVISION`` — :param:`meshpackage` = 3
    Built-in subdivision algorithm.  Supports triangular 2D and tetrahedral 3D
    meshes, and can operate in parallel mode.

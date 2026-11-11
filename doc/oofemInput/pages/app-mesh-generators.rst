.. _meshpackages:

Mesh generator interfaces
=========================

The mesh generator interface provides the link to a specific mesh generator.  It
is selected by the :param:`meshpackage` parameter of an adaptive analysis; see
:ref:`AdaptiveLinearStatic` and :ref:`AdaptiveNonLinearStatic`.  The parameters
of the interface (including :param:`remeshcmd` and the Gmsh parameters below)
are given in the analysis record; they are supported by both adaptive analyses.

All mesher-specific details are hidden behind the interface.  When asked for a
new mesh, the interface either returns the new mesh directly (built-in
subdivision, or an external mesh generator run by the interface), or writes the
required mesh density into a mesher-specific file only, leaving the remeshing to
the user.

``MPT_T3D`` — :param:`meshpackage` = 0
    T3d mesh interface.  The default.  Supports 1D, 2D (triangular) and 3D
    (tetrahedral) meshes.  Writes the background mesh file ``t3d.bmf``.

``MPT_TARGE2`` — :param:`meshpackage` = 1
    Interface to the Targe2 2D mesh generator.  Writes ``targe2.bmf``.

``MPT_FREEM`` — :param:`meshpackage` = 2
    Interface to the Freem 2D mesh generator.  Writes ``freem.bmf``.

``MPT_SUBDIVISION`` — :param:`meshpackage` = 3
    Built-in subdivision algorithm.  Supports triangular 2D and tetrahedral 3D
    meshes, and can operate in parallel mode.  Creates the new mesh directly.

``MPT_GMSH`` — :param:`meshpackage` = 4
    Interface to the `Gmsh <https://gmsh.info>`_ mesh generator.  The required
    element size is written as a Gmsh background mesh (post-processing view,
    file ``<output>.bgm.<n>.pos``), built from the current mesh (triangles,
    quads, tetrahedra, hexahedra and wedges are supported), so the current mesh
    can come from any source.  When the geometry and control files are given,
    the interface runs Gmsh and converts the mesh to the new oofem domain using
    ``tools/unv2oofem/gmsh2oofem.py``; the properties (cross sections,
    materials, boundary conditions, loads and sets) are assigned to the mesh
    based on Gmsh physical groups, as described in the unv2oofem control file.
    The same script generates the initial mesh::

        gmsh2oofem.py [--dim 3] [--elemtype simplex|quad|hex] [--setnumber lc 0.1] model.geo model.ctrl model.in

    (``--setnumber`` sets a geometry parameter declared by ``DefineConstant``
    in the geometry file, e.g. the initial mesh size.)

    Parameters:

    :optparam:`gmshgeo{s}`
        Gmsh geometry file (``.geo``).
    :optparam:`gmshctrl{s}`
        unv2oofem control file.
    :optparam:`gmshopts{s}`
        Additional options of ``gmsh2oofem.py``, for example
        ``"--elemtype hex"`` (hexahedral mesh obtained by subdivision of
        tetrahedra) or ``"--option Mesh.Algorithm=6"``.
    :optparam:`gmshgrading{rn}`
        Limits the gradation of the required element size field written to the
        background mesh: the size at any node is reduced to at most
        :math:`h_j + g\,|x_i - x_j|` for every node :math:`j` of the current mesh,
        so that the transition from fine to coarse mesh is gradual.  Gmsh itself
        provides no such control for background meshes.  Default ``0`` (not
        limited); values of 0.1–0.25 are typical.
    :optparam:`gmshcmd{s}`
        Gmsh executable, default ``gmsh``.
    :optparam:`pythoncmd{s}`
        Python interpreter, default ``python3``.
    :optparam:`gmsh2oofem{s}`
        Path to the conversion script, by default the script in the oofem
        source tree.

External remeshing
------------------

The interfaces writing a mesh density file (T3d, Targe2, Freem, Gmsh) can run an
external mesh generator and converter themselves, when the command is given by
the :param:`remeshcmd` parameter.  The command has to create the new oofem
domain file (the oofem input file starting with the ``domain`` record).  The
following placeholders are substituted before the command is executed:

``%d``
    the mesh density file written by the interface,
``%o``
    the domain file to be created (``<output>.domain.1.<n>.din``),
``%s``
    the serial number of the new domain,
``%n``
    the solution step number,
``%%``
    the percent character.

For example, with T3d::

    remeshcmd "t3d -i job.t3d -m %d -o job.t3d.out && t3d2oofem job.ctrl job.t3d.out job.in && tail -n +4 job.in > %o"

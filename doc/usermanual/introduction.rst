Introduction
============

OOFEM is a general purpose, object-oriented finite element code for solving
mechanical, transport and fluid mechanics problems, primarily focused on
academic and research use. It has been continuously developed since 1997 and at
present consists of more than half a million lines of C++ source code. It has
an active development team and a large user base around the world, and the code
has been successfully applied to the solution of several industrial problems.

OOFEM aims to be one of the best FEM solvers. To make this a reality, we focus
on solving and do not try to provide extensive pre- and post-processing
capabilities; there we rely on third party external tools.

As a FEM solver, OOFEM solves problems described by a set of partial
differential equations. On input it expects the discretized domain and the
corresponding problem parameters, initial and boundary conditions. On output it
provides the primary unknown fields as well as other secondary fields of
interest. The results can be exported in many ways to facilitate
post-processing.

OOFEM is free software, distributed under the GNU LGPL license.


General Features
----------------
* Problem modules: structural mechanics (``sm``), heat and mass transfer
  (``tm``), fluid mechanics (``fm``), lattice models (``lm``), additive
  manufacturing (``am``) and the experimental symbolic multiphysics module
  (``mpm``)
* Support for parallel processing on shared memory (OpenMP) and distributed
  memory (MPI) machines and computer clusters, including dynamic load balancing
* Direct and iterative solvers (built-in sparse solvers, IML++, SPOOLES,
  PETSc, SLEPc, SuperLU_MT, Pardiso/MKL)
* Full restart capability, support for adaptive and staggered analyses
* Post-processing: export to VTK (ParaView), and further export modules for
  Matlab, gnuplot, homogenized quantities and integration point data
* Python bindings allowing to script OOFEM and to implement new components in
  Python


.. _documentation:

Documentation
---------------------
This user manual is only the starting point. The individual records that make
up an input file — the available analyses, elements, materials, boundary
conditions and export modules — are catalogued in the companion manuals. The
full documentation set is available on the
`OOFEM web pages <https://www.oofem.org>`_:

* Input manual, available in [`html <http://www.oofem.org/resources/doc/oofemInput/html/index.html>`__] and [`PDF <http://www.oofem.org/resources/doc/oofemInput/oofeminputmanual.pdf>`__].
* Element Library Manual, available in [`html <http://www.oofem.org/resources/doc/elementlibmanual/html/elementlibmanual.html>`__] and [`PDF <http://www.oofem.org/resources/doc/elementlibmanual/elementlibmanual.pdf>`__].
* Material Model Library Manual, available in [`html <http://www.oofem.org/resources/doc/matlibmanual/html/matlibmanual.html>`__] and [`PDF <http://www.oofem.org/resources/doc/matlibmanual/matlibmanual.pdf>`__].
* Theory manual, available in [`html <http://www.oofem.org/resources/doc/theorymanual/html/index.html>`__].
* Developer guide (formerly the Programmer's manual), available in [`html <http://www.oofem.org/resources/doc/developerguide/html/index.html>`__].
* The [`C++ reference manual <http://www.oofem.org/resources/doc/oofemrefman/index.html>`__], generated from the source code.
* Python bindings documentation, available in [`html <http://www.oofem.org/resources/doc/python/html/index.html>`__].


OOFEM Ecosystem
---------------------

* The `OOFEM forum <https://www.oofem.org/forum/>`_ is the best place to ask
  questions and get support from the developers and the user community.
* The `OOFEM web site <https://www.oofem.org>`_ contains many useful resources,
  a gallery of results and tutorials.
* The `OOFEM GitHub repository <https://github.com/oofem/oofem>`_ holds the
  sources, the issue tracker and the continuous integration builds.

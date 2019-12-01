Introduction
============

OOFEM is a general purpose, object-oriented Finite Element code for solving mechanical. transport and fluid mechanics problems, primarily focused on academic/research use. It has been continuously developed since 1997. At present it consist of more than 230 thousands lines of source C++ code. It has acive development team and large user base around the world. The code has been successfully applied to solution of several industrial problems.

OOFEM is trying to be one of the best FEM solvers. To make this a reality, we focus on this and do not try to provide extensive pre and post processing capabilities. Here we rely on third party external tools.

OOFEM as a FEM solver solves the problems described by a set of partial differential equations. On input it expects the discretized domain and corresponding problem parameters, initial and boundary conditions. On output, it provides the primary unknown fields as well as other secondary fields of interest. The results can be exported in many ways to facilitate post-processing.

OOFEM is free software, distributed under GNU LGPL license.



General Features
----------------
* Structural, heat & mass transfer, and CFD modules
* Support for parallel processing on shared and massively parallel machines, computer clusters, Dynamic load balancing
* Direct & iterative solvers (IML, Spooles, PETSc, SuperLU, Pardisso)
* Full restart capability, support for adaptive and staggered analyses
* Postprocessing: X-Windows post-processing, export to VTK
* Python bindings allowing to script oofem and implement new components in Python


Documentation
---------------------
The extensive documentation is available on `OOFEM web pages <http://www.oofem.org>`_: 

* Input manual, available in [`html <http://www.oofem.org/resources/doc/oofemInput/html/oofemInput.html>`_] and [`PDF <http://www.oofem.org/resources/doc/oofemInput/html/oofemInput.pdf>`_].
* Element Library Manual, available in [`html <http://www.oofem.org/resources/doc/elementlibmanual/html/elementlibmanual.html>`_] and [`PDF <http://www.oofem.org/resources/doc/elementlibmanual/elementlibmanual.pdf>`_] 
* Material Model Library Manual, available in [`html <http://www.oofem.org/resources/doc/matlibmanual/html/matlibmanual.html>`_] and [`PDF <http://www.oofem.org/resources/doc/matlibmanual/matlibmanual.pdf>`_] 

* Programmer's manual, available in [`html <http://www.oofem.org/resources/doc/programmer/html/programmer.html>`_] and [`PDF <http://www.oofem.org/resources/doc/programmer/programmer.pdf>`_].
* The [`Reference manual <http://www.oofem.org/resources/doc/oofemrefman/index.html>`_] 
* Python bindings documentation , available in [`html <http://www.oofem.org/resources/doc/python/html/index.html>`_] and [`PDF <http://www.oofem.org/resources/doc/python/oofempythonbindings.pdf>`_]


OOFEM Ecosystem
---------------------

* `OOFEM forum <http://www.oofem.org/forum>`_ is the best place to ask questions and get support from developpers and user community.
* `OOFEM wiki <http://www.oofem.org/wiki>`_ contains many useful resources, gallery of results and tutorials.
* `OOFEM gitHub repository <https://github.com/oofem/oofem>`_.
.. _installation:

Installation
============

Installation options
---------------------
* Official stable releases (https://www.oofem.org/en/download)

  * Usually a 12 month release cycle
  * Binary packages for Windows (AMD64 version only) and Linux (x86_64
    version, Debian package)
  * Source package (requires compilation)

* Development version

  * Manual installation from the OOFEM GitHub repository
    (https://github.com/oofem/oofem.git). The default branch is ``devel``.

* Python wheel, if you only intend to use OOFEM from Python

  .. code-block:: bash

      $ pip install oofem

In the examples below, ``X.Y`` stands for the release version (3.0 at the time
of writing); substitute the version of the package you actually downloaded.


Binary package installation - Windows
-------------------------------------
* Extract the downloaded zip archive (``oofem_X.Y_AMD64.zip``) into any
  directory
* The extraction should create an ``oofem_X.Y_AMD64`` directory
* Modify the ``PATH`` variable to include the ``oofem_X.Y_AMD64\lib`` directory

  .. code-block:: bat

      set PATH=C:\Users\user\Documents\oofem_3.0_AMD64\lib;%PATH%

* Test run

  .. code-block:: bat

      C:\Users\user\Documents\oofem_3.0_AMD64\bin\oofem -v

  which prints the version, the configured modules and the exact commit the
  binary was built from:

  .. code-block:: text

      OOFEM version 3.0 (AMD64-Windows, fm tm sm sm_elements sm_materials mpm IML++)
      Git RepoURL: https://github.com/oofem/oofem.git
          Branch: devel, Hash: v3.0

      Copyright (C) 1994-2025 Borek Patzak
      This is free software; see the source for copying conditions.  There is NO
      warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Binary installation - Unix/Linux
--------------------------------
You can either choose the generic Linux binary or the Debian package.

Linux - binary
^^^^^^^^^^^^^^
* Extract the downloaded archive (``oofem_X.Y_x86_64.tar.gz``)

  .. code-block:: bash

      $ tar xzvf oofem_3.0_x86_64.tar.gz

* This creates the ``oofem_X.Y_x86_64/bin`` and ``oofem_X.Y_x86_64/lib``
  directories
* Tell the dynamic linker where to find the shared libraries

  .. code-block:: bash

      $ export LD_LIBRARY_PATH=../lib:$LD_LIBRARY_PATH

* Run oofem

  .. code-block:: bash

      $ cd bin; ./oofem -v

Linux - Debian package
^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: bash

    $ sudo apt install ./oofem_3.0_x86_64.deb

(``apt`` needs the leading ``./`` to recognize the argument as a local file
rather than a package name; ``sudo dpkg -i oofem_3.0_x86_64.deb`` works too but
does not pull in dependencies.)


Installation from source
------------------------
Requirements

* A C++ compiler with C++17 support (g++, clang++, Visual Studio, MinGW)
* CMake 3.15 or newer (https://cmake.org/)

CMake controls the software compilation process using a platform independent
configuration and generates native makefiles or projects. To generate the
platform makefile or project configuration, use ``cmake`` (or any of the user
friendly GUIs based on it, such as ``ccmake`` or ``cmake-gui``).

Everything else is optional: the default configuration builds the ``sm``,
``tm``, ``fm`` and ``mpm`` modules together with the built-in solvers and
IML++, and needs no external library beyond a compiler and CMake.

Installation from source - Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the following, we assume the sources are in ``/home/user/oofem.src`` and the
build directory is ``/home/user/oofem.build``. Keep the two separate: an
out-of-source build leaves the sources untouched and lets you keep several
differently configured builds of the same checkout side by side.

.. code-block:: bash

    $ git clone https://github.com/oofem/oofem.git /home/user/oofem.src
    $ mkdir -p /home/user/oofem.build
    $ cd /home/user/oofem.build
    $ ccmake /home/user/oofem.src      # or: cmake /home/user/oofem.src
    $ make -j

.. note::

   When ``CMAKE_BUILD_TYPE`` is not given, CMake configures a **Debug** build,
   which is considerably slower. For production runs configure with
   ``-DCMAKE_BUILD_TYPE=Release``.

Configuration options
^^^^^^^^^^^^^^^^^^^^^
All options are set either interactively in ``ccmake``/``cmake-gui`` or on the
``cmake`` command line as ``-DOPTION=ON``. The most frequently used ones are:

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Option
     - Default
     - Meaning
   * - ``CMAKE_BUILD_TYPE``
     - ``Debug``
     - Build type; use ``Release`` for production runs
   * - ``USE_SM``
     - ``ON``
     - Structural mechanics module
   * - ``USE_TM``
     - ``ON``
     - Heat and mass transfer module
   * - ``USE_FM``
     - ``ON``
     - Fluid mechanics module
   * - ``USE_MPM``
     - ``ON``
     - Experimental symbolic multiphysics module
   * - ``USE_LM``
     - ``OFF``
     - Lattice module (turns on ``USE_SM`` and ``USE_TM``)
   * - ``USE_AM``
     - ``OFF``
     - Additive manufacturing module (turns on ``USE_SM`` and ``USE_TM``)
   * - ``USE_XML``
     - ``ON``
     - XML input format support (requires pugixml)
   * - ``USE_IML``
     - ``ON``
     - IML++ iterative solvers
   * - ``USE_PETSC``
     - ``OFF``
     - PETSc solvers; required for distributed memory parallel runs
   * - ``USE_SLEPC``
     - ``OFF``
     - SLEPc eigenvalue solvers (turns on ``USE_PETSC``)
   * - ``USE_OPENMP_PARALLEL``
     - ``OFF``
     - Shared memory parallelism, enables the ``-t`` command line option
   * - ``USE_MPI_PARALLEL``
     - ``OFF``
     - Distributed memory parallelism, enables the ``-p`` command line option
   * - ``USE_PYBIND_BINDINGS``
     - ``OFF``
     - Build the ``oofempy`` Python module (see :ref:`python`)
   * - ``USE_PYTHON_EXTENSION``
     - ``OFF``
     - Allow input files to call user code written in Python
   * - ``USE_SHARED_LIB``
     - ``ON``
     - Build ``liboofem`` as a shared library
   * - ``USE_VTK``
     - ``OFF``
     - Link VTK, needed for binary (rather than ASCII) VTU export

For the complete list, including the ``*_DIR`` variables used to point CMake at
manually installed third party libraries, inspect the options in ``ccmake`` or
read the ``Options`` section of the top-level ``CMakeLists.txt``.

Verifying the build
^^^^^^^^^^^^^^^^^^^
OOFEM comes with an extensive regression test suite driven by CTest. Each test
is an input file carrying the expected results, which are checked by the
``errorcheck`` export module. After building, run the suite from the build
directory:

.. code-block:: bash

    $ ctest -j4                # run all tests
    $ ctest -L sm              # only the structural mechanics tests
    $ ctest -R beam2d_1        # only tests matching a name
    $ ctest --output-on-failure # show the solver output of failing tests

All tests are expected to pass; a failure usually indicates a misconfigured
third party library. The available labels are ``sm``, ``tm``, ``fm``, ``tmsm``,
``tmfm``, ``mpm``, ``benchmark`` and ``workflow``.

Installation from source - Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Requirements:

* C++ compiler (Visual Studio, MinGW)
* CMake build tools (https://cmake.org/)

Procedure:

* Clone the oofem git repository
* Use cmake to generate the VS solution

  .. image:: figs/Installation_win_cmake.png
      :width: 50%
      :align: center
      :alt: screenshot of cmake gui
      :target: _images/Installation_win_cmake.png

  * Select the compiler
  * Set the source directory
  * Set the build directory
  * Generate the project/solution

* Use the compiler to build the project targets (``oofem`` and ``RUN_TESTS``)

  .. image:: figs/Installation_win_VisualC2.png
      :width: 65%
      :align: center
      :alt: screenshot of the Visual Studio solution explorer
      :target: _images/Installation_win_VisualC2.png

  .. image:: figs/Installation_win_VisualC.png
      :width: 65%
      :align: center
      :alt: screenshot of the Visual Studio build output
      :target: _images/Installation_win_VisualC.png

.. _installation:

Installation
============

Installation options
---------------------
* Official stable releases (http://www.oofem.org/en/download)
    * Usually 12M release cycle
    * Binary packages for Windows (AMD64 version only) and Linux (x86_64 version, Debian package)
    * Source package (requires compilation)

• Development version
    * Manual installation from OOFEM GitHub repository (https://github.com/oofem/oofem.git)

Binary package installation - Windows
-------------------------------------
* Extract downloaded zip archive (oofem_2.5_AMD64.zip) into any directory
* The extraction should create oofem_2.5_AMD64 directory
* Modify PATH variable to include oofem_2.5_AMD64/lib dir

    .. code-block:: bat

        set PATH=C:\Users\user\Documents\oofem_2.5_AMD64\lib;%PATH%

* Test run

    .. code-block:: bat

        C:\Users\user\Documents\oofem_2.5_AMD64\bin\oofem -v
        OOFEM version 2.5 (AMD64-Windows, fm;tm;sm;IML++) of Dec 30 2017 on JAJA
        Copyright (C) 1994-2017 Borek Patzak
        This is free software; see the source for copying conditions. There is NO
        warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


Binary installation – Unix/Linux
--------------------------------
You can either choose generic Linux binary or Debian package. 

Linux – binary
^^^^^^^^^^^^^^
* Extract downloaded zip archive (oofem_2.5_x86.tar.gz)
  .. code-block:: bash

    tr xzvf oofem_2.5_x86.tar.gz

* This creates oofem_2.5_x86_64/bin and oofem_2.5_x86_64/bin directories
* Define where to find libraries
    .. code-block:: bash

        $ export LDD_LIBRARY_PATH=../lib:$LDD_LIBRARY_PATH

* Run oofem
    .. code-block:: bash

        $ cd bin; ./oofem –v

Linux – Debian package
^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: bash

    $ sudo apt install oofem_2.5_x86_64.deb

Installation from source
------------------------
Requirements

* C++ compiler (g++, Visual Studio, MinGW, ...)
* CMake build tools (https://cmake.org/)
* Controls software compilation process by using platform independent configuration and generate native makefiles or projects
* To generate platform makefile or project configuration, use ``cmake`` (or any user friendly based GUI, such as ccmake or cmake-gui)

Instalation from source - Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the following, we assume the sources are in /home/user/oofem.src and build directory in /home/user/oofem.build

.. code-block:: bash

    $ cd /home/user/oofem.build
    $ ccmake /home/user/oofem.src
    $ make 

Installation from source - Windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Requirements:

* C++ compiler (Visual Studio, MinGW)
* CMake build tools (https://cmake.org/)

Procedure:

* Clone oofem git repository
* Use cmake to generate VS solution
    .. image:: figs/Installation_win_cmake.png
        :scale: 10 %
        :alt: screenshot of cmake gui
        :align: right

    * Select compiler
    * Set source directory
    * Set build directory
    * Generate project/solution
* Use compiler to build project targets (oofem and RUN_TESTS)
    .. image:: figs/Installation_win_VisualC2.png
        :scale: 10%
    .. image:: figs/Installation_win_VisualC.png
        :scale: 10%

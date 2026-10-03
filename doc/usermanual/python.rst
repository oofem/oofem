.. _python:

Using OOFEM from Python
=======================

Besides being a standalone solver, OOFEM can be used as a Python module. This
is useful when the analysis is only one step of a larger workflow — a parameter
study, an optimization loop, or postprocessing with numpy and matplotlib — and
when you want to implement a new component (a material model, a boundary
condition, a time function) in Python rather than in C++.

This section is only a short overview. The Python interface is documented in
full in the
`OOFEM Python bindings manual <http://www.oofem.org/resources/doc/python/html/index.html>`_.

Getting the module
------------------
Configure the build with ``USE_PYBIND_BINDINGS`` (see :ref:`installation`),
which produces an ``oofempy`` module in the build directory. Make the build directory visible to Python:

.. code-block:: bash

    $ export PYTHONPATH=/home/user/oofem.build:$PYTHONPATH

Running an existing input file
------------------------------
The shortest useful thing you can do is to run an ordinary input file and drive
the solution loop yourself:

.. code-block:: python

    import oofempy

    dr = oofempy.OOFEMTXTDataReader("beam2d_1.in")
    problem = oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, 0, None, False)
    problem.checkProblemConsistency()
    problem.init()
    problem.solveYourself()
    problem.terminateAnalysis()

This is exactly what the ``oofem`` executable does, and it produces the same
output file. Once the problem object exists, though, you can reach into it
between steps — query nodal unknowns, change boundary condition values, and
solve again — which is what makes the interface worth using.

Building a problem programmatically
-----------------------------------
The whole input file can also be bypassed: every component (domain, nodes,
elements, cross sections, materials, boundary conditions, functions, sets) has
a corresponding Python constructor, so a model can be assembled directly in
Python without any input file at all. See the "Building a problem" chapter of
the Python bindings manual for the pattern, and the
`bindings/python/tests <https://github.com/oofem/oofem/tree/devel/bindings/python/tests>`_
directory for complete working examples, which are part of the test suite and
therefore always up to date with the current interface.

Calling Python from the input file
----------------------------------
The traffic also goes the other way: an ordinary input file can refer to code
written in Python, and OOFEM will call it during the analysis. Two independent
build options are involved.

A build configured with ``USE_PYTHON_EXTENSION`` embeds a Python interpreter
and enables components driven by Python expressions:

``pythonexpression``
    A function whose value (and optionally its time derivatives) is given as a
    Python expression in the ``f``, ``dfdt`` and ``d2fdt2`` fields, or read
    from a file named by ``fpath``. Usable anywhere a function is, so any
    boundary condition can be given an arbitrary time variation without
    recompiling.

``userdefdirichletbc``, ``usrdefboundaryload``
    A Dirichlet boundary condition and a boundary load whose values are
    evaluated by user-supplied Python code, which lets them depend on
    position and on the state of the analysis.

A build configured with ``USE_PYBIND_BINDINGS`` additionally provides
``pythonmaterial``, a material model implemented in Python. The input record
names the Python module to import and the object within it to use (the
``module`` and ``object`` fields). That object must provide
``hasMaterialModeCapability``, ``giveCharacteristicMatrix``,
``giveCharacteristicVector`` and ``giveCharacteristicValue`` — and optionally
``printOutputAt`` — mirroring the C++ material interface; per-integration-point
history is kept for it in a Python dictionary that OOFEM commits or rolls back
along with the rest of the solution state. This is the fastest route to trying
out a new constitutive model, at the cost of the per-call overhead of the
Python interpreter.

Creating OOFEM model
####################

Creating model from input file
==============================
The problem can be instatiated form a stadard oofem input file, as illustrated form the following example

.. code-block:: pycon

    $ python3
    Python 2.7.10 (default, Aug 22 2015, 20:33:39)
    Python 3.6.8 (default, Oct  7 2019, 12:59:55) 
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import oofempy
    >>> dr=oofempy.OOFEMTXTDataReader("patch010.in")
    >>> problem=oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)


After that, the problem and its components can be manipulated from python3

.. code-block:: pycon

    >>> problem.init()
    >>> problem.solveYourself()
    Computing initial guess
    StaticStructural :: solveYourselfAt - Solving step 1, metastep 1, (neq = 3)
    NRSolver: Iteration ForceError
    ----------------------------------------------------------------------------
    NRSolver: 0      D_u:  0.000e+00
    Checking rules...
    EngngModel info: user time consumed by solution step 1: 0.00s
    >>> problem.terminateAnalysis()
    ANALYSIS FINISHED
    Real time consumed: 000h:00m:45s
    User time consumed: 000h:00m:00s

Building model from script
==========================
It is also possible to create individual problem components directly from Python, as illustrated in the following code

.. code-block:: python3

    import oofempy
    import util # some utility functions
    problem = oofempy.linearStatic(nSteps=1, outFile="test2.out") # engngModel
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._2dBeamMode, tstep_all=True, dofman_all=True, element_all=True) # domain aka mesh
    problem.setDomain(1, domain, True) # associate domain to the problem
    

One can create individual components directly from python. The convenience constructors accepting all keyword-value pairs of component records are provided. 

.. code-block:: python3

    # load time function
    ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
    ltfs = (ltf1, )
    # boundary conditions
    # loadTimeFunction parameter can be specified as int value or as LoadTimeFunction itself (valid for all objects with giveNumber() method)
    bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
    nLoad = oofempy.nodalLoad(            2, domain, loadTimeFunction=1,    components=(-18.,24.,0.))
    bcs = (bc1, nLoad)
    
Next we create nodes and other components. See, that reference to other components can be specified by the component number of by passing the object itself.

.. code-block:: python3

    # nodes
    # if one value is passed as parameter where oofem expects array of values, it must be passed as tuple or list (see load in n4)
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(1,1,1))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), bc=(0,0,0), load = (nLoad,))
    nodes = (n1, n2)
    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
    cs  = oofempy.simpleCS(1, domain, area=0.162, Iy=0.0039366, beamShearCoeff=1.e18, thick=0.54)
    # elements
    e1 = oofempy.beam2d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1)
    elems = (e1, )
    # add eveything to domain 
    util.setupDomain (domain, nodes, elems, (mat,), (cs,), bcs, ltfs, ())

After setting up the problem we can solve it.

.. code-block:: pycon

    >>> print("\nSolving problem")
    Solving problem
    >>> problem.checkProblemConsistency()
    >>> problem.init()
    >>> problem.postInitialize()
    >>> problem.setRenumberFlag()
    >>> problem.solveYourself()
    Solving ...
    EngngModel info: user time consumed by solution step 1: 0.00s
    >>> problem.terminateAnalysis()
    ANALYSIS FINISHED
    Real time consumed: 000h:00m:00s
    User time consumed: 000h:00m:00s

The more elaborated example can be found in test_2.py (https://github.com/oofem/oofem/blob/master/bindings/python/tests/test_2.py).
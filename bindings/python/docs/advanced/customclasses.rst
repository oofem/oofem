Extending OOFEM in Python
#########################

The python interface allows user to define custom derived classes in Python and inject them inside oofem.
We first illustate the concept on defining user defined time function (called MyFunc) in Python. The class should be derived from corresponding oofem base class (Function).
It is necessary to implement at least methods to evaluate the function value (evaluateAtTime) and its derivatives (evaluateVelocityAtTime and evaluateAccelerationAtTime).

.. code-block:: python3

    import oofempy
    class MyFunc(oofempy.Function):
        def __init__(self, num, domain):
            oofempy.Function.__init__(self, num, domain)
        def evaluateAtTime (self, t):
            return 2.5*t
        def evaluateVelocityAtTime(self, t):
            return 0.0;
        def evaluateAccelerationAtTime(self, t):
            return 0.0;

After the definition, we can directly use the MyFunc class in the problem setup:

.. code-block:: python3

    problem = oofempy.linearStatic(nSteps=1, outFile="test3.out")
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._1dTrussMode, tstep_all=True, dofman_all=True, element_all=True)
    problem.setDomain(1, domain, True)
    ltf1 = oofempy.peakFunction(1, domain, t=1, f_t=1)
    ltf2 = MyFunc(2, domain) # use custom ltf here
    ltfs = (ltf1, ltf2)
    # boundary conditions
    bc1   = oofempy.boundaryCondition(    1, domain, loadTimeFunction=1,    prescribedValue=0.0)
    n2    = oofempy.nodalLoad(            2, domain, loadTimeFunction=2,    components=(1.,), dofs=(1,))
    bcs = (bc1, n2)
    #nodes
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(bc1,))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), load=(n2,))
    nodes = (n1, n2)
    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., E=30.e6, n=0.2, tAlpha=1.2e-5)
    cs  = oofempy.simpleCS(1, domain, area=0.5, Iy=0.0, beamShearCoeff=1.e18, thick=0.5)
    # elements
    e1 = oofempy.truss1d(1, domain, nodes=(1,n2),  mat=1,   crossSect=1)
    elems = (e1,)
    # add eveything to domain (resize container first to save some time, but it is not necessary 0 see ltfs)
    domain.setup(nodes, elems, (mat,), (cs,), bcs, ltfs, (,))
     
    problem.checkProblemConsistency()
    problem.init()
    problem.postInitialize()
    problem.solveYourself()
    problem.terminateAnalysis()
    print("\nProblem solved")

Follow https://github.com/oofem/oofem/blob/master/bindings/python/tests/test_3.py for a full example.

Custom element
--------------
This section illustrates how to implement custom element in Python. In this case, we are implementing element for structural analysis, so it should be derived from cooresponding oofem base class, which is StructuralElement in this case.
There are various options how to implement an element. The simplest option is to directly define how to evaluate element stifness matrix and internal force vector. Additional functuos to determine degrees of freedom required by the element anre also to be provided.
Other possibilities exists, for example one can rely on default implementation of StructuralElement to evalaute stifness matrix, but then methods delivering strain-displacement and constituve matrix have to be provided together with element integration rules. For details, please consult oofem manuals.

.. code-block:: python3

    import oofempy
    class MyElement(oofempy.StructuralElement):
    def __init__(self, num, domain):
        oofempy.StructuralElement.__init__(self, num, domain)
        self.setNumberOfDofManagers(2)
    def computeStiffnessMatrix(self, answer, rMode, tStep):
        answer.resize(2,2);
        answer[0,0] = 1.0;
        answer[0,1] = -1.0;
        answer[1,0] = 1.0;
        answer[1,1] = -1.0;
    def computeBmatrixAt(self, gp, answer, lowerIndx, upperIndx):
        answer.resize(1,2);
        answer[0,0]=-1.0;
        answer[0,1]=1.0;
    def giveInternalForcesVector(self, answer, tStep, useUpdatedGpRecord):
        u = oofempy.FloatArray()
        k = oofempy.FloatMatrix(2,2)
        self.computeVectorOf(oofempy.ValueModeType.VM_Total, tStep, u)
        self.computeStiffnessMatrix(k, oofempy.CharType.StiffnessMatrix, tStep)
        answer.beProductOf (k, u)
    def giveNumberOfDofs(self):
        return 2;
    def computeNumberOfDofs(self):
        return 2;
    def giveDofManDofIDMask(self, inode, answer):
        print ("giveDofManDofIDMask for %d"%(inode,))
        print (answer)
        answer.resize(1)
        answer[0] = oofempy.DofIDItem.D_u
        #answer.pY()
        print(answer)
    def giveClassName(self):
        return "MyElement"
    def giveInputRecordName(self):
        return "MyElement"

The element can again added into domain and used from oofem

.. code-block:: python3

    # nodes
    n1 = oofempy.node(1, domain, coords=(0.,  0., 0. ), bc=(bc1,))
    n2 = oofempy.node(2, domain, coords=(2.4, 0., 0. ), load=(n2,))
    nodes = (n1, n2)
    # elements
    e1 = MyElement(1, domain) # additinal entries should go to to the custom element constructor
    e1.setDofManagers((1,2))
    ir = oofempy.OOFEMTXTInputRecord()
    ir.setRecordString ("nodes 2 1 2 mat 1 crosssect 1")
    # pass input record to elem
    e1.initializeFrom(ir)
    elems = (e1,)
    ...

You can follow https://github.com/oofem/oofem/blob/master/bindings/python/tests/test_3.py for a complete illustration.

Custom material
---------------
This section illustrates how to implement custom constitutive model in Python.
First we have to define Python class implementing the model derived from corresponding oofem class. In this case, we are going to implement constututive model for structural analysis, so we derive our class from StructuralMaterial.
The presented implementation is mininimalistic one, we overload or define just methods to support 1d stress strain state, by overriding give1dStiffnessMatrix and giveRealStressVector_1d methods, that simply returns constitutive matrix and evaluate 1d stress from given strain. The more elaborate implementation would be necessary to support sevaral stress-strain modes, nonlinear materials (need to create custom status to track history variables). For more details, please refer to oofem programmer's manual.
Each matteril model must define its material status contatining its internal state variables. Even if material models does not need to track internal variables, it should provide status. In such case, it is sufficient to create instance of StructuralMaterialStatus. In oofem, the method responsible for status creation is CreateStatus, that would normally be overriden in Python as well. However, doe to the issue in Pybind11 (https://github.com/pybind/pybind11/issues/1962) this is not yet possible. The workaroud is to override giveStatus, as illustrate in the following example.

.. code-block:: python3

  class MyMaterial(oofempy.StructuralMaterial):
    def __init__(self, num, domain):
        oofempy.StructuralMaterial.__init__(self, num, domain)
        self.k = 1.5;
    def giveClassName(self):
        return "MyMaterial" 
    def giveInputRecordName(self):
        return "MyElement"
    # Overloading this method is not yet possible in pybind11
    # see https://github.com/pybind/pybind11/issues/1962
    # However a workaround is to override giveStatus, see below 
    def CreateStatus(self, gp):
        return oofempy.StructuralMaterialStatus (gp)
    def give1dStressStiffMtrx(self, answer, mode, gp, tStep):
        answer.resize(1,1)
        answer[0,0] = self.k
        return
    def giveRealStressVector_1d (self, answer, gp, reducedStrain, tStep):
        answer.resize(1)
        answer[0] = self.k * reducedStrain[0]
        status = self.giveStatus(gp)
        status.letTempStrainVectorBe(reducedStrain);
        status.letTempStressVectorBe(answer);
        return
    def giveStatus (self, gp):
        print ("getStatus")
        if (gp.giveMaterialStatus() is None):
            print ("getStatus creating")
            status = oofempy.StructuralMaterialStatus (gp)
            gp.setMaterialStatus(status)
        return gp.giveMaterialStatus() 
        

You can follow https://github.com/oofem/oofem/blob/master/bindings/python/tests/test_4.py for a complete illustration.


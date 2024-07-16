#
# Demo code to illustrate new mpm symbolic concept
# Requires oofempy compiled with __MPM_MODULE ON
#
#

import sys
sys.path.extend(['/home/bp/devel/oofem.git/build', '/home/bp/devel/oofem.git/bindings/python'])
import oofempy
import util
import numpy as np
import pyvista as pv
 


def test_8():
    # Requires oofempy compiled with __MPM_MODULE ON
    import sys
    sys.path.extend(['/home/bp/devel/oofem.git/build', '/home/bp/devel/oofem.git/bindings/python'])
    import oofempy
    import util

    

    # Create a new dummy problem (placeholder for our demo) with one domain.
    problem = oofempy.dummyProblem(nSteps=1, outFile='test_7.out')
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._HeatTransferMode, tstep_all=1, dofman_all=0, element_all=0)
    problem.setDomain(1, domain, True)
    
    # Define nodes
    n1 = oofempy.node(1, domain, coords=(0, 0, 0. ))
    n2 = oofempy.node(2, domain, coords=(1., 0.0, 0. ))
    n3 = oofempy.node(3, domain, coords=(1., 1., 0. ))
    n4 = oofempy.node(4, domain, coords=(0, 1, 0. ))
    
    # Defdine elements, note that q1 defines just element geometry.
    q1 = oofempy.q1(1, domain, nodes=(1,2,3,4), mat=1, crossSect=1) # quad element #1
    l1 = oofempy.l1(2, domain, nodes=(2,3), mat=1, crossSect=1)     # boundary element #2

    # Dirichlet Boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, dofs=(1,2), values=(0.,0.), set=1)
    bc2 = oofempy.boundaryCondition(2, domain, loadTimeFunction=1, dofs=(1,),  values=(0.,),   set=2)
    # material and cross section
    mat = oofempy.isoLE(1, domain, d=1., e=1., n=0.3, talpha=1.)
    cs  = oofempy.simpleCS(1, domain, mat=1, thickness=1.0)
    # time functions
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.0)
    # some sets (groups of nodes and elements) for later use
    s1 = oofempy.createSet(1, domain, nodes=(1,))
    s2 = oofempy.createSet(2, domain, nodes=(4,))
    s3 = oofempy.createSet(3, domain, elements=(1,))
    bs1 = oofempy.createSet(4, domain, elements=(2,))
    util.setupDomain(domain, nodes=(n1,n2,n3,n4), elems=(q1,l1), css=(cs,), mats=(mat,), bcs=(bc1,bc2), ics=(), ltfs=(ltf1,), sets=(s1,s2,s3,bs1))


    interpolation = oofempy.linearinterpolation()
    w= u = oofempy.Variable(interpolation, oofempy.VariableQuantity.Displacement, oofempy.VariableType.vector, 2, [1,2], None)
    mt = oofempy.BTSigmaTerm(w, u, oofempy.MaterialMode._PlaneStress)
    lt = oofempy.NTfTerm(w, u, oofempy.MaterialMode._PlaneStress)
    tstep = problem.giveNextStep()


    I1 = oofempy.Integral(domain, s3, mt)
    problem.addIntegral(I1)
    I1.initialize()

    I2 = oofempy.Integral(domain, bs1, lt)
    I2.initialize()

    problem.postInitialize()
    problem.forceEquationNumbering()


    lhs=oofempy.skyline()
    lhs.buildInternalStructure(problem, 1, oofempy.EModelDefaultEquationNumbering());
    I1.assemble_lhs(lhs, oofempy.EModelDefaultEquationNumbering(), tstep)
    lhs.printYourself()


    rhs = oofempy.FloatArray(5)

    I2.assemble_rhs(rhs, oofempy.EModelDefaultEquationNumbering(), tstep)
    rhs.printYourself()

    r = oofempy.FloatArray(5)
    linsolv = oofempy.ldltfactorization(domain, problem)
    linsolv.solve(lhs, rhs, r)
    print ("Displacement vector = ", r)

    assert (round (r[0]-1.0, 4) == 0), "Displacement at node 2, check failed"
    assert (round (r[1]-0.0, 4) == 0), "Displacement at node 2, check failed"
    assert (round (r[2]-1.0, 4) == 0), "Displacement at node 3, check failed"
    assert (round (r[3]+0.3, 4) == 0), "Displacement at node 3, check failed"
    assert (round (r[4]+0.3, 4) == 0), "Displacement at node 4, check failed"


if __name__ == "__main__":
    if oofempy.hasModule('mpm'):
        test_8()
    else:
        print("This example requires mpm module.")
        print("Please recompile oofem with mpm module enabled.")

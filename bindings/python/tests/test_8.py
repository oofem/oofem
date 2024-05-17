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
    # Create a new domain
    problem = oofempy.dummyProblem(nSteps=1, outFile='test_7.out')
    domain = oofempy.domain(1, 1, problem, oofempy.domainType._HeatTransferMode, tstep_all=1, dofman_all=0, element_all=0)
    problem.setDomain(1, domain, True)
    
    #nodes
    n1 = oofempy.node(1, domain, coords=(0, 0, 0. ))
    n2 = oofempy.node(2, domain, coords=(1., 0.0, 0. ))
    n3 = oofempy.node(3, domain, coords=(1., 1., 0. ))
    n4 = oofempy.node(4, domain, coords=(0, 1, 0. ))
    
    # elements
    q1 = oofempy.q1(1, domain, nodes=(1,2,3,4), mat=1, crossSect=1)
    # boundary conditions
    bc1 = oofempy.boundaryCondition(1, domain, loadTimeFunction=1, dofs=(1,2), values=(0.,0.), set=1)
    bc2 = oofempy.boundaryCondition(2, domain, loadTimeFunction=1, dofs=(1,),  values=(0.,),   set=2)
    # material and cross section
    #mat = oofempy.upm(1, domain, d=1., e=1.e7, nu=0.3, k=1.e-15, c=3.33333333e-10, alpha=1.)
    mat = oofempy.isoLE(1, domain, d=1., e=1., n=0.3, talpha=1.)
    cs  = oofempy.simpleCS(1, domain, mat=1, thickness=1.0)
    ltf1 = oofempy.constantFunction(1, domain, f_t=1.0)
    #Initial condition
    ic1 = oofempy.initialCondition(1, domain, Conditions=1, u=50.0, dofs=(10,), set=3)
    s1 = oofempy.createSet(1, domain, nodes=(1,))
    s2 = oofempy.createSet(2, domain, nodes=(4,))
    s3 = oofempy.createSet(3, domain, elements=(1,))
    util.setupDomain(domain, nodes=(n1,n2,n3,n4), elems=(q1,), css=(cs,), mats=(mat,), bcs=(bc1,bc2), ics=(ic1,), ltfs=(ltf1,), sets=(s1,s2,s3))


    interpolation = oofempy.fei2dquadlin(1,2)
    u = oofempy.Variable(interpolation, oofempy.VariableQuantity.Displacement, oofempy.VariableType.vector, 2, [1,2], None)
    mt = oofempy.BTSigmaTerm(u, u, oofempy.MaterialMode._PlaneStress)
    tstep = problem.giveNextStep()

    integral = oofempy.Integral(domain, s3, mt)
    problem.addIntegral(integral)
    integral.initialize()
    problem.postInitialize()

    problem.forceEquationNumbering()
    lhs=oofempy.skyline()
    lhs.buildInternalStructure(problem, 1, oofempy.EModelDefaultEquationNumbering());
    integral.assemble_lhs(lhs, oofempy.EModelDefaultEquationNumbering(), tstep)
    lhs.printYourself()

    rhs = oofempy.FloatArray(5)
    rhs[0]=1.0
    rhs[2]=1.0
    r = oofempy.FloatArray(5)
    linsolv = oofempy.ldltfactorization(domain, problem)
    linsolv.solve(lhs, rhs, r)
    print (r)

if __name__ == "__main__":
    if oofempy.hasModule('mpm'):
        test_8()
    else:
        print("This example requires mpm module.")
        print("Please recompile oofem with mpm module enabled.")

import oofempy
import util
import matplotlib

dr=oofempy.OOFEMTXTDataReader("concrete_3point.in")
problem=oofempy.InstanciateProblem(dr, oofempy.problemMode.processor, False, None, False)
problem.init()
problem.solveYourself()

d=problem.giveDomain(1)
tstep = problem.giveCurrentStep()
    
field = problem.giveField(oofempy.InternalStateType.IST_StressTensor, tstep)
util.plot2D(d, field, 0, tstep).show()

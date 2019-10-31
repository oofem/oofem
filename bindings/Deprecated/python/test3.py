#from __future__ import print_function
import liboofem

dr=liboofem.OOFEMTXTDataReader("inputs/tmpatch42.in")
problem=liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
problem.checkProblemConsistency()
problem.setRenumberFlag()
problem.solveYourself()
problem.terminateAnalysis()
timeStep=problem.giveCurrentStep()

t=liboofem.FieldType.FT_Temperature
f=problem.giveField(liboofem.FieldType(t),timeStep)
print("Time:", timeStep.targetTime, "name:", t,"reference:",str(f))

print('Number of domains:',problem.giveNumberOfDomains())
domain=problem.giveDomain(1)
print('Domain', domain)
elem=domain.giveElement(3)
print('Element type',elem.giveGeometryType())

f=problem.giveField(liboofem.FieldType(t),timeStep)
print('Field at time',timeStep.targetTime,'type',t, f)

for d in range(1,domain.giveNumberOfDofManagers()+1): #Nodes
    dof=domain.giveDofManager(d)
    cc=dof.giveCoordinates()
    valueModeType=liboofem.ValueModeType.VM_Total
    val=f.evaluateAtDman(dof,mode=valueModeType,atTime=timeStep)
    print('DofManager #%d at' % dof.giveNumber(), [cc[i] for i in range(len(cc))], "value", val)

for el in range(1,domain.giveNumberOfElements()+1): #Elements
    elem=domain.giveElement(el)
    gt=elem.giveGeometryType()
    print('Element %d Nodes' % elem.giveNumber(), [elem.giveDofManager(i+1).giveNumber() for i in range(elem.numberOfDofManagers)])

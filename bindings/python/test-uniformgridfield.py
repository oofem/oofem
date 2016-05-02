from __future__ import print_function
import liboofem
import numpy as np
# instantiate uniform field
ugrid=liboofem.UniformGridField()
# set 2d geometry
ugrid.setGeometry(lo=(0,0),hi=(1,1),div=(2,2))
# set data: div is 2x2, hence 3x3=9 node grid with 9 values
ugrid.setValues([0,.5,0, .5,1,.5, 0,.5,.5])

dr = liboofem.OOFEMTXTDataReader("tmpatch42.in")
transpModel = liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
transpModel.checkProblemConsistency()
transpModel.solveYourself()

fieldMan = transpModel.giveContext().giveFieldManager()
timeStep = transpModel.giveCurrentStep()
eField = transpModel.giveField(liboofem.FieldType.FT_Temperature,timeStep)

fieldMan.registerField(eField,liboofem.FieldType.FT_Velocity)
fieldMan.registerField(ugrid,liboofem.FieldType.FT_Unknown)

print(fieldMan.giveRegisteredKeys())
print(fieldMan.giveField(liboofem.FieldType.FT_Velocity))
print(fieldMan.giveField(liboofem.FieldType.FT_Unknown))

transpModel.terminateAnalysis()
timeStep = transpModel.giveCurrentStep()

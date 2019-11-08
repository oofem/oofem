from __future__ import print_function
import liboofem
import numpy as np

# Instantiate uniform field
ugrid=liboofem.UniformGridField()
# Set 2D geometry
ugrid.setGeometry(lo=(0,0),hi=(1,1),div=(2,2))
# Set data: div is 2x2, hence 3x3=9 node grid with 9 values
ugrid.setValues([0,.5,0, .5,1,.5, 0,.5,.5])
#ugrid.setValues([15,15,15,15,15,15,15,15,15])
print (ugrid.evaluateAtPos((0.5,0.2)))

dr = liboofem.OOFEMTXTDataReader("inputs/tmpatch05.in")
transpModel = liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
transpModel.checkProblemConsistency()

fieldMan = transpModel.giveContext().giveFieldManager()
timeStep = transpModel.giveCurrentStep()

#eField = transpModel.giveField(liboofem.FieldType.FT_Temperature,timeStep)
#fieldMan.registerField(eField,liboofem.FieldType.FT_TransportProblemUnknowns)

fieldMan.registerField(ugrid,liboofem.FieldType.FT_TemperatureAmbient)

print(fieldMan.giveRegisteredKeys())
print(fieldMan.giveField(liboofem.FieldType.FT_TemperatureAmbient))

transpModel.solveYourself()
transpModel.terminateAnalysis()


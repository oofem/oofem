import liboofem
dr=liboofem.OOFEMTXTDataReader("tmpatch42.in")
pb=liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0)
pb.checkProblemConsistency()
pb.setRenumberFlag()
pb.solveYourself()
pb.terminateAnalysis()
ts=pb.giveCurrentStep()

for t in liboofem.FieldType.FT_Temperature,liboofem.FieldType.FT_HumidityConcentration:
    f=pb.giveField(liboofem.FieldType(t),ts)
    print ts.targetTime,t,str(f)
#fm=pb.giveContext().giveFieldManager()
#print fm
#print fm.giveRegisteredKeys()
# why is this one None?
#print fm.giveField(liboofem.FieldType.FT_Displacements)
#for i in range(100):
#    print fm.giveField(liboofem.FieldType(i))

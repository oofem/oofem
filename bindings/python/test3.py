from __future__ import print_function
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
    print(ts.targetTime,t,str(f))

#fm=pb.giveContext().giveFieldManager()
#print fm
#print fm.giveRegisteredKeys()
# why is this one None?
#print fm.giveField(liboofem.FieldType.FT_Displacements)
#for i in range(100):
#    print fm.giveField(liboofem.FieldType(i))

print('Number of domains:',pb.giveNumberOfDomains())
# everything starts at 1
do=pb.giveDomain(1)
print('Domain',do)
e=do.giveElement(3)
print('Element type',e.giveGeometryType())
for t in liboofem.FieldType.FT_Temperature,: #,liboofem.FieldType.FT_HumidityConcentration:
    f=pb.giveField(liboofem.FieldType(t),ts)
    print('Field at time',ts.targetTime,'type',t,str(f))
ne=do.giveNumberOfElements()
nd=do.giveNumberOfDofManagers()

# vertices
for d in range(1,nd+1):
    dof=do.giveDofManager(d)
    # print dof
    cc=dof.giveCoordinates()
    print('DofManager #%d at '%(d-1),[cc[i] for i in range(len(cc))])
    # for vmt in liboofem.ValueModeType.VM_Total,liboofem.ValueModeType.VM_Unknown,liboofem.ValueModeType.VM_Velocity,liboofem.ValueModeType.VM_Acceleration,liboofem.ValueModeType.VM_Incremental:
    vmt=liboofem.ValueModeType.VM_Total
    val=f.evaluateAtDman(dof,mode=vmt,atTime=ts)
    if len(val): print('   value',vmt,len(val),val)
    else: print ('   (zero-length array read from this DoF manager)')
# elements
for ei in range(1,ne+1):
    e=do.giveElement(ei)
    gt=e.giveGeometryType()
    nn=e.numberOfDofManagers
    for n in range(1,nn+1):
        dmn=e.giveDofManager(n)
        print('Element %d, DOF #%d -> node %d'%(ei-1,n-1,dmn.giveNumber()))

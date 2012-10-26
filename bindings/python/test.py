import oofemlib

a = oofemlib.FloatArray(2);
b = oofemlib.FloatArray(4);
ind = oofemlib.IntArray(2);

a.zero();
b.zero();
a[1]=15.0;
a[0]=10.0;
print a[1], a[2];
a.printYourself();

c = oofemlib.FloatArray(a)
print c[1], c[2];
c.printYourself();

d = oofemlib.FloatMatrix(2,2);
d.beUnitMatrix();
d[0,1]=1.0;
d.printYourself();

b.beProductOf(d,a);
b.printYourself();

b.resize(4);
b.zero();
ind[0]=1;
ind[1]=3;
b.assemble(a,ind);
b.printYourself();


dr=oofemlib.OOFEMTXTDataReader("patch100.in");
problem=oofemlib.InstanciateProblem(dr,oofemlib.problemMode._processor,0);
problem.checkProblemConsistency();
problem.setRenumberFlag();
problem.solveYourself();
problem.terminateAnalysis();

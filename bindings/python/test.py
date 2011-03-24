import oofemlib

a = oofemlib.FloatArray(2);
b = oofemlib.FloatArray(4);
ind = oofemlib.IntArray(2);

a.zero();
b.zero();
a[2]=15.0;
a[1]=10.0;
print a[1], a[2];
a.printYourself();

c = oofemlib.FloatMatrix(2,2);
c.beUnitMatrix();
c[1,2]=1.0;
c.printYourself();

b.beProductOf(c,a);
b.printYourself();

b.resize(4);
b.zero();
ind[1]=1;
ind[2]=3;
b.assemble(a,ind);
b.printYourself();


dr=oofemlib.OOFEMTXTDataReader("/home/bp/work/patch100.in");
problem=oofemlib.InstanciateProblem(dr,oofemlib.problemMode._processor,0);
problem.checkProblemConsistency();
problem.setRenumberFlag();
problem.solveYourself();
problem.terminateAnalysis();

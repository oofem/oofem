import liboofem

a = liboofem.FloatArray(2);
b = liboofem.FloatArray(4);
ind = liboofem.IntArray(2);

a.zero();
b.zero();
a[1]=15.0;
a[0]=10.0;
print a[0], a[1];
a.printYourself();

c = liboofem.FloatArray(a)
print c[0], c[1];
c.printYourself();

d = liboofem.FloatMatrix(2,2);
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


dr=liboofem.OOFEMTXTDataReader("patch100.in");
problem=liboofem.InstanciateProblem(dr,liboofem.problemMode._processor,0);
problem.checkProblemConsistency();
problem.setRenumberFlag();
problem.solveYourself();
problem.terminateAnalysis();

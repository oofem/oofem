import oofempy

a = oofempy.FloatArray((1.0, 2.0, 3.0))
b = oofempy.FloatArray((0.0, -1.0, 1.0))
x = oofempy.FloatArray(3)

c = a+b
c += a
d = a-b
d -= a
print(a)
print(d)

A = oofempy.FloatMatrix(3,3)

A.beUnitMatrix()

A.solveForRhs(a, x, False);
print(x)
A.giveDeterminant()
Ainv = oofempy.FloatMatrix(3,3);
Ainv.beInverseOf(A)
y = Ainv*x
print(y)
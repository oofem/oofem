import oofempy


def test_1():
    a = oofempy.FloatArray((1.0, 2.0, 3.0))
    b = oofempy.FloatArray((0.0, -1.0, 1.0))
    x = oofempy.FloatArray(3)

    c = a+b
    assert (round(c[0] - 1.0, 6) == 0), "Error in a+b"
    assert (round(c[1] - 1.0, 6) == 0)
    assert (round(c[2] - 4.0, 6) == 0)

    c += a
    assert (round(c[0]-2.0, 6) == 0)
    assert (round(c[1]-3.0, 6) == 0)
    assert (round(c[2]-7.0, 6) == 0)
    
    d = a-b
    assert (round(d[0]-1.0, 6) == 0)
    assert (round(d[1]-3.0, 6) == 0)
    assert (round(d[2]-2.0, 6) == 0)

    d -= a
    assert (round(d[0]-0.0, 6) == 0)
    assert (round(d[1]-1.0, 6) == 0)
    assert (round(d[2]+1.0, 6) == 0)
   
    print(a)
    print(d)

    A = oofempy.FloatMatrix(3,3)
    A.beUnitMatrix()
    
    det = A.giveDeterminant()
    assert (round (det - 1.0, 6) == 0), "Error in geveDeterminant"

    A.solveForRhs(a, x, False);
    assert (round(x[0]-1.0, 6) == 0)
    assert (round(x[1]-2.0, 6) == 0)
    assert (round(x[2]-3.0, 6) == 0)
     
    print(x)
 
    Ainv = oofempy.FloatMatrix(3,3);
    Ainv.beInverseOf(A)
    y = Ainv*x

    assert (round(y[0]-1.0, 6) == 0)
    assert (round(y[1]-2.0, 6) == 0)
    assert (round(y[2]-3.0, 6) == 0)
    print(y)

if __name__ == "__main__":
    test_1()
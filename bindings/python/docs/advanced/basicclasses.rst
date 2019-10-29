Using OOFEM Classes
###################

The binding code is generated for all fundamental classes of OOFEM. 
This section presents the simple use cases, where selected OOFEM classes are used directly in python code.

OOFEm comes with built in representation of vectors and matrices. The full feature Python interface is provided in terms of all methods, but also the convenience constructors (from any sequence) and usual operators are provided.

.. code-block:: pycon

    $ python3
    Python 2.7.10 (default, Aug 22 2015, 20:33:39)
    Python 3.6.8 (default, Oct  7 2019, 12:59:55) 
    [GCC 8.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import oofempy
    >>> a = oofempy.FloatArray((1.0, 2.0, 3.0))
    >>> b = oofempy.FloatArray((0.0, -1.0, 1.0))
    >>> c = a+b
    >>> print (c)
    <oofempy.FloatArray: {1.000000, 1.000000, 4.000000, }>
    >>> c += a

The following example illustrates how to use FloatMatrix class to solve simple linear system

.. code-block:: pycon

    $ python3
    >>> import oofempy
    >>> A = oofempy.FloatMatrix(3,3)
    >>> A.beUnitMtrx()
    >>> A[0,1]=-2.
    >>> A[1,2]=3
    >>> A.printYourself()
    FloatMatrix with dimensions : 3 3
    1.000e+00  -2.000e+00   0.000e+00  
    0.000e+00   1.000e+00   3.000e+00  
    0.000e+00   0.000e+00   1.000e+00
    >>> b=oofempy.FloatArray((-3, 11, 3))
    >>> x = oofempy.FloatArray(3)
    >>> A.solveForRhs(b, x, False)
    >>> print (x)
    <oofempy.FloatArray: {1.000000, 2.000000, 3.000000, }>

 

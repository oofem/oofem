set (dss
    BigMatrix.cpp
    BiSection.cpp
    DenseMatrixArithmeticsNN.cpp
    DSSolver.cpp
    IntArrayList.cpp
    MathTracer.cpp
    SkyLineMtx.cpp
    SkyLineMtxLDL.cpp
    SparseConectivityMtx.cpp
    SparseGridMtx.cpp
    SparseGridMtxLDL.cpp
    SparseGridMtxLL.cpp
    SparseGridMtxLU.cpp
    SparseGridMtxPD.cpp
    SparseMatrixF.cpp
    # OOFEM binding:
    dssmatrix.C
    dsssolver.C
    )

add_library (dss
    ${LIB_TYPE}
    ${dss}
    )

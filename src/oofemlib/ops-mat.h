#pragma once

#ifndef __MATH_OPS
    #error __MATH_OPS must be defined to use this header (only oofem internal headers should include it)
#endif

#include<string>
#include<initializer_list>

#ifndef _MATHOPS_STAGE
    #define _MATHOPS_STAGE 100
#endif

namespace oofem{
    #ifndef _USE_EIGEN
        class FloatMatrix;
        class FloatArray;
    #endif
    class IntArray;
}

namespace oofem::mat{

    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param loc Localization indices.
     */
    void assemble(FloatMatrix& dst, const FloatMatrix &src, const IntArray &loc);
    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param rowind Row localization indices.
     * @param colind Column localization indices.
     */
    void assemble(FloatMatrix& dst, const FloatMatrix &src, const IntArray &rowind, const IntArray &colind);
    /**
     * Assembles the transposed contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param rowind Row localization indices.
     * @param colind Column localization indices.
     */
    void assembleT(FloatMatrix& dst, const FloatMatrix &src, const IntArray &rowind, const IntArray &colind);
    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param rowind Row localization indices.
     * @param colind Column localization indices.
     */
    void assemble(FloatMatrix& dst, const FloatMatrix &src, const int *rowind, const int *colind);

    /**
     * Computes the Frobenius norm of the receiver.
     * The Frobenius norm is defined as the square root of the sum of the absolute squares of its elements.
     * @return Frobenius norm.
     */
    double computeFrobeniusNorm(const FloatMatrix&);
    /**
     * Computes the operator norm of the receiver.
     * @param p Norm type, '1' for 1 norm, '2' for 2 norm.
     * @return Norm of receiver.
     */
    double computeNorm(const FloatMatrix&, char p);
    /**
     * Computes the conditioning of the receiver. From 0 to 1, where 0 is singular and 1 best.
     * The receiver must be square.
     * Works identically as MATLAB/Octaves rcond().
     * @param p Norm type, '1' for 1 norm, '2' for 2 norm.
     * @return Conditioning of receiver.
     */
    double computeReciprocalCondition(const FloatMatrix&, char p = '1');
    /*
     * Computes the eigenvalues of a symmetric matrix.
     * The receiver must be square and symmetric.
     * @param lambda Eigenvalues.
     * @param v Eigenvectors (stored column wise).
     * @param neigs If set, only neigs largest eigenvalues are computed.
     * @return True if successful.
     */
    //bool computeEigenValuesSymmetric(FloatArray &lambda, FloatMatrix &v, int neigs = 0) const;
    /**
     * Modifies receiver to be a diagonal matrix with the components specified in diag.
     * @return Determinant of receiver.
     */

    /**
     * Makes receiver the local coordinate for the given normal.
     * Implemented for 2D and 3D.
     * @param normal Normal (normalized).
     */
    FloatMatrix beLocalCoordSys(const FloatArray &normal);


    /**
     * Modifies receiver to become inverse of given parameter. Size of receiver will be adjusted.
     * @param src Matrix to be inverted.
     * @return False if K is singular, otherwise true.
     */
    bool beInverseOf(FloatMatrix& inv, const FloatMatrix &src);
    FloatMatrix beInverseOf(const FloatMatrix& src);
    /**
     * Solves the  system of linear equations @f$ K\cdot a = b @f$ . Uses Gaussian elimination with pivoting directly on receiver.
     * @param b RHS of linear system.
     * @param answer Solution of linear equations.
     * @param transpose Solves for the transpose of K.
     * @return False if K is singular, otherwise true.
     */
    bool solveForRhs(FloatMatrix& m, const FloatArray &b, FloatArray &answer, bool transpose = false);
    /**
     * Solves the  system of linear equations @f$ K\cdot A = B @f$ . Uses Gaussian elimination with pivoting directly on receiver.
     * @param B RHS of linear system.
     * @param answer Solution of linear equations, each column corresponding to columns in B.
     * @param transpose Solves for the transpose of K.
     */
    bool solveForRhs(FloatMatrix& m, const FloatMatrix &B, FloatMatrix &answer, bool transpose = false);

    /**
     * Adds to the receiver the product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Assumes that receiver and product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$ are symmetric matrices. Computes only the
     * upper half of receiver.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    void plusProductSymmUpper(FloatMatrix& recv, const FloatMatrix &a, const FloatMatrix &b, double dV);
    /**
     * Adds to the receiver the dyadic product @f$ a \otimes a \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Computes only the upper half of receiver.
     * @param a Array a in equation.
     * @param dV Scaling factor.
     */
    void plusDyadSymmUpper(FloatMatrix& recv, const FloatArray &a, double dV);
    /**
     * Adds to the receiver the product @f$a^{\mathrm{T}} \cdot b \mathrm{d}V@f$. If the receiver has zero size, it is expanded.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    void plusProductUnsym(FloatMatrix& recv, const FloatMatrix &a, const FloatMatrix &b, double dV);
    /**
     * Adds to the receiver the product @f$a \otimes b \mathrm{d}V@f$. If the receiver has zero size, it is expanded.
     * @param a Array a in equation.
     * @param b Array b in equation.
     * @param dV Scaling factor.
     */
    void plusDyadUnsym(FloatMatrix& recv, const FloatArray &a, const FloatArray &b, double dV);



    /**
     * Reciever will be a 3x3 matrix formed from a vector with either 9 or 6 components.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12, 32, 31, 21
     * If size(aArray) = 6, a symmetric matrix will be created.
     * @param aArray Array to transform.
     */
    FloatMatrix beMatrixFormOfStress(const FloatArray &aArray);
    FloatMatrix beMatrixForm(const FloatArray &aArray);

    /**
     * Swaps the indices in the 6x6 matrix such that it converts between OOFEM's
     * and Abaqus' way of writing matrices. Currently used to convert the 6x6 Jacobian
     * from Abaqus UMAT to OOFEM.
     */
    void changeComponentOrder(FloatMatrix&);


    /**
     * Computes eigenvalues and eigenvectors of receiver (must be symmetric)
     * The receiver is preserved.
     * @param eval Requested eigenvalues.
     * @param v Requested eigenvectors (stored colum wise).
     * @param nf Number of significant figures.
     * @return True if ok,otherwise false.
     */
    bool jaco_(FloatMatrix& M, FloatArray &eval, FloatMatrix &v, int nf);

    /// Prints matrix to stdout. Useful for debugging.
    void printYourself(const FloatMatrix& m);
    /**
     * Print receiver on stdout with custom name.
     * @param name Display name of reciever.
     */
    void printYourself(const FloatMatrix& m, const std::string &name);
    
    /**
     * Print matrix to file.
     * @param filename Output filename
     * @parap showDimensions Determines if dimensions should be included in output. Default is true.
     */
    void printYourselfToFile(const FloatMatrix& m, const std::string filename, const bool showDimensions=true);
    
    /// Higher accuracy than printYourself.
    void pY(const FloatMatrix& m);

    /**
     * Writes receiver as CSV (comma seperated values)
     * @param name Filename
     */
    void writeCSV(const FloatMatrix& m, const std :: string &name);

    /**
     * Returns the receiver 'a' transformed using give transformation matrix r.
     * The method performs the operation  @f$ a = r^{\mathrm{T}} \cdot a \cdot r@f$ .
     * @param r Transformation matrix.
     * @param mode If set to 't' then the transpose of the rotation matrix is used instead.
     */
    void rotatedWith(FloatMatrix& a, const FloatMatrix &r, char mode = 'n');
} /* namespace oofem::mat */

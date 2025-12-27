/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef numerics_h
#define numerics_h

#ifdef _USE_EIGEN
    #define EIGEN_INITIALIZE_MATRICES_BY_ZERO
    #include<Eigen/Core>
    #include <Eigen/Dense>
#endif

namespace oofem{
    #ifdef _USE_EIGEN
        #define OOFEM_EIGEN_DERIVED(MyKlass,EigenBase) \
            MyKlass(void): EigenBase() {} \
            /* This constructor allows you to construct MyKlass from Eigen expressions */ \
            template<typename OtherDerived> MyKlass(const Eigen::MatrixBase<OtherDerived>& other): EigenBase(other) { } \
            /* This method allows you to assign Eigen expressions to MyKlass */ \
            template<typename OtherDerived> MyKlass& operator=(const Eigen::MatrixBase<OtherDerived>& other) { this->EigenBase::operator=(other); return *this; }

        typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> MatrixXXd;
        typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;

        /* typedef for fixed-size base classes */
        // columns vectors are required to be col-major, even though it makes not difference for storage
        template<std::size_t N>
        using VectorNd = Eigen::Matrix<double,N,1,Eigen::ColMajor>; 

        // FloatMatrixF with N==1 (one row) must be RowMajor
        template<std::size_t R, std::size_t C>
        using MatrixRCd = Eigen::Matrix<double,R,C,(R==1?Eigen::RowMajor:Eigen::ColMajor)>;

        class FloatArray;
        class FloatMatrix;
        template<std::size_t N> class FloatArrayF;
        template<std::size_t R, std::size_t C> class FloatMatrixF;
        const double NaN(std::numeric_limits<double>::signaling_NaN());
        // using Eigen::Index;
        typedef int Index;
    #else
        typedef int Index;
    #endif
} // end namespace oofem



#ifdef __LAPACK_MODULE
// Some forward declarations for LAPACK. Remember to append the underscore to the function name.
// consumed in floatmatrix.C and floatarray.C
extern "C" {
    extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda, const double *x,
                       const int *incx, const double *beta, double *y, const int *incy, int aColumns, int xSize, int ySize);
    // Y = Y + alpha * X
    extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy, int xsize, int ysize);

    /// Computes the reciprocal condition number for a LU decomposed function.
    extern void dgecon_(const char *norm, const int *n, const double *a, const int *lda,
                        const double *anorm, double *rcond, double *work, int *iwork, int *info, int norm_len);
    /// Replaces a with the LU-decomposition.
    extern int dgetrf_(const int *m, const int *n, double *a, const int *lda, int *lpiv, int *info);
    /// Replaces a with its inverse.
    extern int dgetri_(const int *n, double *a, const int *lda, int *ipiv, double *work, const int *lwork, int *info);
    /// Solves a system of equations.
    extern int dgesv_(const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, const double *b, const int *ldb, int *info);
    /// Computes the norm.
    extern double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work, int norm_len);
    /// Computes eigenvalues and vectors.
    extern int dsyevx_(const char *jobz,  const char *range, const char *uplo, const int *n, double *a, const int *lda,
                       const double *vl, const double *vu, const int *il, const int *iu,
                       const double *abstol, int *m, double *w, double *z, const int *ldz,
                       double *work, int *lwork, int *iwork, int *ifail, int *info,
                       int jobz_len, int range_len, int uplo_len);
    /// Solves system which has been LU-factorized.
    extern void dgetrs_(const char *trans, const int *n, const int *nrhs, double *a, const int *lda, int *ipiv, const double *b, const int *ldb, int *info);
    /// General matrix multiplication
    extern void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
                       const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc,
                       int a_columns, int b_columns, int c_columns);
    /// General dyad product of vectors
    extern void dger_(const int *m, const int *n, const double *alpha, const double *x, const int *incx,
                      const double *y, const int *incy, double *a, const int *lda,
                      int x_len, int y_len, int a_columns);
    /// Symmetric dyad product of vector
    extern void dsyr_(const char *uplo, const int *n, const double *alpha, const double *x, const int *incx,
                      double *a, const int *lda, int x_len, int a_columns);
    /// Y = Y + alpha * X
    extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy, int xsize, int ysize);
    /// X = alpha * X
    extern void dscal_(const int *n, const double *alpha, const double *x, const int *incx, int size);
}
#endif


#endif // numerics_h


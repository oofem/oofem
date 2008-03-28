/* $Header: /home/cvs/bp/oofem/oofemlib/src/flotmtrx.h,v 1.17.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/
/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/


//   **************************
//   *** CLASS FLOAT MATRIX ***
//   **************************


#ifndef flotmtrx_h
#define flotmtrx_h

#include "compiler.h"
#include "matrix.h"
#include "debug.h"
#include "freestor.h"
#include "cltypes.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#endif

class FloatArray ; class DiagonalMatrix ; class IntArray;
class DataStream;
#ifdef __PARALLEL_MODE
class CommunicationBuffer;
#endif

/**
 Implementation of matrix containing floatin point numbers. FloatMatrix can grow and shring
 to requested size if required. Implementation includes many computing and manipulation methods.
 Rows and Columns are indexed starting from 1.
*/
class FloatMatrix : public Matrix
{
/*
   This class implements a matrix containing floating-point numbers.
 DESCRIPTION :
   The matrix stores its nRows*nColumns coefficients in 'values'. The coef-
   ficients are stored column by column, like in Fortran.
 TASKS :
   - storing and retrieveing a coefficient (method 'at') ;
   - performing standard operations : addition, multiplication, transposition,
     inversion, lumping, rotation, etc ;
   - introduced allocatedSize variable to allow dynamic rescaling of matrix
     size possibly without memory realocation. At startup matrix occupies space
     given by allocatedSpace = nRows*nColums. Then there can be 
     1) further request for resizeing matrix to smaller dimension
        then we only change nRows and nColumns wariables , but allocatedSize
 variable remain untouched - expecting possible matrix grow and then re-using
 previously allocated space.
     2) if further request for growing then is necessary memory realocation.
     This process is controlled in resize member function. 
*/

 protected:
  /// values of matrix stored columwise
  double*  values ;
   /// allocated size for values
   int      allocatedSize;
  
   public:
 /** Constructor. Creates zeroed FloatMatrix of given size.
  @param n number of rows
  @param m requested number of columns
  */
 FloatMatrix (int n,int m):Matrix(n,m) {allocatedSize = n*m; values = allocDouble(n*m);}
 /// Constructor, creates zero sized matrix.
 FloatMatrix () : Matrix(0,0)          {allocatedSize = 0; values = NULL ;}
 /** Constructor. Creates float matrix from float vector. Vector may be stored rowwise 
  or columnwise, depending on second parameter.
  @param vector float vector from which matrix is constructed
  @param transpose if transpose==0 (default) then matrix of size (vector->giveSize(),1)
  will be created and initialized. if transpose!=0 then (1,vector->giveSize()) matrix
  will be created.
  */
 FloatMatrix (FloatArray *vector ,int transpose =0 ) ;
 /// copy constructor
 FloatMatrix (const FloatMatrix& );    // copy constructor
 /// destructor
 virtual ~FloatMatrix ()               {if(values) freeDouble(values);}
 /// assingnment operator, adjusts size of the receiver if necessary
 FloatMatrix&  operator=  (const FloatMatrix&);  // assignment: cleanup and copy


#ifdef DEBUG
 /** Coefficient access function. Returns value of  coeffiicient at given 
  position of the receiver. Implements 1-based indexing.
  @param i row position of coefficient
  @param j column coefficient possition
  */
 double    at (int i,int j) const;

 /** Coefficient access function. Returns l-value of coeffiicient at given 
  position of the receiver. Implements 1-based indexing.
  @param i row position of coefficient
  @param j column coefficient possition
  */
 double&   at (int i,int j) ;
#else
 double    at (int i,int j) const   {return values[(j-1)*nRows+i-1];}
 double&   at (int i,int j)         {return values[(j-1)*nRows+i-1];}
#endif

 /** Coefficient access function. Returns l-value of coeffiicient at given 
  position of the receiver. Implements 0-based indexing.
  @param i row position of coefficient
  @param j column coefficient possition
  */
  double&       operator()(int, int) ;
 /** Coefficient access function. Implements 0-based indexing.
  @param i row position of coefficient
  @param j column coefficient possition
  */
 double  operator()(int, int) const;          


/* new reference based interface */
 /**
  Assembles the contribution using localization array into receiver. The receiver must
  have dimensions large enough to localize contribution.
  @param src source to be assembled.
  @param loc localization numbers.
  */
 void          assemble (const FloatMatrix& src, const IntArray& loc);
 /**
  Assembles the contribution using localization array into receiver. The receiver must
  have dimensions large enough to localize contribution.
  @param src source to be assembled.
  @param rowind row localization numbers
  @param colind column localization numbers
  */
 void          assemble (const FloatMatrix& src, const IntArray& rowind, const IntArray& colind);
 /**
  Assembles the contribution using localization array into receiver. The receiver must
  have dimensions large enough to localize contribution.
  @param src source to be assembled.
  @param rowind row localization numbers
  @param colind column localization numbers
  */
 void          assemble (const FloatMatrix& src, const int* rowind, const int* colind);
 

 /** Computes the Frobenius norm of the receiver 
    (defined as the square root of the sum of the absolute squares of its elements)
 */
 double        computeFrobeniusNorm();

 /** Returns determinant of the receiver. Receiver should be square matrix. 
  Current implementation works for (3,3) and smaller matrices.
  @return determinat of receiver
  */
 double        giveDeterminant () ;
 /// Zeroes all coeficint of receiver.
 void          zero () const;
 /// sets receiver to unity matrix.
 void beUnitMatrix  ();

 /** Assigns to the receiver the tharansposition of parameter.
    Grows or shrinks if necessary */
 void          beTranspositionOf (const FloatMatrix& src) ;
 /** Assigns to the receiver product of aMatrix * bMatrix.
    Grows or shrinks if necessary */
 void          beProductOf (const FloatMatrix& aMatrix, const FloatMatrix& bMatrix);
 /** Assigns to the receiver product of aMatrix^T * bMatrix.
    Grows or shrinks if necessary */
 void          beTProductOf (const FloatMatrix& aMatrix, const FloatMatrix& bMatrix);
 /** Assigns to the receiver product of aMatrix * bMatrix^T.
    Grows or shrinks if necessary */
 void          beProductTOf (const FloatMatrix& aMatrix, const FloatMatrix& bMatrix);
 /**
  Adds the given matrix as sub-matrix to receiver. The sub-matrix values will be added to receivers
  corresponding receiver's values at positions (ri...ri+src.nrows, ci....ci+src.ncolumns). 
  The size of receiver will be adjusted, if necesary.
  @param src the sub-matrix to be added
  @param sr determines the row position (receiver's 1-based index) of first src value to be added.
  @param sc determines the column position (receiver's 1-based index) of first src value to be added.
  */
 void addSubMatrix (const FloatMatrix& src, int sr, int sc);
 /** Assigns to the receiver the sub matrix of another matrix.
    @param src matrix from which sub matrix is taken
   @param topRow index of top row, where sub matrix row index starts
   @param bottomRow index of bottom row, where sub matrix ends (including this row)
   @param topCol index of top column of submatrix
   @param bottomCol index of bottom column of submatrix
  */
 void          beSubMatrixOf(const FloatMatrix& src, int topRow, 
                             int bottomRow, int topCol, int bottomCol);
 /** Assigns to the receiver the sub matrix of another matrix.
    Works only for square matrices. Should produce bigger matrix than source matrix.
   @param src source matrix
   @param indx receiver (i.e. submatrix) will resize to square matrix of size given by
   maximum value found in indx parameter. On its position (indx->at(i), indx->at(j)) the
   value of src->at(i,j) will be stored. If one of indx->at(i), indx->at(j) indexes is 
   less or equal to zero, then on this position zero will be stored.
  */
 void          beSubMatrixOf (const FloatMatrix& src, const IntArray &indx);
 /** Assigns to the receiver the sub matrix of another matrix.
    Works only for square matrices. Should produce bigger matrix than source matrix.
   @param src source matrix
   @param indx Describes submatrix extaraction. On receiver  position 
   (indx->at(i), indx->at(j)) the
   value of src->at(i,j) will be stored. If one of indx->at(i), indx->at(j) indexes is 
   less or equal to zero, then on this position zero will be stored.
   @param size Receiver becomes square (size, size) matrix.
  */
 void          beSubMatrixOfSizeOf (const FloatMatrix& src, const IntArray &indx, int size);
 /** Adds given vector to receiver row sr, starting at column sc */
 void addSubVectorRow (const FloatArray& src, int sr, int sc);
 /** Copy (set)  given vector to receiver row sr, starting at column sc */
 void copySubVectorRow (const FloatArray& src, int sr, int sc);


 /**
  Modifies receiver to become inverse of given parameter. Size of receiver will be adjusted.
  */
 void          beInverseOf (const FloatMatrix& src);
 /**
  Solves the  system of linear equations. Uses Gaussian elimination with pivoting.
  @param b RHS of linear system.
  @param answer solution of linear equations.
  */
 void          solveForRhs (const FloatArray& b, FloatArray& answer) ;
 /**
  Solves the  system of linear equations. Uses Gaussian elimination with pivoting.
  @param b RHS of linear system, stored columnwise.
  @param answer solution of linear equations, each column corresponding to b column.
  */
 void          solveForRhs (const FloatMatrix& b, FloatMatrix& answer);
 /**
  Adds to the receiver the product  \f$a^T b dV\f$. If the receiver has zero size, it is expanded.
  Assumes that receiver and product \f$a^T b dV\f$ are symmetric matrices. Computes only the 
  upper half of receiver.
  @param a float matrix 
  @param b float matrix
  @param dV double value
  */
 void          plusProductSymmUpper (const FloatMatrix& a, const FloatMatrix& b,double dV) ;
 /**
  Adds to the receiver the product  \f$a^T b dV\f$. If the receiver has zero size, it is expanded.
  @param a float matrix 
  @param b float matrix
  @param dV double value
  */
 void          plusProductUnsym (const FloatMatrix& a, const FloatMatrix& b,double dV) ;
 /**
  Adds matrix to the receiver. If receiver has zero size, size is accordingly adjusted.
  @param aMatrix matrix to be added.
  */
 void          plus        (const FloatMatrix& aMatrix);
 /**
  Assigns to receiver one column or one row matrix containing vector. 
  @param transposed if (transposed == 0) then (vector->giveSize(),1) FloatMatrix is assigned
  else (1,vector->giveSize()) FloatMatrix is assigned.
  */
 void          initFromVector (const FloatArray& vector, int transposed);
 // void          add (const FloatMatrix& src);
 /**
  Initializes the lower half of the receiver according to the upper half.
  */
 void symmetrized () ;
 /**
    Returns the receiver 'a' transformed using givem transformation matrix r.
    The method performs the operation  a = r(transp) * a * r .
  Warning : this works only for square matrices and  the receiver and r must have the same size!
  @param r transformation matrix
  */
 void rotatedWith (const FloatMatrix& r) ;
 /**
  Sets size of receiver to be an empty matrix. It will have zero rows and zero columns size.
  */
 void          beEmptyMtrx () {this->resize(0,0);}
 /** Checks size of receiver towards requested bounds.
  If dimension mismatch, size is adjusted accordingly.
  Warning: after this operation array values are in undefined state, programmer should
  zero receiver 
  @param allocChunk if reallocation needed, an aditional space for allocChunk values 
  */
 void          resize (int, int, int allocChunk = 0); 
 /** Checks size of receiver towards requested bounds.
  If dimension mismatch, size is adjusted accordingly.
  Note: New coefficients are initialized to zero, old are kept.
  */
 void          resizeWithData (int, int); 

/* old pointer like functions */

 virtual FloatMatrix*  GiveCopy () ;
 //virtual FloatMatrix*  beCopyOf (FloatMatrix*);
 //FloatArray*   SolveForRhs (FloatArray *b) ;
 //FloatMatrix*  GiveInverse () ;
 FloatMatrix*  GiveTransposition () ;
 //FloatMatrix*  Lumped () ;
 //void            beLumpedOf (const FloatMatrix& src);
 //FloatMatrix*  Over (double f)         {return this->Times(1./f) ;}
 //FloatMatrix*  plus (FloatMatrix*) ;
 //FloatMatrix*  plusDiagonalMatrix (DiagonalMatrix*) ;
 void          assemble (FloatMatrix* src, IntArray* loc) {this-> assemble (*src, *loc);}
 //void          plusProductSymmUpper (FloatMatrix*,FloatMatrix*,double) ;
 //FloatMatrix*  RotatedWith (FloatMatrix*) ;
 /**
  Multiplies the receiver by given number.
  */
 FloatMatrix*  times (double f) ;
 FloatMatrix*  Times (double f)     {return this->GiveCopy()->times(f);}
 FloatArray*   Times (FloatArray*) ;
 FloatMatrix*  Times (FloatMatrix*) ;
 double*      givePointer ()  const     { return values ;} 
 //FloatArray*   timesTo (FloatArray*, FloatArray*);
 //FloatMatrix*  timesTo (FloatMatrix*, FloatMatrix*);
 //FloatMatrix*  GiveSubMatrix ( int, int, int, int);
 //FloatMatrix*  GiveSubMatrix ( IntArray *);
 //FloatMatrix*  GiveSubMatrixOfSize ( IntArray *, int );

 //void          Jacobi(FloatArray *d, FloatMatrix *v, int *nrot) ;
  /**
   Computes eigenvalues and eigenvectors of receiver (must be symmetric)
   The receiver is preserved. 
   @param eval eigenvalues
   @param v eigen vectors (stored columwise)
   @param nf number of significant figures
   @return 0 if o.k.
   */
  int           jaco_(FloatArray &eval, FloatMatrix &v, int nf);

 FloatMatrix*  hardResize (int, int); // hard memory realoocation 
 void          printYourself () const;
 contextIOResultType           storeYourself(DataStream* stream, ContextMode mode);
 contextIOResultType           restoreYourself(DataStream* stream, ContextMode mode);

#ifdef __PARALLEL_MODE
 /**@name Methods for  packing/unpacking to/from communication buffer */
  //@{
  /**
   Packs receiver into communication buffer.
   @param buff buffer to pack itself into.
   @return nonzero if succesfull
  */
  int packToCommBuffer (CommunicationBuffer& buff) const;
  /**
   Unpacks receiver from communication buffer.
   @param buff buffer from which unpack itself.
   @return nonzero if succesfull
  */
  int unpackFromCommBuffer (CommunicationBuffer& buff);
  //@}
  /**
   Returns how much space is needed to pack receivers message.
   @param buff buffer used for packing
   */
  int givePackSize (CommunicationBuffer& buff);
#endif

 protected:
} ;

#endif // flotmtrx_h


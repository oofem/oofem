/* $Header: /home/cvs/bp/oofem/oofemlib/src/sparsemtrx.h,v 1.12.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *************************************************************************
//   *** CLASS SPARSE MATRIX BASE CLASS                           ************
//   *************************************************************************

#ifndef sparsemtrx_h 
#define sparsemtrx_h 

//#ifdef __GNUG__
//#pragma interface
//#endif

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "cltypes.h"
#include "matrix.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "iml.h"


class EngngModel;
class TimeStep;
/**
 Base class for all matrices stored in sparse format. Basicaly sparce matrix 
 contains contribution of local element matrices. Localization of local element
 matrix into global (structural) matrix is determined using element code numbers.
 Basic methods include building internal structure of sparse matrix 
 (according to code numbers of elements), assembling of local element matrices,
 multiplication by array, and possible factorization and back substitution.
*/
class SparseMtrx : public Matrix
{
/*
  This class implem
ents a base class for sparse matrices containing 
  floating point numbers.
 DESCRIPTION :
  The sparse matrix contains contributions (local matrices) from FE element.
  Mapping between local element matrix values and SparseMatrix values is determined
  by code number array of each particular element.
 TASKS :
   - building its internal storage structure (method 'buildInternalStructure')
   - store and localize local mtrices (method 'localize')
  - performing standard operations : multiplication by array (method 'times')
  - possible factorization and backSubstitution (recognized by nonzero result of 
    canBeFactorized) (methods 'factorize' and 'bacSobstitution')
  - setting all coefficients to zero (method 'zero')
 */ 
 public:
 typedef long SparseMtrxVersionType;
 protected:
  /** 
   Allows to track if receiver changes.
   Any change on receiver should increment this variable.
   The factorization should not change version (nrsolver direct controll relies on that).
   This version info is used for example by preconditioners associated to
   particular matrix; the preconditioner initialization can be demanding
   and this versioning allows to reuse initilized preconditioner for same
   matrix, if there is no change;
   */
  SparseMtrxVersionType version;

 public:
 /** Constructor, creates (n,m) sparse matrix. Due to sparsity character of matrix, 
  not all coefficient are physicaly stored (in general, zero members are ommited).
  */
 SparseMtrx (int n,int m) : Matrix (n,m) {version = 0;}
 /// Constructor
 SparseMtrx () : Matrix () {version = 0;}
 // plus copy and assignment operators defined by derived classes
 //SparseMtrx& SparseMtrx::operator=(const SparseMtrx &C)  ;
 //SparseMtrx::SparseMtrx(const SparseMtrx &S) ;
 

 /** Return receiver version */
 SparseMtrxVersionType giveVersion() {return this->version;}

 /** Returns {\bf newly allocated} copy of receiver. Programmer must take 
  care about proper deallocation of allocated space.
  @return newly allocated copy of receiver */ 
  virtual SparseMtrx* GiveCopy () const = 0;

 /** Evaluates a product of receiver with vector. 
  @param x array to be multiplied with receiver
  @param answer result of product of receiver and x parameter
  */
 virtual void times (const FloatArray& x, FloatArray & answer) const = 0;
 /** Multiplies receiver by scalar value.
  @param x value to multiply receiver
  */
 virtual void times (double x) = 0;
 /**
  Builds internal structure of receiver. This method determines the internal profile
  of sparse matrix, allocates necessary space for storing nonzero coefficients and
  initializes receiver. In general, the profile of sparse matrix is determined 
  using one (or more) loop over local code numbers of elements.
  This method must be called before any operation, like assembly, zeroing, 
  or multiplication.
  The UnknownType parameter allows to distinguis between several possible governing equations, that 
  can be numbered separately.
  @param eModel pointer to corresponding engineering model
  @param di domain index specify which domain to use 
  @param ut unknown type
  */
 virtual int buildInternalStructure (EngngModel* eModel, int di, EquationID ut) = 0; 
 // virtual int assemble (FloatMatrix*, IntArray*) = 0;
 /** 
  Assembles sparse matrix from contribution of local elements. This method for
  each element adds its contribution to itself. Mapping between local element 
  contribution and its global position is given by local code numbers of element.
  @param loc location array. The values corresponding to zero loc array value are not assembled.
  @param mat contribution to be assembled using loc array.
  */
 virtual int assemble (const IntArray& loc, const FloatMatrix& mat) = 0;
 /** 
  Assembles sparse matrix from contribution of local elements. This method for
  each element adds its contribution to itself. Mapping between local element 
  contribution and its global position is given by row and column local code numbers.
  @param rloc row location array. The values corresponding to zero loc array value are not assembled.
  @param cloc column location array. The values corresponding to zero loc array value are not assembled.
  @param mat contribution to be assembled using rloc and cloc arrays. The rloc position determines the row, the 
  cloc determines the corresponding column.
  */
 virtual int assemble (const IntArray& rloc, const IntArray& cloc, const FloatMatrix& mat) = 0;

 virtual int assembleBegin () {return 1;}
 virtual int assembleEnd   () {return 1;}
 
 /// Determines, whether receiver can be factorized.
 virtual int canBeFactorized () const = 0;
 /**
  Returns the receiver factorized. \f$L^T D L\f$ form is used.
  @return pointer to the receiver
  */
 virtual SparseMtrx* factorized () {return NULL;}
 /**
  Computes the solution of linear system \f$A x = y\f$. A is receiver. 
  solution vector x overwrites the right hand side vector y.
  Receiver must be in factorized form.
  @param y right hand side on input, solution on output.
  @return pointer to y array
  @see factorized method
  */
 virtual FloatArray* backSubstitutionWith (FloatArray& y) const {return NULL;}
 /// Zeroes the receiver.
 virtual SparseMtrx* zero () = 0;

 /// Returns coefficient at position (i,j).
 virtual double& at (int i, int j) = 0;
 /// Returns coefficient at position (i,j).
 virtual double at (int i, int j) const = 0;
 virtual void toFloatMatrix (FloatMatrix& answer) const = 0;
 /// Prints the receiver statistics (one-line) to stdout.
 virtual void printStatistics () const {}
 /// Prints receiver to stdout. Works only for relatively small matrices.
 virtual void printYourself () const = 0;

 /// Sparse matrix type identification
 virtual SparseMtrxType  giveType() const = 0;
 /// Returns nonzero if anti-symmetric
 virtual int isAntisymmetric () const = 0;


#ifdef IML_COMPAT
// /***********************************/
//  /*  Matrix/Vector multiply         */
//  /***********************************/

 virtual FloatArray operator*(const FloatArray &x) const
  {FloatArray answer; this->times (x, answer) ; return answer;}
 virtual FloatArray trans_mult(const FloatArray &x) const = 0;

#endif



};

#endif // sparsemtrx_h 


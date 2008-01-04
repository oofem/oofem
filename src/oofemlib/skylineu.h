/* $Header: /home/cvs/bp/oofem/oofemlib/src/skylineu.h,v 1.2 2003/04/06 14:08:25 bp Exp $ */
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

//   ***********************************
//   *** CLASS SKYLINE (UNSYMMETRIC) ***
//   ***********************************

#ifndef skylineu_h

#include "sparsemtrx.h"
#include "rowcol.h"


/// "zero" pivot for SkylineUnsym class
#define SkylineUnsym_TINY_PIVOT 1.e-30

/**
 This class implements a nonsymmetric matrix stored in a compacted
 (skyline) form. Its shape is symmetric, but not its coefficients.
 The Skyline contains 'size' row-column segments, each of them of any size;
 these are stored in the attribute 'rowColumns'.
*/
class SkylineUnsym : public SparseMtrx
{

protected:
 /// Row column segments
 RowColumn**  rowColumns ;
 /// size of receiver
 int size;
 /// factorizationFlag
 int isFactorized;

public:
   /** Constructor. Before any operation an internal profile must be built.
    @see builInternalStructure
   */
      SkylineUnsym (int n) ;   
    /** Constructor. Before any operation an internal profile must be built.
     @see builInternalStructure
   */
    SkylineUnsym () ;
    /// Destructor
      ~SkylineUnsym () ;                           // destructor

    /** Returns {\bf newly allocated} copy of receiver. Programmer must take 
     care about proper deallocation of allocated space.
     @return newly allocated copy of receiver */ 
    SparseMtrx* GiveCopy () const ;
      /** Evaluates a product of receiver with vector. 
        @param x array to be multiplied with receiver
        @param answer result of product of receiver and x parameter
     */
     void times (const FloatArray& x, FloatArray& answer) const;
     /** Multiplies receiver by scalar value.
        @param x value to multiply receiver
     */
     virtual void times (double x) ;
     /// Builds internal structure of receiver
     int buildInternalStructure (EngngModel*, int, EquationID) ; 
     /**
        allocates and built structure according to given
        array of maximal column heights
     */
     int setInternalStructure (IntArray* a) ;

 /** Assembles receiver from local element contributions.
    @param loc location array. The values corresponding to zero loc array value are not assembled.
    @param mat contribution to be assembled using loc array.
   */
 int assemble (const IntArray& loc, const FloatMatrix& mat) ;
 /** Assembles receiver from local element contributions.
    @param rloc row location array. The values corresponding to zero loc array value are not assembled.
    @param cloc column location array. The values corresponding to zero loc array value are not assembled.
    @param mat contribution to be assembled using loc array.
   */
 int assemble (const IntArray& rloc, const IntArray& cloc, const FloatMatrix& mat) ;
 
 /// Determines, whether receiver can be factorized.
 int canBeFactorized () const {return 1;}
 /// Performs factorization of receiver
 SparseMtrx* factorized () ;
 /**
  Computes the solution of linear system \f$A x = y\f$. A is receiver. 
  solution vector x overwrites the right hand side vector y.
  Receiver must be in factorized form.
  @param y right hand side on input, solution on output.
  @return pointer to y array
  @see factorized method
  */
 FloatArray* backSubstitutionWith (FloatArray& ) const;
 /// Zeroes the receiver.
 virtual SparseMtrx* zero ();

 /// Returns coefficient at position (i,j).
 virtual double& at (int i, int j) ;
 virtual double  at (int i, int j) const;
 virtual void toFloatMatrix (FloatMatrix& answer) const;
 /// Prints receiver to stdout. Works only for relatively small matrices.
 virtual void printYourself () const;
 virtual void printStatistics () const;
 SparseMtrxType  giveType() const {return SMT_SkylineU;}
 int isAntisymmetric () const {return 1;}

protected:
/*
      FloatMatrix*  AsFloatMatrix () ;
      double&       at (int i,int j) ;
      void          assemble (FloatMatrix*,IntArray*) ;
      void          assemble (FloatMatrix*,IntArray*,IntArray*) ;
      FloatArray*   backSubstitutionWith (FloatArray*) ;
      void          carveYourselfFor (Domain*) ;
      int           computeNumberNegativeEigenvalue();
      void          checkSizeTowards (IntArray*) ;
      Skyline*      diagonalScalingWith (FloatArray*) ;
      Skyline*      factorized () ;
      Skyline*      forwardReductionWith (FloatArray*) ;
*/
 void          checkSizeTowards (const IntArray&) ;
 void          checkSizeTowards (const IntArray& rloc, const IntArray& cloc);
 RowColumn*    giveRowColumn (int j) const ;
 void          growTo (int) ;
/*
  void          printYourself () ;
  Skyline*      reinitialized () ;
  SkylineSym*   giveSymmetricPart();
  int       computeNumberNegativePivotsOfSymPart();
*/
 
 SkylineUnsym (RowColumn**, int, int);

#ifdef IML_COMPAT
 /***********************************/
/*  Matrix/Vector multiply         */
/***********************************/

 virtual FloatArray trans_mult(const FloatArray &x) const;

#endif

} ;


#define skylineu_h
#endif


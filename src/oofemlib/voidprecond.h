/* $Header: /home/cvs/bp/oofem/oofemlib/src/voidprecond.h,v 1.6 2003/04/06 14:08:26 bp Exp $ */
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


#ifndef voidprecond_h
#define voidprecond_h

#include "flotarry.h"
#include "sparsemtrx.h"
#include "precond.h"


/**
 Class implementing void preconditioner.
*/
class VoidPreconditioner : public Preconditioner {

 private:

 public:
  /// Constructor. Creates the empty precontitioner.
  VoidPreconditioner (const SparseMtrx &a, InputRecord& attributes)  ;
  /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
  VoidPreconditioner () ;
  /// Destructor
  ~VoidPreconditioner (void) { };

  /// Construct preconditioning matrix of given matrix
  void init (const SparseMtrx &a) {};
  /// Solves the linear system
  void solve (const FloatArray &x, FloatArray&y) const {y=x;}
  /// Solves transposed system
  void trans_solve (const FloatArray &x, FloatArray&y) const {y=x;};
  /// returns the preconditioner name
  virtual const char*  giveClassName () const {return "VoidPreconditioner";}

};

#endif

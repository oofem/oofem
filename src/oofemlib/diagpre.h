/* $Header: /home/cvs/bp/oofem/oofemlib/src/diagpre.h,v 1.6 2003/04/06 14:08:23 bp Exp $ */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifndef diagpre_h
#define diagpre_h

#include "flotarry.h"
#include "sparsemtrx.h"
#include "precond.h"


/**
 Implementation of diagonal preconditioner
*/
class DiagPreconditioner : public Preconditioner {

 private:
  FloatArray diag_;

public:
 /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
 DiagPreconditioner (const SparseMtrx &, InputRecord& attributes) ;
 /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
 DiagPreconditioner () : Preconditioner() {}
 /// Destructor
  ~DiagPreconditioner (void) { };
 
 /**
  Initializes the receiver (constructs the precontioning matrix M) of given matrix.
 */
 virtual void init (const SparseMtrx& a);
 
 /// Solves the linear system
  void solve (const FloatArray &x, FloatArray&y) const;
 /// Solves transposed system
  void trans_solve (const FloatArray &x, FloatArray&y) const;
 /// returns the preconditioner name
 virtual const char*  giveClassName () const {return "DiagPre";}
  
  const double&     diag(int i) const { return diag_(i); }
  double&           diag(int i) { return diag_(i); }
};

#endif // diagpre_h

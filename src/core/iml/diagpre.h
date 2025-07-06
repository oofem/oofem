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

#include "floatarray.h"
#include "sparsemtrx.h"
#include "precond.h"

namespace oofem {
/**
 * Implementation of diagonal preconditioner
 */
class OOFEM_EXPORT DiagPreconditioner : public Preconditioner
{
private:
    FloatArray diag;

public:
    /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
    DiagPreconditioner(const SparseMtrx &, InputRecord & attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    DiagPreconditioner() : Preconditioner() { }
    /// Destructor
    virtual ~DiagPreconditioner(void) { }

    void init(const SparseMtrx &a) override;

    void solve(const FloatArray &rhs, FloatArray &solution) const override;
    void trans_solve(const FloatArray &rhs, FloatArray &solution) const override;

    const char *giveClassName() const override { return "DiagPre"; }
};
} // end namespace oofem
#endif // diagpre_h

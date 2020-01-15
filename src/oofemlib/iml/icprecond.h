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

#ifndef icprecond_h
#define icprecond_h

#include "floatarray.h"
#include "intarray.h"
#include "symcompcol.h"
#include "precond.h"

namespace oofem {
/**
 * Incomplete Cholesky IC(0) (no fill - up) preconditioner for
 * symmetric, positive definite matrices.
 */
class OOFEM_EXPORT CompCol_ICPreconditioner : public Preconditioner
{
private:
    FloatArray val;
    IntArray pntr;
    IntArray indx;
    int nz;
    int dim [ 2 ];

public:
    /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
    CompCol_ICPreconditioner(const SparseMtrx & A, InputRecord & attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    CompCol_ICPreconditioner() : Preconditioner() { }
    /// Destructor.
    virtual ~CompCol_ICPreconditioner(void) { }

    void init(const SparseMtrx &a) override;

    void initialize(const CompCol &A);

    void solve(const FloatArray &rhs, FloatArray &solution) const override;
    void trans_solve(const FloatArray &rhs, FloatArray &solution) const override;

    const char *giveClassName() const override { return "ICP"; }
    void initializeFrom(InputRecord &ir) override;

protected:
    void qsortRow(IntArray &, FloatArray &, int l, int r);
    int  qsortRowPartition(IntArray &, FloatArray &, int l, int r);

    void ICSolve(FloatArray &dest) const;
    void ICFactor();
};
} // end namespace oofem
#endif // icprecond_h

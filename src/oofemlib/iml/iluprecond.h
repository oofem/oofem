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

#ifndef iluprecond_h
#define iluprecond_h

#include "floatarray.h"
#include "intarray.h"
#include "compcol.h"
#include "dyncompcol.h"
#include "precond.h"

namespace oofem {
/**
 * Implemantation of ILU (Incomplete LU) Preconditioner.
 * No fill-up - ILU(0).
 */
class OOFEM_EXPORT CompCol_ILUPreconditioner : public Preconditioner
{
private:
    FloatArray l_val;
    IntArray l_colptr;
    IntArray l_rowind;
    int l_nz;

    FloatArray u_val;
    IntArray u_colptr;
    IntArray u_rowind;
    int u_nz;

    int dim [ 2 ];

public:
    /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
    CompCol_ILUPreconditioner(const SparseMtrx & A, InputRecord & attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    CompCol_ILUPreconditioner() : Preconditioner() { }
    /// Destructor
    virtual ~CompCol_ILUPreconditioner(void) { }

    void init(const SparseMtrx &) override;

    void initialize(const CompCol &A);
    void initialize(const DynCompCol &A);

    void solve(const FloatArray &x, FloatArray &y) const override;
    void trans_solve(const FloatArray &x, FloatArray &y) const override;

    const char *giveClassName() const override { return "ILU"; }
    void initializeFrom(InputRecord &ir) override;

protected:
    void qsortRow(IntArray &, FloatArray &, int l, int r);
    int  qsortRowPartition(IntArray &, FloatArray &, int l, int r);
};
} // end namespace oofem
#endif // iluprecond_h

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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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


#ifndef ilucomprowprecond_h
#define ilucomprowprecond_h

#include "floatarray.h"
#include "intarray.h"
#include "dyncomprow.h"
#include "precond.h"

///@name Input fields for CompRowPrecond
//@{
#define _IFT_CompRow_ILUPrecond_droptol "droptol"
#define _IFT_CompRow_ILUPrecond_partfill "partfill"
//@}

namespace oofem {
/**
 * Implemantation of ILU (Incomplete LU) Preconditioner for compressed row sparse matrices.
 * Fill - up supported.
 */
class OOFEM_EXPORT CompRow_ILUPreconditioner : public Preconditioner
{
private:
    DynCompRow A;

    double drop_tol;
    int part_fill;

public:
    /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
    CompRow_ILUPreconditioner(const SparseMtrx & A, InputRecord & attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    CompRow_ILUPreconditioner() : Preconditioner() { }
    /// Destructor
    virtual ~CompRow_ILUPreconditioner(void) { }

    void init(const SparseMtrx &a) override;

    //void initialize (const CompCol &A);
    void initialize(const DynCompRow &A);

    void solve(const FloatArray &x, FloatArray &y) const override;
    void trans_solve(const FloatArray &x, FloatArray &y) const override;

    const char *giveClassName() const override { return "ILUT"; }
    void initializeFrom(InputRecord &ir) override;

protected:
    void qsortCol(IntArray &, FloatArray &, int l, int r);
    int  qsortColPartition(IntArray &, FloatArray &, int l, int r);
};
} // end namespace oofem
#endif // ilucomprowprecond_h

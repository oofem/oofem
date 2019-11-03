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
 *               Copyright (C) 1993 - 2017   Borek Patzak
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

#ifndef superlusolver_h
#define superlusolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "SUPERLU_MT/include/slu_mt_ddefs.h"

#define _IFT_SuperLUSolver_Name "superlu"

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;
/**
 * Class implementig interface to SuperLU_MT solver.
 * Dependencies: iml module required (CompCol)
 * Note: There is a conflict in include file (set.h) as defined by libstdc and superlu_mt 3.1.
 *       Therefore it is assumed that superlu include files are included using #include "SUPERLU_MT/include".
 *       The ${SUPERLU_MT_DIR} is added into compiler include path.
 *       The SUPERLU_MT/include has to be manually added pointing to ${SUPERLU_MT__DIR}/src directory
 */
class OOFEM_EXPORT SuperLUSolver : public SparseLinearSystemNM
{
private:

public:
    SuperLUSolver(Domain * d, EngngModel * m);
    /// Destructor
    virtual ~SuperLUSolver();
    /**
     * Solves the given linear system iteratively by method described by SuperLU SolverType.
     * @param A Coefficient matrix.
     * @param b Right hand side.
     * @param x Solution array.
     * @return Status value.
     */
    NM_Status solve(SparseMtrx &A, FloatArray &b, FloatArray &x) override;
    SparseMtrxType giveRecommendedMatrix(bool symmetric) const override { return SMT_CompCol; }


    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "SuperLUSolver"; }
    LinSystSolverType giveLinSystSolverType() const override { return ST_SuperLU_MT; }

private:

    int_t cholnzcnt(int_t neqns, int_t *xadj, int_t *adjncy, int_t *perm, int_t *invp, int_t *etpar, int_t *colcnt, int_t *nlnz, int_t *part_super_L);
    void convertRhs(SuperMatrix *A, FloatArray &x);
    int_t dPrint_CompCol_Matrix(SuperMatrix *A);
    int_t dPrint_Dense_Matrix(SuperMatrix *A);
};
} // end namespace oofem
#endif // imlsolver_h

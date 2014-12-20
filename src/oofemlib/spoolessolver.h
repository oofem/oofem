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
#ifndef spoolessolver_h
#define spoolessolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
extern "C" {
#include <spooles/misc.h>
#include <spooles/FrontMtx.h>
#include <spooles/SymbFac.h>
};


///@name Input fields for SpoolesSolver
//@{
#define _IFT_SpoolesSolver_msglvl "msglvl"
#define _IFT_SpoolesSolver_msgfile "msgfile"
//@}

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x = b @f$ using solvers
 * from SPOOLES library. Can work with only SPOOLES sparse matrix implementation.
 */
class OOFEM_EXPORT SpoolesSolver : public SparseLinearSystemNM
{
private:
    /// Last mapped LHS matrix
    SparseMtrx *Lhs;
    /// Last mapped matrix version
    SparseMtrx :: SparseMtrxVersionType lhsVersion;
    int msglvl;
    FILE *msgFile;
    int msgFileCloseFlag;

    FrontMtx *frontmtx;
    IV *oldToNewIV, *newToOldIV;
    ETree *frontETree;
    IVL *adjIVL, *symbfacIVL;
    SubMtxManager *mtxmanager;
    Graph *graph;

public:
    /**
     * Constructor.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    SpoolesSolver(Domain * d, EngngModel * m);

    ///Destructor
    virtual ~SpoolesSolver();

    /**
     * Solves the given linear system by LDL^T factorization.
     */
    virtual NM_Status solve(SparseMtrx &A, FloatArray &b, FloatArray &x);

    /// Initializes receiver from given record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "SpoolesSolver"; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_Spooles; }
    virtual SparseMtrxType giveRecommendedMatrix(bool symmetric) const { return SMT_SpoolesMtrx; }
};
} // end namespace oofem
#endif // spoolessolver_h

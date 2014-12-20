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

//   ****************************************
//   *** CLASS SUBSPACE ITERATION  SOLVER ***
//   ****************************************


#ifndef subspaceit_h
#define subspaceit_h

#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements the class NumericalMethod instance Subspace Iteration Eigen Value Problem Solver
 *
 * DESCRIPTION :
 * Perform solution of eigen value problem in the form
 * K y = (omega)^2 M y
 *
 * TASKS :
 *
 * - solving problem
 *   solveYourselfAt.
 * - returning results (eigen values and associated eigen vectors).
 *
 * Variable description  :
 *
 *      - - - INPUT  DATA - - -
 *
 *      A(NWK)  -  STIFFNESS MATRIX
 *      B(NWM)  -  MASS MARTRIX
 *      NN  -      SIZE OF PROBLEM
 *      NNM  -  NN+1
 *      NROOT  -  REQUIRED NUMBER OF
 *      RTOL  -  KRITERIUM KONVERGENCE VLASTNICH CISEL
 *      NC  -  POCET VEKTORU SIMULTANNI ITERACE, DOPORUCUJE SE VOLIT
 *      NC = MIN (2*NROOT , NROOT+8 )
 *      NITEM  -  MAXIMALNI POCET ITERACI (OBYC. 16)
 *
 *      - - - PRACOVNI POLE - - -
 *
 *      TT(NN),W(NN),D(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC)
 *      AR(NC,NC)  -  PRACOVNI MATICE - PROJEKCE MATICE  A
 *      BR(NC,NC)  -  PROJEKCE MATICE  B
 *
 *      - - - VYSTUPNI DATA - - -
 *
 *      EIGV(NROOT)  -  VLASTNI CISLA
 *      R(NN,NROOT)  -  VLASTNI VEKTORY
 */

class OOFEM_EXPORT SubspaceIteration : public SparseGeneralEigenValueSystemNM
{
private:
    int n, nc, nitem;
    int solved;

public:
    SubspaceIteration(Domain * d, EngngModel * m);
    virtual ~SubspaceIteration();

    virtual NM_Status solve(SparseMtrx &A, SparseMtrx &B, FloatArray &x, FloatMatrix &v, double rtol, int nroot);
    virtual const char *giveClassName() const { return "SubspaceIterationSolver"; }
};
} // end namespace oofem
#endif // subspaceit_h

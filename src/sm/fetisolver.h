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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef fetisolver_h
#define fetisolver_h

#ifdef __PARALLEL_MODE

 #include "sparselinsystemnm.h"
 #include "sparsemtrx.h"
 #include "flotarry.h"
 #include "flotmtrx.h"
 #include "processcomm.h"
 #include "feticommunicator.h"
 #include "engngm.h"

namespace oofem {
class Domain;
class EngngModel;


 #define FETISOLVER_MAX_RBM 6
 /// computer zero
 #define FETISOLVER_ZERONUM 1.e-40

/**
 * This class implements the class NumericalMethod instance FETI
 * linear algebraic equation parallel solver.
 */
class FETISolver : public SparseLinearSystemNM
{
private:
    int nse;
    /// Max number of iterations.
    int ni;
    /// Max allowed error.
    double err;
    /// Linear dep./indep. trigger.
    double limit;
    /// Rigid body motions.
    FloatMatrix rbm;
    /// Rigid body motions of all partitions. On master only.
    FloatMatrix l;
    /// Indices of singular equations.
    IntArray se;
    /// Addresses of initial partition contribution to rbm matrix.
    IntArray rbmAddr;
    /// Primary unknowns.
    FloatArray w;
    //
    FloatArray qq, q, dd, g, d, p, pp, gamma, localGammas;
    //
    IntArray nsem;
    ///
    ProcessCommunicatorBuff pcbuff;
    ProcessCommunicator processCommunicator;
    /// Common Communicator buffer.
    CommunicatorBuff *commBuff;
    FETICommunicator *masterCommunicator;
    /// List of local nodes (at master) participating in communication (list of boundary dof managers).
    IntArray masterCommMap;
    /// Flag indicating computation of energy norm.
    int energyNorm_comput_flag;
public:
    FETISolver(int i, Domain *d, EngngModel *m);
    virtual ~FETISolver();

    /**
     * Solves the given linear system by LDL^T factorization.
     * Implementation rely on factorization support provided by mapped sparse matrix.
     * It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
     * @param A Coefficient matrix
     * @param b Right hand side
     * @param x Solution array
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);

    int estimateMaxPackSize(IntArray &, CommunicationBuffer &, int &);
    /// Sets up the communication maps
    void setUpCommunicationMaps();

    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "FETISolver"; }
    virtual classType giveClassID() const { return NumericalMethodClass; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_Feti; }

    void projection(FloatArray &v, FloatMatrix &l, FloatMatrix &l1);

    int packRBM(ProcessCommunicator &processComm);
    int masterUnpackRBM(ProcessCommunicator &processComm);
    int packQQProducts(ProcessCommunicator &processComm);
    int masterUnpackQQProduct(ProcessCommunicator &processComm);
    int packSolution(ProcessCommunicator &processComm);
    int unpackSolution(ProcessCommunicator &processComm);
    int packResiduals(ProcessCommunicator &processComm);
    int unpackResiduals(ProcessCommunicator &processComm);
    int packDirectionVector(ProcessCommunicator &processComm);
    int unpackDirectionVector(ProcessCommunicator &processComm);
    int packPPVector(ProcessCommunicator &processComm);
    int unpackPPVector(ProcessCommunicator &processComm);
    int packGammas(ProcessCommunicator &processComm);
    int unpackGammas(ProcessCommunicator &processComm);

    int masterMapRBM();
    int masterMapQQProduct();
    int masterMapSolution();
    int masterMapResiduals();
    int masterMapDirectionVector();
    int masterMapPPVector();
    int masterMapGammas();

    enum { FETISolverZeroTag, NumberOfRBMMsg, RBMMessage, QQMessage, SolutionMessage, ResidualMessage, DirectionVectorMessage, PPVectorMessage, GammasMessage, FETISolverIterationContinue, FETISolverIterationBreak };
};
} // end namespace oofem
#endif
#endif // fetisolver_h

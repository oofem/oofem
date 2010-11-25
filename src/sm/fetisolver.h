/* $Header: /home/cvs/bp/oofem/sm/src/fetisolver.h,v 1.6.4.1 2004/04/05 15:19:46 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   **********************************
//   *** CLASS FETI PARALLEL SOLVER ***
//   **********************************

#ifndef fetisolver_h
#define fetisolver_h

#ifdef __PARALLEL_MODE

 #include "sparselinsystemnm.h"
 #include "sparsemtrx.h"
 #include "flotarry.h"

 #include "skyline.h"
 #include "flotmtrx.h"
 #include "processcomm.h"
 #include "feticommunicator.h"
 #include "engngm.h"

 #ifndef __MAKEDEPEND
  #include <stdio.h>
 #endif

namespace oofem {
class Domain;
class EngngModel;


 #define FETISOLVER_MAX_RBM 6
/// computer zero
 #define FETISOLVER_ZERONUM 1.e-40

class FETISolver : public SparseLinearSystemNM
{
    /*
     * This class implements the class NumericalMethod instance FETI
     * linear algebraic equation parallel solver
     */
private:
    int nse;
    /// max number of iterations
    int ni;
    /// max allowed error
    double err;
    /// linear dep/indep. trigger
    double limit;
    //Skyline* partitionStiffness;
    //FloatArray* partitionLoad;
    //FloatArray* partitionSolution;
    /// rigid body motions
    FloatMatrix rbm;
    /// rigid body motiona of all partitions. On master only
    FloatMatrix l;
    /// indexes of singular equations
    IntArray se;
    // List of local nodes participating in communication (list of boundary dof managers).
    // IntArray    commMap;
    /// Adresses of initial partititon contribution to rbm matrix
    IntArray rbmAddr;
    /*  vektor neznamych  */
    FloatArray w;
    //
    FloatArray qq, q, dd, g, d, p, pp, gamma, localGammas;
    //
    IntArray nsem;
    ///
    ProcessCommunicatorBuff pcbuff;
    ProcessCommunicator processCommunicator;
    ///
    /// Common Communicator buffer
    CommunicatorBuff *commBuff;
    FETICommunicator *masterCommunicator;
    /// List of local nodes (at master) participating in communication (list of boundary dof managers).
    IntArray masterCommMap;
    /// flag indicating computation of energy norm
    int energyNorm_comput_flag;
public:
    FETISolver(int i, Domain *d, EngngModel *m);
    // constructor
    ~FETISolver();                    // destructor

    /**
     * Solves the given linear system by LDL^T factorization.
     * Implementation rely on factorization support provided by mapped sparse matrix.
     * It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
     * @param A coefficient matrix
     * @param b right hand side
     * @param x solution array
     * @return NM_Status value
     * @param tNow time step
     */
    NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);

    int                estimateMaxPackSize(IntArray &, CommunicationBuffer &, int &);
    /// Sets up the communication maps
    void               setUpCommunicationMaps();
    // management  components

    //      void               instanciateFromString (char* initString) {}
    IRResultType initializeFrom(InputRecord *ir);

    // identification
    const char *giveClassName() const { return "FETISolver"; }
    classType giveClassID() const { return NumericalMethodClass; }
    LinSystSolverType giveLinSystSolverType() const { return ST_Feti; }


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

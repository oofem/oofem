/* $Header: /home/cvs/bp/oofem/sm/src/nlinearstatic.h,v 1.12.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//
// Class NonLinearStatic
//

#ifndef nlinearstatic_h
#define nlinearstatic_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "linearstatic.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"

#ifdef __PARALLEL_MODE
#include "problemcomm.h"
#include "processcomm.h"
#endif

/**
 * Type determing the stiffness mode
 * nls_tangentStiffness - the tangent stiffness is used and updated whenever requsted
 * nls_secantStiffness  - the secant stiffness is used and updated whenever requsted
 * nls_elasticStiffness - the initial elastic stiffness is used in the whole solution
 * nls_secantInitialStiffness - the secant stiffness is used and updated only at the benining of new load step
 */
enum    NonLinearStatic_stifnessMode {
    nls_tangentStiffness = 0, nls_secantStiffness = 1, nls_elasticStiffness = 2,
    nls_secantInitialStiffness = 3
};
/**
 * Type determining type of loading controll. This type determines the solver to be used.
 * - the nls_indirectControll describes the indirect controll - a generalized norm of displacement and loading vectors
 * is controlled. In current implementation, the CALM solver is used, the reference load vector is FIXED.
 * - nls_directControll describes the direct controll where load or displacement (or both) are controlled.
 */
enum    NonLinearStatic_controllType { nls_indirectControll = 0, nls_directControll = 1, nls_directControll2 = 2 };

class NonLinearStatic : public LinearStatic
{
    /*
     * This class implements Non - LinearStatic Engineering problem.
     * DESCRIPTION:
     * Solution of this problem is performed  as a series of increments (loading or displacement).
     * At start of Each increment we assemble new tangent stiffness, and iteratively trying
     * to fullfill balance of external and real internal forces
     * at end of load step (see numerical method ).
     * The loading applied can bo of two types:
     * - proportional incremental  loading
     * - non-proportional fixed loading, reflecting the previous history,
     *  butt could not be scaled (like dead weight).
     *
     * TASK:
     * Creating Numerical method for solving nonlinear problem.
     * Assembling tangent stiffness matrix
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */


protected:
    double prevStepLength, currentStepLength;
    FloatArray totalDisplacement,  incrementOfDisplacement;
    FloatArray internalForces;
    /// initialLoadVector is a load vector already applied, which does not scales
    FloatArray initialLoadVector;
    /**
     * incremental Load Vector which is applied in loading step (as a whole for direct controll
     * or proportionally for indirect controll).
     */
    FloatArray incrementalLoadVector;
    /// initialLoadVectorOfPrescribed is a load vector which does not scale for prescribed DOFs
    FloatArray initialLoadVectorOfPrescribed;
    /// incremental Load Vector for prescribed DOFs
    FloatArray incrementalLoadVectorOfPrescribed;

    FloatArray incrementalBCLoadVector; // for direct controll
    double loadLevel, cumulatedLoadLevel, rtolv;
    bool mstepCumulateLoadLevelFlag;
    int currentIterations;
    NonLinearStatic_stifnessMode stiffMode;
    int loadInitFlag;
    int nonlocalStiffnessFlag;
    NM_Status numMetStatus;
    /// Numerical method used to solve the problem
    SparseNonLinearSystemNM *nMethod;
    /**
     * Characterizes the type of controll used.
     */
    NonLinearStatic_controllType controllMode;
    /**
     * Intrinsic time increment
     */
    double deltaT;
    /**
     * The following parameter allows to specify how the reference load vector
     * is obtained from given totalLoadVector and initialLoadVector.
     * The initialLoadVector desribes the part of loading which does not scale.
     * If refLoadInputMode is rlm_total (default) then the reference incremental load vector is defined as
     * totalLoadVector assembled at given time.
     * If refLoadInputMode is rlm_inceremental then the reference load vector is
     * obtained as incremental load vector at given time.
     */
    SparseNonLinearSystemNM :: referenceLoadInputModeType refLoadInputMode;

#ifdef __PARALLEL_MODE
    /// Message tags
    enum { InternalForcesExchangeTag, MassExchangeTag, LoadExchangeTag, RemoteElementsExchangeTag };
#endif

public:
    NonLinearStatic(int i, EngngModel *_master = NULL);
    ~NonLinearStatic();
    // solving
    //void solveYourself ();
    /**
     * Starts solution process. Invokes  for each time step
     * solveYourselfAt function with time step as parameter. Time steps are created
     * using giveNextStep function (this will set current time step to newly created,
     * and updates previous step).
     */
    void solveYourself();
    void solveYourselfAt(TimeStep *);
    void terminate(TimeStep *);

    /**
     * Prints output of receiver to ouput domain stream, for given time step.
     * Corresponding function for element gauss points is invoked
     * (gaussPoint::printOutputAt).
     */
    void                  printOutputAt(FILE *, TimeStep *);

    //int requiresNewLhs () {return 1;}
    virtual void               updateYourself(TimeStep *);
    virtual void updateComponent(TimeStep *, NumericalCmpn, Domain *);
    void                       updateAttributes(TimeStep *);

    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    double giveUnknownComponent(UnknownType, ValueModeType, TimeStep *, Domain *, Dof *);
    IRResultType initializeFrom(InputRecord *ir);
    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // identification
    const char *giveClassName() const { return "NonLinearStatic"; }
    classType giveClassID() const { return NonLinearStaticClass; }
    int isIncremental() { return 1; }
    fMode giveFormulation() { return TL; }
    /// Returns nonzero if nonlocal stiffness option activated.
    virtual int useNonlocalStiffnessOption() { return this->nonlocalStiffnessFlag; }
    /// For load balancing purposes we store all values with same EquationID; so hash is computed from mode value only
    virtual int       giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN)
    { return ( int ) mode; }

#ifdef __OOFEG
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    void               showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime);
#endif

#ifdef __PARALLEL_MODE
    /**
     * Exchanges necessary remote element data with remote partitions. The receiver's nonlocalExt flag must be set.
     * Uses receiver nonlocCommunicator to perform the task using packRemoteElementData and unpackRemoteElementData
     * receiver's services.
     * @return nonzero if success.
     */
    int exchangeRemoteElementData();
    /**
     * Determines the space necessary for send/receive buffer.
     * It uses related communication map pattern to determine the maximum size needed.
     * @param commMap communication map used to send/receive messages
     * @param buff communication buffer
     * @return upper bound of space needed
     */
    int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);

    /**
     * Initializes communication maps of the receiver
     */
    void initializeCommMaps(bool forceInit = false);

    /** returns reference to receiver's load balancer*/
    virtual LoadBalancer *giveLoadBalancer();
    /** returns reference to receiver's load balancer monitor*/
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor();


#endif

protected:
    void       assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type, Domain *domain);
    void giveInternalForces(FloatArray &answer, const FloatArray &DeltaR, Domain *d, TimeStep *);
    void proceedStep(int di, TimeStep *);
    void updateLoadVectors(TimeStep *tStep);
    // void        updateInternalStepState (const FloatArray &, TimeStep* );
    // void assembleInitialLoadVector (FloatArray& answer, TimeStep* atTime);
    /**
     * Computes the contribution external loading to reaction forces in given domain.
     * @param reactions contains the comuted contributions
     * @param tStep solution step
     * @param domain number
     */
    void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
    void assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                 FloatArray &_incrementalLoadVectorOfPrescribed,
                                                 SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                 Domain *sourceDomain, EquationID ut, TimeStep *tStep);
#ifdef __PARALLEL_MODE
    /** Packs receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renubering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void packMigratingData(TimeStep *);
    /** Unpacks receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renubering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void unpackMigratingData(TimeStep *);

#endif
};

#endif // nlinearstatic_h

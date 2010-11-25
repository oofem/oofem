/* $Header: /home/cvs/bp/oofem/sm/src/pnldeidynamic.h,v 1.6.4.2 2004/05/14 13:45:45 bp Exp $ */
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
// Class NlDEIDynamic - DirectExplicitIntegrationDynamic
//

#ifndef pnldeidynamic_h
#define pnldeidynamic_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif
#include "structengngmodel.h"
#include "skyline.h"

#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"


class PNlDEIDynamic;
class ProblemCommunicator;

/**
 * Pointer to packing function type.  Functions of this type are used to pack
 * various nodal data into domainCommunicators send buffer. Only nodes in toSend
 * map are considered.
 */
//typedef int (PNlDEIDynamic::*PNlDEIDynamic_Pack_func) (PNlDEIDynamicDomainCommunicator&);
/**
 * Pointer to unpacking function type.  Functions of this type are used to unpack
 * various nodal data from domainCommunicators receive buffer. Only nodes in toRecv
 * map are considered.
 */
//typedef int (PNlDEIDynamic::*PNlDEIDynamic_Unpack_func) (PNlDEIDynamicDomainCommunicator&);
/**
 * NlDEIDynamicCommunicatorMode determines the valid mode.
 * The mode is used to set up communication pattern, which differ for
 * node and element cut algorithms.
 * Additional remote element mode hass been added to capture the case, when CommunicatorM is intended to
 * support remote element data exchange (for example when nonlocal material models are present).
 */

#endif


/**
 * This class implements NonLinear (- may be changed) solution of dynamic
 * problems using Direct Explicit Integration scheme - Central Difference
 * Method. For efficiency reasons it uses diagonal mass matrix. It is formulated
 * in increments of displacements rather than in total variables.
 *
 * Solution of this problem is series of loading cases, maintained as sequence of
 * time-steps. For obtaining diagonal mass matrix from possibly non-diagonal one
 * returned from Element::giveMassMatrix() a FloatMatrix::Lumped() is called
 * to obtain diagonal form.
 *
 * Analysis starts assembling governing equations at time step 0 ( 0 given by boundary and initial cond.)
 * they result in response at time step 1. For time step 0 we need special start code.
 * Because this method is explicit, when solving equations for step t, we obtain
 * solution in step t+dt. But printing is performed for step t.
 * So, when analyst specifies initial conditions, then he/she specifies them in time step 0.
 *
 * WARNING - FloatMatrix::Lumped() works only for elements with Linear displacement filed !
 *
 * Current implementation supports parallel processing. Both node- and element cut strategies can
 * be used.
 * \begin{itemize}
 * \item
 * In node cut strategy, partitions are divided using cut, which goes through nodes.
 * These cutted nodes are called "shared" nodes. Generally, unknown values in shared nodes are
 * composed from local partition contributions as well as from contributions from remote partitions
 * sharing this node. Particullary, masses and real nodal forces have to be exchaneged for shared
 * nodes.
 * \item
 * In element cut strategy, partitions are divided using cut running through elements. The cutted elements are
 * replicated on neighbouring partitions. The nodes belonging to  replicated elements belonging to
 * remote partitions are called remote nodes. The are mirrors or remote copies of cooresponding
 * nodes on neighbouring partition.
 * \item additional mode has been introduced remote element mode. It introduces the "remote" elements, the
 * exact local mirrors of remote counterparts. Introduced to support general nonlocal constitutive models,
 * in order to provide efficient way, how to average local data without need of fine grain communication.
 * \end{itemize}
 */
class PNlDEIDynamic : public StructuralEngngModel
{
    /*
     * This class implements NonLinear (- may be changed) solution of dynamic
     * problems using Direct Explicit Integration scheme - Central Difference
     * Method. For efficiency reasons it uses diagonal mass matrix. It is formulated
     * in increments of displacements rather than in total variables.
     *
     * DESCRIPTION:
     * Solution of this problem is series of loading cases, maintained as sequence of
     * time-steps. For obtaining diagonal mass matrix from possibly non-diagonal one
     * returned from Element::giveMassMatrix() a FloatMatrix::Lumped() is called
     * to obtain diagonal form.
     *
     * we start assemble governing equations at time step 0 ( 0 given by boundary and initial cond.)
     * they result in response at time step 1.
     * for time step 0 we need special start code.
     * so we obtain solution for time step 1 and next.
     * because this method is explicit, when solving equations for step t, we obtain
     * solution in step t+dt. But printing is performed for step t.
     * see diidynamic.h for difference.
     * So, when You specify initial conditions, you specify them in time step 0.
     *
     * WARNING - FloatMatrix::Lumped() works only for elements with Linear displacement filed !
     *
     * TASK:
     * Creating Numerical method for solving Ax=b
     * Interfacing Numerical method to Elements
     * Managing time  steps
     */

protected:
    /// Mass matrix
    FloatArray massMatrix;
    /// Load vector
    FloatArray loadVector;
    /// Vector storing  displacement increnents
    FloatArray previousIncrementOfDisplacementVector;
    /// Displacement, velocity and acceleration vectors
    FloatArray displacementVector, velocityVector, accelerationVector;
    /// vector of real nodal forces
    FloatArray internalForces;
    /// dumping coefficient (C = dumpingCoef * MassMtrx)
    double dumpingCoef;
    /// Time step
    double deltaT;
    /// Flag indicating the need for initialization
    int initFlag;

    // dynamic relaxation specifiv vars
    /// flag indicating whether dynamic relaxation takes place
    int drFlag;
    /// reference load vector
    FloatArray loadRefVector;
    /// parameter determining rate of the loading process
    double c;
    /// load level
    double pt;
    /// end of time interval
    double Tau;
    /// estimate of loadRefVector^T*displacementVector(Tau)
    double pyEstimate;
    /// product of p^tM^(-1)p; where p is reference load vector
    double pMp;

#ifdef __PARALLEL_MODE
    // public:
    /// Message tags
    enum { InternalForcesExchangeTag, MassExchangeTag, LoadExchangeTag, RemoteElementsExchangeTag };
#endif

public:
    /// Constructor.
    PNlDEIDynamic(int i, EngngModel *_master = NULL);
    /// Destructor.
    ~PNlDEIDynamic();
    // solving
    /**
     * Starts solution process. Invokes  for each time step
     * solveYourselfAt function with time step as parameter. Time steps are created
     * using giveNextStep function (this will set current time step to newly created,
     * and updates previous step).
     */
    void solveYourself();
    /**
     * Solves problem for given time step.
     */
    void solveYourselfAt(TimeStep *);
    //int requiresNewLhs () {return 0;}
    /**
     * Updates internal state after finishing time step.
     */
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    virtual void               updateYourself(TimeStep *);
    /**
     * Returns requested unknown. Unknown at give time step is characterized by its type and mode
     * and by its equation number. This function is used by Dofs, when they are requsted for
     * their associated unknowns. Supports DisplacementVector type with following modes:
     * TotalMode, IncrementalMode, VelocityMode, AccelerationMode.
     * @see Dof::giveUnknown method
     */
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    /**
     * Reads receiver description from record stored in initString.
     */
    IRResultType initializeFrom(InputRecord *ir);
    /// Returns next time step (next to current step) of receiver.
    TimeStep *giveNextStep();
    NumericalMethod *giveNumericalMethod(TimeStep *);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Terminates the solution of time step.
     */
    void    terminate(TimeStep *);
    /**
     *  Prints output of receiver to ouput stream, for given time step.
     */
    void                  printOutputAt(FILE *, TimeStep *);
    /**
     * Assembles the nodal internal forces vector. It assembles the contribution from all elements in
     * particular domain. If runs in parallel mode, the nodal forces for shared nodes are exchanged and
     * updated acordingly.
     */
    void    giveInternalForces(FloatArray &answer, TimeStep *stepN);

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);



    // identification
    const char *giveClassName() const { return "PNlDEIDynamic"; }
    classType giveClassID() const { return PNlDEIDynamicClass; }
    fMode giveFormulation() { return TL; }

    /// Returns first step number
    virtual int        giveNumberOfFirstStep() { return 0; }
    /// Returns time step number, for which initial conditions apply.
    virtual int        giveNumberOfTimeStepWhenIcApply() { return 0; }
protected:
    /**
     * Assembles the load vector.
     * If in parallel mode, the loads of shared/remote nodes are exchanged and remote contributions are taken into account.
     * @param answer load vector
     * @param mode value type mode of load vector
     * @param stepN solution step
     */
    void computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *stepN);
    /**
     * Assembles the diagonal mass matrix of receiver.
     * Local or Global variant of zero mass elements replacement is performed.
     * If runs in parallel, the masses of shared nodes are exchanged and
     * remote contributions are added accordingly.
     * @param mass assembled mass matrix
     * @param maximum estimate of eigen frequency
     * @param solution step
     */
    void    computeMassMtrx(FloatArray &mass, double &maxOm, TimeStep *tStep);

#ifdef __PARALLEL_MODE
    /**
     * Packing function for masses. Pascks mass of shared DofManagers
     * into send communication buffer of given process communicator.
     * @param processComm task communicator for which to pack masses
     * @return nonzero if successfull.
     */
    int packMasses(ProcessCommunicator &processComm);
    /**
     * Unpacking function for masses. Unpacks mass of shared DofManagers
     * from  receive communication buffer of given process communicator.
     * @param processComm task communicator for which to unpack masses
     * @return nonzero if successfull.
     */
    int unpackMasses(ProcessCommunicator &processComm);
    /**
     * Exchanges necessary remote element data with remote partitions. The receiver's nonlocalExt flag must be set.
     * Uses receiver nonlocCommunicator to perform the task using packRemoteElementData and unpackRemoteElementData
     * receiver's services.
     * @return nonzero if success.
     */
    int exchangeRemoteElementData();
public:
    /**
     * Determines the space necessary for send/receive buffer.
     * It uses related communication map pattern to determine the maximum size needed.
     * @param commMap communication map used to send/receive messages
     * @param buff communication buffer
     * @return upper bound of space needed
     */
    int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);
    /**
     * Initializes remote dof managers and their dofs acording to initial and boundary conditions.
     */
    //void initializeRemoteDofs ();
    /**
     * Updates remote dofs velocity, accelaretaions and displacement incerents values.
     */
    //void updateRemoteDofs ();
    /**
     * Updates displacement of remote dof managers dofs.
     */
    //void updateRemoteDofDisplacement ();
    /**
     * Initailizes the list of remote dof managers in current partition.
     * @return nonzero if success
     */
    //int initRemoteDofManList ();
#endif
};
} // end namespace oofem
#endif // pnldeidynamic_h

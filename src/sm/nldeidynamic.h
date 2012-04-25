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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef nldeidynamic_h
#define nldeidynamic_h

#include "structengngmodel.h"
#include "skyline.h"

#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
#define LOCAL_ZERO_MASS_REPLACEMENT 1
#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"

class NlDEIDynamic;
class ProblemCommunicator;

/*
 * Pointer to packing function type.  Functions of this type are used to pack
 * various nodal data into domainCommunicators send buffer. Only nodes in toSend
 * map are considered.
 */
//typedef int (NlDEIDynamic::*NlDEIDynamic_Pack_func) (NlDEIDynamicDomainCommunicator&);
/*
 * Pointer to unpacking function type.  Functions of this type are used to unpack
 * various nodal data from domainCommunicators receive buffer. Only nodes in toRecv
 * map are considered.
 */
//typedef int (NlDEIDynamic::*NlDEIDynamic_Unpack_func) (NlDEIDynamicDomainCommunicator&);
/*
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
 * @note FloatMatrix::Lumped() works only for elements with Linear displacement filed !
 *
 * Current implementation supports parallel processing. Both node- and element cut strategies can
 * be used.
 * - In node cut strategy, partitions are divided using cut, which goes through nodes.
 *   These cut nodes are called "shared" nodes. Generally, unknown values in shared nodes are
 *   composed from local partition contributions as well as from contributions from remote partitions
 *   sharing this node. Particularly, masses and real nodal forces have to be exchanged for shared
 *   nodes.
 * - In element cut strategy, partitions are divided using cut running through elements. The cut elements are
 *   replicated on neighbouring partitions. The nodes belonging to  replicated elements belonging to
 *   remote partitions are called remote nodes. The are mirrors or remote copies of corresponding
 *   nodes on neighbouring partition.
 * - Additional mode has been introduced remote element mode. It introduces the "remote" elements, the
 *   exact local mirrors of remote counterparts. Introduced to support general nonlocal constitutive models,
 *   in order to provide efficient way, how to average local data without need of fine grain communication.
 */
class NlDEIDynamic : public StructuralEngngModel
{
protected:
    /// Mass matrix.
    FloatArray massMatrix;
    /// Load vector.
    FloatArray loadVector;
    /// Vector storing  displacement increments.
    FloatArray previousIncrementOfDisplacementVector;
    /// Displacement, velocity and acceleration vectors.
    FloatArray displacementVector, velocityVector, accelerationVector;
    /// Vector of real nodal forces.
    FloatArray internalForces;
    /// Dumping coefficient (C = dumpingCoef * MassMtrx).
    double dumpingCoef;
    /// Time step.
    double deltaT;
    /// Flag indicating the need for initialization.
    int initFlag;

    // dynamic relaxation specific vars
    /// Flag indicating whether dynamic relaxation takes place.
    int drFlag;
    /// Reference load vector.
    FloatArray loadRefVector;
    /// Parameter determining rate of the loading process.
    double c;
    /// Load level.
    double pt;
    /// End of time interval.
    double Tau;
    /// Estimate of loadRefVector^T*displacementVector(Tau).
    double pyEstimate;
    /// Product of p^tM^(-1)p; where p is reference load vector.
    double pMp;

#ifdef __PARALLEL_MODE
    // public:
    /// Message tags
    enum { InternalForcesExchangeTag, MassExchangeTag, LoadExchangeTag, RemoteElementsExchangeTag };
#endif

public:
    /// Constructor.
    NlDEIDynamic(int i, EngngModel *_master = NULL);
    /// Destructor.
    virtual ~NlDEIDynamic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void terminate(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /**
     * Assembles the nodal internal forces vector. It assembles the contribution from all elements in
     * particular domain. If runs in parallel mode, the nodal forces for shared nodes are exchanged and
     * updated accordingly.
     */
    void giveInternalForces(FloatArray &answer, TimeStep *stepN);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    virtual const char *giveClassName() const { return "NlDEIDynamic"; }
    virtual classType giveClassID() const { return NlDEIDynamicClass; }
    virtual fMode giveFormulation() { return TL; }

    virtual int giveNumberOfFirstStep() { return 0; }
    virtual int giveNumberOfTimeStepWhenIcApply() { return 0; }

protected:
    /**
     * Assembles the load vector.
     * If in parallel mode, the loads of shared/remote nodes are exchanged and remote contributions are taken into account.
     * @param answer Load vector.
     * @param mode Value type mode of load vector.
     * @param stepN Solution step.
     */
    void computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *stepN);
    /**
     * Assembles the diagonal mass matrix of receiver.
     * Local or Global variant of zero mass elements replacement is performed.
     * If runs in parallel, the masses of shared nodes are exchanged and
     * remote contributions are added accordingly.
     * @param mass Assembled mass matrix.
     * @param maximum Estimate of eigenfrequency.
     * @param solution Step.
     */
    void computeMassMtrx(FloatArray &mass, double &maxOm, TimeStep *tStep);

#ifdef __PARALLEL_MODE
    /**
     * Packing function for masses. Pascks mass of shared DofManagers
     * into send communication buffer of given process communicator.
     * @param processComm Task communicator for which to pack masses.
     * @return Nonzero if successful.
     */
    int packMasses(ProcessCommunicator &processComm);
    /**
     * Unpacking function for masses. Unpacks mass of shared DofManagers
     * from  receive communication buffer of given process communicator.
     * @param processComm Task communicator for which to unpack masses.
     * @return Nonzero if successful.
     */
    int unpackMasses(ProcessCommunicator &processComm);
    /**
     * Exchanges necessary remote element data with remote partitions. The receiver's nonlocalExt flag must be set.
     * Uses receiver nonlocCommunicator to perform the task using packRemoteElementData and unpackRemoteElementData
     * receiver's services.
     * @return Nonzero if success.
     */
    int exchangeRemoteElementData();

public:
    /**
     * Determines the space necessary for send/receive buffer.
     * It uses related communication map pattern to determine the maximum size needed.
     * @param commMap Communication map used to send/receive messages.
     * @param buff Communication buffer.
     * @return Upper bound of space needed.
     */
    int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);
    /*
     * Initializes remote dof managers and their dofs acording to initial and boundary conditions.
     */
    //void initializeRemoteDofs ();
    /*
     * Updates remote dofs velocity, accelaretaions and displacement incerents values.
     */
    //void updateRemoteDofs ();
    /*
     * Updates displacement of remote dof managers dofs.
     */
    //void updateRemoteDofDisplacement ();
    /*
     * Initailizes the list of remote dof managers in current partition.
     * @return nonzero if success
     */
    //int initRemoteDofManList ();
#endif
};
} // end namespace oofem
#endif // nldeidynamic_h

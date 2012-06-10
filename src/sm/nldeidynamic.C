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

#include "nldeidynamic.h"
#include "timestep.h"
#include "dofmanager.h"
#include "element.h"
#include "dof.h"
#include "verbose.h"
#include "outputmanager.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
#endif

namespace oofem {
#define ZERO_REL_MASS  1.E-6

NlDEIDynamic ::  NlDEIDynamic(int i, EngngModel *_master) : StructuralEngngModel(i, _master), massMatrix(), loadVector(),
    previousIncrementOfDisplacementVector(), displacementVector(),
    velocityVector(), accelerationVector(), internalForces()
{
#ifdef __PARALLEL_MODE
    commMode = ProblemCommMode__UNKNOWN_MODE;
    nonlocalExt = 0;
    communicator = nonlocCommunicator = NULL;
    commBuff = NULL;
#endif
    ndomains = 1;
    initFlag = 1;
}


NlDEIDynamic :: ~NlDEIDynamic()
{
}

NumericalMethod *NlDEIDynamic :: giveNumericalMethod(MetaStep *mStep)
// Only one has reason for NlDEIDynamic
//     - SolutionOfLinearEquations

{
    return NULL;  // Not necessary here - Diagonal matrix and simple inversion is used.
}

IRResultType
NlDEIDynamic :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    StructuralEngngModel :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, dumpingCoef, IFT_NlDEIDynamic_dumpcoef, "dumpcoef"); // C = dumpingCoef * M // Macro
    IR_GIVE_FIELD(ir, deltaT, IFT_NlDEIDynamic_deltat, "deltat"); // Macro

    drFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, drFlag, IFT_NlDEIDynamic_drflag, "drflag"); // Macro
    if ( drFlag ) {
        IR_GIVE_FIELD(ir, Tau, IFT_NlDEIDynamic_tau, "tau");
        IR_GIVE_FIELD(ir, pyEstimate, IFT_NlDEIDynamic_py, "py");
    }

#ifdef __PARALLEL_MODE
    if ( ir->hasField(IFT_NlDEIDynamic_nodecutmode, "nodecutmode") ) {
        commMode = ProblemCommMode__NODE_CUT;
    } else if ( ir->hasField(IFT_NlDEIDynamic_elementcutmode, "elementcutmode") ) {
        commMode = ProblemCommMode__ELEMENT_CUT;
    } else {
        _error("instanciateFrom: NlDEIDynamicCommunicatorMode not specified");
    }

    commBuff = new CommunicatorBuff( this->giveNumberOfProcesses() );
    communicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                           this->giveNumberOfProcesses(),
                                           this->commMode);

    if ( ir->hasField(IFT_NlDEIDynamic_nonlocalext, "nonlocalext") ) {
        nonlocalExt = 1;
        nonlocCommunicator = new ProblemCommunicator(this, commBuff, this->giveRank(),
                                                     this->giveNumberOfProcesses(),
                                                     ProblemCommMode__REMOTE_ELEMENT_MODE);
    }

#endif

    return IRRT_OK;
}





double NlDEIDynamic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                              TimeStep *tStep, Domain *d, Dof *dof)
// Returns unknown quantity like displacement, velocity of equation eq.
// This function translates this request to numerical method language.
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }


    if ( chc != EID_MomentumBalance ) {
        _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
        return 0.;
    }

    switch ( mode ) {
    case VM_Total:
        return displacementVector.at(eq);

    case VM_Incremental:
        return previousIncrementOfDisplacementVector.at(eq);

    case VM_Velocity:
        return velocityVector.at(eq);

    case VM_Acceleration:
        return accelerationVector.at(eq);

    default:
        _error("giveUnknownComponent: Unknown is of undefined type for this problem");
    }

    return 0.;
}

TimeStep *NlDEIDynamic :: giveNextStep()
{
    int istep = 0;
    double totalTime = 0;
    StateCounterType counter = 1;

    if (previousStep != NULL){
        delete previousStep;
    }

    if ( currentStep != NULL ) {
        totalTime = currentStep->giveTargetTime() + deltaT;
        istep     = currentStep->giveNumber() + 1;
        counter   = currentStep->giveSolutionStateCounter() + 1;
    }

    previousStep = currentStep;
    currentStep  = new TimeStep(istep, this, 1, totalTime, deltaT, counter);

    return currentStep;
}



void NlDEIDynamic :: solveYourself()
{
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
    // Force equation numbering before setting up comm maps.
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    OOFEM_LOG_INFO("[process rank %d] neq is %d\n", this->giveRank(), neq);
 #endif

    // Set up communication patterns,
    communicator->setUpCommunicationMaps(this, true);
    if ( nonlocalExt ) {
        nonlocCommunicator->setUpCommunicationMaps(this, true);
    }
#endif

    StructuralEngngModel :: solveYourself();
}

void NlDEIDynamic :: solveYourselfAt(TimeStep *tStep)
{
    //
    // Creates system of governing eq's and solves them at given time step.
    //

    Domain *domain = this->giveDomain(1);
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    int nman  = domain->giveNumberOfDofManagers();

    DofManager *node;
    Dof *iDof;

    int nDofs;
    int i, k, j, jj;
    double coeff, maxDt, maxOm = 0.;
    double prevIncrOfDisplacement, incrOfDisplacement;

    if ( initFlag ) {
#ifdef VERBOSE
        OOFEM_LOG_DEBUG("Assembling mass matrix\n");
#endif

        //
        // Assemble mass matrix.
        //
        this->computeMassMtrx(massMatrix, maxOm, tStep);

        if ( drFlag ) {
            // If dynamic relaxation: Assemble amplitude load vector.
            loadRefVector.resize(neq);
            loadRefVector.zero();

            this->computeLoadVector(loadRefVector, VM_Total, tStep);

#ifdef __PARALLEL_MODE
            // Compute the processor part of load vector norm pMp
            this->pMp = 0.0;
            double my_pMp = 0.0, coeff = 1.0;
            int eqNum, ndofs, ndofman = domain->giveNumberOfDofManagers();
            dofManagerParallelMode dofmanmode;
            DofManager *dman;
            Dof *jdof;
            for ( int dm = 1; dm <= ndofman; dm++ ) {
                dman = domain->giveDofManager(dm);
                ndofs = dman->giveNumberOfDofs();
                dofmanmode = dman->giveParallelMode();

                // Skip all remote and null dofmanagers
                coeff = 1.0;
                if ( ( dofmanmode == DofManager_remote ) || ( ( dofmanmode == DofManager_null ) ) ) {
                    continue;
                } else if ( dofmanmode == DofManager_shared ) {
                    coeff = 1. / dman->givePartitionsConnectivitySize();
                }

                // For shared nodes we add locally an average = 1/givePartitionsConnectivitySize()*contribution,
                for ( j = 1; j <= ndofs; j++ ) {
                    jdof = dman->giveDof(j);
                    if ( jdof->isPrimaryDof() && ( eqNum = jdof->__giveEquationNumber() ) ) {
                        my_pMp += coeff * loadRefVector.at(eqNum) * loadRefVector.at(eqNum) / massMatrix.at(eqNum);
                    }
                }
            }

            // Sum up the contributions from processors.
            MPI_Allreduce(& my_pMp, & pMp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
            this->pMp = 0.0;
            for ( i = 1; i <= neq; i++ ) {
                pMp += loadRefVector.at(i) * loadRefVector.at(i) / massMatrix.at(i);
            }
#endif
            // Solve for rate of loading process (parameter "c") (undamped system assumed),
            if ( dumpingCoef < 1.e-3 ) {
                c = 3.0 * this->pyEstimate / pMp / Tau / Tau;
            } else {
                c = this->pyEstimate * Tau * dumpingCoef * dumpingCoef * dumpingCoef / pMp /
                    ( -3.0 / 2.0 + dumpingCoef * Tau + 2.0 * exp(-dumpingCoef * Tau) - 0.5 * exp(-2.0 * dumpingCoef * Tau) );
            }
        }

        initFlag = 0;
    }


    if ( tStep->giveNumber() == giveNumberOfFirstStep() ) {
        //
        // Special init step - Compute displacements at tstep 0.
        //
        displacementVector.resize(neq);
        displacementVector.zero();
        previousIncrementOfDisplacementVector.resize(neq);
        previousIncrementOfDisplacementVector.zero();
        velocityVector.resize(neq);
        velocityVector.zero();
        accelerationVector.resize(neq);
        accelerationVector.zero();

        for ( j = 1; j <= nman; j++ ) {
            node = domain->giveDofManager(j);
            nDofs = node->giveNumberOfDofs();

            for ( k = 1; k <= nDofs; k++ ) {
                // Ask for initial values obtained from
                // bc (boundary conditions) and ic (initial conditions)
                // all dofs are expected to be  DisplacementVector type.
                iDof  =  node->giveDof(k);
                if ( !iDof->isPrimaryDof() ) {
                    continue;
                }

                jj = iDof->__giveEquationNumber();
                if ( jj ) {
                    displacementVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
                    velocityVector.at(jj)     = iDof->giveUnknown(EID_MomentumBalance, VM_Velocity, tStep);
                    accelerationVector.at(jj)    = iDof->giveUnknown(EID_MomentumBalance, VM_Acceleration, tStep) ;
                }
            }
        }

        //
        // Set-up numerical model.
        //

        // Try to determine the best deltaT,
        maxDt = 2.0 / sqrt(maxOm);
        if ( deltaT > maxDt ) {
            // Print reduced time step increment and minimum period Tmin
            OOFEM_LOG_RELEVANT("deltaT reduced to %e, Tmin is %e\n", maxDt, maxDt * M_PI);
            deltaT = maxDt;
            tStep->setTimeIncrement(deltaT);
        }

        for ( j = 1; j <= neq; j++ ) {
            previousIncrementOfDisplacementVector.at(j) =  velocityVector.at(j) * ( deltaT );
            displacementVector.at(j) -= previousIncrementOfDisplacementVector.at(j);
        }
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT( "\n\nSolving [Step number %8d, Time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif
        return;
    } // end of init step

#ifdef VERBOSE
    OOFEM_LOG_DEBUG("Assembling right hand side\n");
#endif

    for ( i = 1; i <= neq; i++ ) {
        displacementVector.at(i) += previousIncrementOfDisplacementVector.at(i);
    }

    // Update solution state counter
    tStep->incrementStateCounter();

    // Compute internal forces.
    this->giveInternalForces( internalForces, false, 1, tStep );

    if ( !drFlag ) {
        //
        // Assembling the element part of load vector.
        //
        this->computeLoadVector(loadVector, VM_Total, tStep);
        //
        // Assembling additional parts of right hand side.
        //
        for ( k = 1; k <= neq; k++ ) {
            loadVector.at(k) -= internalForces.at(k);
        }
    } else {
        // Dynamic relaxation
        // compute load factor
        pt = 0.0;

#ifdef __PARALLEL_MODE
        double my_pt = 0.0, coeff = 1.0;
        int eqNum, ndofs, ndofman = domain->giveNumberOfDofManagers();
        dofManagerParallelMode dofmanmode;
        DofManager *dman;
        Dof *jdof;
        for ( int dm = 1; dm <= ndofman; dm++ ) {
            dman = domain->giveDofManager(dm);
            ndofs = dman->giveNumberOfDofs();
            dofmanmode = dman->giveParallelMode();
            // skip all remote and null dofmanagers
            coeff = 1.0;
            if ( ( dofmanmode == DofManager_remote ) || ( dofmanmode == DofManager_null ) ) {
                continue;
            } else if ( dofmanmode == DofManager_shared ) {
                coeff = 1. / dman->givePartitionsConnectivitySize();
            }

            // For shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution.
            for ( j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( eqNum = jdof->__giveEquationNumber() ) ) {
                    my_pt += coeff * internalForces.at(eqNum) * loadRefVector.at(eqNum) / massMatrix.at(eqNum);
                }
            }
        }

        // Sum up the contributions from processors.
        MPI_Allreduce(& my_pt, & pt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        for ( k = 1; k <= neq; k++ ) {
            pt += internalForces.at(k) * loadRefVector.at(k) / massMatrix.at(k);
        }

#endif
        pt = pt / pMp;
        if ( dumpingCoef < 1.e-3 ) {
            pt += c * ( Tau - tStep->giveTargetTime() ) / Tau;
        } else {
            pt += c * ( 1.0 - exp( dumpingCoef * ( tStep->giveTargetTime() - Tau ) ) ) / dumpingCoef / Tau;
        }

        loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
        for ( k = 1; k <= neq; k++ ) {
            loadVector.at(k) = pt * loadRefVector.at(k) - internalForces.at(k);
        }


        // Compute relative error.
        double err = 0.0;
#ifdef __PARALLEL_MODE
        double my_err = 0.0;

        for ( int dm = 1; dm <= ndofman; dm++ ) {
            dman = domain->giveDofManager(dm);
            ndofs = dman->giveNumberOfDofs();
            dofmanmode = dman->giveParallelMode();
            // Skip all remote and null dofmanagers.
            coeff = 1.0;
            if ( ( dofmanmode == DofManager_remote ) || ( dofmanmode == DofManager_null ) ) {
                continue;
            } else if ( dofmanmode == DofManager_shared ) {
                coeff = 1. / dman->givePartitionsConnectivitySize();
            }

            // For shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution.
            for ( j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( eqNum = jdof->__giveEquationNumber() ) ) {
                    my_err += coeff * loadVector.at(eqNum) * loadVector.at(eqNum) / massMatrix.at(eqNum);
                }
            }
        }

        // Sum up the contributions from processors.
        MPI_Allreduce(& my_err, & err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
        for ( k = 1; k <= neq; k++ ) {
            err = loadVector.at(k) * loadVector.at(k) / massMatrix.at(k);
        }

#endif
        err = err / ( pMp * pt * pt );
        OOFEM_LOG_RELEVANT("Relative error is %e, loadlevel is %e\n", err, pt);
    }

    for ( j = 1; j <= neq; j++ ) {
        coeff =  massMatrix.at(j);
        loadVector.at(j) +=
            coeff * ( ( 1. / ( deltaT * deltaT ) ) - dumpingCoef * 1. / ( 2. * deltaT ) ) *
            previousIncrementOfDisplacementVector.at(j);
    }

    //
    // Set-up numerical model
    //
    /* it is not necesary to call numerical method
     * approach used here is not good, but effective enough
     * inverse of diagonal mass matrix is done here
     */
    //
    // call numerical model to solve arised problem - done localy here
    //
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "\n\nSolving [Step number %8d, Time %15e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    for ( i = 1; i <= neq; i++ ) {
        prevIncrOfDisplacement = previousIncrementOfDisplacementVector.at(i);
        incrOfDisplacement = loadVector.at(i) /
                             ( massMatrix.at(i) * ( 1. / ( deltaT * deltaT ) + dumpingCoef / ( 2. * deltaT ) ) );
        accelerationVector.at(i) = ( incrOfDisplacement - prevIncrOfDisplacement ) / ( deltaT * deltaT );
        velocityVector.at(i)     = ( incrOfDisplacement + prevIncrOfDisplacement ) / ( 2. * deltaT );
        previousIncrementOfDisplacementVector.at(i) = incrOfDisplacement;
    }
}


void NlDEIDynamic :: updateYourself(TimeStep *stepN)
{
    // Updates internal state to reached one
    // all internal variables are directly updated by
    // numerical method - void function here

    // this->updateInternalState(stepN);
    StructuralEngngModel :: updateYourself(stepN);
}


void
NlDEIDynamic :: computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *stepN)
{
    answer.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    answer.zero();

    //
    // Assemble the nodal part of load vector.
    //
    this->assembleVector( answer, stepN, EID_MomentumBalance, ExternalForcesVector, mode,
                          EModelDefaultEquationNumbering(), this->giveDomain(1));

    //
    // Exchange contributions.
    //
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers( answer, LoadExchangeTag );
#endif
}


void
NlDEIDynamic :: computeMassMtrx(FloatArray &massMatrix, double &maxOm, TimeStep *tStep)
{
    Domain *domain = this->giveDomain(1);
    int nelem = domain->giveNumberOfElements();
    int neq = this->giveNumberOfEquations(EID_MomentumBalance);
    int i, j, jj, n;
    double maxOmi, maxOmEl;
    FloatMatrix charMtrx, charMtrx2;
    IntArray loc;
    Element *element;
    EModelDefaultEquationNumbering en;
#ifdef __PARALLEL_MODE
    int result;
#endif

#ifndef LOCAL_ZERO_MASS_REPLACEMENT
    FloatArray diagonalStiffMtrx;
#endif

    maxOm = 0.;
    massMatrix.resize(neq);
    massMatrix.zero();
    for ( i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);

#ifdef __PARALLEL_MODE
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }
#endif

        element->giveLocationArray(loc, EID_MomentumBalance, en);
        element->giveCharacteristicMatrix(charMtrx, LumpedMassMatrix, tStep);

#ifdef LOCAL_ZERO_MASS_REPLACEMENT
        element->giveCharacteristicMatrix(charMtrx2, StiffnessMatrix, tStep);
#endif

#ifdef DEBUG
        if ( ( n = loc.giveSize() ) != charMtrx.giveNumberOfRows() ) {
            _error("solveYourselfAt : dimension mismatch");
        }
#endif

        n = loc.giveSize();

#ifdef LOCAL_ZERO_MASS_REPLACEMENT
        maxOmEl = 0.;

        double maxElmass = -1.0;
        for ( j = 1; j <= n; j++ ) {
            maxElmass = max( maxElmass, charMtrx.at(j, j) );
        }

        if ( maxElmass <= 0.0 ) {
            _warning2("solveYourselfAt: Element (%d) with zero (or negative) lumped mass encountered\n", i);
        }

        for ( j = 1; j <= n; j++ ) {
            if ( charMtrx.at(j, j) > maxElmass * ZERO_REL_MASS ) {
                maxOmi =  charMtrx2.at(j, j) / charMtrx.at(j, j);
                maxOmEl = ( maxOmEl > maxOmi ) ? ( maxOmEl ) : ( maxOmi );
            }
        }

        maxOm = ( maxOm > maxOmEl ) ? ( maxOm ) : ( maxOmEl );

        for ( j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( ( jj ) && ( charMtrx.at(j, j) <= maxElmass * ZERO_REL_MASS ) ) {
                charMtrx.at(j, j) = charMtrx2.at(j, j) / maxOmEl;
            }
        }
#endif

        for ( j = 1; j <= n; j++ ) {
            jj = loc.at(j);
            if ( jj ) {
                massMatrix.at(jj) += charMtrx.at(j, j);
            }
        }
    }

#ifndef LOCAL_ZERO_MASS_REPLACEMENT
    // If init step - find minimun period of vibration in order to
    // determine maximal admisible time step
    // global variant
    for (i=1; i<=nelem; i++)
    {
        element = domain -> giveElement(i);
        element->giveLocationArray(loc, EID_MomentumBalance, en);
        element->giveCharacteristicMatrix( charMtrx, StiffnessMatrix, tStep );
        n = loc.giveSize() ;
        for (j=1; j<=n; j++) {
            jj = loc.at(j);
            if (jj) {
                diagonalStiffMtrx.at(jj) += charMtrx.at(j,j);
            }
        }
    }

    // Find find global minimun period of vibration
    double maxElmass = -1.0;
    for (j=1 ; j<=n; j++) {
        maxElmass = max(maxElmass,charMtrx.at(j,j));
    }

    if ( maxElmass <= 0.0 ) {
        _error("solveYourselfAt: Element with zero (or negative) lumped mass encountered\n");
    }

    for (j=1; j<= neq; j++) {
        if (massMatrix.at(j) > maxElmass * ZERO_REL_MASS ) {
            maxOmi =  diagonalStiffMtrx.at(j) / massMatrix.at(j);
            maxOm  = (maxOm > maxOmi) ? (maxOm) : (maxOmi);
        }
    }

    // Set ZERO MASS members in massMatrix to value which corresponds to global maxOm.
    for (i=1; i<= neq; i++) {
        if (massMatrix.at(i) <= maxElmass*ZERO_REL_MASS) {
            massMatrix.at(i) = diagonalStiffMtrx.at(i) / maxOm;
        }
    }
#endif

#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers( massMatrix, MassExchangeTag );

    // Determine maxOm over all processes.
 #ifdef __USE_MPI
    double globalMaxOm;

  #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT( "NlDEIDynamic :: computeMassMtrx", "Reduce of maxOm started", this->giveRank() );
  #endif

    result = MPI_Allreduce(& maxOm, & globalMaxOm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  #ifdef __VERBOSE_PARALLEL
    VERBOSEPARALLEL_PRINT( "NlDEIDynamic :: computeMassMtrx", "Reduce of maxOm finished", this->giveRank() );
  #endif

    if ( result != MPI_SUCCESS ) {
        _error("setUpCommunicationMaps: MPI_Allreduce failed");
    }

    maxOm = globalMaxOm;
 #else
    WARNING: NOT SUPPORTED MESSAGE PARSING LIBRARY
 #endif

#endif
}

#ifdef __PARALLEL_MODE
int
NlDEIDynamic :: estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType)
{
    int mapSize = commMap.giveSize();
    int i, j, ndofs, count = 0, pcount = 0;
    IntArray locationArray;
    Domain *domain = this->giveDomain(1);
    DofManager *dman;
    Dof *jdof;

    if ( packUnpackType == ProblemCommMode__ELEMENT_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveDofManager( commMap.at(i) )->giveNumberOfDofs();
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * count );
    } else if ( packUnpackType == ProblemCommMode__NODE_CUT ) {
        for ( i = 1; i <= mapSize; i++ ) {
            ndofs = ( dman = domain->giveDofManager( commMap.at(i) ) )->giveNumberOfDofs();
            for ( j = 1; j <= ndofs; j++ ) {
                jdof = dman->giveDof(j);
                if ( jdof->isPrimaryDof() && ( jdof->__giveEquationNumber() ) ) {
                    count++;
                } else {
                    pcount++;
                }
            }
        }

        return ( buff.givePackSize(MPI_DOUBLE, 1) * max(count, pcount) );
    } else if ( packUnpackType == ProblemCommMode__REMOTE_ELEMENT_MODE ) {
        for ( i = 1; i <= mapSize; i++ ) {
            count += domain->giveElement( commMap.at(i) )->estimatePackSize(buff);
        }

        return count;
    }

    return 0;
}

#endif



contextIOResultType NlDEIDynamic :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int closeFlag = 0;
    FILE *file;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = StructuralEngngModel :: saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // Ensure consistent records

    return CIO_OK;
}



contextIOResultType NlDEIDynamic :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int closeFlag = 0;
    int istep, iversion;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // Save element context.
    if ( ( iores = StructuralEngngModel :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = previousIncrementOfDisplacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = displacementVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = velocityVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = accelerationVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& deltaT, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    } // ensure consistent records

    return CIO_OK;
}



void
NlDEIDynamic :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime)
{
    static char dofchar[] = "dva";
    static ValueModeType dofmodes[] = {
        VM_Total, VM_Velocity, VM_Acceleration
    };

    iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}

void
NlDEIDynamic :: terminate(TimeStep *tStep)
{
    StructuralEngngModel :: terminate(tStep);
    this->printReactionForces(tStep, 1);
}


void
NlDEIDynamic :: printOutputAt(FILE *File, TimeStep *stepN)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(stepN) ) {
        return;                                                                      // do not print even Solution step header
    }

    fprintf( File, "\n\nOutput for time % .3e, solution step number %d\n", stepN->giveTargetTime(), stepN->giveNumber() );
    if ( drFlag ) {
        fprintf(File, "Reached load level : %e\n\n", this->pt);
    }

    this->giveDomain(1)->giveOutputManager()->doDofManOutput(File, stepN);
    this->giveDomain(1)->giveOutputManager()->doElementOutput(File, stepN);
}
} // end namespace oofem

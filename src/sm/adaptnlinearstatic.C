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

#include "mathfem.h"
#include "adaptnlinearstatic.h"
#include "verbose.h"
#include "clock.h"
#include "metastep.h"
#include "timestep.h"
#include "nummet.h"
#include "element.h"
#include "node.h"
#include "domain.h"
#include "datareader.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "dof.h"
#include "eleminterpunknownmapper.h"
#include "errorestimator.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "contextioerr.h"
#include "oofem_terminate.h"
#include "util.h"

#ifdef __PETSC_MODULE
 #include "petsccontext.h"
#endif

namespace oofem {
AdaptiveNonLinearStatic :: AdaptiveNonLinearStatic(int i, EngngModel *_master) : NonLinearStatic(i, _master),
    d2_totalDisplacement(), d2_incrementOfDisplacement(), timeStepLoadLevels()
{
    meshPackage = MPT_T3D;
    equilibrateMappedConfigurationFlag = 0;

#ifdef __PARALLEL_MODE
    this->preMappingLoadBalancingFlag = false;
#endif
}


AdaptiveNonLinearStatic :: ~AdaptiveNonLinearStatic()
{
}


IRResultType
AdaptiveNonLinearStatic :: initializeFrom(InputRecord *ir)
// input from inputString
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int _val;

    NonLinearStatic :: initializeFrom(ir);
    int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, IFT_AdaptiveNonLinearStatic_meshpackage, "meshpackage"); // Macro
    meshPackage = ( MeshPackageType ) meshPackageId;

    equilibrateMappedConfigurationFlag =  0;
    IR_GIVE_OPTIONAL_FIELD(ir, equilibrateMappedConfigurationFlag, IFT_AdaptiveNonLinearStatic_equilmc, "equilmc"); // Macro
    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_AdaptiveNonLinearStatic_preMappingLoadBalancingFlag, "premaplbflag"); // Macro
    preMappingLoadBalancingFlag = _val;

    return IRRT_OK;
}

void
AdaptiveNonLinearStatic :: solveYourselfAt(TimeStep *tStep) {
    proceedStep(1, tStep);
    this->updateYourself(tStep);

#ifdef __OOFEG
    ESIEventLoop( YES, oofem_tmpstr("AdaptiveNonLinearStatic: Solution finished; Press Ctrl-p to continue") );
#endif

    this->terminate( this->giveCurrentStep() );

#ifdef __PARALLEL_MODE
    if ( preMappingLoadBalancingFlag ) {
        this->balanceLoad( this->giveCurrentStep() );
    }

#endif

    // evaluate error of the reached solution
    this->defaultErrEstimator->estimateError( equilibratedEM, this->giveCurrentStep() );
    //this->defaultErrEstimator->estimateError( temporaryEM, this->giveCurrentStep() );
    this->defaultErrEstimator->giveRemeshingCrit()->estimateMeshDensities( this->giveCurrentStep() );
    RemeshingStrategy strategy = this->defaultErrEstimator->giveRemeshingCrit()->giveRemeshingStrategy( this->giveCurrentStep() );

    // if ((strategy == RemeshingFromCurrentState_RS) && (this->giveDomain(1)->giveSerialNumber() == 0))
    //  strategy = RemeshingFromPreviousState_RS;

    if ( strategy == NoRemeshing_RS ) {
        //
    } else if ( ( strategy == RemeshingFromCurrentState_RS ) || ( strategy == RemeshingFromPreviousState_RS ) ) {
        // do remeshing
        MesherInterface *mesher = CreateUsrDefMesherInterface( meshPackage, this->giveDomain(1) );

        Domain *newDomain;
        MesherInterface :: returnCode result = mesher->createMesh(this->giveCurrentStep(), 1,
                                                                  this->giveDomain(1)->giveSerialNumber() + 1, & newDomain);

        delete mesher;

        if ( result == MesherInterface :: MI_OK ) {
            this->initFlag = 1;
            this->adaptiveRemap(newDomain);
        } else if ( result == MesherInterface :: MI_NEEDS_EXTERNAL_ACTION ) {
            if ( strategy == RemeshingFromCurrentState_RS ) {
                // ensure the updating the step
                this->setContextOutputMode(COM_Always);
                //this->terminate (this->giveCurrentStep());
            } else {
                // save previous step (because update not called)
            }

            this->terminateAnalysis();
            throw OOFEM_Terminate();
        } else {
            _error("solveYourselfAt: MesherInterface::createMesh failed");
        }
    }
}

void
AdaptiveNonLinearStatic :: updateYourself(TimeStep *atTime)
{
    if ( timeStepLoadLevels.isEmpty() ) {
        timeStepLoadLevels.resize( this->giveNumberOfSteps() );
    }

    // in case of adaptive restart from given timestep
    // a time step with same number and incremented version will be generated.
    // Then load level reached on new discretization will overwrite the old one obtained for
    // old discretization for this step. But this is consistent, since when initialLoadVector
    // is requested to be recovered the reference load vectors are assembled
    // on actual discretization.
    timeStepLoadLevels.at( atTime->giveNumber() ) = loadLevel;

    NonLinearStatic :: updateYourself(atTime);
}





double AdaptiveNonLinearStatic ::  giveUnknownComponent(EquationID chc, ValueModeType mode,
                                                        TimeStep *tStep, Domain *d, Dof *dof)
// returns unknown quantity like displacement, velocity of equation eq
// This function translates this request to numerical method language
{
    int eq = dof->__giveEquationNumber();
    if ( eq == 0 ) {
        _error("giveUnknownComponent: invalid equation number");
    }

    if ( tStep != this->giveCurrentStep() ) {
        _error("giveUnknownComponent: unknown time step encountered");
        return 0.;
    }

    if ( d->giveNumber() == 2 ) {
        if ( chc == EID_MomentumBalance ) {
            switch ( mode ) {
            case VM_Incremental:
                if ( d2_incrementOfDisplacement.isNotEmpty() ) {
                    return d2_incrementOfDisplacement.at(eq);
                } else {
                    return 0.;
                }

            case VM_Total:
                if ( d2_totalDisplacement.isNotEmpty() ) {
                    return d2_totalDisplacement.at(eq);
                } else {
                    return 0.;
                }

            default:
                _error("giveUnknownComponent: Unknown is of undefined ValueModeType for this problem");
            }
        } else {
            _error("giveUnknownComponent: Unknown is of undefined CharType for this problem");
            return 0.;
        }
    } else {
        return NonLinearStatic :: giveUnknownComponent(chc, mode, tStep, d, dof);
    }

    return 0.0;
}



int
AdaptiveNonLinearStatic :: initializeAdaptiveFrom(EngngModel *sourceProblem)
{
    int ielem, nelem, result = 1;

    // measure time consumed by mapping
    //clock_t sc = this->getClock();
    //clock_t mc, mc2, ec ;
    oofem_timeval st, st1, st2, mc1, mc2, mc3, ec;
    getUtime(st);

    if ( sourceProblem->giveClassID() != AdaptiveNonLinearStaticClass ) {
        _error("sory");
    }

    this->currentStep = new TimeStep( * ( sourceProblem->giveCurrentStep() ) );
    if ( sourceProblem->givePreviousStep() ) {
        this->previousStep = new TimeStep( * ( sourceProblem->givePreviousStep() ) );
    }

    // map primary unknowns
    EIPrimaryUnknownMapper mapper;

    totalDisplacement.resize( this->giveNumberOfDomainEquations(1, EID_MomentumBalance) );
    incrementOfDisplacement.resize( this->giveNumberOfDomainEquations(1, EID_MomentumBalance) );
    totalDisplacement.zero();
    incrementOfDisplacement.zero();

    result &= mapper.mapAndUpdate( totalDisplacement, VM_Total, EID_MomentumBalance,
                                  sourceProblem->giveDomain(1), this->giveDomain(1), sourceProblem->giveCurrentStep() );

    result &= mapper.mapAndUpdate( incrementOfDisplacement, VM_Incremental, EID_MomentumBalance,
                                  sourceProblem->giveDomain(1), this->giveDomain(1), sourceProblem->giveCurrentStep() );

    //mc = this->getClock();
    getRelativeUtime(mc1, st);
    getUtime(st1);

    // map internal ip state
    nelem = this->giveDomain(1)->giveNumberOfElements();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        result &= this->giveDomain(1)->giveElement(ielem)->adaptiveMap( sourceProblem->giveDomain(1),
                                                                       sourceProblem->giveCurrentStep() );
    }

    //mc2 = this->getClock();
    getRelativeUtime(mc2, st1);
    getUtime(st2);

    // computes the stresses and calls updateYourself to mapped state
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        result &= this->giveDomain(1)->giveElement(ielem)->adaptiveUpdate(currentStep);
    }

    // finish mapping process
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        result &= this->giveDomain(1)->giveElement(ielem)->adaptiveFinish(currentStep);
    }


    // increment time step if mapped state will be considered as new solution stepL
    /*
     * this->giveNextStep();
     * if (equilibrateMappedConfigurationFlag) {
     * // we need to  assemble the load vector in same time as the restarted step,
     * // so new time step is generated with same intrincic time as has the
     * // previous step if equilibrateMappedConfigurationFlag is set.
     * // this allows to equlibrate the previously reached state
     * TimeStep* cts = this->giveCurrentStep();
     * cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
     * cts = this->givePreviousStep();
     * cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
     * }
     *
     *
     * if (this->giveCurrentStep()->giveNumber() ==
     * this->giveMetaStep(this->giveCurrentStep()->giveMetaStepNumber())->giveFirstStepNumber()) {
     * this->updateAttributes (this->giveCurrentStep());
     * }
     */
    this->updateAttributes( this->giveCurrentMetaStep() );

    // increment solution state counter - not needed, IPs are updated by adaptiveUpdate previously called
    // and there is no change in primary vars.
    // this->giveCurrentStep()->incrementStateCounter();

    // assemble new initial load for new discretization
    this->assembleInitialLoadVector( initialLoadVector, initialLoadVectorOfPrescribed,
                                    ( AdaptiveNonLinearStatic * ) sourceProblem, 1, this->giveCurrentStep() );
    // assemble new total load for new discretization
    // this->assembleCurrentTotalLoadVector (totalLoadVector, totalLoadVectorOfPrescribed, this->giveCurrentStep());
    // set bcloadVector to zero (no increment within same step)

    //ec = this->getClock();
    getRelativeUtime(mc3, st2);
    getRelativeUtime(ec, st);

    // compute processor time used by the program
    OOFEM_LOG_INFO( "user time consumed by primary mapping: %.2fs\n", ( double ) ( mc1.tv_sec + mc1.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by ip mapping:      %.2fs\n", ( double ) ( mc2.tv_sec + mc2.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by ip update:       %.2fs\n", ( double ) ( mc3.tv_sec + mc3.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by mapping:         %.2fs\n", ( double ) ( ec.tv_sec + ec.tv_usec / ( double ) OOFEM_USEC_LIM ) );

    //

    //
    // bring mapped configuration into equilibrium
    //
    if ( equilibrateMappedConfigurationFlag ) {
        // use secant stiffness to resrore equilibrium
        NonLinearStatic_stifnessMode oldStiffMode = this->stiffMode;
        stiffMode = nls_secantStiffness;



        if ( initFlag ) {
            stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
            if ( stiffnessMatrix == NULL ) {
                _error("proceedStep: sparse matrix creation failed");
            }

            if ( nonlocalStiffnessFlag ) {
                if ( !stiffnessMatrix->isAsymmetric() ) {
                    _error("proceedStep: stiffnessMatrix does not support asymmetric storage");
                }
            }

            stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble( stiffnessMatrix, this->giveCurrentStep(), EID_MomentumBalance, SecantStiffnessMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            initFlag = 0;
        }

        // updateYourself() not necessary - the adaptiveUpdate previously called does the job
        //this->updateYourself(this->giveCurrentStep());
#ifdef VERBOSE
        OOFEM_LOG_INFO( "Equilibrating mapped configuration [step number %5d.%d]\n",
                       this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //double deltaL = nMethod->giveUnknownComponent (StepLength, 0);
        double deltaL = nMethod->giveCurrentStepLength();
        this->assembleIncrementalReferenceLoadVectors( incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                      refLoadInputMode, this->giveDomain(1), EID_MomentumBalance, this->giveCurrentStep() );
        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT( "Solving [step number %5d.%d]\n",
                           this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //nMethod -> solveYourselfAt(this->giveCurrentStep()) ;
        nMethod->setStepLength(deltaL / 5.0);
        if ( initialLoadVector.isNotEmpty() ) {
            numMetStatus = nMethod->solve( stiffnessMatrix, & incrementalLoadVector, & initialLoadVector,
                                          & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        } else {
            numMetStatus = nMethod->solve( stiffnessMatrix, & incrementalLoadVector, NULL,
                                          & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        }


        loadVector.zero();

        this->updateYourself( this->giveCurrentStep() );
        this->terminate( this->giveCurrentStep() );
        this->updateLoadVectors( this->giveCurrentStep() );

        // restore old step length
        nMethod->setStepLength(deltaL);

        stiffMode = oldStiffMode;
    } else {
        // comment this, if output for mapped configuration (not equilibrated) not wanted
        this->printOutputAt( this->giveOutputStream(), this->giveCurrentStep() );
    }

    return result;
}


int
AdaptiveNonLinearStatic :: initializeAdaptive(int stepNumber)
{
    int stepinfo [ 2 ];

    stepinfo [ 0 ] = stepNumber;
    stepinfo [ 1 ] = 0;

    try {
        this->restoreContext(NULL, CM_State, ( void * ) stepinfo);
    } catch(ContextIOERR & c) {
        c.print();
        exit(1);
    }

    this->initStepIncrements();

    int sernum = this->giveDomain(1)->giveSerialNumber();
    OOFEM_LOG_INFO("restoring domain %d.%d\n", 1, sernum + 1);
    Domain *dNew = new Domain(2, sernum + 1, this);
    DataReader *domainDr = this->GiveDomainDataReader(1, sernum + 1, contextMode_read);
    if ( !dNew->instanciateYourself(domainDr) ) {
        _error("initializeAdaptive: domain Instanciation failed");
    }

    delete domainDr;

    // remap solution to new domain
    return this->adaptiveRemap(dNew);
}


int
AdaptiveNonLinearStatic :: adaptiveRemap(Domain *dNew)
{
    int ielem, nelem, result = 1;

    this->initStepIncrements();

    this->ndomains = 2;
    this->domainNeqs.resize(2);
    this->domainPrescribedNeqs.resize(2);
    this->domainNeqs.at(2) = 0;
    this->domainPrescribedNeqs.at(2) = 0;
    this->domainList->put(2, dNew);

#ifdef __PETSC_MODULE
    PetscContext *pcNew = new PetscContext(this, EID_MomentumBalance);

    this->petscContextList->put(2, pcNew);
#endif

    // init equation numbering
    //this->forceEquationNumbering(2);
    this->forceEquationNumbering();

    // measure time consumed by mapping
    //clock_t sc = this->getClock();
    //clock_t mc, mc2, ec ;
    oofem_timeval st, st1, st2, mc1, mc2, mc3, ec;
    getUtime(st);


    // map primary unknowns
    EIPrimaryUnknownMapper mapper;

    d2_totalDisplacement.resize( this->giveNumberOfDomainEquations(2, EID_MomentumBalance) );
    d2_incrementOfDisplacement.resize( this->giveNumberOfDomainEquations(2, EID_MomentumBalance) );
    d2_totalDisplacement.zero();
    d2_incrementOfDisplacement.zero();

    result &= mapper.mapAndUpdate( d2_totalDisplacement, VM_Total, EID_MomentumBalance,
                                  this->giveDomain(1), this->giveDomain(2), this->giveCurrentStep() );

    result &= mapper.mapAndUpdate( d2_incrementOfDisplacement, VM_Incremental, EID_MomentumBalance,
                                  this->giveDomain(1), this->giveDomain(2), this->giveCurrentStep() );

    //mc = this->getClock();
    getRelativeUtime(mc1, st);
    getUtime(st1);

    // map internal ip state
    nelem = this->giveDomain(2)->giveNumberOfElements();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( this->giveDomain(2)->giveElement(ielem)->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif


        result &= this->giveDomain(2)->giveElement(ielem)->adaptiveMap( this->giveDomain(1), this->giveCurrentStep() );
    }

    /* replace domains */
    OOFEM_LOG_DEBUG("deleting old domain\n");
    //delete domainList->at(1);
    //domainList->put(1, dNew);
    //dNew->setNumber(1);
    //domainList->put(2, NULL);
    domainList->put( 1, domainList->unlink(2) );
    domainList->at(1)->setNumber(1);

#ifdef __PETSC_MODULE
    petscContextList->put( 1, petscContextList->unlink(2) );
    petscContextList->growTo(1);
#endif

    // keep equation numbering of new domain
    this->numberOfEquations = this->domainNeqs.at(1) = this->domainNeqs.at(2);
    this->numberOfPrescribedEquations = this->domainPrescribedNeqs.at(1) = this->domainPrescribedNeqs.at(2);
    this->equationNumberingCompleted = 1;

    // update solution
    totalDisplacement = d2_totalDisplacement;
    incrementOfDisplacement = d2_incrementOfDisplacement;


    this->ndomains = 1;
    // init equation numbering
    // this->forceEquationNumbering();
    this->updateDomainLinks();

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        // set up communication patterns
        this->initializeCommMaps(true);
        this->exchangeRemoteElementData( RemoteElementExchangeTag );
    }

#endif

    //mc2 = this->getClock();
    getRelativeUtime(mc2, st1);
    getUtime(st2);

    // computes the stresses and calls updateYourself to mapped state
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( this->giveDomain(1)->giveElement(ielem)->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        result &= this->giveDomain(1)->giveElement(ielem)->adaptiveUpdate( this->giveCurrentStep() );
    }

    // finish mapping process
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( this->giveDomain(1)->giveElement(ielem)->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        result &= this->giveDomain(1)->giveElement(ielem)->adaptiveFinish( this->giveCurrentStep() );
    }

    nMethod->reinitialize();


    // increment time step if mapped state will be considered as new solution stepL
    // this->giveNextStep();
    if ( equilibrateMappedConfigurationFlag ) {
        // we need to  assemble the load vector in same time as the restarted step,
        // so new time step is generated with same intrincic time as has the
        // previous step if equilibrateMappedConfigurationFlag is set.
        // this allows to equlibrate the previously reached state
        TimeStep *cts = this->giveCurrentStep();
        // increment version of solution step
        cts->incrementVersion();

        //cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
        //cts = this->givePreviousStep();
        //cts->setTime(cts->giveTime()-cts->giveTimeIncrement());
    }

    if ( this->giveCurrentStep()->giveNumber() == this->giveCurrentMetaStep()->giveFirstStepNumber() ) {
        this->updateAttributes( this->giveCurrentMetaStep() );
    }

    // increment solution state counter - not needed, IPs are updated by adaptiveUpdate previously called
    // and there is no change in primary vars.
    // this->giveCurrentStep()->incrementStateCounter();

    // assemble new initial load for new discretization
    this->assembleInitialLoadVector( initialLoadVector, initialLoadVectorOfPrescribed,
                                    this, 1, this->giveCurrentStep() );
    this->assembleIncrementalReferenceLoadVectors( incrementalLoadVector, incrementalLoadVectorOfPrescribed,
                                                  refLoadInputMode, this->giveDomain(1), EID_MomentumBalance,
                                                  this->giveCurrentStep() );

    // assemble new total load for new discretization
    // this->assembleCurrentTotalLoadVector (totalLoadVector, totalLoadVectorOfPrescribed, this->giveCurrentStep());
    // set bcloadVector to zero (no increment within same step)

    //ec = this->getClock();
    getRelativeUtime(mc3, st2);
    getRelativeUtime(ec, st);

    // compute processor time used by the program
    OOFEM_LOG_INFO( "user time consumed by primary mapping: %.2fs\n", ( double ) ( mc1.tv_sec + mc1.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by ip mapping:      %.2fs\n", ( double ) ( mc2.tv_sec + mc2.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by ip update:       %.2fs\n", ( double ) ( mc3.tv_sec + mc3.tv_usec / ( double ) OOFEM_USEC_LIM ) );
    OOFEM_LOG_INFO( "user time consumed by mapping:         %.2fs\n", ( double ) ( ec.tv_sec + ec.tv_usec / ( double ) OOFEM_USEC_LIM ) );

    //

    /********
     * #if 0
     * {
     *
     * // evaluate the force error of mapped configuration
     * this->updateComponent (this->giveCurrentStep(), InternalRhs);
     * FloatArray rhs;
     *
     * loadVector.resize (this->numberOfEquations);
     * loadVector.zero();
     * this->assemble (loadVector, this->giveCurrentStep(), ExternalForcesVector_Total, this->giveDomain(1)) ;
     * this->assemble (loadVector, this->giveCurrentStep(), ExternalForcesVector_Total, this->giveDomain(1));
     *
     * rhs =  loadVector;
     * rhs.times(loadLevel);
     * if (initialLoadVector.isNotEmpty()) rhs.add(initialLoadVector);
     * rhs.subtract(internalForces);
     *
     * //
     * // compute forceError
     * //
     * // err is relative error of unbalanced forces
     * double RR, RR0, forceErr = dotProduct (rhs.givePointer(),rhs.givePointer(),rhs.giveSize());
     * if (initialLoadVector.isNotEmpty())
     * RR0 = dotProduct (initialLoadVector.givePointer(), initialLoadVector.givePointer(), initialLoadVector.giveSize());
     * else
     * RR0 = 0.0;
     * RR = dotProduct(loadVector.givePointer(),loadVector.givePointer(),loadVector.giveSize());
     * // we compute a relative error norm
     * if ((RR0 + RR * loadLevel * loadLevel) < calm_SMALL_NUM) forceErr = 0.;
     * else forceErr = sqrt (forceErr / (RR0+RR * loadLevel * loadLevel));
     *
     * printf ("Relative Force Error of Mapped Configuration is %-15e\n", forceErr);
     *
     * }
     *#endif
     *************/

#ifdef __OOFEG
    ESIEventLoop( YES, oofem_tmpstr("AdaptiveRemap: Press Ctrl-p to continue") );
#endif

    //
    // bring mapped configuration into equilibrium
    //
    if ( equilibrateMappedConfigurationFlag ) {
        // use secant stiffness to resrore equilibrium
        NonLinearStatic_stifnessMode oldStiffMode = this->stiffMode;
        stiffMode = nls_secantStiffness;



        if ( initFlag ) {
            if ( !stiffnessMatrix ) {
                stiffnessMatrix = CreateUsrDefSparseMtrx(sparseMtrxType);
                if ( stiffnessMatrix == NULL ) {
                    _error("proceedStep: sparse matrix creation failed");
                }
            }

            if ( nonlocalStiffnessFlag ) {
                if ( !stiffnessMatrix->isAsymmetric() ) {
                    _error("proceedStep: stiffnessMatrix does not support asymmetric storage");
                }
            }

            stiffnessMatrix->buildInternalStructure( this, 1, EID_MomentumBalance, EModelDefaultEquationNumbering() );
            stiffnessMatrix->zero(); // zero stiffness matrix
            this->assemble( stiffnessMatrix, this->giveCurrentStep(), EID_MomentumBalance, SecantStiffnessMatrix,
                           EModelDefaultEquationNumbering(), this->giveDomain(1) );
            initFlag = 0;
        }

        // updateYourself() not necessary - the adaptiveUpdate previously called does the job
        //this->updateYourself(this->giveCurrentStep());
#ifdef VERBOSE
        OOFEM_LOG_INFO( "Equilibrating mapped configuration [step number %5d.%d]\n",
                       this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //double deltaL = nMethod->giveUnknownComponent (StepLength, 0);
        double deltaL = nMethod->giveCurrentStepLength();
        //
        // call numerical model to solve arised problem
        //
#ifdef VERBOSE
        OOFEM_LOG_RELEVANT( "Solving [step number %5d.%d]\n",
                           this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion() );
#endif

        //nMethod -> solveYourselfAt(this->giveCurrentStep()) ;
        nMethod->setStepLength(deltaL / 5.0);
        if ( initialLoadVector.isNotEmpty() ) {
            numMetStatus = nMethod->solve( stiffnessMatrix, & incrementalLoadVector, & initialLoadVector,
                                          & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        } else {
            numMetStatus = nMethod->solve( stiffnessMatrix, & incrementalLoadVector, NULL,
                                          & totalDisplacement, & incrementOfDisplacement, & internalForces,
                                          internalForcesEBENorm, loadLevel, refLoadInputMode, currentIterations, this->giveCurrentStep() );
        }


        loadVector.zero();

        this->updateYourself( this->giveCurrentStep() );
        this->terminate( this->giveCurrentStep() );
        // this->updateLoadVectors (this->giveCurrentStep()); // already in terminate

        // restore old step length
        nMethod->setStepLength(deltaL);

        stiffMode = oldStiffMode;
    } else {
        // comment this, if output for mapped configuration (not equilibrated) not wanted
        this->printOutputAt( this->giveOutputStream(), this->giveCurrentStep() );
    }

    return result;
}


contextIOResultType
AdaptiveNonLinearStatic :: saveContext(DataStream *stream, ContextMode mode, void *obj) {
    int closeFlag = 0;
    contextIOResultType iores;
    FILE *file;

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = NonLinearStatic :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = timeStepLoadLevels.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                       // ensure consistent records

    return CIO_OK;
}

contextIOResultType
AdaptiveNonLinearStatic :: restoreContext(DataStream *stream, ContextMode mode, void *obj) {
    int closeFlag = 0;
    int istep, iversion;
    contextIOResultType iores;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    if ( ( iores = NonLinearStatic :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = timeStepLoadLevels.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                       // ensure consistent records

    return CIO_OK;
}


void
AdaptiveNonLinearStatic :: updateDomainLinks() {
    NonLinearStatic :: updateDomainLinks();
    this->defaultErrEstimator->setDomain( this->giveDomain(1) );
}


void
AdaptiveNonLinearStatic :: assembleInitialLoadVector(FloatArray &loadVector, FloatArray &loadVectorOfPrescribed,
                                                     AdaptiveNonLinearStatic *sourceProblem, int domainIndx,
                                                     TimeStep *atTime)
{
    const char *__proc = "assembleInitialLoadVector"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                           // Required by IR_GIVE_FIELD macro

    int mstepNum = atTime->giveMetaStepNumber();
    int imstep;
    int hasfixed, mode;
    InputRecord *ir;
    MetaStep *iMStep;
    FloatArray _incrementalLoadVector, _incrementalLoadVectorOfPrescribed;
    SparseNonLinearSystemNM :: referenceLoadInputModeType rlm;
    //Domain* sourceDomain = sourceProblem->giveDomain(domainIndx);

    loadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    loadVectorOfPrescribed.resize( this->giveNumberOfPrescribedEquations(EID_MomentumBalance) );
    loadVector.zero();
    loadVectorOfPrescribed.zero();
    _incrementalLoadVector.resize( this->giveNumberOfEquations(EID_MomentumBalance) );
    _incrementalLoadVectorOfPrescribed.resize( this->giveNumberOfPrescribedEquations(EID_MomentumBalance) );
    _incrementalLoadVector.zero();
    _incrementalLoadVectorOfPrescribed.zero();

    for ( imstep = 1; imstep < mstepNum; imstep++ ) {
        iMStep = this->giveMetaStep(imstep);
        ir = iMStep->giveAttributesRecord();
        //hasfixed = ir->hasField("fixload");
        hasfixed = 1;
        if ( hasfixed ) {
            // test for control mode
            // here the algorithm works only for direct load control.
            // Direct displacement control requires to know the quasi-rections, and the controlled nodes
            // should have corresponding node on new mesh -> not supported
            // Indirect control -> the load level from prevous steps is required, currently nt supported.

            // additional problem: direct load control supports the reduction of step legth if convergence fails
            // if this happens, this implementation does not work correctly.
            // But there is NO WAY HOW TO TEST IF THIS HAPPEN

            mode = 0;
            IR_GIVE_OPTIONAL_FIELD(ir, mode, IFT_AdaptiveNonLinearStatic_controlmode, "controlmode"); // Macro

            // check if displacement control takes place
            if ( ir->hasField(IFT_AdaptiveNonLinearStatic_ddm, "ddm") ) {
                _error("assembleInitialLoadVector: fixload recovery not supported for direct displacement control");
            }

            int firststep = iMStep->giveFirstStepNumber();
            int laststep  = iMStep->giveLastStepNumber();

            int _val = 0;
            IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_AdaptiveNonLinearStatic_refloadmode, "refloadmode"); // Macro
            rlm = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

            if ( mode == ( int ) nls_directControl ) { // and only load control
                for ( int istep = firststep; istep <= laststep; istep++ ) {
                    // bad practise here
                    TimeStep *old = new TimeStep(istep, this, imstep, istep - 1.0, deltaT, 0);
                    this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                                  rlm, this->giveDomain(domainIndx), EID_MomentumBalance, old);

                    _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(istep) );
                    loadVector.add(_incrementalLoadVector);
                    loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
                }
            } else if ( mode == ( int ) nls_indirectControl ) {
                // bad practise here
                if ( !ir->hasField(IFT_NonLinearStatic_donotfixload, "donotfixload") ) {
                    TimeStep *old = new TimeStep(firststep, this, imstep, firststep - 1.0, deltaT, 0);
                    this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                                  rlm, this->giveDomain(domainIndx), EID_MomentumBalance, old);

                    _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(laststep) );
                    loadVector.add(_incrementalLoadVector);
                    loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
                }
            } else {
                _error("assembleInitialLoadVector: fixload recovery not supported");
            }
        }
    } // end loop over meta-steps

    /* if direct control; add to initial load also previous steps in same metestep */
    iMStep = this->giveMetaStep(mstepNum);
    ir = iMStep->giveAttributesRecord();
    mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, IFT_AdaptiveNonLinearStatic_controlmode, "controlmode"); // Macro
    int firststep = iMStep->giveFirstStepNumber();
    int laststep  = atTime->giveNumber();
    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_AdaptiveNonLinearStatic_refloadmode, "refloadmode"); // Macro
    rlm = ( SparseNonLinearSystemNM :: referenceLoadInputModeType ) _val;

    if ( mode == ( int ) nls_directControl ) { // and only load control
        iMStep = this->giveMetaStep(imstep);
        for ( int istep = firststep; istep <= laststep; istep++ ) {
            // bad practise here
            TimeStep *old = new TimeStep(istep, this, imstep, istep - 1.0, deltaT, 0);
            this->assembleIncrementalReferenceLoadVectors(_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
                                                          rlm, this->giveDomain(domainIndx), EID_MomentumBalance, old);

            _incrementalLoadVector.times( sourceProblem->giveTimeStepLoadLevel(istep) );
            loadVector.add(_incrementalLoadVector);
            loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
        }
    }
}

/*
 * void
 * AdaptiveNonLinearStatic::assembleCurrentTotalLoadVector (FloatArray& loadVector,
 *                           FloatArray& loadVectorOfPrescribed,
 *                           AdaptiveNonLinearStatic* sourceProblem, int domainIndx,
 *                           TimeStep* atTime)
 * {
 * const char *__proc = "assembleInitialLoadVector"; // Required by IR_GIVE_FIELD macro
 * IRResultType result;                              // Required by IR_GIVE_FIELD macro
 *
 * int mstepNum = atTime->giveMetaStepNumber() ;
 * int mode;
 * InputRecord* ir;
 * MetaStep* mStep = sourceProblem->giveMetaStep(mstepNum);
 * FloatArray _incrementalLoadVector, _incrementalLoadVectorOfPrescribed;
 * SparseNonLinearSystemNM::referenceLoadInputModeType rlm;
 * //Domain* sourceDomain = sourceProblem->giveDomain(domainIndx);
 *
 * loadVector.resize(this->giveNumberOfEquations(EID_MomentumBalance));
 * loadVectorOfPrescribed.resize(this->giveNumberOfPrescribedEquations(EID_MomentumBalance));
 * loadVector.zero(); loadVectorOfPrescribed.zero();
 * _incrementalLoadVector.resize(this->giveNumberOfEquations(EID_MomentumBalance));
 * _incrementalLoadVectorOfPrescribed.resize(this->giveNumberOfPrescribedEquations(EID_MomentumBalance));
 * _incrementalLoadVector.zero(); _incrementalLoadVectorOfPrescribed.zero();
 *
 * // ASK CURRENT MSTEP FOR ITS RECORD
 * ir = mStep->giveAttributesRecord();
 *
 * mode = 0;
 * IR_GIVE_OPTIONAL_FIELD (ir, mode, IFT_AdaptiveNonLinearStatic_controlmode, "controlmode"); // Macro
 *
 * // check if displacement control takes place
 * if (ir->hasField(IFT_AdaptiveNonLinearStatic_ddm, "ddm"))
 * _error ("assembleCurrentTotalLoadVector: fixload recovery not supported for direct displacement control");
 * int _val = 0;
 * IR_GIVE_OPTIONAL_FIELD (ir, _val, IFT_AdaptiveNonLinearStatic_refloadmode, "refloadmode"); // Macro
 *
 * int firststep = mStep->giveFirstStepNumber();
 * int laststep  = atTime->giveNumber()-1;
 *
 * rlm = (SparseNonLinearSystemNM::referenceLoadInputModeType) _val;
 *
 * if (mode == (int)nls_directControl) { // and only load control
 * for (int istep = firststep; istep<=laststep; istep++) {
 * // bad practise here
 * TimeStep* old = new TimeStep (istep, this, mstepNum, istep-1.0, deltaT, 0);
 * this->assembleIncrementalReferenceLoadVectors (_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
 *                         rlm, this->giveDomain(domainIndx), EID_MomentumBalance, old);
 *
 * _incrementalLoadVector.times(sourceProblem->giveTimeStepLoadLevel(istep));
 * loadVector.add(_incrementalLoadVector);
 * loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
 * }
 * } else if (mode == (int)nls_indirectControl) {
 * // bad practise here
 * TimeStep* old = new TimeStep (firststep, this, mstepNum, firststep-1.0, deltaT, 0);
 * this->assembleIncrementalReferenceLoadVectors (_incrementalLoadVector, _incrementalLoadVectorOfPrescribed,
 *                        rlm, this->giveDomain(domainIndx), EID_MomentumBalance, old);
 *
 * _incrementalLoadVector.times(sourceProblem->giveTimeStepLoadLevel(laststep));
 * loadVector.add(_incrementalLoadVector);
 * loadVectorOfPrescribed.add(_incrementalLoadVectorOfPrescribed);
 * } else {
 * _error ("assembleCurrentTotalLoadVector: totalload recovery not supported");
 * }
 *
 *
 * }
 */


double
AdaptiveNonLinearStatic :: giveTimeStepLoadLevel(int istep)
{
    if ( ( istep >= this->giveNumberOfFirstStep() ) && ( istep <= this->giveNumberOfSteps() ) ) {
        return timeStepLoadLevels.at(istep);
    } else {
        _error("solution step out of range");
    }

    return 0.0; // to make compiler happy
}


#ifdef __PARALLEL_MODE
LoadBalancer *
AdaptiveNonLinearStatic :: giveLoadBalancer()
{
    if ( lb ) {
        return lb;
    }

    if ( loadBalancingFlag || preMappingLoadBalancingFlag ) {
        lb = CreateUsrDefLoadBalancerOfType( ParmetisLoadBalancerClass, this->giveDomain(1) );
        return lb;
    } else {
        return NULL;
    }
}
LoadBalancerMonitor *
AdaptiveNonLinearStatic :: giveLoadBalancerMonitor()
{
    if ( lbm ) {
        return lbm;
    }

    if ( loadBalancingFlag || preMappingLoadBalancingFlag ) {
        lbm = CreateUsrDefLoadBalancerMonitorOfType(WallClockLoadBalancerMonitorClass, this);
        return lbm;
    } else {
        return NULL;
    }
}
#endif
} // end namespace oofem

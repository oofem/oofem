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

#include "huertaerrorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "load.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "loadtime.h"
#include "timestep.h"
#include "metastep.h"
#include "integrationrule.h"
#include "conTable.h"
#include "crosssection.h"
#include "dof.h"
#include "util.h"
#include "eleminterpunknownmapper.h"
#include "buffereddatareader.h"
#include "adaptnlinearstatic.h"
#include "verbose.h"
#include "datastream.h"
#ifndef __MAKEDEPEND
 #include <vector>
 #include <string>
#endif
#include "contextioerr.h"

namespace oofem {
//#define STIFFNESS_TYPE       TangentStiffnessMatrix
#define STIFFNESS_TYPE       ElasticStiffnessMatrix
//#define STIFFNESS_TYPE       SecantStiffnessMatrix

#define ERROR_EXCESS                0.1   // 10 %


// debug defines

//#define USE_FILE

#define PRINT_ERROR

#ifdef VERBOSE
 #define INFO
//#define TIME_INFO
 #define EXACT_ERROR
 #define PRINT_ERROR
#endif

#ifdef USE_FILE
 #define USE_INPUT_FILE
//#define USE_OUTPUT_FILE
//#define USE_CONTEXT_FILE
#endif

#ifdef PRINT_ERROR
//#define PRINT_FINE_ERROR
 #define PRINT_COARSE_ERROR
#endif

static bool masterRun = true;

static bool exactFlag = false;

#ifdef EXACT_ERROR
static bool wholeFlag = false, huertaFlag = false;

double exactENorm, coarseUNorm, fineUNorm, mixedNorm;

 #ifdef PRINT_ERROR
static int finePos;
static FloatArray exactFineError;
static FloatArray exactCoarseError;
 #endif
#endif

//static FloatArray uNormArray;

static BufferedDataReader refinedReader;

static int impMat, perMat;
static FloatArray impPos;

static int globalNelems;


int
HuertaErrorEstimator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
    Domain *d = this->domain;
    EngngModel *model = d->giveEngngModel();
    int ielem, nelems = d->giveNumberOfElements();
    int inode, nnodes = d->giveNumberOfDofManagers();
    double pe;
    IntArray localNodeIdArray, globalNodeIdArray;    // these arrays are declared here to
                                                     // prevent their repeated creation for
                                                     // each element and patch problem

    if ( this->stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

    this->globalENorm = 0.0;
    this->globalUNorm = 0.0;
    this->globalWENorm = 0.0;

    // initiate to prevent oofeg failure (which loads eNorms from context)
    this->eNorms.resize(nelems);
    this->eNorms.zero();

    if ( initialSkipSteps != 0 ) {
        initialSkipSteps--;

        skippedNelems = nelems;
        this->stateCounter = tStep->giveSolutionStateCounter();
        OOFEM_LOG_RELEVANT( "Relative error estimate [step number %5d]: skipped\n",
                           d->giveEngngModel()->giveCurrentStep()->giveNumber() );
        return 1;
    }

    if ( stepsToSkip != 0 ) {
        stepsToSkip--;
        skippedSteps++;

        skippedNelems = nelems;
        this->stateCounter = tStep->giveSolutionStateCounter();
        OOFEM_LOG_RELEVANT( "\nRelative error estimate [step number %5d]: skipped\n",
                           d->giveEngngModel()->giveCurrentStep()->giveNumber() );
        return 1;
    }

#ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Estimating error [step number %5d]\n",
                   d->giveEngngModel()->giveRank(),
                   d->giveEngngModel()->giveCurrentStep()->giveNumber() );
#else
    OOFEM_LOG_INFO( "Estimating error [step number %5d]\n",
                   d->giveEngngModel()->giveCurrentStep()->giveNumber() );
#endif

    switch ( d->giveEngngModel()->giveClassID() ) {
    case AdaptiveLinearStaticClass:
        this->mode = HEE_linear;
        break;
    case AdaptiveNonLinearStaticClass:
        this->mode = HEE_nlinear;
        break;
    default:
        _error("estimateError: Unsupported analysis type");
    }

    // check if each node has default number of dofs
    // check if first load time function is constant
    // check if only constant edge or surface load is used

    this->buildRefinedMesh();

    localNodeIdArray.resize(this->refinedMesh.nodes);
    localNodeIdArray.zero();

    this->skippedNelems = 0;

#ifdef EXACT_ERROR
    if ( exactFlag == true ) {
 #ifdef PRINT_ERROR
        finePos = 0;
        exactFineError.resize(this->refinedMesh.elems);
        exactCoarseError.resize(nelems);
 #endif

        wholeFlag = true;
        this->solveRefinedWholeProblem(localNodeIdArray, globalNodeIdArray, tStep);
        wholeFlag = false;

        if ( huertaFlag == false ) {
            OOFEM_LOG_INFO("\n\n");
            OOFEM_LOG_INFO("Exact eNorm2:        %15.8e\n", exactENorm);
            OOFEM_LOG_INFO("Exact coarse uNorm2: %15.8e\n", coarseUNorm);
            //   fprintf(stdout, "Exact mixed euNorm:  %15.8e\n", mixedNorm);
            OOFEM_LOG_INFO("Exact fine uNorm2:   %15.8e\n", fineUNorm);

            pe = sqrt( exactENorm / ( exactENorm + coarseUNorm ) );
            if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
                OOFEM_LOG_INFO("Exact relative error (coarse): %6.3f%% (L2 norm)\n", pe * 100.0);
            }

            if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
                OOFEM_LOG_INFO("Exact relative error (coarse): %6.3f%% (energy norm)\n", pe * 100.0);
            }

            pe = sqrt(exactENorm / fineUNorm);
            if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
                OOFEM_LOG_INFO("Exact relative error (fine):   %6.3f%% (L2 norm)\n", pe * 100.0);
            }

            if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
                OOFEM_LOG_INFO("Exact relative error (fine):   %6.3f%% (energy norm)\n", pe * 100.0);
            }

            exit(1);
        }
    }

#endif

    primaryUnknownError.resize( this->refinedMesh.nodes * d->giveDefaultNodeDofIDArry().giveSize() );
    primaryUnknownError.zero();

    // uNormArray.resize(nelems);

    // first loop over patches and then over elements;
    // this is advantageous, because when orthogonalizing element and patch error I have
    // to assemble again the stifness matrix that can be then used (as whole) for evaluation
    // of error in the energy norm;
    // when using opposite way, I have to evaluate error on fine elements and than sum
    // contributions to coarse elements;
    // in this case it would be better to evaluate the error on the patch as a whole
    // and assoociating it with the corresponding node;
    // but this would complicate remeshing criterion

    // however there is problem that in nonlinear analysis stiffness matrix of the same fine
    // element is different for element and for patch problem because each may be at slightly
    // different point of the solution

    // freopen("/dev/null", "w", stdout);

    for ( inode = 1; inode <= nnodes; inode++ ) {
        this->solveRefinedPatchProblem(inode, localNodeIdArray, globalNodeIdArray, tStep);
    }

    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        this->solveRefinedElementProblem(ielem, localNodeIdArray, globalNodeIdArray, tStep);
    }

#ifdef __PARALLEL_MODE
 #ifdef __USE_MPI
    double buffer_out [ 4 ], buffer_in [ 4 ];

    buffer_out [ 0 ] = globalENorm;
    buffer_out [ 1 ] = globalUNorm;
    buffer_out [ 2 ] = ( double ) nelems;
    buffer_out [ 3 ] = ( double ) skippedNelems;

    MPI_Allreduce(buffer_out, buffer_in, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    globalENorm = buffer_in [ 0 ];
    globalUNorm = buffer_in [ 1 ];
    globalNelems = ( int ) ( buffer_in [ 2 ] + 0.1 );
    skippedNelems = ( int ) ( buffer_in [ 3 ] + 0.1 );
 #endif
#else
    globalNelems = nelems;
#endif

#ifdef PRINT_COARSE_ERROR
    OOFEM_LOG_DEBUG("\n");
    if ( exactFlag == true ) {
        OOFEM_LOG_DEBUG("  elem        a_Error2        x_Error2         a/x rate2 \n");
        OOFEM_LOG_DEBUG("---------------------------------------------------------\n");
    } else {
        OOFEM_LOG_DEBUG("  elem        a_Error2 \n");
        OOFEM_LOG_DEBUG("-----------------------\n");
    }

    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        if ( exactFlag == false ) {
 #ifdef __PARALLEL_MODE
            OOFEM_LOG_DEBUG("[%d] %5d: %15.8e %s\n", d->giveEngngModel()->giveRank(), ielem,
                            this->eNorms.at(ielem) * this->eNorms.at(ielem),
                            ( this->skipRegion( d->giveElement(ielem)->giveRegionNumber() ) != 0 ) ? "(skipped)" :
                            ( d->giveElement(ielem)->giveParallelMode() == Element_remote ) ? "(remote)" : "");
 #else
            OOFEM_LOG_DEBUG("%5d: %15.8e %s\n", ielem,
                            this->eNorms.at(ielem) * this->eNorms.at(ielem),
                            ( this->skipRegion( d->giveElement(ielem)->giveRegionNumber() ) != 0 ) ? "(skipped)" : "");
 #endif
        }

 #ifdef EXACT_ERROR
        else {
            if ( fabs( exactCoarseError.at(ielem) ) > 1.0e-30 && this->eNorms.at(ielem) != 0.0 ) {
                OOFEM_LOG_DEBUG("%5d: %15.8e %15.8e   %15.8e %s\n", ielem,
                                this->eNorms.at(ielem) * this->eNorms.at(ielem), exactCoarseError.at(ielem),
                                this->eNorms.at(ielem) / sqrt( exactCoarseError.at(ielem) ),
                                ( this->skipRegion( d->giveElement(ielem)->giveRegionNumber() ) != 0 ) ? "(skipped)" : "");
            } else {
                OOFEM_LOG_DEBUG("%5d: %15.8e %15.8e          N/A %s\n", ielem,
                                this->eNorms.at(ielem) * this->eNorms.at(ielem), exactCoarseError.at(ielem),
                                ( this->skipRegion( d->giveElement(ielem)->giveRegionNumber() ) != 0 ) ? "(skipped)" : "");
            }
        }
 #endif
    }

#endif

    double pwe = 0.0;
    if ( wError == true ) {
        // calculate weighted error
        double elemErrLimit;

        elemErrLimit = sqrt( ( this->globalENorm + this->globalUNorm ) / ( globalNelems - skippedNelems ) ) * requiredError;
        if ( elemErrLimit != 0.0 ) {
            double eerror, iratio;

            for ( ielem = 1; ielem <= nelems; ielem++ ) {
                eerror = this->eNorms.at(ielem);
                iratio = eerror / elemErrLimit;
                globalWENorm += eerror * eerror * iratio;
            }

#ifdef __PARALLEL_MODE
 #ifdef __USE_MPI
            double myGlobalWENorm = globalWENorm;
            MPI_Allreduce(& myGlobalWENorm, & globalWENorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 #endif
#endif

            pwe = sqrt( globalWENorm / ( this->globalENorm + this->globalUNorm ) );
        }
    }

    OOFEM_LOG_INFO("\n");
    OOFEM_LOG_INFO("Global eNorm2: %15.8e\n", this->globalENorm);
    if ( wError == true ) {
        OOFEM_LOG_INFO("Global wNorm2: %15.8e\n", globalWENorm);
    }

    OOFEM_LOG_INFO("Global uNorm2: %15.8e\n", this->globalUNorm);

    // report the error estimate
    pe = sqrt( this->globalENorm / ( this->globalENorm + this->globalUNorm ) );
    if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
        OOFEM_LOG_INFO("Relative error estimate [step number %5d]: %6.3f%% (L2 norm)\n",
                       d->giveEngngModel()->giveCurrentStep()->giveNumber(), pe * 100.0);
    }

    if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
        OOFEM_LOG_INFO("Relative error estimate [step number %5d]: %6.3f%% (energy norm)\n",
                       d->giveEngngModel()->giveCurrentStep()->giveNumber(), pe * 100.0);
    }

    if ( wError == true ) {
        OOFEM_LOG_INFO("Relative error estimate [step number %5d]: %6.3f%% (weighted)\n",
                       d->giveEngngModel()->giveCurrentStep()->giveNumber(), pwe * 100.0);
    }

    //fflush(stdout);

    if ( maxSkipSteps != 0 && this->mode == HEE_nlinear ) {
        if ( lastError < 0.0 ) {
            lastError = pe;
            stepsToSkip = 0;
        } else {
            if ( requiredError > pe ) {
                if ( pe <= lastError ) {
                    stepsToSkip = maxSkipSteps;
                } else {
                    // estimate number of steps to skip using linear extrapolation
                    stepsToSkip = ( int ) ( ( requiredError - pe ) / ( pe - lastError ) * ( skippedSteps + 1 ) );
                    // make the number of steps to skip safe with respect to how many steps have been skipped last time
                    stepsToSkip = stepsToSkip / ( skippedSteps + 1 );
                }

                // make the number of steps to skip safe with respect to the absolute value of the error
                // (decrease by 0.5 for pe equal requiredError)
                stepsToSkip = ( int ) ( stepsToSkip * ( ( 0.5 - 1.0 ) / requiredError * pe + 1 ) + 0.5 );
                if ( skippedSteps == 0 ) {
                    stepsToSkip--;
                }

                if ( stepsToSkip > maxSkipSteps ) {
                    stepsToSkip = maxSkipSteps;
                }

                if ( stepsToSkip < 0 ) {
                    stepsToSkip = 0;
                }
            } else {
                // allow error excess and prevent recursion
                if ( requiredError * ( 1.00 + ERROR_EXCESS ) < pe && skippedSteps != 0 ) {
                    /* HUHU CHEATING do not return just by one step */
                    if ( maxSkipSteps == 1 ) {
                        return 1;
                    }

                    int stepNumber, curNumber = tStep->giveNumber();

                    stepNumber = curNumber - skippedSteps;

                    OOFEM_LOG_INFO("Returning to step %5d\n", stepNumber);

                    skippedSteps = 0;                // prevent recursion when returning
                    maxSkipSteps /= 2;               // prevent oscillation

                    // I trying to avoid repeated solution of (all, not only the one to which I am returning) previous steps
                    // IMPORTANT: I must restore context only from steps corresponding to this session !!!
                    //            this means I cannot go in previous adaptive run, because there was a different domain
                    // it would be much cleaner to call restore from engng model
                    while ( stepNumber < curNumber ) {
                        try {
                            model->restoreContext(NULL, CM_State, ( void * ) & stepNumber);
                        } catch(ContextIOERR & c) {
                            c.print();
                            exit(1);
                        }

                        stepsToSkip = 0;
                        model->giveCurrentStep()->incrementStateCounter();
                        this->estimateError( temporaryEM, model->giveCurrentStep() );

                        if ( lastError > requiredError ) {
                            return 1;
                        }

                        stepNumber += ( skippedSteps = stepsToSkip ) + 1;
                    }

                    return 1;
                }
            }

            lastError = pe;
        }

        skippedSteps = 0;
    }

#ifdef EXACT_ERROR
    if ( exactFlag == true ) {
 #ifdef __PARALLEL_MODE
  #ifdef __USE_MPI
        buffer_out [ 0 ] = exactENorm;
        buffer_out [ 1 ] = coarseUNorm;
        buffer_out [ 2 ] = mixedNorm;
        buffer_out [ 3 ] = fineUNorm;

        MPI_Allreduce(buffer_out, buffer_in, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        exactENorm = buffer_in [ 0 ];
        coarseUNorm = buffer_in [ 1 ];
        mixedNorm = buffer_in [ 2 ];
        fineUNorm = buffer_in [ 3 ];
  #endif
 #endif

        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("Exact eNorm2:        %15.8e\n", exactENorm);
        OOFEM_LOG_INFO("Exact coarse uNorm2: %15.8e\n", coarseUNorm);
        //fprintf(stdout, "Exact mixed euNorm:  %15.8e\n", mixedNorm);
        OOFEM_LOG_INFO("Exact fine uNorm2:   %15.8e\n", fineUNorm);

        pe = sqrt( exactENorm / ( exactENorm + coarseUNorm ) );
        if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
            OOFEM_LOG_INFO("Exact relative error (coarse): %6.3f%% (L2 norm)\n", pe * 100.0);
        }

        if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
            OOFEM_LOG_INFO("Exact relative error (coarse): %6.3f%% (energy norm)\n", pe * 100.0);
        }

        pe = sqrt(exactENorm / fineUNorm);
        if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
            OOFEM_LOG_INFO("Exact relative error (fine):   %6.3f%% (L2 norm)\n", pe * 100.0);
        }

        if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
            OOFEM_LOG_INFO("Exact relative error (fine):   %6.3f%% (energy norm)\n", pe * 100.0);
        }
    }

#endif

    this->globalENorm = sqrt(this->globalENorm);
    this->globalUNorm = sqrt(this->globalUNorm);
    if ( wError == true ) {
        this->globalWENorm = sqrt(this->globalWENorm);
    }

    this->stateCounter = tStep->giveSolutionStateCounter();

    return 1;
}



double
HuertaErrorEstimator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    this->estimateError(equilibratedEM, tStep);
    if ( type == primaryUnknownET ) {
        return this->eNorms.at( elem->giveNumber() );
    }

    return 0.0;
}



double
HuertaErrorEstimator :: giveValue(EE_ValueType type, TimeStep *tStep)
{
    this->estimateError(equilibratedEM, tStep);
    if ( type == globalErrorEEV ) {
        return this->globalENorm;
    } else if ( type == globalNormEEV ) {
        return this->globalUNorm;
    } else if ( type == globalWeightedErrorEEV ) {
        return this->globalWENorm;
    }

    return 0.0;
}


RemeshingCriteria *
HuertaErrorEstimator :: giveRemeshingCrit()
{
    if ( this->rc ) {
        return this->rc;
    }

    return ( this->rc = new HuertaRemeshingCriteria(1, this) );
}


IRResultType
HuertaErrorEstimator :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";    // Required by IR_GIVE_FIELD macro
    IRResultType result;                       // Required by IR_GIVE_FIELD macro
    int n, level, wErrorFlag = 0;

    ErrorEstimator :: initializeFrom(ir);
    n = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, n, IFT_HuertaErrorEstimator_normtype, "normtype"); // Macro
    if ( n == 0 ) {
        this->normType = L2Norm;
    } else {
        this->normType = EnergyNorm; // default
    }

    level = this->refineLevel;
    IR_GIVE_OPTIONAL_FIELD(ir, level, IFT_HuertaErrorEstimator_refinelevel, "refinelevel"); // Macro
    if ( level >= 0 ) {
        this->refineLevel = level;
    }

    IR_GIVE_FIELD(ir, this->requiredError, IFT_HuertaErrorEstimator_requirederror, "requirederror"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, maxSkipSteps, IFT_HuertaErrorEstimator_skipsteps, "skipsteps"); // Macro
    if ( maxSkipSteps < 0 ) {
        maxSkipSteps = 0;
    }

    if ( maxSkipSteps > 5 ) {
        maxSkipSteps = 5;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, initialSkipSteps, IFT_HuertaErrorEstimator_initialskipsteps, "initialskipsteps"); // Macro
    if ( initialSkipSteps < 0 ) {
        initialSkipSteps = 0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, wErrorFlag, IFT_HuertaErrorEstimator_werror, "werror"); // Macro
    if ( wErrorFlag != 0 ) {
        wError = true;
    }

    if ( masterRun == true ) {  // prevent overwriting of static variables
        masterRun = false;

        perMat = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, perMat, IFT_HuertaErrorEstimator_permat, "permat"); // Macro

        impMat = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, impMat, IFT_HuertaErrorEstimator_impmat, "impmat"); // Macro
        IR_GIVE_OPTIONAL_FIELD(ir, impPos, IFT_HuertaErrorEstimator_imppos, "imppos"); // Macro

        if ( impMat != 0 && perMat == 0 ) {
            _error("initializeFrom: Missing perfect material specification");
        }

#ifdef EXACT_ERROR
        n = 0;
        IR_GIVE_OPTIONAL_FIELD(ir, n, IFT_HuertaErrorEstimator_exact, "exact"); // Macro
        if ( n > 0 ) {
            exactFlag = true;
            if ( n != 1 ) {
                huertaFlag = true;         // run also error estimate
            }
        } else {
            exactFlag = false;
        }

#endif
    }

    return this->giveRemeshingCrit()->initializeFrom(ir);
}


contextIOResultType
HuertaErrorEstimator :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    if ( this->stateCounter != tStep->giveSolutionStateCounter() ) {
        this->estimateError(equilibratedEM, tStep);
    }

    // save parent class status
    if ( ( iores = ErrorEstimator :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = this->eNorms.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( !stream->write(& stateCounter, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
HuertaErrorEstimator :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = ErrorEstimator :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = eNorms.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( !stream->read(& stateCounter, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


HuertaRemeshingCriteria :: HuertaRemeshingCriteria(int n, ErrorEstimator *e) : RemeshingCriteria(n, e)
{
    this->mode = primaryUnknownBased;
    this->stateCounter = 0;
    this->refineCoeff = 1.0;
    this->noRemesh = false;
    this->wError = false;
}


double
HuertaRemeshingCriteria :: giveRequiredDofManDensity(int num, TimeStep *tStep, int relative)
{
    double size;

    this->estimateMeshDensities(tStep);
    size = this->nodalDensities.at(num);
    size = max(minElemSize, size);
    if ( relative ) {
        return size / this->giveDofManDensity(num);
    }

    return size;
}


RemeshingStrategy
HuertaRemeshingCriteria :: giveRemeshingStrategy(TimeStep *tStep)
{
    this->estimateMeshDensities(tStep);
    return this->remeshingStrategy;
}


#define MAX_COARSE_RATE     2.0
#define MAX_REFINE_RATE     5.0

int
HuertaRemeshingCriteria :: estimateMeshDensities(TimeStep *tStep)
{
    int i, j, k, nelem, nnode, jnode, elemPolyOrder, ielemNodes, skipped;
    double globValNorm = 0.0, globValErrorNorm = 0.0, globValWErrorNorm = 0.0;
    double globValNorm2, globValErrorNorm2, globValWErrorNorm2;
    double eerror, iratio, currDensity, elemSize, elemErrLimit, percentError;
    Element *ielem;
    EE_ErrorType errorType = unknownET;
    HuertaRemeshingCriteriaInterface *interface;
    bool refine = false;
    IntegrationRule *iRule;
    int nip, result;
    double sval, maxVal;
    FloatArray val;

    static int run = 0;

    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

    this->remeshingStrategy = NoRemeshing_RS;
    if ( this->noRemesh == true ) {
        // remember time stamp
        stateCounter = tStep->giveSolutionStateCounter();
        return 1;
    }

    nelem = this->domain->giveNumberOfElements();
    nnode = this->domain->giveNumberOfDofManagers();

    skipped = this->ee->giveNumberOfSkippedElements();

#ifndef __PARALLEL_MODE
    globalNelems = nelem;
#endif

    if ( skipped == globalNelems ) {
        // remember time stamp
        stateCounter = tStep->giveSolutionStateCounter();
        return 1;
    }

    // compute element error limit based on equally distribution error principle
    if ( mode == primaryUnknownBased ) {
        globValNorm      = this->ee->giveValue(globalNormEEV, tStep);
        globValErrorNorm = this->ee->giveValue(globalErrorEEV, tStep);
        globValWErrorNorm = this->ee->giveValue(globalWeightedErrorEEV, tStep);
        errorType = primaryUnknownET;
    } else {
        _error("estimateMeshDensities: unsupported mode");
    }

    globValNorm2 = globValNorm * globValNorm;
    globValErrorNorm2 = globValErrorNorm * globValErrorNorm;
    globValWErrorNorm2 = globValWErrorNorm * globValWErrorNorm;

    elemErrLimit = sqrt( ( globValNorm2 + globValErrorNorm2 ) / ( globalNelems - skipped ) ) * this->requiredError;
    if ( elemErrLimit == 0.0 ) {
        return 0;
    }

    run++;

    this->nodalDensities.resize(nnode);

    std :: vector< int >connectedElems(nnode);
    for ( i = 0; i < nnode; i++ ) {
        connectedElems [ i ] = 0;
    }

    percentError = sqrt( globValErrorNorm2 / ( globValErrorNorm2 + globValNorm2 ) );
    // test if solution is allowable
    if ( percentError > this->requiredError ) {
        this->remeshingStrategy = RemeshingFromPreviousState_RS;
    } else {
        if ( run > 100 ) {
            OOFEM_LOG_INFO("hehe\n");
            for ( i = 1; i <= nelem; i++ ) {
                ielem = domain->giveElement(i);

                // skipping these elements will result in reusing their current size;
                // in adaptive analysis this leads to monotonical decrease of their size !!!
                // it is better to either try to enlarge these elements (as they have zero error)
                // or to presribe zero nodal mesh density (but this can be handled only by t3d)
                /*
                 * if(skipped != 0){
                 * if(this -> ee -> skipRegion(ielem -> giveRegionNumber()) != 0)continue;
                 * }
                 */
                interface = ( HuertaRemeshingCriteriaInterface * ) ielem->giveInterface(HuertaRemeshingCriteriaInterfaceType);
                if ( !interface ) {
                    _error("estimateMeshDensities: element does not support HuertaRemeshingCriteriaInterface");
                }

                currDensity = interface->HuertaRemeshingCriteriaI_giveCharacteristicSize();

                //#ifdef HUHU
                // toto je treba udelat obecne
                // zjistit, zda material na prvku je nelokalni, po te z prvku vytahnout danou velicinu
                // a pokud presahne danou mez pouzit mensi z elemSize a dane size
                // chtelo by tez zaridit spusteni remeshingu, kdyz tato situace nastane -> combined
                // huerta + error indicator
                maxVal = 0.0;
                iRule = ielem->giveDefaultIntegrationRulePtr();
                nip = iRule->getNumberOfIntegrationPoints();
                for ( k = 0; k < nip; k++ ) {
                    result = ielem->giveIPValue(val, iRule->getIntegrationPoint(k), IST_PrincipalDamageTensor, tStep);
                    if ( result ) {
                        sval = val.computeNorm();
                        maxVal = max(maxVal, sval);
                    }
                }

                if ( maxVal > 0.75 ) {
                    if ( currDensity > 1.0 ) {
                        OOFEM_LOG_INFO("step %d elem %d dam %e size %e\n", tStep->giveNumber(), i, maxVal, currDensity);
                        this->remeshingStrategy = RemeshingFromPreviousState_RS;
                    }
                }
            }
        }
    }

    if ( this->remeshingStrategy == NoRemeshing_RS ) {
        // check whether to remesh because of low level of error
        if ( this->ee->giveErrorEstimatorType() == EET_HEE ) {
            /* HUHU toto je divne !!! pak se neuplatni refine viz dale */
            if ( ( ( HuertaErrorEstimator * ) ( this->ee ) )->giveAnalysisMode() == HuertaErrorEstimator :: HEE_linear ) {
                stateCounter = tStep->giveSolutionStateCounter();
                return 1;
            }

            if ( wError == false ) {
                if ( percentError <= this->requiredError ) {
                    if ( percentError >= this->requiredError * 0.5 * this->refineCoeff || run <= 5 ) {
                        for ( i = 1; i <= nnode; i++ ) {
                            nodalDensities.at(i) = this->giveDofManDensity(i);
                        }

                        stateCounter = tStep->giveSolutionStateCounter();
                        OOFEM_LOG_INFO("huhu\n");
                        return 1;
                    }
                }
            } else {
                double pwe = percentError;

                pwe = sqrt( globValWErrorNorm2 / ( globValErrorNorm2 + globValNorm2 ) );
                if ( pwe <= this->requiredError * 1.1 ) {
                    if ( pwe >= this->requiredError * 0.5 * this->refineCoeff || run <= 5 ) {
                        for ( i = 1; i <= nnode; i++ ) {
                            nodalDensities.at(i) = this->giveDofManDensity(i);
                        }

                        stateCounter = tStep->giveSolutionStateCounter();
                        return 1;
                    }
                }
            }

            OOFEM_LOG_INFO("haha\n");
            this->remeshingStrategy = RemeshingFromPreviousState_RS;
        }
    }

    /* HUHU ma zamezit aby se opakovane restartovalo ze stejneho kroku, na coz davky nejsou pripravene */
    if ( run == 1 && tStep->giveNumber() != 1 ) {
        this->remeshingStrategy = NoRemeshing_RS;
        for ( i = 1; i <= nnode; i++ ) {
            nodalDensities.at(i) = this->giveDofManDensity(i);
        }

        stateCounter = tStep->giveSolutionStateCounter();
        return 1;
    }

    for ( i = 1; i <= nelem; i++ ) {
        ielem = domain->giveElement(i);

        // skipping these elements will result in reusing their current size;
        // in adaptive analysis this leads to monotonical decrease of their size !!!
        // it is better to either try to enlarge these elements (as they have zero error)
        // or to presribe zero nodal mesh density (but this can be handled only by t3d)
        /*
         * if(skipped != 0){
         * if(this -> ee -> skipRegion(ielem -> giveRegionNumber()) != 0)continue;
         * }
         */
        interface = ( HuertaRemeshingCriteriaInterface * ) ielem->giveInterface(HuertaRemeshingCriteriaInterfaceType);
        if ( !interface ) {
            _error("estimateMeshDensities: element does not support HuertaRemeshingCriteriaInterface");
        }

        eerror = this->ee->giveElementError(errorType, ielem, tStep);
        iratio = eerror / elemErrLimit;
        if ( iratio > 1.0 ) {
            refine = true;

            // checking the step number avoids application of the refine coefficient for linear problems
            // or linear stage in nonlinear problems
            //   if(tStep -> giveNumber() != 1)iratio /= this -> refineCoeff;
            iratio /= this->refineCoeff;
            // limit the refinement (note there is also minElemSize)
            //  if (iratio > MAX_REFINE_RATE)iratio = MAX_REFINE_RATE;
        } else {
            // limit the coarsening
            if ( iratio < 1.0 / MAX_COARSE_RATE ) {
                iratio = 1.0 / MAX_COARSE_RATE;
            }
        }

        currDensity = interface->HuertaRemeshingCriteriaI_giveCharacteristicSize();
        elemPolyOrder = interface->HuertaRemeshingCriteriaI_givePolynOrder();
        elemSize = currDensity / pow(iratio, 1.0 / elemPolyOrder);
        //#ifdef HUHU
        // toto je treba udelat obecne
        // zjistit, zda material na prvku je nelokalni, po te z prvku vytahnout danou velicinu
        // a pokud presahne danou mez pouzit mensi z elemSize a dane size
        // chtelo by tez zaridit spusteni remeshingu, kdyz tato situace nastane -> combined
        // huerta + error indicator
        maxVal = 0.0;
        iRule = ielem->giveDefaultIntegrationRulePtr();
        nip = iRule->getNumberOfIntegrationPoints();
        for ( k = 0; k < nip; k++ ) {
            result = ielem->giveIPValue(val, iRule->getIntegrationPoint(k), IST_PrincipalDamageTensor, tStep);
            if ( result ) {
                sval = val.computeNorm();
                maxVal = max(maxVal, sval);
            }
        }

        if ( maxVal > 0.75 ) {
            elemSize = min(elemSize, 1.0);
        }

        //#endif
        ielemNodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= ielemNodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            if ( connectedElems [ jnode - 1 ] != 0 ) {
                //this -> nodalDensities.at(jnode) = min(this -> nodalDensities.at(jnode), elemSize);
                this->nodalDensities.at(jnode) += elemSize;
            } else {
                this->nodalDensities.at(jnode) = elemSize;
            }

            connectedElems [ jnode - 1 ]++;
        }
    }

    // init non-initialized nodes -> those in skip regions
    for ( i = 0; i < nnode; i++ ) {
        if ( connectedElems [ i ] == 0 ) {
            this->nodalDensities.at(i + 1) = this->giveDofManDensity(i + 1);
        } else {
            this->nodalDensities.at(i + 1) /= connectedElems [ i ];
        }
    }

    if ( refine == true ) {
        if ( this->ee->giveErrorEstimatorType() == EET_HEE ) {
            if ( ( ( HuertaErrorEstimator * ) ( this->ee ) )->giveAnalysisMode() == HuertaErrorEstimator :: HEE_linear ) {
                this->remeshingStrategy = RemeshingFromPreviousState_RS;
            }
        }
    }

    // remember time stamp
    stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}


IRResultType
HuertaRemeshingCriteria :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double coeff;
    int noRemeshFlag = 0, wErrorFlag = 0;

    IR_GIVE_FIELD(ir, this->requiredError, IFT_HuertaRemeshingCriteria_requirederror, "requirederror"); // Macro
    IR_GIVE_FIELD(ir, this->minElemSize, IFT_HuertaRemeshingCriteria_minelemsize, "minelemsize"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, noRemeshFlag, IFT_HuertaRemeshingCriteria_noremesh, "noremesh"); // Macro
    if ( noRemeshFlag != 0 ) {
        this->noRemesh = true;
    }


    IR_GIVE_OPTIONAL_FIELD(ir, wErrorFlag, IFT_HuertaRemeshingCriteria_werror, "werror"); // Macro
    if ( wErrorFlag != 0 ) {
        this->wError = true;
    }

    coeff = this->refineCoeff;
    IR_GIVE_OPTIONAL_FIELD(ir, coeff, IFT_HuertaRemeshingCriteria_refinecoeff, "refinecoeff"); // Macro
    if ( coeff > 0.0 && coeff <= 1.0 ) {
        this->refineCoeff = coeff;
    }

    return IRRT_OK;
}


double
HuertaRemeshingCriteria :: giveDofManDensity(int num)
{
    int i, isize;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    const IntArray *con;
    HuertaRemeshingCriteriaInterface *interface;
    double density;

    con = ct->giveDofManConnectivityArray(num);
    isize = con->giveSize();

    /*
     * // Minimum density
     *
     * for(i = 1; i <= isize; i++) {
     * interface = (HuertaRemeshingCriteriaInterface*)
     * domain -> giveElement(con -> at(i)) -> giveInterface(HuertaRemeshingCriteriaInterfaceType);
     * if (!interface) {
     * _error ("giveDofManDensity: element does not support HuertaRemeshingCriteriaInterface");
     * }
     * if (i==1) density = interface -> HuertaRemeshingCriteriaI_giveCharacteristicSize();
     * else density = min (density, interface -> HuertaRemeshingCriteriaI_giveCharacteristicSize());
     * }
     */

    // Average density

    density = 0.0;
    for ( i = 1; i <= isize; i++ ) {
        interface = ( HuertaRemeshingCriteriaInterface * )
                    domain->giveElement( con->at(i) )->giveInterface(HuertaRemeshingCriteriaInterfaceType);
        if ( !interface ) {
            _error("giveDofManDensity: element does not support HuertaRemeshingCriteriaInterface");
        }

        density += interface->HuertaRemeshingCriteriaI_giveCharacteristicSize();
    }

    density /= isize;

    return density;
}


void
HuertaErrorEstimator :: buildRefinedMesh(void)
{
    // Element *element;
    RefinedElement *refinedElement;
    int ielem, nelems;
    Domain *d = this->domain;

    if ( this->refinedMesh.completed == 1 ) {
        return;
    }

    nelems = d->giveNumberOfElements();
    this->refinedElementList.growTo(nelems);

    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        //  element = d -> giveElement(ielem);
        refinedElement = new RefinedElement(d, ielem, this->refineLevel);
        this->refinedElementList.put(ielem, refinedElement);
    }

    if ( refinedMesh.refineMeshGlobally(d, this->refineLevel, this->refinedElementList) != 0 ) {
        _error("buildRefinedMesh: refineMeshGlobally failed");
    }

    this->refinedMesh.completed = 1;
}



void
HuertaErrorEstimatorInterface :: setupRefinedElementProblem1D(Element *element, RefinedElement *refinedElement,
                                                              int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                              HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                                              FloatArray **corner, FloatArray &midNode,
                                                              int &localNodeId, int &localElemId, int &localBcId,
                                                              IntArray &controlNode, IntArray &controlDof,
                                                              HuertaErrorEstimator :: AnalysisMode aMode, const char *edgetype)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    IntArray *connectivity, boundary(1);
    int startNode, endNode, inode, m, pos, nd, bc, dofs;
    Node *node;
    std :: string str;

    if ( nodeId != 0 ) {
        startNode = endNode = nodeId;
    } else {
        startNode = 1;
        endNode = nodes;
    }

    dofs = element->giveDomain()->giveDefaultNodeDofIDArry().giveSize();

    if ( mode == HuertaErrorEstimatorInterface :: CountMode ) {
        for ( inode = startNode; inode <= endNode; inode++ ) {
            localElemId += ( level + 1 );

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);

            pos = 1;
            for ( m = 0; m < level + 2; m++, pos++ ) {
                nd = connectivity->at(pos);
                if ( localNodeIdArray.at(nd) == 0 ) {
                    localNodeIdArray.at(nd) = ++localNodeId;

                    // count boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))

                    bc = 0;
                    if ( nodeId == 0 ) {
                        if ( m == 0 && boundary.at(1) == 0 ) {
                            bc = 1;
                        }
                    } else {
                        if ( m == level + 1 ) {
                            bc = 1;
                        }
                    }

#ifdef EXACT_ERROR
                    if ( wholeFlag == true ) {
                        bc = 0;
                    }

#endif

                    if ( bc == 1 ) {
                        localBcId += dofs;
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: NodeMode ) {
        double x, y, z, u, du = 1.0 / ( level + 1 );
        double xc, yc, zc, xm, ym, zm;
        int idof, iload, loads, sideNumBc, bcId, bcDofId;
        IntArray sideBcDofId, dofIdArray, *loadArray;
        FloatMatrix *lcs;
        bool hasBc;
        Dof *nodeDof;

        dofIdArray = element->giveDomain()->giveDefaultNodeDofIDArry();

        for ( inode = startNode; inode <= endNode; inode++ ) {
            xc = corner [ inode - 1 ]->at(1);
            yc = corner [ inode - 1 ]->at(2);
            zc = corner [ inode - 1 ]->at(3);

            xm = midNode.at(1);
            ym = midNode.at(2);
            zm = midNode.at(3);

            node = element->giveNode(inode);

            if ( node->giveNumberOfDofs() != dofs ) {
                OOFEM_ERROR("HuertaErrorEstimatorInterface::setupRefinedElementProblem1D : Dof mismatch");
            }

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasBc = refinedElement->giveBcDofArray1D(inode, element, & sideBcDofId, sideNumBc, tStep);

            pos = 1;
            u = 0.0;
            for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                nd = connectivity->at(pos);
                if ( localNodeIdArray.at(nd) == 0 ) {
                    localNodeIdArray.at(nd) = ++localNodeId;
                    globalNodeIdArray.at(localNodeId) = nd;

                    x = xc * ( 1.0 - u ) + xm * u;
                    y = yc * ( 1.0 - u ) + ym * u;
                    z = zc * ( 1.0 - u ) + zm * u;

                    sprintf(line, "node %d coords 3 %f %f %f", localNodeId, x, y, z);
                    str = line;

                    if ( ( lcs = node->giveLocalCoordinateTriplet() ) != NULL ) {
                        sprintf( line, " lcs 6 %f %f %f %f %f %f",
                                lcs->at(1, 1), lcs->at(1, 2), lcs->at(1, 3),
                                lcs->at(2, 1), lcs->at(2, 2), lcs->at(2, 3) );
                        str += line;
                    }

                    // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                    bc = 0;
                    if ( nodeId == 0 ) {
                        if ( m == 0 && boundary.at(1) == 0 ) {
                            bc = 1;
                        }
                    } else {
                        if ( m == level + 1 ) {
                            bc = 1;
                        }
                    }

#ifdef EXACT_ERROR
                    if ( wholeFlag == true ) {
                        bc = 0;
                    }

#endif

                    // jak jsou razeny stupne volnosti v bc a ic (je-li jich jiny pocet, nez default dofs)
                    // razeni je zrejme shodne s razenim dof v dofmanageru

                    if ( bc == 1 ) {
                        if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
                            sprintf(line, " bc %d", dofs);
                            str += line;
                            for ( idof = 0; idof < dofs; idof++ ) {
                                sprintf(line, " %d", ++localBcId);
                                str += line;
                            }
                        }
                    } else {
                        if ( hasBc == true ) {
                            // it is necessary to reproduce bc from coarse mesh

                            if ( m == 0 ) {       // at node
                                sprintf(line, " bc %d", dofs);
                                str += line;

                                for ( idof = 1; idof <= dofs; idof++ ) {
                                    bcDofId = 0;
                                    nodeDof = node->giveDof(idof);
                                    if ( nodeDof->hasBc(tStep) != 0 ) {
                                        bcDofId = nodeDof->giveBcId();
                                    }

                                    sprintf(line, " %d", bcDofId);
                                    str += line;
                                }

                                // copy node load

                                if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                    sprintf(line, " load %d", loads);
                                    str += line;
                                    for ( iload = 1; iload <= loads; iload++ ) {
                                        sprintf( line, " %d", loadArray->at(iload) );
                                        str += line;
                                    }
                                }
                            } else {
                                if ( sideNumBc != 0 ) {
                                    sprintf(line, " bc %d", dofs);
                                    str += line;

                                    // I rely on the fact that bc dofs to be reproduced are ordered with respect to the dof ordering of the corner node

                                    bcId = 1;
                                    for ( idof = 1; idof <= dofs; idof++ ) {
                                        bcDofId = 0;
                                        if ( bcId <= sideNumBc ) {
                                            nodeDof = node->giveDof( sideBcDofId.at(bcId) );
                                            if ( nodeDof->giveDofID() == dofIdArray.at(idof) ) {
                                                bcDofId = nodeDof->giveBcId();
                                                bcId++;
                                            }
                                        }

                                        sprintf(line, " %d", bcDofId);
                                        str += line;
                                    }
                                }
                            }
                        } else {
                            // copy node load

                            if ( m == 0 ) {
                                if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                    sprintf(line, " load %d", loads);
                                    str += line;
                                    for ( iload = 1; iload <= loads; iload++ ) {
                                        sprintf( line, " %d", loadArray->at(iload) );
                                        str += line;
                                    }
                                }
                            }
                        }
                    }

                    refinedReader.appendInputString(str);
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: ElemMode ) {
        int mat, csect, iload, loads, material;
        int nd1, nd2;
        IntArray *loadArray, boundaryLoadArray;
        bool hasLoad;

        mat = element->giveMaterial()->giveNumber();
        csect = element->giveCrossSection()->giveNumber();

        for ( inode = startNode; inode <= endNode; inode++ ) {
            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasLoad = refinedElement->giveBoundaryLoadArray1D(inode, element, boundaryLoadArray);

            for ( m = 0; m < level + 1; m++ ) {
                localElemId++;

                nd = m + 1;

                nd1 = localNodeIdArray.at( connectivity->at(nd) );
                nd2 = localNodeIdArray.at( connectivity->at(nd + 1) );

                material = mat;
                if ( impMat != 0 && impMat == mat ) {
                    int i;
                    //     InputRecord &ir();
                    FloatArray coordinates1, coordinates2;
                    const char *__proc = "sesetupRefinedElementProblem1D"; // Required by IR_GIVE_FIELD macro
                    IRResultType result;                    // Required by IR_GIVE_FIELD macro

                    refinedReader.seek(nd1 + 6);
                    InputRecord *ir = refinedReader.giveInputRecord(DataReader :: IR_dofmanRec, nd1);
                    IR_GIVE_FIELD(ir, coordinates1, IFT_HuertaErrorEstimatorInterface_coords, "coords"); // Macro

                    refinedReader.seek(nd2 + 6);
                    ir = refinedReader.giveInputRecord(DataReader :: IR_dofmanRec, nd2);
                    IR_GIVE_FIELD(ir, coordinates2, IFT_HuertaErrorEstimatorInterface_coords, "coords"); // Macro

                    material = impMat;
                    for ( i = 1; i <= impPos.giveSize(); i++ ) {
                        if ( impPos.at(i) < min( coordinates1.at(i), coordinates2.at(i) ) ) {
                            material = perMat;
                            break;
                        }

                        if ( impPos.at(i) > max( coordinates1.at(i), coordinates2.at(i) ) ) {
                            material = perMat;
                            break;
                        }
                    }
                }

                sprintf(line, "%s %d nodes 2 %d %d mat %d crossSect %d",
                        edgetype, localElemId, nd1, nd2, material, csect);
                str = line;

                // copy body and boundary loads

                if ( ( loads = ( loadArray = element->giveBodyLoadArray() )->giveSize() ) != 0 ) {
                    sprintf(line, " bodyLoads %d", loads);
                    str += line;
                    for ( iload = 1; iload <= loads; iload++ ) {
                        sprintf( line, " %d", loadArray->at(iload) );
                        str += line;
                    }
                }

                if ( hasLoad == true ) {
                    if ( ( loads = boundaryLoadArray.giveSize() ) != 0 ) {
                        sprintf(line, " boundaryLoads %d", loads);
                        str += line;
                        for ( iload = 1; iload <= loads; iload++ ) {
                            sprintf( line, " %d", boundaryLoadArray.at(iload) );
                            str += line;
                        }
                    }
                }

                refinedReader.appendInputString(str);
            }
        }
    }

    if ( mode == HuertaErrorEstimatorInterface :: BCMode ) {
        if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
            int idof;
            double u, du = 1.0 / ( level + 1 );
            double xc, yc, zc, xm, ym, zm;
            MaterialMode mode;
            GaussPoint *gp;
            FloatArray globCoord(3), * locCoord;
            FloatMatrix Nmatrix;
            FloatArray uCoarse, uFine;

            mode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();

            // create a fictitious integration point
            locCoord = new FloatArray;
            IntegrationRule ir(1, element);
            //gp = new GaussPoint(element, 1, locCoord, 1.0, mode);
            gp = new GaussPoint( &ir, 1, locCoord, 1.0, mode);

            for ( inode = startNode; inode <= endNode; inode++ ) {
                xc = corner [ inode - 1 ]->at(1);
                yc = corner [ inode - 1 ]->at(2);
                zc = corner [ inode - 1 ]->at(3);

                xm = midNode.at(1);
                ym = midNode.at(2);
                zm = midNode.at(3);

                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                // get corner displacements
                element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, uCoarse);

                pos = 1;
                u = 0.0;
                for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                    nd = connectivity->at(pos);
                    if ( localNodeIdArray.at(nd) > 0 ) {
                        localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                        // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                        bc = 0;
                        if ( nodeId == 0 ) {
                            if ( m == 0 && boundary.at(1) == 0 ) {
                                bc = 1;
                            }
                        } else {
                            if ( m == level + 1 ) {
                                bc = 1;
                            }
                        }

#ifdef EXACT_ERROR
                        if ( wholeFlag == true ) {
                            bc = 0;
                        }

#endif

                        if ( bc == 1 ) {
                            globCoord.at(1) = xc * ( 1.0 - u ) + xm * u;
                            globCoord.at(2) = yc * ( 1.0 - u ) + ym * u;
                            globCoord.at(3) = zc * ( 1.0 - u ) + zm * u;

                            // this effectively rewrites the local coordinates of the fictitious integration point
                            this->HuertaErrorEstimatorI_computeLocalCoords(* locCoord, globCoord);
                            // get N matrix at the fictitious integration point
                            this->HuertaErrorEstimatorI_computeNmatrixAt(gp, Nmatrix);
                            // get displacement at the fictitious integration point
                            uFine.beProductOf(Nmatrix, uCoarse);

                            // first loadtime function must be constant 1.0
                            for ( idof = 1; idof <= dofs; idof++ ) {
                                sprintf( line, "BoundaryCondition %d loadTimeFunction %d prescribedvalue %e",
                                        ++localBcId, 1, uFine.at(idof) );
                                refinedReader.appendInputString(line);
                            }
                        }
                    }
                }
            }

            delete gp;
        } else {
            int idof;

            for ( inode = startNode; inode <= endNode; inode++ ) {
                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                pos = 1;
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    nd = connectivity->at(pos);
                    if ( localNodeIdArray.at(nd) > 0 ) {
                        localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                        // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                        bc = 0;
                        if ( nodeId == 0 ) {
                            if ( m == 0 && boundary.at(1) == 0 ) {
                                bc = 1;
                            }
                        } else {
                            if ( m == level + 1 ) {
                                bc = 1;
                            }
                        }

#ifdef EXACT_ERROR
                        if ( wholeFlag == true ) {
                            bc = 0;
                        }

#endif

                        if ( bc == 1 ) {
                            for ( idof = 1; idof <= dofs; idof++ ) {
                                localBcId++;
                                controlNode.at(localBcId) = -localNodeIdArray.at(nd);
                                controlDof.at(localBcId) = idof;
                            }
                        }
                    }
                }
            }
        }

        return;
    }
}



void
HuertaErrorEstimatorInterface :: setupRefinedElementProblem2D(Element *element, RefinedElement *refinedElement,
                                                              int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                              HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                                              FloatArray **corner, FloatArray *midSide, FloatArray &midNode,
                                                              int &localNodeId, int &localElemId, int &localBcId,
                                                              IntArray &controlNode, IntArray &controlDof,
                                                              HuertaErrorEstimator :: AnalysisMode aMode, const char *quadtype)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    IntArray *connectivity, boundary(2);
    int startNode, endNode, inode, n, m, pos, nd, bc, dofs;
    Node *node;
    std :: string str;

    if ( nodeId != 0 ) {
        startNode = endNode = nodeId;
    } else {
        startNode = 1;
        endNode = nodes;
    }

    dofs = element->giveDomain()->giveDefaultNodeDofIDArry().giveSize();

    if ( mode == HuertaErrorEstimatorInterface :: CountMode ) {
        for ( inode = startNode; inode <= endNode; inode++ ) {
            localElemId += ( level + 1 ) * ( level + 1 );

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);

            pos = 1;
            for ( n = 0; n < level + 2; n++ ) {
                for ( m = 0; m < level + 2; m++, pos++ ) {
                    nd = connectivity->at(pos);
                    if ( localNodeIdArray.at(nd) == 0 ) {
                        localNodeIdArray.at(nd) = ++localNodeId;

                        // count boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                        bc = 0;
                        if ( nodeId == 0 ) {
                            if ( m == 0 && boundary.at(1) == 0 ) {
                                bc = 1;
                            }

                            if ( n == 0 && boundary.at(2) == 0 ) {
                                bc = 1;
                            }
                        } else {
                            if ( n == level + 1 || m == level + 1 ) {
                                bc = 1;
                            }

                            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                            if ( bc == 0 ) {
                                if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                    if ( m == 0 && boundary.at(1) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( n == 0 && boundary.at(2) == 0 ) {
                                        bc = 1;
                                    }
                                }
                            }

#endif
                        }

#ifdef EXACT_ERROR
                        if ( wholeFlag == true ) {
                            bc = 0;
                        }

#endif

                        if ( bc == 1 ) {
                            localBcId += dofs;
                        }
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: NodeMode ) {
        int iside, s1, s2, idof, index, iload, loads;
        double x, y, z, u, v, du = 1.0 / ( level + 1 ), dv = 1.0 / ( level + 1 );
        double xc, yc, zc, xs1, ys1, zs1, xs2, ys2, zs2, xm, ym, zm;
        int bcId, bcDofId;
        AList< IntArray >sideBcDofIdList;
        IntArray *sideBcDofId, sideNumBc(2), dofIdArray, * loadArray;
        bool hasBc;
        Dof *nodeDof;
        FloatMatrix *lcs;

        dofIdArray = element->giveDomain()->giveDefaultNodeDofIDArry();

        sideBcDofIdList.growTo(2);
        for ( iside = 1; iside <= 2; iside++ ) {
            sideBcDofId = new IntArray;
            sideBcDofIdList.put(iside, sideBcDofId);
        }

        for ( inode = startNode; inode <= endNode; inode++ ) {
            s1 = inode;
            if ( ( s2 = inode - 1 ) == 0 ) {
                s2 = nodes;
            }

            xc = corner [ inode - 1 ]->at(1);
            yc = corner [ inode - 1 ]->at(2);
            zc = corner [ inode - 1 ]->at(3);

            xs1 = midSide [ s1 - 1 ].at(1);
            ys1 = midSide [ s1 - 1 ].at(2);
            zs1 = midSide [ s1 - 1 ].at(3);

            xs2 = midSide [ s2 - 1 ].at(1);
            ys2 = midSide [ s2 - 1 ].at(2);
            zs2 = midSide [ s2 - 1 ].at(3);

            xm = midNode.at(1);
            ym = midNode.at(2);
            zm = midNode.at(3);

            node = element->giveNode(inode);

            if ( node->giveNumberOfDofs() != dofs ) {
                OOFEM_ERROR("HuertaErrorEstimatorInterface::setupRefinedElementProblem2D : Dof mismatch");
            }

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasBc = refinedElement->giveBcDofArray2D(inode, element, sideBcDofIdList, sideNumBc, tStep);

            pos = 1;
            v = 0.0;
            for ( n = 0; n < level + 2; n++, v += dv ) {
                u = 0.0;
                for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                    nd = connectivity->at(pos);
                    if ( localNodeIdArray.at(nd) == 0 ) {
                        localNodeIdArray.at(nd) = ++localNodeId;
                        globalNodeIdArray.at(localNodeId) = nd;

                        x = ( xc * ( 1.0 - u ) + xs1 * u ) * ( 1.0 - v ) + ( xs2 * ( 1.0 - u ) + xm * u ) * v;
                        y = ( yc * ( 1.0 - u ) + ys1 * u ) * ( 1.0 - v ) + ( ys2 * ( 1.0 - u ) + ym * u ) * v;
                        z = ( zc * ( 1.0 - u ) + zs1 * u ) * ( 1.0 - v ) + ( zs2 * ( 1.0 - u ) + zm * u ) * v;

                        sprintf(line, "node %d coords 3 %f %f %f", localNodeId, x, y, z);
                        str = line;

                        if ( ( lcs = node->giveLocalCoordinateTriplet() ) != NULL ) {
                            sprintf( line, " lcs 6 %f %f %f %f %f %f",
                                    lcs->at(1, 1), lcs->at(1, 2), lcs->at(1, 3),
                                    lcs->at(2, 1), lcs->at(2, 2), lcs->at(2, 3) );
                            str += line;
                        }

                        // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                        bc = 0;
                        if ( nodeId == 0 ) {
                            if ( m == 0 && boundary.at(1) == 0 ) {
                                bc = 1;
                            }

                            if ( n == 0 && boundary.at(2) == 0 ) {
                                bc = 1;
                            }
                        } else {
                            if ( n == level + 1 || m == level + 1 ) {
                                bc = 1;
                            }

                            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                            if ( bc == 0 ) {
                                if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                    if ( m == 0 && boundary.at(1) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( n == 0 && boundary.at(2) == 0 ) {
                                        bc = 1;
                                    }
                                }
                            }

#endif
                        }

#ifdef EXACT_ERROR
                        if ( wholeFlag == true ) {
                            bc = 0;
                        }

#endif

                        // jak jsou razeny stupne volnosti v bc a ic (je-li jich jiny pocet, nez default dofs)
                        // razeni je zrejme shodne s razenim dof v dofmanageru

                        if ( bc == 1 ) {
                            if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
                                sprintf(line, " bc %d", dofs);
                                str += line;
                                for ( idof = 0; idof < dofs; idof++ ) {
                                    sprintf(line, " %d", ++localBcId);
                                    str += line;
                                }
                            }
                        } else {
                            if ( hasBc == true && ( m == 0 || n == 0 ) ) {
                                // it is necessary to reproduce bc from coarse mesh

                                if ( m == 0 && n == 0 ) {    // at node
                                    sprintf(line, " bc %d", dofs);
                                    str += line;

                                    for ( idof = 1; idof <= dofs; idof++ ) {
                                        bcDofId = 0;
                                        nodeDof = node->giveDof(idof);
                                        if ( nodeDof->hasBc(tStep) != 0 ) {
                                            bcDofId = nodeDof->giveBcId();
                                        }

                                        sprintf(line, " %d", bcDofId);
                                        str += line;
                                    }

                                    // copy node load

                                    if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                        sprintf(line, " load %d", loads);
                                        str += line;
                                        for ( iload = 1; iload <= loads; iload++ ) {
                                            sprintf( line, " %d", loadArray->at(iload) );
                                            str += line;
                                        }
                                    }
                                } else {
                                    index = 0;
                                    if ( m == 0 && sideNumBc.at(1) != 0 ) {
                                        index = 1;          // along next side
                                    }

                                    if ( n == 0 && sideNumBc.at(2) != 0 ) {
                                        index = 2;          // along prev side
                                    }

                                    if ( index != 0 ) {
                                        sprintf(line, " bc %d", dofs);
                                        str += line;

                                        // I rely on the fact that bc dofs to be reproduced are ordered with respect to the dof ordering of the corner node

                                        sideBcDofId = sideBcDofIdList.at(index);
                                        bcId = 1;
                                        for ( idof = 1; idof <= dofs; idof++ ) {
                                            bcDofId = 0;
                                            if ( bcId <= sideNumBc.at(index) ) {
                                                nodeDof = node->giveDof( sideBcDofId->at(bcId) );
                                                if ( nodeDof->giveDofID() == dofIdArray.at(idof) ) {
                                                    bcDofId = nodeDof->giveBcId();
                                                    bcId++;
                                                }
                                            }

                                            sprintf(line, " %d", bcDofId);
                                            str += line;
                                        }
                                    }
                                }
                            } else {
                                // copy node load

                                if ( m == 0 && n == 0 ) {
                                    if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                        sprintf(line, " load %d", loads);
                                        str += line;
                                        for ( iload = 1; iload <= loads; iload++ ) {
                                            sprintf( line, " %d", loadArray->at(iload) );
                                            str += line;
                                        }
                                    }
                                }
                            }
                        }

                        refinedReader.appendInputString(str);
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: ElemMode ) {
        int mat, csect, iside, iload, loads, index;
        int nd1, nd2, nd3, nd4;
        IntArray *loadArray;
        AList< IntArray >boundaryLoadList;
        bool hasLoad;

        mat = element->giveMaterial()->giveNumber();
        csect = element->giveCrossSection()->giveNumber();

        boundaryLoadList.growTo(2);
        for ( iside = 1; iside <= 2; iside++ ) {
            loadArray = new IntArray;
            boundaryLoadList.put(iside, loadArray);
        }

        for ( inode = startNode; inode <= endNode; inode++ ) {
            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasLoad = refinedElement->giveBoundaryLoadArray2D(inode, element, boundaryLoadList);

            for ( n = 0; n < level + 1; n++ ) {
                for ( m = 0; m < level + 1; m++ ) {
                    localElemId++;

                    nd = n * ( level + 2 ) + m + 1;

                    nd1 = localNodeIdArray.at( connectivity->at(nd) );
                    nd2 = localNodeIdArray.at( connectivity->at(nd + 1) );
                    nd3 = localNodeIdArray.at( connectivity->at(nd + level + 3) );
                    nd4 = localNodeIdArray.at( connectivity->at(nd + level + 2) );

                    sprintf(line, "%s %d nodes 4 %d %d %d %d mat %d crossSect %d",
                            quadtype, localElemId, nd1, nd2, nd3, nd4, mat, csect);
                    str = line;

                    // copy body and boundary loads

                    if ( ( loads = ( loadArray = element->giveBodyLoadArray() )->giveSize() ) != 0 ) {
                        sprintf(line, " bodyLoads %d", loads);
                        str += line;
                        for ( iload = 1; iload <= loads; iload++ ) {
                            sprintf( line, " %d", loadArray->at(iload) );
                            str += line;
                        }
                    }

                    // boundary load is not copied on non-boundary sides

                    if ( hasLoad == true && ( m == 0 || n == 0 ) ) {
                        if ( m == 0 && n == 0 ) {
                            loads = 0;
                            for ( iside = 1; iside <= 2; iside++ ) {
                                if ( boundary.at(iside) == 0 ) {
                                    continue;
                                }

                                loads += boundaryLoadList.at(iside)->giveSize();
                            }

                            if ( loads != 0 ) {
                                sprintf(line, " boundaryLoads %d", loads);
                                str += line;
                                for ( iside = 1; iside <= 2; iside++ ) {
                                    if ( boundary.at(iside) == 0 ) {
                                        continue;
                                    }

                                    if ( ( loads = ( loadArray = boundaryLoadList.at(iside) )->giveSize() ) == 0 ) {
                                        continue;
                                    }

                                    for ( iload = 1; iload <= loads; iload++ ) {
                                        sprintf( line, " %d", loadArray->at(iload) );
                                        str += line;
                                    }
                                }
                            }
                        } else {
                            index = 0;
                            if ( m == 0 && boundary.at(1) != 0 && boundaryLoadList.at(1)->giveSize() != 0 ) {
                                index = 1;
                            }

                            if ( n == 0 && boundary.at(2) != 0 && boundaryLoadList.at(2)->giveSize() != 0 ) {
                                index = 2;
                            }

                            if ( index != 0 ) {
                                loads = ( loadArray = boundaryLoadList.at(index) )->giveSize();
                                sprintf(line, " boundaryLoads %d", loads);
                                str += line;
                                for ( iload = 1; iload <= loads; iload++ ) {
                                    sprintf( line, " %d", loadArray->at(iload) );
                                    str += line;
                                }
                            }
                        }
                    }

                    refinedReader.appendInputString(str);
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: BCMode ) {
        if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
            int s1, s2, idof;
            double u, v, du = 1.0 / ( level + 1 ), dv = 1.0 / ( level + 1 );
            double xc, yc, zc, xs1, ys1, zs1, xs2, ys2, zs2, xm, ym, zm;
            MaterialMode mode;
            GaussPoint *gp;
            FloatArray globCoord(3), * locCoord;
            FloatMatrix Nmatrix;
            FloatArray uCoarse, uFine;

            mode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();

            // create a fictitious integration point
            locCoord = new FloatArray;
            IntegrationRule ir(0, element);
            gp = new GaussPoint( &ir, 1, locCoord, 1.0, mode);

            for ( inode = startNode; inode <= endNode; inode++ ) {
                s1 = inode;
                if ( ( s2 = inode - 1 ) == 0 ) {
                    s2 = nodes;
                }

                xc = corner [ inode - 1 ]->at(1);
                yc = corner [ inode - 1 ]->at(2);
                zc = corner [ inode - 1 ]->at(3);

                xs1 = midSide [ s1 - 1 ].at(1);
                ys1 = midSide [ s1 - 1 ].at(2);
                zs1 = midSide [ s1 - 1 ].at(3);

                xs2 = midSide [ s2 - 1 ].at(1);
                ys2 = midSide [ s2 - 1 ].at(2);
                zs2 = midSide [ s2 - 1 ].at(3);

                xm = midNode.at(1);
                ym = midNode.at(2);
                zm = midNode.at(3);

                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                // get corner displacements
                element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, uCoarse);

                pos = 1;
                v = 0.0;
                for ( n = 0; n < level + 2; n++, v += dv ) {
                    u = 0.0;
                    for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                        nd = connectivity->at(pos);
                        if ( localNodeIdArray.at(nd) > 0 ) {
                            localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                            // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                            bc = 0;
                            if ( nodeId == 0 ) {
                                if ( m == 0 && boundary.at(1) == 0 ) {
                                    bc = 1;
                                }

                                if ( n == 0 && boundary.at(2) == 0 ) {
                                    bc = 1;
                                }
                            } else {
                                if ( n == level + 1 || m == level + 1 ) {
                                    bc = 1;
                                }

                                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                if ( bc == 0 ) {
                                    if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                        if ( m == 0 && boundary.at(1) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( n == 0 && boundary.at(2) == 0 ) {
                                            bc = 1;
                                        }
                                    }
                                }

#endif
                            }

#ifdef EXACT_ERROR
                            if ( wholeFlag == true ) {
                                bc = 0;
                            }

#endif

                            if ( bc == 1 ) {
                                globCoord.at(1) = ( xc * ( 1.0 - u ) + xs1 * u ) * ( 1.0 - v ) + ( xs2 * ( 1.0 - u ) + xm * u ) * v;
                                globCoord.at(2) = ( yc * ( 1.0 - u ) + ys1 * u ) * ( 1.0 - v ) + ( ys2 * ( 1.0 - u ) + ym * u ) * v;
                                globCoord.at(3) = ( zc * ( 1.0 - u ) + zs1 * u ) * ( 1.0 - v ) + ( zs2 * ( 1.0 - u ) + zm * u ) * v;

                                // this effectively rewrites the local coordinates of the fictitious integration point
                                this->HuertaErrorEstimatorI_computeLocalCoords(* locCoord, globCoord);
                                // get N matrix at the fictitious integration point
                                this->HuertaErrorEstimatorI_computeNmatrixAt(gp, Nmatrix);
                                // get displacement at the fictitious integration point
                                uFine.beProductOf(Nmatrix, uCoarse);

                                // first loadtime function must be constant 1.0
                                for ( idof = 1; idof <= dofs; idof++ ) {
                                    sprintf( line, "BoundaryCondition %d loadTimeFunction %d prescribedvalue %e",
                                            ++localBcId, 1, uFine.at(idof) );
                                    refinedReader.appendInputString(line);
                                }
                            }
                        }
                    }
                }
            }

            delete gp;
        } else {
            int idof;

            for ( inode = startNode; inode <= endNode; inode++ ) {
                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                pos = 1;
                for ( n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        nd = connectivity->at(pos);
                        if ( localNodeIdArray.at(nd) > 0 ) {
                            localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                            // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                            bc = 0;
                            if ( nodeId == 0 ) {
                                if ( m == 0 && boundary.at(1) == 0 ) {
                                    bc = 1;
                                }

                                if ( n == 0 && boundary.at(2) == 0 ) {
                                    bc = 1;
                                }
                            } else {
                                if ( n == level + 1 || m == level + 1 ) {
                                    bc = 1;
                                }

                                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                if ( bc == 0 ) {
                                    if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                        if ( m == 0 && boundary.at(1) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( n == 0 && boundary.at(2) == 0 ) {
                                            bc = 1;
                                        }
                                    }
                                }

#endif
                            }

#ifdef EXACT_ERROR
                            if ( wholeFlag == true ) {
                                bc = 0;
                            }

#endif

                            if ( bc == 1 ) {
                                for ( idof = 1; idof <= dofs; idof++ ) {
                                    localBcId++;
                                    controlNode.at(localBcId) = -localNodeIdArray.at(nd);
                                    controlDof.at(localBcId) = idof;
                                }
                            }
                        }
                    }
                }
            }
        }

        return;
    }
}



void
HuertaErrorEstimatorInterface :: setupRefinedElementProblem3D(Element *element, RefinedElement *refinedElement,
                                                              int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                              HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                                              FloatArray **corner, FloatArray *midSide, FloatArray *midFace, FloatArray &midNode,
                                                              int &localNodeId, int &localElemId, int &localBcId,
                                                              int hexaSideNode [ 1 ] [ 3 ], int hexaFaceNode [ 1 ] [ 3 ],
                                                              IntArray &controlNode, IntArray &controlDof,
                                                              HuertaErrorEstimator :: AnalysisMode aMode, const char *hexatype)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    IntArray *connectivity, boundary(3);
    int startNode, endNode, inode, k, n, m, pos, nd, bc, dofs;
    Node *node;
    std :: string str;

    if ( nodeId != 0 ) {
        startNode = endNode = nodeId;
    } else {
        startNode = 1;
        endNode = nodes;
    }

    dofs = element->giveDomain()->giveDefaultNodeDofIDArry().giveSize();

    if ( mode == HuertaErrorEstimatorInterface :: CountMode ) {
#ifdef DEBUG
        if ( nodes / 2 * 2 != nodes ) {
            abort();
        }

        //   _error ("setupRefinedElementProblem3D: unexpected situation");

        /* number of internal quad faces = nodes * 3 / 2;
         *
         * at each node there are 3 internal quad faces and each internal quad face is shared by two internal hexas;
         * it is clear that this does not work for element with odd number of nodes such as pyramid;
         * however, pyramid is not decomposable into hexas and therefore should not provide this interface */
#endif

        for ( inode = startNode; inode <= endNode; inode++ ) {
            localElemId += ( level + 1 ) * ( level + 1 ) * ( level + 1 );

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);

            pos = 1;
            for ( k = 0; k < level + 2; k++ ) {
                for ( n = 0; n < level + 2; n++ ) {
                    for ( m = 0; m < level + 2; m++, pos++ ) {
                        nd = connectivity->at(pos);
                        if ( localNodeIdArray.at(nd) == 0 ) {
                            localNodeIdArray.at(nd) = ++localNodeId;

                            // count boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                            bc = 0;
                            if ( nodeId == 0 ) {
                                if ( m == 0 && boundary.at(1) == 0 ) {
                                    bc = 1;
                                }

                                if ( n == 0 && boundary.at(2) == 0 ) {
                                    bc = 1;
                                }

                                if ( k == 0 && boundary.at(3) == 0 ) {
                                    bc = 1;
                                }
                            } else {
                                if ( k == level + 1 || n == level + 1 || m == level + 1 ) {
                                    bc = 1;
                                }

                                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                if ( bc == 0 ) {
                                    if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                        if ( m == 0 && boundary.at(1) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( n == 0 && boundary.at(2) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( k == 0 && boundary.at(3) == 0 ) {
                                            bc = 1;
                                        }
                                    }
                                }

#endif
                            }

#ifdef EXACT_ERROR
                            if ( wholeFlag == true ) {
                                bc = 0;
                            }

#endif

                            if ( bc == 1 ) {
                                localBcId += dofs;
                            }
                        }
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: NodeMode ) {
        int iside, iface, s1, s2, s3, f1, f2, f3, idof, index, iload, loads;
        double x, y, z, u, v, w, du = 1.0 / ( level + 1 ), dv = 1.0 / ( level + 1 ), dw = 1.0 / ( level + 1 );
        double xc, yc, zc, xm, ym, zm;
        double xs1, ys1, zs1, xs2, ys2, zs2, xs3, ys3, zs3;
        double xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3;
        int bcId, bcDofId;
        AList< IntArray >sideBcDofIdList, faceBcDofIdList;
        IntArray *sideBcDofId, *faceBcDofId, sideNumBc(3), faceNumBc(3), dofIdArray, * loadArray;
        bool hasBc;
        Dof *nodeDof;
        FloatMatrix *lcs;

        dofIdArray = element->giveDomain()->giveDefaultNodeDofIDArry();

        sideBcDofIdList.growTo(3);
        for ( iside = 1; iside <= 3; iside++ ) {
            sideBcDofId = new IntArray;
            sideBcDofIdList.put(iside, sideBcDofId);
        }

        faceBcDofIdList.growTo(3);
        for ( iface = 1; iface <= 3; iface++ ) {
            faceBcDofId = new IntArray;
            faceBcDofIdList.put(iface, faceBcDofId);
        }

        for ( inode = startNode; inode <= endNode; inode++ ) {
            s1 = hexaSideNode [ inode - 1 ] [ 0 ];
            s2 = hexaSideNode [ inode - 1 ] [ 1 ];
            s3 = hexaSideNode [ inode - 1 ] [ 2 ];
            f1 = hexaFaceNode [ inode - 1 ] [ 0 ];
            f2 = hexaFaceNode [ inode - 1 ] [ 1 ];
            f3 = hexaFaceNode [ inode - 1 ] [ 2 ];

            xc = corner [ inode - 1 ]->at(1);
            yc = corner [ inode - 1 ]->at(2);
            zc = corner [ inode - 1 ]->at(3);

            xs1 = midSide [ s1 - 1 ].at(1);
            ys1 = midSide [ s1 - 1 ].at(2);
            zs1 = midSide [ s1 - 1 ].at(3);

            xs2 = midSide [ s2 - 1 ].at(1);
            ys2 = midSide [ s2 - 1 ].at(2);
            zs2 = midSide [ s2 - 1 ].at(3);

            xs3 = midSide [ s3 - 1 ].at(1);
            ys3 = midSide [ s3 - 1 ].at(2);
            zs3 = midSide [ s3 - 1 ].at(3);

            xf1 = midFace [ f1 - 1 ].at(1);
            yf1 = midFace [ f1 - 1 ].at(2);
            zf1 = midFace [ f1 - 1 ].at(3);

            xf2 = midFace [ f2 - 1 ].at(1);
            yf2 = midFace [ f2 - 1 ].at(2);
            zf2 = midFace [ f2 - 1 ].at(3);

            xf3 = midFace [ f3 - 1 ].at(1);
            yf3 = midFace [ f3 - 1 ].at(2);
            zf3 = midFace [ f3 - 1 ].at(3);

            xm = midNode.at(1);
            ym = midNode.at(2);
            zm = midNode.at(3);

            node = element->giveNode(inode);

            if ( node->giveNumberOfDofs() != dofs ) {
                OOFEM_ERROR("HuertaErrorEstimatorInterface::setupRefinedElementProblem3D : Dof mismatch");
            }

            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasBc = refinedElement->giveBcDofArray3D(inode, element, sideBcDofIdList, sideNumBc,
                                                     faceBcDofIdList, faceNumBc, tStep);
            pos = 1;
            w = 0.0;
            for ( k = 0; k < level + 2; k++, w += dw ) {
                v = 0.0;
                for ( n = 0; n < level + 2; n++, v += dv ) {
                    u = 0.0;
                    for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                        nd = connectivity->at(pos);
                        if ( localNodeIdArray.at(nd) == 0 ) {
                            localNodeIdArray.at(nd) = ++localNodeId;
                            globalNodeIdArray.at(localNodeId) = nd;

                            x = ( ( xc * ( 1.0 - u ) + xs1 * u ) * ( 1.0 - v ) + ( xs2 * ( 1.0 - u ) + xf1 * u ) * v ) * ( 1.0 - w )
                                + ( ( xs3 * ( 1.0 - u ) + xf2 * u ) * ( 1.0 - v ) + ( xf3 * ( 1.0 - u ) + xm * u ) * v ) * w;
                            y = ( ( yc * ( 1.0 - u ) + ys1 * u ) * ( 1.0 - v ) + ( ys2 * ( 1.0 - u ) + yf1 * u ) * v ) * ( 1.0 - w )
                                + ( ( ys3 * ( 1.0 - u ) + yf2 * u ) * ( 1.0 - v ) + ( yf3 * ( 1.0 - u ) + ym * u ) * v ) * w;
                            z = ( ( zc * ( 1.0 - u ) + zs1 * u ) * ( 1.0 - v ) + ( zs2 * ( 1.0 - u ) + zf1 * u ) * v ) * ( 1.0 - w )
                                + ( ( zs3 * ( 1.0 - u ) + zf2 * u ) * ( 1.0 - v ) + ( zf3 * ( 1.0 - u ) + zm * u ) * v ) * w;

                            sprintf(line, "node %d coords 3 %f %f %f", localNodeId, x, y, z);
                            str = line;

                            if ( ( lcs = node->giveLocalCoordinateTriplet() ) != NULL ) {
                                sprintf( line, " lcs 6 %f %f %f %f %f %f",
                                        lcs->at(1, 1), lcs->at(1, 2), lcs->at(1, 3),
                                        lcs->at(2, 1), lcs->at(2, 2), lcs->at(2, 3) );
                                str += line;
                            }

                            // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                            bc = 0;
                            if ( nodeId == 0 ) {
                                if ( m == 0 && boundary.at(1) == 0 ) {
                                    bc = 1;
                                }

                                if ( n == 0 && boundary.at(2) == 0 ) {
                                    bc = 1;
                                }

                                if ( k == 0 && boundary.at(3) == 0 ) {
                                    bc = 1;
                                }
                            } else {
                                if ( k == level + 1 || n == level + 1 || m == level + 1 ) {
                                    bc = 1;
                                }

                                /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                if ( bc == 0 ) {
                                    if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                        if ( m == 0 && boundary.at(1) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( n == 0 && boundary.at(2) == 0 ) {
                                            bc = 1;
                                        }

                                        if ( k == 0 && boundary.at(3) == 0 ) {
                                            bc = 1;
                                        }
                                    }
                                }

#endif
                            }

#ifdef EXACT_ERROR
                            if ( wholeFlag == true ) {
                                bc = 0;
                            }

#endif

                            // jak jsou razeny stupne volnosti v bc a ic (je-li jich jiny pocet, nez default dofs)
                            // razeni je zrejme shodne s razenim dof v dofmanageru

                            if ( bc == 1 ) {
                                if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
                                    sprintf(line, " bc %d", dofs);
                                    str += line;
                                    for ( idof = 0; idof < dofs; idof++ ) {
                                        sprintf(line, " %d", ++localBcId);
                                        str += line;
                                    }
                                }
                            } else {
                                if ( hasBc == true && ( m == 0 || n == 0 || k == 0 ) ) {
                                    // it is necessary to reproduce bc from coarse mesh

                                    if ( m == 0 && n == 0 && k == 0 ) { // at node
                                        sprintf(line, " bc %d", dofs);
                                        str += line;

                                        for ( idof = 1; idof <= dofs; idof++ ) {
                                            bcDofId = 0;
                                            nodeDof = node->giveDof(idof);
                                            if ( nodeDof->hasBc(tStep) != 0 ) {
                                                bcDofId = nodeDof->giveBcId();
                                            }

                                            sprintf(line, " %d", bcDofId);
                                            str += line;
                                        }

                                        // copy node load

                                        if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                            sprintf(line, " load %d", loads);
                                            str += line;
                                            for ( iload = 1; iload <= loads; iload++ ) {
                                                sprintf( line, " %d", loadArray->at(iload) );
                                                str += line;
                                            }
                                        }
                                    } else {
                                        if ( ( m == 0 && n == 0 ) || ( m == 0 && k == 0 ) || ( n == 0 && k == 0 ) ) {
                                            index = 0;
                                            if ( n == 0 && k == 0 && sideNumBc.at(1) != 0 ) {
                                                index = 1;
                                            }

                                            if ( m == 0 && k == 0 && sideNumBc.at(2) != 0 ) {
                                                index = 2;
                                            }

                                            if ( m == 0 && n == 0 && sideNumBc.at(3) != 0 ) {
                                                index = 3;
                                            }

                                            if ( index != 0 ) {
                                                sprintf(line, " bc %d", dofs);
                                                str += line;

                                                // I rely on the fact that bc dofs to be reproduced are ordered with respect to the dof ordering of the corner node

                                                sideBcDofId = sideBcDofIdList.at(index);
                                                bcId = 1;
                                                for ( idof = 1; idof <= dofs; idof++ ) {
                                                    bcDofId = 0;
                                                    if ( bcId <= sideNumBc.at(index) ) {
                                                        nodeDof = node->giveDof( sideBcDofId->at(bcId) );
                                                        if ( nodeDof->giveDofID() == dofIdArray.at(idof) ) {
                                                            bcDofId = nodeDof->giveBcId();
                                                            bcId++;
                                                        }
                                                    }

                                                    sprintf(line, " %d", bcDofId);
                                                    str += line;
                                                }
                                            }
                                        } else {
                                            index = 0;
                                            if ( m == 0 && faceNumBc.at(1) != 0 ) {
                                                index = 1;
                                            }

                                            if ( n == 0 && faceNumBc.at(2) != 0 ) {
                                                index = 2;
                                            }

                                            if ( k == 0 && faceNumBc.at(3) != 0 ) {
                                                index = 3;
                                            }

                                            if ( index != 0 ) {
                                                sprintf(line, " bc %d", dofs);
                                                str += line;

                                                // I rely on the fact that bc dofs to be reproduced are ordered with respect to the dof ordering of the corner node

                                                faceBcDofId = faceBcDofIdList.at(index);
                                                bcId = 1;
                                                for ( idof = 1; idof <= dofs; idof++ ) {
                                                    bcDofId = 0;
                                                    if ( bcId <= faceNumBc.at(index) ) {
                                                        nodeDof = node->giveDof( faceBcDofId->at(bcId) );
                                                        if ( nodeDof->giveDofID() == dofIdArray.at(idof) ) {
                                                            bcDofId = nodeDof->giveBcId();
                                                            bcId++;
                                                        }
                                                    }

                                                    sprintf(line, " %d", bcDofId);
                                                    str += line;
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    // copy node load

                                    if ( m == 0 && n == 0 ) {
                                        if ( ( loads = ( loadArray = node->giveLoadArray() )->giveSize() ) != 0 ) {
                                            sprintf(line, " load %d", loads);
                                            str += line;
                                            for ( iload = 1; iload <= loads; iload++ ) {
                                                sprintf( line, " %d", loadArray->at(iload) );
                                                str += line;
                                            }
                                        }
                                    }
                                }
                            }

                            refinedReader.appendInputString(str);
                        }
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: ElemMode ) {
        int mat, csect, iside, iload, loads;
        int nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8;
        IntArray *loadArray;
        AList< IntArray >boundaryLoadList;
        bool hasLoad;

        mat = element->giveMaterial()->giveNumber();
        csect = element->giveCrossSection()->giveNumber();

        boundaryLoadList.growTo(3);
        for ( iside = 1; iside <= 3; iside++ ) {
            loadArray = new IntArray;
            boundaryLoadList.put(iside, loadArray);
        }

        for ( inode = startNode; inode <= endNode; inode++ ) {
            connectivity = refinedElement->giveFineNodeArray(inode);
            refinedElement->giveBoundaryFlagArray(inode, element, boundary);
            hasLoad = refinedElement->giveBoundaryLoadArray3D(inode, element, boundaryLoadList);

            for ( k = 0; k < level + 1; k++ ) {
                for ( n = 0; n < level + 1; n++ ) {
                    for ( m = 0; m < level + 1; m++ ) {
                        localElemId++;

                        nd = k * ( level + 2 ) * ( level + 2 ) + n * ( level + 2 ) + m + 1;

                        nd1 = localNodeIdArray.at( connectivity->at(nd) );
                        nd2 = localNodeIdArray.at( connectivity->at(nd + 1) );
                        nd3 = localNodeIdArray.at( connectivity->at(nd + level + 3) );
                        nd4 = localNodeIdArray.at( connectivity->at(nd + level + 2) );

                        nd += ( level + 2 ) * ( level + 2 );

                        nd5 = localNodeIdArray.at( connectivity->at(nd) );
                        nd6 = localNodeIdArray.at( connectivity->at(nd + 1) );
                        nd7 = localNodeIdArray.at( connectivity->at(nd + level + 3) );
                        nd8 = localNodeIdArray.at( connectivity->at(nd + level + 2) );

                        sprintf(line, "%s %d nodes 8 %d %d %d %d %d %d %d %d mat %d crossSect %d",
                                hexatype, localElemId, nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, mat, csect);
                        str = line;

                        // copy body and boundary loads

                        if ( ( loads = ( loadArray = element->giveBodyLoadArray() )->giveSize() ) != 0 ) {
                            sprintf(line, " bodyLoads %d", loads);
                            str += line;
                            for ( iload = 1; iload <= loads; iload++ ) {
                                sprintf( line, " %d", loadArray->at(iload) );
                                str += line;
                            }
                        }

                        // boundary load is not copied on non-boundary sides

                        if ( hasLoad == true && ( m == 0 || n == 0 || k == 0 ) ) {
                            loads = 0;
                            for ( iside = 1; iside <= 3; iside++ ) {
                                if ( boundary.at(iside) == 0 ) {
                                    continue;
                                }

                                if ( m != 0 && iside == 1 ) {
                                    continue;
                                }

                                if ( n != 0 && iside == 2 ) {
                                    continue;
                                }

                                if ( k != 0 && iside == 3 ) {
                                    continue;
                                }

                                loads += boundaryLoadList.at(iside)->giveSize();
                            }

                            if ( loads != 0 ) {
                                sprintf(line, " boundaryLoads %d", loads);
                                str += line;
                                for ( iside = 1; iside <= 3; iside++ ) {
                                    if ( boundary.at(iside) == 0 ) {
                                        continue;
                                    }

                                    if ( m != 0 && iside == 1 ) {
                                        continue;
                                    }

                                    if ( n != 0 && iside == 2 ) {
                                        continue;
                                    }

                                    if ( k != 0 && iside == 3 ) {
                                        continue;
                                    }

                                    if ( ( loads = ( loadArray = boundaryLoadList.at(iside) )->giveSize() ) == 0 ) {
                                        continue;
                                    }

                                    for ( iload = 1; iload <= loads; iload++ ) {
                                        sprintf( line, " %d", loadArray->at(iload) );
                                        str += line;
                                    }
                                }
                            }
                        }

                        refinedReader.appendInputString(str);
                    }
                }
            }
        }

        return;
    }

    if ( mode == HuertaErrorEstimatorInterface :: BCMode ) {
        if ( aMode == HuertaErrorEstimator :: HEE_linear ) {
            int s1, s2, s3, f1, f2, f3, idof;
            double u, v, w, du = 1.0 / ( level + 1 ), dv = 1.0 / ( level + 1 ), dw = 1.0 / ( level + 1 );
            double xc, yc, zc, xm, ym, zm;
            double xs1, ys1, zs1, xs2, ys2, zs2, xs3, ys3, zs3;
            double xf1, yf1, zf1, xf2, yf2, zf2, xf3, yf3, zf3;
            MaterialMode mode;
            GaussPoint *gp;
            FloatArray globCoord(3), * locCoord;
            FloatMatrix Nmatrix;
            FloatArray uCoarse, uFine;

            mode = element->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0)->giveMaterialMode();

            // create a fictitious integration point
            locCoord = new FloatArray;
            IntegrationRule ir(0, element);
            gp = new GaussPoint( &ir, 1, locCoord, 1.0, mode);

            for ( inode = startNode; inode <= endNode; inode++ ) {
                s1 = hexaSideNode [ inode - 1 ] [ 0 ];
                s2 = hexaSideNode [ inode - 1 ] [ 1 ];
                s3 = hexaSideNode [ inode - 1 ] [ 2 ];
                f1 = hexaFaceNode [ inode - 1 ] [ 0 ];
                f2 = hexaFaceNode [ inode - 1 ] [ 1 ];
                f3 = hexaFaceNode [ inode - 1 ] [ 2 ];

                xc = corner [ inode - 1 ]->at(1);
                yc = corner [ inode - 1 ]->at(2);
                zc = corner [ inode - 1 ]->at(3);

                xs1 = midSide [ s1 - 1 ].at(1);
                ys1 = midSide [ s1 - 1 ].at(2);
                zs1 = midSide [ s1 - 1 ].at(3);

                xs2 = midSide [ s2 - 1 ].at(1);
                ys2 = midSide [ s2 - 1 ].at(2);
                zs2 = midSide [ s2 - 1 ].at(3);

                xs3 = midSide [ s3 - 1 ].at(1);
                ys3 = midSide [ s3 - 1 ].at(2);
                zs3 = midSide [ s3 - 1 ].at(3);

                xf1 = midFace [ f1 - 1 ].at(1);
                yf1 = midFace [ f1 - 1 ].at(2);
                zf1 = midFace [ f1 - 1 ].at(3);

                xf2 = midFace [ f2 - 1 ].at(1);
                yf2 = midFace [ f2 - 1 ].at(2);
                zf2 = midFace [ f2 - 1 ].at(3);

                xf3 = midFace [ f3 - 1 ].at(1);
                yf3 = midFace [ f3 - 1 ].at(2);
                zf3 = midFace [ f3 - 1 ].at(3);

                xm = midNode.at(1);
                ym = midNode.at(2);
                zm = midNode.at(3);

                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                // get corner displacements
                element->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, uCoarse);

                pos = 1;
                w = 0.0;
                for ( k = 0; k < level + 2; k++, w += dw ) {
                    v = 0.0;
                    for ( n = 0; n < level + 2; n++, v += dv ) {
                        u = 0.0;
                        for ( m = 0; m < level + 2; m++, u += du, pos++ ) {
                            nd = connectivity->at(pos);
                            if ( localNodeIdArray.at(nd) > 0 ) {
                                localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                                // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                                bc = 0;
                                if ( nodeId == 0 ) {
                                    if ( m == 0 && boundary.at(1) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( n == 0 && boundary.at(2) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( k == 0 && boundary.at(3) == 0 ) {
                                        bc = 1;
                                    }
                                } else {
                                    if ( k == level + 1 || n == level + 1 || m == level + 1 ) {
                                        bc = 1;
                                    }

                                    /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                    if ( bc == 0 ) {
                                        if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                            if ( m == 0 && boundary.at(1) == 0 ) {
                                                bc = 1;
                                            }

                                            if ( n == 0 && boundary.at(2) == 0 ) {
                                                bc = 1;
                                            }

                                            if ( k == 0 && boundary.at(3) == 0 ) {
                                                bc = 1;
                                            }
                                        }
                                    }

#endif
                                }

#ifdef EXACT_ERROR
                                if ( wholeFlag == true ) {
                                    bc = 0;
                                }

#endif

                                if ( bc == 1 ) {
                                    globCoord.at(1) = ( ( xc * ( 1.0 - u ) + xs1 * u ) * ( 1.0 - v ) + ( xs2 * ( 1.0 - u ) + xf1 * u ) * v ) * ( 1.0 - w )
                                                      + ( ( xs3 * ( 1.0 - u ) + xf2 * u ) * ( 1.0 - v ) + ( xf3 * ( 1.0 - u ) + xm * u ) * v ) * w;
                                    globCoord.at(2) = ( ( yc * ( 1.0 - u ) + ys1 * u ) * ( 1.0 - v ) + ( ys2 * ( 1.0 - u ) + yf1 * u ) * v ) * ( 1.0 - w )
                                                      + ( ( ys3 * ( 1.0 - u ) + yf2 * u ) * ( 1.0 - v ) + ( yf3 * ( 1.0 - u ) + ym * u ) * v ) * w;
                                    globCoord.at(3) = ( ( zc * ( 1.0 - u ) + zs1 * u ) * ( 1.0 - v ) + ( zs2 * ( 1.0 - u ) + zf1 * u ) * v ) * ( 1.0 - w )
                                                      + ( ( zs3 * ( 1.0 - u ) + zf2 * u ) * ( 1.0 - v ) + ( zf3 * ( 1.0 - u ) + zm * u ) * v ) * w;

                                    // this effectively rewrites the local coordinates of the fictitious integration point
                                    this->HuertaErrorEstimatorI_computeLocalCoords(* locCoord, globCoord);
                                    // get N matrix at the fictitious integration point
                                    this->HuertaErrorEstimatorI_computeNmatrixAt(gp, Nmatrix);
                                    // get displacement at the fictitious integration point
                                    uFine.beProductOf(Nmatrix, uCoarse);

                                    // first loadtime function must be constant 1.0
                                    for ( idof = 1; idof <= dofs; idof++ ) {
                                        sprintf( line, "BoundaryCondition %d loadTimeFunction %d prescribedvalue %e",
                                                ++localBcId, 1, uFine.at(idof) );
                                        refinedReader.appendInputString(line);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            delete gp;
        } else {
            int idof;

            for ( inode = startNode; inode <= endNode; inode++ ) {
                connectivity = refinedElement->giveFineNodeArray(inode);
                refinedElement->giveBoundaryFlagArray(inode, element, boundary);

                pos = 1;
                for ( k = 0; k < level + 2; k++ ) {
                    for ( n = 0; n < level + 2; n++ ) {
                        for ( m = 0; m < level + 2; m++, pos++ ) {
                            nd = connectivity->at(pos);
                            if ( localNodeIdArray.at(nd) > 0 ) {
                                localNodeIdArray.at(nd) = -localNodeIdArray.at(nd); // prevent repeated bc specification

                                // setup boundary conditions (for element (nodeId == 0), for patch (nodeId != 0))
                                bc = 0;
                                if ( nodeId == 0 ) {
                                    if ( m == 0 && boundary.at(1) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( n == 0 && boundary.at(2) == 0 ) {
                                        bc = 1;
                                    }

                                    if ( k == 0 && boundary.at(3) == 0 ) {
                                        bc = 1;
                                    }
                                } else {
                                    if ( k == level + 1 || n == level + 1 || m == level + 1 ) {
                                        bc = 1;
                                    }

                                    /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
                                    if ( bc == 0 ) {
                                        if ( element->giveNode(nodeId)->giveParallelMode() == DofManager_shared ) {
                                            if ( m == 0 && boundary.at(1) == 0 ) {
                                                bc = 1;
                                            }

                                            if ( n == 0 && boundary.at(2) == 0 ) {
                                                bc = 1;
                                            }

                                            if ( k == 0 && boundary.at(3) == 0 ) {
                                                bc = 1;
                                            }
                                        }
                                    }

#endif
                                }

#ifdef EXACT_ERROR
                                if ( wholeFlag == true ) {
                                    bc = 0;
                                }

#endif

                                if ( bc == 1 ) {
                                    for ( idof = 1; idof <= dofs; idof++ ) {
                                        localBcId++;
                                        controlNode.at(localBcId) = -localNodeIdArray.at(nd);
                                        controlDof.at(localBcId) = idof;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return;
    }
}





void
HuertaErrorEstimator :: solveRefinedElementProblem(int elemId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                   TimeStep *tStep)
{
    int contextFlag = 0;
    Element *element;
    RefinedElement *refinedElement;
    HuertaErrorEstimatorInterface *interface;
    EngngModel *problem, *refinedProblem;
    int localNodeId, localElemId, localBcId, localLtf;
    int mats, csects, loads, ltfuncs, nlbarriers;
    int inode, idof, dofs, pos, ielem, size;
    IntArray dofIdArray;
    FloatArray nodeSolution, uCoarse, elementVector, patchVector, coarseVector;
    FloatArray elementError, patchError, coarseSolution;
    Domain *domain = this->domain, *refinedDomain;
    Node *node;
    TimeStep *refinedTStep;
    double coarseSol;
    FloatMatrix mat;
    EIPrimaryUnknownMapper mapper;
    Dof *nodeDof;
    std :: string irec;
    double coeff, elementNorm, patchNorm, mixedNorm, eNorm = 0.0, uNorm = 0.0;
    IntArray controlNode, controlDof;

#ifdef TIME_INFO
    oofem_timeval st_total, et_total, st_setup, et_setup, st_init, et_init, st_solve, et_solve, st_error, et_error;

    :: getUtime(st_total);
    :: getUtime(st_setup);
#endif

    element = domain->giveElement(elemId);

#ifdef __PARALLEL_MODE
    if ( element->giveParallelMode() == Element_remote ) {
        this->skippedNelems++;
        this->eNorms.at(elemId) = 0.0;
        //  uNormArray.at(elemId) = 0.0;
        return;
    }

#endif

    if ( this->skipRegion( element->giveRegionNumber() ) != 0 ) {
        this->skippedNelems++;
        this->eNorms.at(elemId) = 0.0;
        //  uNormArray.at(elemId) = 0.0;

#ifdef INFO
        //  printf("\nElement no %d: skipped          [step number %5d]\n", elemId, tStep -> giveNumber());
#endif

        return;
    }

#ifdef INFO
 #ifdef __PARALLEL_MODE
    OOFEM_LOG_DEBUG( "[%d] Element no %d: estimating error [step number %5d]\n",
                    domain->giveEngngModel()->giveRank(), elemId, tStep->giveNumber() );
 #else
    OOFEM_LOG_DEBUG( "Element no %d: estimating error [step number %5d]\n", elemId, tStep->giveNumber() );
 #endif
#endif

    refinedElement = this->refinedElementList.at(elemId);
    interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
    if ( interface == NULL ) {
        _error("solveRefinedElementProblem: Element has no Huerta error estimator interface defined");
    }

    problem = domain->giveEngngModel();

    mats = domain->giveNumberOfMaterialModels();
    csects = domain->giveNumberOfCrossSectionModels();
    loads = domain->giveNumberOfBoundaryConditions();
    ltfuncs = domain->giveNumberOfLoadTimeFunctions();
    nlbarriers = domain->giveNumberOfNonlocalBarriers();

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = 0;
    localLtf = 0;

    interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                localNodeIdArray, globalNodeIdArray,
                                                                HuertaErrorEstimatorInterface :: CountMode, tStep,
                                                                localNodeId, localElemId, localBcId,
                                                                controlNode, controlDof,
                                                                this->mode);

    if ( this->mode == HEE_nlinear ) {
        controlDof.resize(localBcId);
        controlNode.resize(localBcId);

        localBcId = 0;
        interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                    localNodeIdArray, globalNodeIdArray,
                                                                    HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                    localNodeId, localElemId, localBcId,
                                                                    controlNode, controlDof,
                                                                    this->mode);

        localBcId = 0;
        localLtf = 1;
    }

    setupRefinedProblemProlog("element", elemId, localNodeIdArray, localNodeId, localElemId,
                              mats, csects, loads + localBcId, ltfuncs + localLtf,
                              controlNode, controlDof, tStep);

    globalNodeIdArray.resize(localNodeId);

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = loads;

    interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                localNodeIdArray, globalNodeIdArray,
                                                                HuertaErrorEstimatorInterface :: NodeMode, tStep,
                                                                localNodeId, localElemId, localBcId,
                                                                controlNode, controlDof,
                                                                this->mode);
    interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                localNodeIdArray, globalNodeIdArray,
                                                                HuertaErrorEstimatorInterface :: ElemMode, tStep,
                                                                localNodeId, localElemId, localBcId,
                                                                controlNode, controlDof,
                                                                this->mode);


    setupRefinedProblemEpilog1(csects, mats, loads, nlbarriers);

    if ( this->mode == HEE_linear ) {
        localBcId = loads;
        interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                    localNodeIdArray, globalNodeIdArray,
                                                                    HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                    localNodeId, localElemId, localBcId,
                                                                    controlNode, controlDof,
                                                                    this->mode);
    }

    setupRefinedProblemEpilog2(ltfuncs);

#ifdef TIME_INFO
    :: getRelativeUtime(et_setup, st_setup);
#endif

    dofs = domain->giveDefaultNodeDofIDArry().giveSize();
    dofIdArray = domain->giveDefaultNodeDofIDArry();

#ifdef USE_INPUT_FILE
    std::ostringstream fileName;
    fileName << "/home/dr/Huerta/element_" << elemId << ".in";
    refinedReader.writeToFile(fileName.str().c_str());
#endif

#ifdef TIME_INFO
    :: getUtime(st_init);
#endif
    refinedReader.rewind();
    refinedProblem = InstanciateProblem(& refinedReader, _processor, contextFlag);
    refinedReader.finish();
#ifdef TIME_INFO
    :: getRelativeUtime(et_init, st_init);
#endif

#ifdef DEBUG
    refinedProblem->checkConsistency();
#endif

    refinedDomain = refinedProblem->giveDomain(1);

    // solve the problem first and then map the coarse solution;
    // when mapping the coarse solution at first, tstep is NULL and it cannot be accessed;
    // this makes some overhead for nonlinear problems because the coarse solution is mapped twice
    // when initiating the solution and after the solution

#ifdef TIME_INFO
    :: getUtime(st_solve);
#endif
    if ( this->mode == HEE_linear ) {
        refinedProblem->solveYourself();
        refinedProblem->terminateAnalysis();
    } else {
        if ( refinedProblem->giveClassID() == AdaptiveNonLinearStaticClass ) {
            ( ( AdaptiveNonLinearStatic * ) refinedProblem )->initializeAdaptiveFrom(problem);
        } else {
            _error("sorry");
        }
    }

#ifdef TIME_INFO
    :: getRelativeUtime(et_solve, st_solve);
#endif

    //fprintf(stdout, "\n");

#ifdef TIME_INFO
    :: getUtime(st_error);
#endif
    refinedTStep = refinedProblem->giveCurrentStep();

    size = refinedDomain->giveNumberOfDofManagers() * dofs;
    coarseSolution.resize(size);
    elementError.resize(size);
    patchError.resize(size);

    // map coarse solution
    uCoarse.resize( refinedProblem->giveNumberOfDomainEquations(1, EID_MomentumBalance) );
    uCoarse.zero();
    mapper.mapAndUpdate(uCoarse, VM_Total, EID_MomentumBalance, domain, refinedDomain, tStep);

    // get coarse solution and element and patch error (including BC !!!)
    pos = 1;
    for ( inode = 1; inode <= localNodeId; inode++ ) {
        node = refinedDomain->giveNode(inode);
        node->giveUnknownVector(nodeSolution, dofIdArray,
                                EID_MomentumBalance, VM_Total, refinedTStep);
        for ( idof = 1; idof <= dofs; idof++, pos++ ) {
            nodeDof = node->giveDof(idof);
            if ( nodeDof->hasBc(refinedTStep) == 0 ) {
                coarseSol = uCoarse.at( nodeDof->__giveEquationNumber() );
            } else {
                // coarse solution is identical with fine solution at BC
                coarseSol = nodeSolution.at(idof);
            }

            //    coarseSol = nodeDof -> giveBcValue(VM_Total, refinedTStep);

            coarseSolution.at(pos) = coarseSol;
            elementError.at(pos) = nodeSolution.at(idof) - coarseSol;
            patchError.at(pos) = primaryUnknownError.at( ( globalNodeIdArray.at(inode) - 1 ) * dofs + idof ) - coarseSol;
        }
    }

    /* enforce zero error on element and patch boundary (this is just for 1D !!!)
     * (for nonlinear problem there may be a nonzero error due to the tolerance) */
    /*
     * elementError.at(1) = 0.0;
     * elementError.at(this -> refineLevel + 3) = 0.0;
     *
     * patchError.at(this -> refineLevel + 2) = 0.0;
     */
    if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
        FloatMatrix Nmatrix;
        FloatArray elementVectorGp, patchVectorGp, coarseVectorGp;
        IntegrationRule *iRule;
        GaussPoint *gp;
        int igp;
        double dV;

        eNorm = uNorm = 0.0;
        for ( ielem = 1; ielem <= localElemId; ielem++ ) {
            element = refinedDomain->giveElement(ielem);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
            iRule = element->giveDefaultIntegrationRulePtr();

            elementNorm = patchNorm = mixedNorm = 0.0;
            for ( igp = 0; igp < iRule->getNumberOfIntegrationPoints(); igp++ ) {
                gp = iRule->getIntegrationPoint(igp);
                dV = element->computeVolumeAround(gp);

                interface->HuertaErrorEstimatorI_computeNmatrixAt(gp, Nmatrix);

                this->extractVectorFrom(element, elementError, elementVector, dofs, refinedTStep);
                elementVectorGp.beProductOf(Nmatrix, elementVector);
                this->extractVectorFrom(element, patchError, patchVector, dofs, refinedTStep);
                patchVectorGp.beProductOf(Nmatrix, patchVector);
                this->extractVectorFrom(element, coarseSolution, coarseVector, dofs, refinedTStep);
                coarseVectorGp.beProductOf(Nmatrix, coarseVector);

                elementNorm += elementVectorGp.computeSquaredNorm() * dV;
                mixedNorm += elementVectorGp.dotProduct(patchVectorGp) * dV;
                patchNorm += patchVectorGp.computeSquaredNorm() * dV;
                uNorm += coarseVectorGp.computeSquaredNorm() * dV;
            }

            if ( fabs(elementNorm) < 1.0e-30 ) {
                if ( elementNorm == 0.0 ) {
                    coeff = 0.0;
                } else {
                    if ( fabs(mixedNorm) > 1.0e6 * fabs(elementNorm) ) {
                        _error("solveRefinedElementProblem: division by zero");
                    }

                    coeff = mixedNorm / elementNorm;
                }
            } else {
                coeff = mixedNorm / elementNorm;
            }

            eNorm += ( 1.0 + coeff * coeff ) * elementNorm + patchNorm - 2.0 * coeff * mixedNorm;
        }
    } else if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
        FloatArray tmpVector;
        double eEnorm, pEnorm;

#ifdef PRINT_FINE_ERROR
        OOFEM_LOG_DEBUG("\n");
        if ( exactFlag == true ) {
            OOFEM_LOG_DEBUG(" elem  sub         e_Error2        p_Error2         a_Error2         x_Error2 \n");
            OOFEM_LOG_DEBUG("------------------------------------------------------------------------------\n");
        } else {
            OOFEM_LOG_DEBUG(" elem  sub         e_Error2        p_Error2         a_Error2 \n");
            OOFEM_LOG_DEBUG("-------------------------------------------------------------\n");
        }

#endif

        eNorm = uNorm = 0.0;
        for ( ielem = 1; ielem <= localElemId; ielem++ ) {
            element = refinedDomain->giveElement(ielem);
            refinedProblem->giveElementCharacteristicMatrix(mat, ielem, STIFFNESS_TYPE, refinedTStep, refinedDomain);

            this->extractVectorFrom(element, elementError, elementVector, dofs, refinedTStep);
            this->extractVectorFrom(element, patchError, patchVector, dofs, refinedTStep);
            this->extractVectorFrom(element, coarseSolution, coarseVector, dofs, refinedTStep);

            tmpVector.beProductOf(mat, elementVector);
            elementNorm = tmpVector.dotProduct(elementVector);
            mixedNorm = tmpVector.dotProduct(patchVector);

            tmpVector.beProductOf(mat, patchVector);
            patchNorm = tmpVector.dotProduct(patchVector);

            if ( fabs(elementNorm) < 1.0e-30 ) {
                if ( elementNorm == 0.0 ) {
                    coeff = 0.0;
                } else {
                    if ( fabs(mixedNorm) > 1.0e6 * fabs(elementNorm) ) {
                        _error("solveRefinedElementProblem: division by zero");
                    }

                    coeff = mixedNorm / elementNorm;
                }
            } else {
                coeff = mixedNorm / elementNorm;
            }

            eNorm += ( 1.0 + coeff * coeff ) * elementNorm + patchNorm - 2.0 * coeff * mixedNorm;

            eEnorm = elementNorm;
            pEnorm = coeff * coeff * elementNorm + patchNorm - 2.0 * coeff * mixedNorm;
            /*
             * elementVector.times(coeff);
             * patchVector.subtract(elementVector);
             * elementVector.times(1.0/coeff);
             *
             * tmpVector.beProductOf(mat, patchVector);
             * pEnorm = dotProduct(tmpVector.givePointer(), patchVector.givePointer(), patchVector.giveSize());
             * eEnorm = dotProduct(tmpVector.givePointer(), elementVector.givePointer(), elementVector.giveSize());
             */
            tmpVector.beProductOf(mat, coarseVector);
            uNorm += tmpVector.dotProduct(coarseVector);

#ifdef PRINT_FINE_ERROR
            if ( exactFlag == false ) {
                OOFEM_LOG_DEBUG("%5d: %3d  %15.8e %15.8e  %15.8e\n",
                                elemId, ielem, eEnorm, pEnorm, eEnorm + pEnorm);
            }

 #ifdef EXACT_ERROR
            else {
                OOFEM_LOG_DEBUG( "%5d: %3d  %15.8e %15.8e  %15.8e  %15.8e\n",
                                elemId, ielem, eEnorm, pEnorm, eEnorm + pEnorm, exactFineError.at(++finePos) );
            }
 #endif
#endif
        }

#ifdef PRINT_FINE_ERROR
        OOFEM_LOG_DEBUG("\n");
#endif
    } else {
        _error("solveRefinedElementProblem: Unsupported norm type");
    }

    // update primaryUnknownError
    // (this is usefull for postprocessing only, but how to draw fine elements ?)
    // be carefull to not overwrite data needed in further calculations
    // pos = 1;
    // for(inode = 1; inode <= localNodeId; inode++){
    //  for(idof = 1; idof <= dofs; idof++, pos++)primaryUnknownError.at((globalNodeIdArray.at(inode) - 1) * dofs + idof) =
    //   elementError.at(pos) * (1.0 - coeff) + patchError.at(pos);
    // }

    element = domain->giveElement(elemId);

    double eeeNorm = 0.0;

#ifdef HUHU
    IntegrationRule *iRule;
    int nip, result;
    double sval, maxVal, eeeNorm = 0.0;
    int k;

    maxVal = 0.0;
    iRule = element->giveDefaultIntegrationRulePtr();
    nip = iRule->getNumberOfIntegrationPoints();
    for ( k = 0; k < nip; k++ ) {
        result = element->giveIPValue(val, iRule->getIntegrationPoint(k), IST_PrincipalDamageTensor, tStep);
        if ( result ) {
            sval = sqrt( dotProduct( val, val, val.giveSize() ) );
            maxVal = max(maxVal, sval);
        }
    }

    if ( maxVal > 0.25 ) {
        double rate, size = 1.0, currDensity;

        HuertaRemeshingCriteriaInterface *remeshInterface;
        remeshInterface = ( HuertaRemeshingCriteriaInterface * ) element->giveInterface(HuertaRemeshingCriteriaInterfaceType);
        if ( !remeshInterface ) {
            _error("estimateMeshDensities: element does not support HuertaRemeshingCriteriaInterface");
        }

        currDensity = remeshInterface->HuertaRemeshingCriteriaI_giveCharacteristicSize();

        rate = currDensity / size;
        if ( rate < 1.0 ) {
            rate = 1.0;
        }

        OOFEM_LOG_DEBUG("koko %d dam %e curr %e rate %e\n", elemId, maxVal, currDensity, rate);

        // eNorm *= rate * rate;
        eeeNorm = eNorm * rate;
    }

#endif

    this->eNorms.at(elemId) = sqrt(eNorm);
    // uNormArray.at(elemId) = sqrt(uNorm);

    this->globalENorm += eNorm + eeeNorm;
    this->globalUNorm += uNorm;

#ifdef TIME_INFO
    :: getRelativeUtime(et_error, st_error);
    :: getRelativeUtime(et_total, st_total);

    OOFEM_LOG_DEBUG( "HEE info: element %d: user time total %.2f s (setup %.2f s, init %.2f s, solve %.2f s, error %.2f s)\n",
                    elemId,
                    ( double ) ( et_total.tv_sec + et_total.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_setup.tv_sec + et_setup.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_init.tv_sec + et_init.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_solve.tv_sec + et_solve.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_error.tv_sec + et_error.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif

    delete refinedProblem;
}



void
HuertaErrorEstimator :: solveRefinedPatchProblem(int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                 TimeStep *tStep)
{
    int contextFlag = 0;
    Element *element;
    RefinedElement *refinedElement;
    HuertaErrorEstimatorInterface *interface;
    EngngModel *problem, *refinedProblem;
    int localNodeId, localElemId, localBcId, localLtf;
    int mats, csects, loads, ltfuncs, nlbarriers;
    int inode, elemId, ielem, elems, skipped = 0;
    const IntArray *con;
    int idof, dofs, pos;
    IntArray dofIdArray;
    FloatArray nodeSolution;
    Domain *domain = this->domain, *refinedDomain;
    TimeStep *refinedTStep;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    IntArray controlNode, controlDof;

#ifdef TIME_INFO
    oofem_timeval st_total, et_total, st_setup, et_setup, st_init, et_init, st_solve, et_solve, st_error, et_error;

    :: getUtime(st_total);
    :: getUtime(st_setup);
#endif

#ifdef __PARALLEL_MODE
    dofManagerParallelMode parMode;

    parMode = domain->giveDofManager(nodeId)->giveParallelMode();
    if ( parMode == DofManager_remote || parMode == DofManager_null ) {
        return;
    }

#endif

    con = ct->giveDofManConnectivityArray(nodeId);
    elems = con->giveSize();

    for ( ielem = 1; ielem <= elems; ielem++ ) {
        elemId = con->at(ielem);
        element = domain->giveElement(elemId);

        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() == Element_remote ) {
            skipped++;
            continue;
        }

#endif

        if ( this->skipRegion( element->giveRegionNumber() ) != 0 ) {
            skipped++;
        }
    }

    if ( skipped == elems ) {
#ifdef INFO
        //  printf("\nPatch no %d: skipped          [step number %5d]\n", nodeId, tStep -> giveNumber());
#endif

        return;
    }

#ifdef INFO
 #ifdef __PARALLEL_MODE
    OOFEM_LOG_INFO( "[%d] Patch no %d: estimating error [step number %5d]\n",
                   domain->giveEngngModel()->giveRank(), nodeId, tStep->giveNumber() );
 #else
    OOFEM_LOG_INFO( "Patch no %d: estimating error [step number %5d]\n", nodeId, tStep->giveNumber() );
 #endif
#endif

    problem = domain->giveEngngModel();

    mats = domain->giveNumberOfMaterialModels();
    csects = domain->giveNumberOfCrossSectionModels();
    loads = domain->giveNumberOfBoundaryConditions();
    ltfuncs = domain->giveNumberOfLoadTimeFunctions();
    nlbarriers = domain->giveNumberOfNonlocalBarriers();

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = 0;
    localLtf = 0;

    for ( ielem = 1; ielem <= elems; ielem++ ) {
        elemId = con->at(ielem);
        element = domain->giveElement(elemId);

        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
        if ( interface == NULL ) {
            _error("solveRefinedPatchProblem: Element has no Huerta error estimator interface defined");
        }

        for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
            if ( element->giveNode(inode)->giveNumber() != nodeId ) {
                continue;
            }

            interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, inode,
                                                                        localNodeIdArray, globalNodeIdArray,
                                                                        HuertaErrorEstimatorInterface :: CountMode, tStep,
                                                                        localNodeId, localElemId, localBcId,
                                                                        controlNode, controlDof,
                                                                        this->mode);
            break;
        }
    }

    if ( this->mode == HEE_nlinear ) {
        controlDof.resize(localBcId);
        controlNode.resize(localBcId);

        localBcId = 0;
        for ( ielem = 1; ielem <= elems; ielem++ ) {
            elemId = con->at(ielem);
            element = domain->giveElement(elemId);

            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
            if ( element->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif

            refinedElement = this->refinedElementList.at(elemId);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);

            for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
                if ( element->giveNode(inode)->giveNumber() != nodeId ) {
                    continue;
                }

                interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, inode,
                                                                            localNodeIdArray, globalNodeIdArray,
                                                                            HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                            localNodeId, localElemId, localBcId,
                                                                            controlNode, controlDof,
                                                                            this->mode);
                break;
            }
        }

        localBcId = 0;
        localLtf = 1;
    }

    setupRefinedProblemProlog("patch", nodeId, localNodeIdArray, localNodeId, localElemId,
                              mats, csects, loads + localBcId, ltfuncs + localLtf,
                              controlNode, controlDof, tStep);

    globalNodeIdArray.resize(localNodeId);

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = loads;

    for ( ielem = 1; ielem <= elems; ielem++ ) {
        elemId = con->at(ielem);
        element = domain->giveElement(elemId);

        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);

        for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
            if ( element->giveNode(inode)->giveNumber() != nodeId ) {
                continue;
            }

            interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, inode,
                                                                        localNodeIdArray, globalNodeIdArray,
                                                                        HuertaErrorEstimatorInterface :: NodeMode, tStep,
                                                                        localNodeId, localElemId, localBcId,
                                                                        controlNode, controlDof,
                                                                        this->mode);
            break;
        }
    }

    for ( ielem = 1; ielem <= elems; ielem++ ) {
        elemId = con->at(ielem);
        element = domain->giveElement(elemId);

        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif

        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);

        for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
            if ( element->giveNode(inode)->giveNumber() != nodeId ) {
                continue;
            }

            interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, inode,
                                                                        localNodeIdArray, globalNodeIdArray,
                                                                        HuertaErrorEstimatorInterface :: ElemMode, tStep,
                                                                        localNodeId, localElemId, localBcId,
                                                                        controlNode, controlDof,
                                                                        this->mode);
            break;
        }
    }

    setupRefinedProblemEpilog1(csects, mats, loads, nlbarriers);

    if ( this->mode == HEE_linear ) {
        localBcId = loads;
        for ( ielem = 1; ielem <= elems; ielem++ ) {
            elemId = con->at(ielem);
            element = domain->giveElement(elemId);

            /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
            if ( element->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif

            refinedElement = this->refinedElementList.at(elemId);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);

            for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
                if ( element->giveNode(inode)->giveNumber() != nodeId ) {
                    continue;
                }

                interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, inode,
                                                                            localNodeIdArray, globalNodeIdArray,
                                                                            HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                            localNodeId, localElemId, localBcId,
                                                                            controlNode, controlDof,
                                                                            this->mode);
                break;
            }
        }
    }

    setupRefinedProblemEpilog2(ltfuncs);

#ifdef TIME_INFO
    :: getRelativeUtime(et_setup, st_setup);
#endif

    dofs = domain->giveDefaultNodeDofIDArry().giveSize();
    dofIdArray = domain->giveDefaultNodeDofIDArry();

#ifdef USE_INPUT_FILE
    std::ostringstream fileName;
    fileName << "/home/dr/Huerta/patch_" << nodeId << ".in";
    refinedReader.writeToFile(fileName.str().c_str());
#endif

#ifdef TIME_INFO
    :: getUtime(st_init);
#endif
    refinedReader.rewind();
    refinedProblem = InstanciateProblem(& refinedReader, _processor, contextFlag);
    refinedReader.finish();
#ifdef TIME_INFO
    :: getRelativeUtime(et_init, st_init);
#endif

#ifdef DEBUG
    refinedProblem->checkConsistency();
#endif

    refinedDomain = refinedProblem->giveDomain(1);

#ifdef TIME_INFO
    :: getUtime(st_solve);
#endif
    if ( this->mode == HEE_linear ) {
        refinedProblem->solveYourself();
        refinedProblem->terminateAnalysis();
    } else {
        if ( refinedProblem->giveClassID() == AdaptiveNonLinearStaticClass ) {
            ( ( AdaptiveNonLinearStatic * ) refinedProblem )->initializeAdaptiveFrom(problem);
        } else {
            _error("sorry");
        }
    }

#ifdef TIME_INFO
    :: getRelativeUtime(et_solve, st_solve);
#endif

    //fprintf(stdout, "\n");

#ifdef TIME_INFO
    :: getUtime(st_error);
#endif
    refinedTStep = refinedProblem->giveCurrentStep();

    // store fine solution in primaryUnknownError
    for ( inode = 1; inode <= localNodeId; inode++ ) {
        refinedDomain->giveNode(inode)->giveUnknownVector(nodeSolution, dofIdArray,
                                                          EID_MomentumBalance, VM_Total,
                                                          refinedTStep);
        pos = globalNodeIdArray.at(inode);
        for ( idof = 1; idof <= dofs; idof++ ) {
            primaryUnknownError.at( ( pos - 1 ) * dofs + idof ) = nodeSolution.at(idof);
        }
    }

#ifdef TIME_INFO
    :: getRelativeUtime(et_error, st_error);
    :: getRelativeUtime(et_total, st_total);

    OOFEM_LOG_DEBUG( "HEE info: patch %d: user time total %.2f s (setup %.2f s, init %.2f s, solve %.2f s, error %.2f s)\n",
                    nodeId,
                    ( double ) ( et_total.tv_sec + et_total.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_setup.tv_sec + et_setup.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_init.tv_sec + et_init.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_solve.tv_sec + et_solve.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_error.tv_sec + et_error.tv_usec / ( double ) OOFEM_USEC_LIM ) );
#endif

    delete refinedProblem;
}


#ifndef EXACT_ERROR
void
HuertaErrorEstimator :: solveRefinedWholeProblem(IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                 TimeStep *tStep) { }
#else
void
HuertaErrorEstimator :: solveRefinedWholeProblem(IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                 TimeStep *tStep)
{
    int contextFlag = 0;
    Element *element;
    RefinedElement *refinedElement;
    HuertaErrorEstimatorInterface *interface;
    EngngModel *problem, *refinedProblem;
    int localNodeId, localElemId, localBcId, localLtf;
    int mats, csects, loads, ltfuncs, nlbarriers;
    int inode, idof, dofs, pos, elemId, ielem, elems, size;
    IntArray dofIdArray;
    FloatArray nodeSolution, uCoarse, errorVector, coarseVector, fineVector;
    FloatArray fineSolution, coarseSolution, errorSolution;
    Domain *domain = this->domain, *refinedDomain;
    Node *node;
    TimeStep *refinedTStep;
    FloatMatrix mat;
    EIPrimaryUnknownMapper mapper;
    Dof *nodeDof;
    IntArray controlNode, controlDof;

 #ifdef TIME_INFO
    oofem_timeval st_total, et_total, st_setup, et_setup, st_init, et_init, st_solve, et_solve, st_error, et_error;

    :: getUtime(st_total);
    :: getUtime(st_setup);
 #endif
 #ifdef INFO
    OOFEM_LOG_INFO( "Whole 0: estimating error [step number %5d]\n", tStep->giveNumber() );
 #endif

    elems = domain->giveNumberOfElements();

    problem = domain->giveEngngModel();

    mats = domain->giveNumberOfMaterialModels();
    csects = domain->giveNumberOfCrossSectionModels();
    loads = domain->giveNumberOfBoundaryConditions();
    ltfuncs = domain->giveNumberOfLoadTimeFunctions();
    nlbarriers = domain->giveNumberOfNonlocalBarriers();

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = 0;
    localLtf = 0;

    for ( elemId = 1; elemId <= elems; elemId++ ) {
        element = domain->giveElement(elemId);
        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
        if ( interface == NULL ) {
            _error("solveRefinedWholeProblem: Element has no Huerta error estimator interface defined");
        }

        interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                    localNodeIdArray, globalNodeIdArray,
                                                                    HuertaErrorEstimatorInterface :: CountMode, tStep,
                                                                    localNodeId, localElemId, localBcId,
                                                                    controlNode, controlDof,
                                                                    this->mode);
    }

    if ( this->mode == HEE_nlinear ) {
        controlDof.resize(localBcId);
        controlNode.resize(localBcId);

        localBcId = 0;
        for ( elemId = 1; elemId <= elems; elemId++ ) {
            element = domain->giveElement(elemId);
            refinedElement = this->refinedElementList.at(elemId);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);

            interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                        localNodeIdArray, globalNodeIdArray,
                                                                        HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                        localNodeId, localElemId, localBcId,
                                                                        controlNode, controlDof,
                                                                        this->mode);
        }

        localBcId = 0;
        localLtf = 1;
    }

    setupRefinedProblemProlog("whole", 0, localNodeIdArray, localNodeId, localElemId,
                              mats, csects, loads + localBcId, ltfuncs + localLtf,
                              controlNode, controlDof, tStep);

    globalNodeIdArray.resize(localNodeId);

    localNodeIdArray.zero();
    localNodeId = 0;
    localElemId = 0;
    localBcId = loads;

    for ( elemId = 1; elemId <= elems; elemId++ ) {
        element = domain->giveElement(elemId);
        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
        interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                    localNodeIdArray, globalNodeIdArray,
                                                                    HuertaErrorEstimatorInterface :: NodeMode, tStep,
                                                                    localNodeId, localElemId, localBcId,
                                                                    controlNode, controlDof,
                                                                    this->mode);
    }

    for ( elemId = 1; elemId <= elems; elemId++ ) {
        element = domain->giveElement(elemId);
        refinedElement = this->refinedElementList.at(elemId);
        interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
        interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                    localNodeIdArray, globalNodeIdArray,
                                                                    HuertaErrorEstimatorInterface :: ElemMode, tStep,
                                                                    localNodeId, localElemId, localBcId,
                                                                    controlNode, controlDof,
                                                                    this->mode);
    }

    setupRefinedProblemEpilog1(csects, mats, loads, nlbarriers);

    if ( this->mode == HEE_linear ) {
        localBcId = loads;
        for ( elemId = 1; elemId <= elems; elemId++ ) {
            element = domain->giveElement(elemId);
            refinedElement = this->refinedElementList.at(elemId);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
            interface->HuertaErrorEstimatorI_setupRefinedElementProblem(refinedElement, this->refineLevel, 0,
                                                                        localNodeIdArray, globalNodeIdArray,
                                                                        HuertaErrorEstimatorInterface :: BCMode, tStep,
                                                                        localNodeId, localElemId, localBcId,
                                                                        controlNode, controlDof,
                                                                        this->mode);
        }
    }

    setupRefinedProblemEpilog2(ltfuncs);

 #ifdef TIME_INFO
    :: getRelativeUtime(et_setup, st_setup);
 #endif

    dofs = domain->giveDefaultNodeDofIDArry().giveSize();
    dofIdArray = domain->giveDefaultNodeDofIDArry();

 #ifdef USE_INPUT_FILE
    std::ostringstream fileName;
    fileName << "/home/dr/Huerta/whole_" << 0 << ".in";
    refinedReader.writeToFile(fileName.str().c_str());
 #endif

 #ifdef TIME_INFO
    :: getUtime(st_init);
 #endif
    refinedReader.rewind();
    refinedProblem = InstanciateProblem(& refinedReader, _processor, contextFlag);
    refinedReader.finish();
 #ifdef TIME_INFO
    :: getRelativeUtime(et_init, st_init);
 #endif

 #ifdef DEBUG
    refinedProblem->checkConsistency();
 #endif

    refinedDomain = refinedProblem->giveDomain(1);

    // solve the problem first and then map the coarse solution;
    // when mapping the coarse solution at first, tstep is NULL and it cannot be accessed;

 #ifdef TIME_INFO
    :: getUtime(st_solve);
 #endif
    refinedProblem->solveYourself();
    refinedProblem->terminateAnalysis();
 #ifdef TIME_INFO
    :: getRelativeUtime(et_solve, st_solve);
 #endif

    //fprintf(stdout, "\n");

 #ifdef TIME_INFO
    :: getUtime(st_error);
 #endif
    refinedTStep = refinedProblem->giveCurrentStep();

    size = refinedDomain->giveNumberOfDofManagers() * dofs;
    fineSolution.resize(size);
    coarseSolution.resize(size);
    errorSolution.resize(size);

    // map coarse solution
    uCoarse.resize( refinedProblem->giveNumberOfDomainEquations(1, EID_MomentumBalance) );
    uCoarse.zero();
    mapper.mapAndUpdate(uCoarse, VM_Total, EID_MomentumBalance, domain, refinedDomain, tStep);

    // get exact and coarse solution (including BC !!!)
    pos = 1;
    for ( inode = 1; inode <= localNodeId; inode++ ) {
        node = refinedDomain->giveNode(inode);
        node->giveUnknownVector(nodeSolution, dofIdArray,
                                EID_MomentumBalance, VM_Total, refinedTStep);
        for ( idof = 1; idof <= dofs; idof++, pos++ ) {
            fineSolution.at(pos) = nodeSolution.at(idof);
            nodeDof = node->giveDof(idof);
            if ( nodeDof->hasBc(refinedTStep) == 0 ) {
                coarseSolution.at(pos) = uCoarse.at( nodeDof->__giveEquationNumber() );
            } else {
                // coarse solution is identical with fine solution at BC
                coarseSolution.at(pos) = nodeSolution.at(idof);
            }

            //    coarseSolution.at(pos) = nodeDof -> giveBcValue(VM_Total, refinedTStep);
        }
    }

    errorSolution = fineSolution;
    errorSolution.subtract(coarseSolution);

    if ( this->normType == HuertaErrorEstimator :: L2Norm ) {
        FloatMatrix Nmatrix;
        FloatArray errorVectorGp, coarseVectorGp, fineVectorGp;
        IntegrationRule *iRule;
        GaussPoint *gp;
        int igp;
        double dV;
        /*
         * exactENorm = dotProduct(errorSolution.givePointer(), errorSolution.givePointer(), size);
         * coarseUNorm = dotProduct(coarseSolution.givePointer(), coarseSolution.givePointer(), size);
         * fineUNorm = dotProduct(fineSolution.givePointer(), fineSolution.givePointer(), size);
         */
        exactENorm = coarseUNorm = fineUNorm = mixedNorm = 0.0;
        for ( ielem = 1; ielem <= localElemId; ielem++ ) {
            element = refinedDomain->giveElement(ielem);
            interface = ( HuertaErrorEstimatorInterface * ) element->giveInterface(HuertaErrorEstimatorInterfaceType);
            iRule = element->giveDefaultIntegrationRulePtr();

            for ( igp = 0; igp < iRule->getNumberOfIntegrationPoints(); igp++ ) {
                gp = iRule->getIntegrationPoint(igp);
                dV = element->computeVolumeAround(gp);

                interface->HuertaErrorEstimatorI_computeNmatrixAt(gp, Nmatrix);

                this->extractVectorFrom(element, errorSolution, errorVector, dofs, refinedTStep);
                errorVectorGp.beProductOf(Nmatrix, errorVector);
                exactENorm += errorVectorGp.computeSquaredNorm() * dV;

                this->extractVectorFrom(element, coarseSolution, coarseVector, dofs, refinedTStep);
                coarseVectorGp.beProductOf(Nmatrix, coarseVector);
                coarseUNorm += coarseVectorGp.computeSquaredNorm() * dV;

                mixedNorm += coarseVectorGp.dotProduct(errorVectorGp) * dV;

                this->extractVectorFrom(element, fineSolution, fineVector, dofs, refinedTStep);
                fineVectorGp.beProductOf(Nmatrix, fineVector);
                fineUNorm += fineVectorGp.computeSquaredNorm() * dV;
            }
        }
    } else if ( this->normType == HuertaErrorEstimator :: EnergyNorm ) {
        FloatArray tmpVector;

 #ifdef PRINT_ERROR
        double fineENorm, coarseENorm;
        int dim, pos, nelems;

        elemId = 1;
        element = domain->giveElement(elemId);
        dim = element->giveSpatialDimension();
        nelems = this->refineLevel + 1;
        while ( --dim ) {
            nelems *= this->refineLevel + 1;
        }

        nelems *= element->giveNumberOfNodes();
        coarseENorm = 0.0;
        pos = 0;
 #endif

        exactENorm = coarseUNorm = fineUNorm = mixedNorm = 0.0;
        for ( ielem = 1; ielem <= localElemId; ielem++ ) {
            if ( this->skipRegion( element->giveRegionNumber() ) != 0 ) {
 #ifdef PRINT_ERROR
                exactFineError.at(ielem) = 0.0;
                if ( ++pos == nelems ) {
                    exactCoarseError.at(elemId) = coarseENorm;
                    if ( ielem < localElemId ) {
                        elemId++;
                        element = domain->giveElement(elemId);
                        dim = element->giveSpatialDimension();
                        nelems = this->refineLevel + 1;
                        while ( --dim ) {
                            nelems *= this->refineLevel + 1;
                        }

                        nelems *= element->giveNumberOfNodes();
                        coarseENorm = 0.0;
                        pos = 0;
                    }
                }

 #endif

                continue;
            }

            element = refinedDomain->giveElement(ielem);
            refinedProblem->giveElementCharacteristicMatrix(mat, ielem, STIFFNESS_TYPE, refinedTStep, refinedDomain);

            this->extractVectorFrom(element, errorSolution, errorVector, dofs, refinedTStep);
            tmpVector.beProductOf(mat, errorVector);
            exactENorm += tmpVector.dotProduct(errorVector); // Coarse and fine are identical? Also, unused.
            fineENorm = tmpVector.dotProduct(errorVector);
            coarseENorm += fineENorm;

            this->extractVectorFrom(element, coarseSolution, coarseVector, dofs, refinedTStep);
            tmpVector.beProductOf(mat, coarseVector);
            coarseUNorm += tmpVector.dotProduct(coarseVector);

            mixedNorm += tmpVector.dotProduct(errorVector);

            this->extractVectorFrom(element, fineSolution, fineVector, dofs, refinedTStep);
            tmpVector.beProductOf(mat, fineVector);
            fineUNorm += tmpVector.dotProduct(fineVector);

 #ifdef PRINT_ERROR
            exactFineError.at(ielem) = fineENorm;
            if ( ++pos == nelems ) {
                exactCoarseError.at(elemId) = coarseENorm;
                if ( ielem < localElemId ) {
                    elemId++;
                    element = domain->giveElement(elemId);
                    dim = element->giveSpatialDimension();
                    nelems = this->refineLevel + 1;
                    while ( --dim ) {
                        nelems *= this->refineLevel + 1;
                    }

                    nelems *= element->giveNumberOfNodes();
                    coarseENorm = 0.0;
                    pos = 0;
                }
            }

 #endif
        }
    } else {
        _error("solveRefinedWholeProblem: Unsupported norm type");
    }

 #ifdef TIME_INFO
    :: getRelativeUtime(et_error, st_error);
    :: getRelativeUtime(et_total, st_total);

    OOFEM_LOG_DEBUG( "HEE info: whole 0: user time total %.2f s (setup %.2f s, init %.2f s, solve %.2f s, error %.2f s)\n",
                    ( double ) ( et_total.tv_sec + et_total.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_setup.tv_sec + et_setup.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_init.tv_sec + et_init.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_solve.tv_sec + et_solve.tv_usec / ( double ) OOFEM_USEC_LIM ),
                    ( double ) ( et_error.tv_sec + et_error.tv_usec / ( double ) OOFEM_USEC_LIM ) );
 #endif

    delete refinedProblem;
}
#endif


// pokud toto neni dostatecne obecne, da se funkce extractVectorFrom
// do HuertaErrorEstimatorInterface a kazdy prvek si ji muze predefinovat

void
HuertaErrorEstimator :: extractVectorFrom(Element *element, FloatArray &vector, FloatArray &answer,
                                          int dofs, TimeStep *tStep)
{
    int inode, idof, pos, p = 0;

    answer.resize(element->giveNumberOfDofManagers() * dofs);
    for ( inode = 1; inode <= element->giveNumberOfNodes(); inode++ ) {
        pos = ( element->giveDofManager(inode)->giveNumber() - 1 ) * dofs;
        for ( idof = 1; idof <= dofs; idof++ ) {
            answer.at(++p) = vector.at(pos + idof);
        }
    }
}



void
HuertaErrorEstimator :: setupRefinedProblemProlog(const char *problemName, int problemId, IntArray &localNodeIdArray,
                                                  int nodes, int elems, int csects, int mats, int loads, int ltfuncs,
                                                  IntArray &controlNode, IntArray &controlDof, TimeStep *tStep)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    EngngModel *problem = this->domain->giveEngngModel();
    InputRecord *ir;
    std :: string str;
    int i, nmstep, nsteps = 0;
    int ddltf = 0, ddmSize = 0, ddvSize = 0, hpcSize = 0, hpcwSize = 0, skipUpdate = 0, renumber = 1;
    int controlMode = 0, hpcMode = 0, stiffMode = 0, maxIter = 30, reqIter = 3, manrmsteps = 0;
    double rtolv, minStepLength = 0.0, initialStepLength, stepLength, psi = 1.0;
    IntArray ddm, hpc;
    FloatArray ddv, hpcw;
    const char *__proc = "setupRefinedProblemProlog"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                           // Required by IR_GIVE_FIELD macro

    char skipUpdateString [ 32 ] = "", parallelFlagString [ 32 ] = "", useContextString [ 32 ] = "";

#if defined ( USE_OUTPUT_FILE ) || defined ( USE_CONTEXT_FILE )
    sprintf(line, "/home/dr/Huerta/%s_%d.out", problemName, problemId);
    skipUpdate = 0;
#else
    sprintf(line, "/dev/null");
    skipUpdate = 1;
#endif

    /* sprintf(skipUpdateString, "skipUpdate %d ", skipUpdate); */

#ifdef USE_CONTEXT_FILE
    int contextOutputStep = 1;
    sprintf(useContextString, "contextOutputStep %d ", contextOutputStep);
#endif
#ifdef __PARALLEL_MODE
    sprintf(parallelFlagString, "parallelFlag 0 ");
#endif

    refinedReader.appendInputString(line);

    sprintf(line, "Refined problem on %s %d", problemName, problemId);
    refinedReader.appendInputString(line);

    switch ( problem->giveClassID() ) {
    case AdaptiveLinearStaticClass:
        sprintf(line, "LinearStatic nsteps 1 renumber %d %s %s %s",
                renumber, skipUpdateString, useContextString, parallelFlagString);
        refinedReader.appendInputString(line);
        break;
    case AdaptiveNonLinearStaticClass:
        nmstep = tStep->giveMetaStepNumber();
        ir = problem->giveMetaStep(nmstep)->giveAttributesRecord();

        IR_GIVE_OPTIONAL_FIELD(ir, stiffMode, IFT_NonLinearStatic_stiffmode, "stiffmode"); //macro
        IR_GIVE_OPTIONAL_FIELD(ir, controlMode, IFT_NonLinearStatic_controlmode, "controlmode"); //macro

        switch ( controlMode ) {
        case 0:
            IR_GIVE_OPTIONAL_FIELD(ir, maxIter, IFT_CylindricalALM_maxiter, "maxiter"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, reqIter, IFT_CylindricalALM_reqiterations, "reqiterations"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, IFT_CylindricalALM_minsteplength, "minsteplength"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, manrmsteps, IFT_CylindricalALM_manrmsteps, "manrmsteps"); // Macro
            IR_GIVE_FIELD(ir, stepLength, IFT_CylindricalALM_steplength, "steplength"); // Macro

            initialStepLength = stepLength;
            IR_GIVE_OPTIONAL_FIELD(ir, initialStepLength, IFT_CylindricalALM_initialsteplength, "initialsteplength"); // Macro

            IR_GIVE_OPTIONAL_FIELD(ir, psi, IFT_CylindricalALM_psi, "psi"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, hpcMode, IFT_CylindricalALM_hpcmode, "hpcmode"); // Macro

            IR_GIVE_OPTIONAL_FIELD(ir, hpc, IFT_CylindricalALM_hpc, "hpc"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, hpcw, IFT_CylindricalALM_hpcw, "hpcw"); // Macro
            IR_GIVE_FIELD(ir, rtolv, IFT_CylindricalALM_rtolv, "rtolv"); //macro

            hpcSize = hpc.giveSize();
            hpcwSize = hpcw.giveSize();
            break;
        case 1:
            IR_GIVE_OPTIONAL_FIELD(ir, maxIter, IFT_NRSolver_maxiter, "maxiter"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, IFT_NRSolver_minsteplength, "minsteplength"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, manrmsteps, IFT_NRSolver_manrmsteps, "manrmsteps"); // Macro

            IR_GIVE_OPTIONAL_FIELD(ir, ddm, IFT_NRSolver_ddm, "ddm"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, ddv, IFT_NRSolver_ddv, "ddv"); // Macro
            IR_GIVE_OPTIONAL_FIELD(ir, ddltf, IFT_NRSolver_ddltf, "ddltf"); // Macro
            IR_GIVE_FIELD(ir, rtolv, IFT_NRSolver_rtolv, "rtolv"); //macro

            ddmSize = ddm.giveSize();
            ddvSize = ddv.giveSize();
            break;
        default:
            _error("setupRefinedProblemProlog: Unsupported control mode");
        }

        if ( problemId != 0 ) {
            int bcSize = controlNode.giveSize(), ddActiveSize = 0;

            if ( controlMode == 1 ) {
                // count original dd active in refined problem
                for ( i = 1; i <= ddmSize; i += 2 ) {
                    if ( localNodeIdArray.at( ddm.at(i) ) == 0 ) {
                        continue;
                    }

                    ddActiveSize++;
                }
            }

            // note: it is impossible to check the solution using the external file;
            // first of all adaptnlinearstatic must be changed to nonlinearstatic to prevent adaptivity;
            // however the most severe problem is that in ddv only zeros are specified !!!
            // it would be desirable to print there coarse mesh displacement

            if ( nmstep == 1 ) {
                // do not specify context because that would result in recursive call to HEE
                // because adaptnlinearstatic analysis type is used

                sprintf(line, "adaptnlinearstatic nsteps %d renumber %d rtolv %e maxiter %d reqiterations %d minsteplength %e stiffmode %d manrmsteps %d equilmc 1 controlmode 1 %s %s ",
                        tStep->giveNumber(), renumber, rtolv, maxIter, reqIter, minStepLength, stiffMode, manrmsteps,
                        skipUpdateString, parallelFlagString);
                str = line;

                // this is not relevant but it is required
                // the refined problem is made adaptive only to enable call to initializeAdaptiveFrom
                sprintf(line, "eetype 3 meshpackage 0 requiredError 0.10 minelemsize 0.01 ");
                str += line;
            } else {
                // if there are more than just a single metastep, it is necessary to produce input with at least the
                // same number of metasteps as is the number of active metastep,
                // because initializeAdaptiveFrom uses a copy constructor for timestep;
                // only the active metastep should be filled by real values

                // do not specify context because that would result in recursive call to HEE
                // because adaptnlinearstatic analysis type is used

                sprintf(line, "adaptnlinearstatic nsteps 1 nmsteps %d renumber %d equilmc 1 controlmode 1 %s %s ",
                        nmstep, renumber, skipUpdateString, parallelFlagString);
                str = line;

                // this is not relevant but it is required
                // the refined problem is made adaptive only to enable call to initializeAdaptiveFrom
                sprintf(line, "eetype 3 meshpackage 0 requiredError 0.10 minelemsize 0.01 ");
                str += line;

                refinedReader.appendInputString(str);

                // IMPORTANT: the total number of steps should be equal to the current step number
                //            (to enable skipUpdate)
                //            therefore the number of steps in the current metastep is modified

                // create fictitious metasteps
                for ( i = 1; i < nmstep; i++ ) {
                    sprintf(line, "metastep %d nsteps %d rtolv %e",
                            i, problem->giveMetaStep(i)->giveNumberOfSteps(), rtolv);
                    refinedReader.appendInputString(line);
                    nsteps += problem->giveMetaStep(i)->giveNumberOfSteps();
                }

                // create active metastep
                sprintf(line, "metastep %d nsteps %d rtolv %e maxiter %d reqiterations %d minsteplength %e stiffmode %d manrmsteps %d equilmc 1 controlmode 1 ",
                        nmstep, tStep->giveNumber() - nsteps,
                        rtolv, maxIter, reqIter, minStepLength, stiffMode, manrmsteps);
                str = line;
            }

            sprintf( line, "ddm %d ", 2 * ( bcSize + ddActiveSize ) );
            str += line;

            // apply refine problem bc as dd
            for ( i = 1; i <= bcSize; i++ ) {
                sprintf( line, "%d %d ", controlNode.at(i), controlDof.at(i) );
                str += line;
            }

            if ( controlMode == 1 ) {
                // apply active original dd
                for ( i = 1; i <= ddmSize; i += 2 ) {
                    if ( localNodeIdArray.at( ddm.at(i) ) == 0 ) {
                        continue;
                    }

                    sprintf( line, "%d %d ", -localNodeIdArray.at( ddm.at(i) ), ddm.at(i + 1) );
                    str += line;
                }
            }

            sprintf(line, "ddv %d ", bcSize + ddActiveSize);
            str += line;

            // apply refined problem bc values
            for ( i = 1; i <= bcSize; i++ ) {
                sprintf(line, "%3.1f ", 0.0); // no further increment required, just to recover equilibrium of mapped state
                str += line;
            }

            if ( controlMode == 1 ) {
                // apply active original dd values
                for ( i = 1; i <= ddvSize; i++ ) {
                    if ( localNodeIdArray.at( ddm.at(2 * i - 1) ) == 0 ) {
                        continue;
                    }

                    sprintf( line, "%e ", ddv.at(i) );
                    str += line;
                }
            }

            // use last ltf to control dd
            sprintf(line, "ddltf %d ", ltfuncs);
            str += line;

            refinedReader.appendInputString(str);
        } else {
            // it makes not too much sense to solve exact problem from beginning if adaptive restart is used
            // because the mesh may be already derefined in regions of no interest (but anyway ...)
            // however if adaptive restart is applied, number of current step does not correspond to the time
            // (step number = time + 1) because step number was encreased when recovering equilibrium at the last time;
            // therefore problem -> giveCurrentStep() -> giveNumber() is not used and the number of steps is
            // recovered from the current time

            // do not prescribe skipUpdate here !!!

            // HUHU dodelat metasteps, dodelat paralelni zpracovani
            sprintf(line, "nonlinearstatic nsteps %d renumber %d rtolv %e maxiter %d reqiterations %d minsteplength %e stiffmode %d manrmsteps %d equilmc 1 controlmode %d %s ",
                    ( int ) ( problem->giveCurrentStep()->giveTargetTime() + 1.5 ), renumber,
                    rtolv, maxIter, reqIter, minStepLength, stiffMode, manrmsteps, controlMode, useContextString);
            str = line;

            switch ( controlMode ) {
            case 0:
                sprintf(line, "stepLength %e initialStepLength %e psi %e hpcmode %d ",
                        stepLength, initialStepLength, psi, hpcMode);
                str += line;

                if ( hpcSize != 0 ) {
                    // apply all original hpc
                    sprintf(line, "hpc %d ", hpcSize);
                    str += line;
                    for ( i = 1; i <= hpcSize; i += 2 ) {
                        sprintf( line, "%d %d ", -localNodeIdArray.at( hpc.at(i) ), hpc.at(i + 1) );
                        str += line;
                    }
                }

                if ( hpcwSize != 0 ) {
                    sprintf(line, "hpcw %d ", hpcwSize);
                    str += line;
                    for ( i = 1; i <= hpcwSize; i++ ) {
                        sprintf( line, "%e ", hpcw.at(i) );
                        str += line;
                    }
                }

                break;
            case 1:
                if ( ddmSize != 0 ) {
                    // apply all original dd
                    sprintf(line, "ddm %d ", ddmSize);
                    str += line;
                    for ( i = 1; i <= ddmSize; i += 2 ) {
                        sprintf( line, "%d %d ", -localNodeIdArray.at( ddm.at(i) ), ddm.at(i + 1) );
                        str += line;
                    }
                }

                if ( ddvSize != 0 ) {
                    sprintf(line, "ddv %d ", ddvSize);
                    str += line;
                    for ( i = 1; i <= ddvSize; i++ ) {
                        sprintf( line, "%e ", ddv.at(i) );
                        str += line;
                    }
                }

                // use the original ltf to control dd
                sprintf(line, "ddltf %d ", ddltf);
                str += line;
                break;
            }

            refinedReader.appendInputString(str);
        }

        break;
    default:
        _error("setupRefinedProblemProlog: Unsupported analysis type");
    }

    switch ( this->domain->giveDomainType() ) {
    case _1dTrussMode:
        sprintf(line, "domain 1dtruss");
        break;
    case _PlaneStrainMode:
        sprintf(line, "domain planestrain");
        break;
    case _2dPlaneStressMode:
        sprintf(line, "domain 2dPlaneStress");
        break;
    case _3dMode:
        sprintf(line, "domain 3d");
        break;
    default:
        _error("solveRefinedElementProblem: Unsupported domain type");
    }

    refinedReader.appendInputString(line);

#ifdef USE_OUTPUT_FILE
    sprintf(line, "OutputManager tstep_all dofman_all element_all");
#else
    sprintf(line, "OutputManager");
#endif

    refinedReader.appendInputString(line);

    sprintf(line, "%d %d %d %d %d %d", nodes, elems, mats, csects, loads, ltfuncs);
    refinedReader.appendInputString(line);
}



void
HuertaErrorEstimator :: setupRefinedProblemEpilog1(int csects, int mats, int loads, int nlbarriers)
{
    Domain *domain = this->domain;
    std :: string str;
    int i;

    /* copy csects, mats, loads */

    for ( i = 1; i <= csects; i++ ) {
        domain->giveCrossSection(i)->giveInputRecordString(str);
        refinedReader.appendInputString(str);
    }

    for ( i = 1; i <= mats; i++ ) {
        domain->giveMaterial(i)->giveInputRecordString(str);
        refinedReader.appendInputString(str);
    }

    for ( i = 1; i <= loads; i++ ) {
        domain->giveLoad(i)->giveInputRecordString(str);
        refinedReader.appendInputString(str);
    }
}



void
HuertaErrorEstimator :: setupRefinedProblemEpilog2(int ltfuncs)
{
    char line [ OOFEM_MAX_LINE_LENGTH + 1 ];
    Domain *domain = this->domain;
    std :: string str;
    int i;

    /* copy tfuncs */

    for ( i = 1; i <= ltfuncs; i++ ) {
        domain->giveLoadTimeFunction(i)->giveInputRecordString(str);
        refinedReader.appendInputString(str);
    }

    if ( this->mode == HEE_nlinear ) {
        sprintf(line, "heavisideltf %d origin %e value 1.0", ltfuncs + 1,
                this->domain->giveEngngModel()->giveCurrentStep()->giveTargetTime() - 0.1);
        refinedReader.appendInputString(line);
    }
}
} // end namespace oofem

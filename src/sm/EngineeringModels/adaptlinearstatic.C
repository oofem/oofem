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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "sm/EngineeringModels/adaptlinearstatic.h"
#include "remeshingcrit.h"
#include "mesherinterface.h"
#include "errorestimator.h"
#include "domain.h"
#include "timestep.h"
#include "metastep.h"
#include "timestepcontroller.h"
#include "sparsemtrx.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "outputmanager.h"
#include "logger.h"

namespace oofem {
REGISTER_EngngModel(AdaptiveLinearStatic);

AdaptiveLinearStatic :: AdaptiveLinearStatic(int i, EngngModel *master) : LinearStatic(i, master),
    meshPackage(MPT_T3D), adaptFlag(false), maxAdaptSteps(5), remeshingStrategy(NoRemeshing_RS)
{ }


AdaptiveLinearStatic :: ~AdaptiveLinearStatic() { }


TimeStep *
AdaptiveLinearStatic :: giveNextStep()
{
    // single load case: all (adaptive) solution steps have the same target time
    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std :: make_unique< TimeStep >(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., 1., 0);
    }
    previousStep = std :: move(currentStep);
    currentStep = std :: make_unique< TimeStep >(previousStep->giveNumber() + 1, this, 1, 1.0, 1.0,
                                                 previousStep->giveSolutionStateCounter() + 1);
    return currentStep.get();
}


void
AdaptiveLinearStatic :: solveYourself()
{
    if ( this->isParallel() ) {
        OOFEM_ERROR("parallel mode not supported");
    }

    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    auto activeMStep = this->giveMetaStep(1);
    timeStepController->setCurrentMetaStepNumber(0);
    timeStepController->initMetaStepAttributes(activeMStep);

    for ( int iadapt = 0; ; iadapt++ ) {
        this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

        this->preInitializeNextStep();
        TimeStep *tStep = this->giveNextStep();
        tStep->setMetaStepNumber(1);

        if ( this->requiresEquationRenumbering(tStep) ) {
            this->forceEquationNumbering();
        }

        OOFEM_LOG_INFO("\nAdaptive step %d: %d nodes, %d elements, %d equations\n", iadapt,
                       this->giveDomain(1)->giveNumberOfDofManagers(), this->giveDomain(1)->giveNumberOfElements(),
                       this->giveNumberOfDomainEquations( 1, this->giveEquationNumbering() ) );

        this->initializeYourself(tStep);
        this->solveYourselfAt(tStep);
        // update + error estimation
        this->updateYourself(tStep);

        this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        double steptime = this->giveSolutionStepTime();
        tStep->solutionTime = steptime;

        this->terminate(tStep);

        OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n", tStep->giveNumber(), steptime);
        if ( !suppressOutput ) {
            fprintf(this->giveOutputStream(), "\nUser time consumed by solution step %d: %.3f [s]\n\n", tStep->giveNumber(), steptime);
        }

        if ( remeshingStrategy == NoRemeshing_RS ) {
            OOFEM_LOG_INFO("Error estimate acceptable, no remeshing required\n");
            break;
        }

        // ask mesher for new mesh (writes at least the mesh density suggestion)
        Domain *dNew = nullptr;
        MesherInterface :: returnCode result =
            mesher->createMesh(tStep, 1, this->giveDomain(1)->giveSerialNumber() + 1, & dNew);

        if ( result == MesherInterface :: MI_FAILED ) {
            OOFEM_ERROR("createMesh failed");
        } else if ( result == MesherInterface :: MI_NEEDS_EXTERNAL_ACTION ) {
            OOFEM_LOG_INFO("Remeshing required, mesh density file created\n");
            break;
        }

        if ( !adaptFlag || iadapt >= maxAdaptSteps ) {
            // new mesh not used
            delete dNew;
            if ( adaptFlag ) {
                OOFEM_WARNING("Maximum number of adaptive steps (%d) reached", maxAdaptSteps);
            }
            break;
        }

        this->replaceDomain(dNew);
    }
}


void
AdaptiveLinearStatic :: replaceDomain(Domain *dNew)
{
    dNew->setNumber(1);
    // deletes the old domain
    this->setDomain(1, dNew);

    // new discretization -> renumber equations, reassemble stiffness matrix
    this->equationNumberingCompleted = 0;
    this->forceEquationNumbering();
    this->stiffnessMatrix.reset();
    this->initFlag = 1;

    // relink numerical method, error estimator, mesher, export modules to new domain
    this->updateDomainLinks();
}


void
AdaptiveLinearStatic :: updateYourself(TimeStep *tStep)
{
    LinearStatic :: updateYourself(tStep);
    // evaluate error of the reached solution
    this->defaultErrEstimator->estimateError(temporaryEM, tStep);
    this->remeshingStrategy = this->defaultErrEstimator->giveRemeshingCrit()->giveRemeshingStrategy(tStep);
}


void
AdaptiveLinearStatic :: printOutputAt(FILE *file, TimeStep *tStep)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return;
    }

    LinearStatic :: printOutputAt(file, tStep);
    fprintf(file, "\nAdaptive step %d: %d nodes, %d elements\n", tStep->giveNumber() - 1,
            this->giveDomain(1)->giveNumberOfDofManagers(), this->giveDomain(1)->giveNumberOfElements() );
    this->defaultErrEstimator->printOutputAt(file, tStep);
}


int
AdaptiveLinearStatic :: initializeAdaptive(int tStepNumber)
{
    // The whole adaptive loop is performed within single run (see solveYourself),
    // restart is not needed due to linear character of the problem.
    return 1;
}


void AdaptiveLinearStatic :: restoreContext(DataStream &stream, ContextMode mode)
{
    LinearStatic :: restoreContext(stream, mode);
}

void
AdaptiveLinearStatic :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    LinearStatic :: initializeFrom(ir);

    int meshPackageId = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, meshPackageId, _IFT_AdaptiveLinearStatic_meshpackage);
    // same numbering as AdaptiveNonLinearStatic (0-T3D, 1-Targe2, 2-Freem, 3-Subdivision, 4-Gmsh)
    meshPackage = ( MeshPackageType ) meshPackageId;

    int adapt = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, adapt, _IFT_AdaptiveLinearStatic_adapt);
    adaptFlag = adapt != 0;
    IR_GIVE_OPTIONAL_FIELD(ir, maxAdaptSteps, _IFT_AdaptiveLinearStatic_maxadaptsteps);

    if ( !this->defaultErrEstimator ) {
        throw ValueInputException(ir, "eetype", "error estimator not defined");
    }

    // mesher specific parameters are part of this record
    mesher = classFactory.createMesherInterface( meshPackage, this->giveDomain(1) );
    if ( !mesher ) {
        throw ValueInputException(ir, _IFT_AdaptiveLinearStatic_meshpackage, "unknown mesh package");
    }
    mesher->initializeFrom(ir);
}


void
AdaptiveLinearStatic :: updateDomainLinks()
{
    LinearStatic :: updateDomainLinks();
    // associate ee to possibly newly restored mesh
    this->defaultErrEstimator->setDomain( this->giveDomain(1) );
    if ( mesher ) {
        mesher->setDomain( this->giveDomain(1) );
    }
}
} // end namespace oofem

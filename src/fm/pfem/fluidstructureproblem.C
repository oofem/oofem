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

#include "fluidstructureproblem.h"
#include "engngm.h"
#include "timestep.h"
#include "function.h"
#include "metastep.h"
#include "exportmodulemanager.h"
#include "mathfem.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "verbose.h"
#include "classfactory.h"

#include "pfem.h"
#include "interactionpfemparticle.h"
#include "../../sm/EngineeringModels/diidynamic.h"
#include "../../sm/EngineeringModels/nlineardynamic.h"

#include <stdlib.h>

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel(FluidStructureProblem);

FluidStructureProblem :: FluidStructureProblem(int i, EngngModel *_master) : StaggeredProblem(i, _master)
{
    ndomains = 1; // domain is needed to store the time step ltf

    iterationNumber = 0;
    this->setRenumberFlag();
}

FluidStructureProblem :: ~FluidStructureProblem()
{}


IRResultType
FluidStructureProblem :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StaggeredProblem :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    maxiter = 50;
    IR_GIVE_OPTIONAL_FIELD(ir, maxiter, _IFT_FluidStructureProblem_maxiter);

    rtolv = 1.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, rtolv, _IFT_FluidStructureProblem_rtolv);

    rtolp = 1.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, rtolp, _IFT_FluidStructureProblem_rtolp);

    return IRRT_OK;
}

void
FluidStructureProblem :: initializeYourself(TimeStep *tStep)
{
    for ( auto &emodel: emodelList ) {
        emodel->initializeYourself(tStep);

        DIIDynamic *dynamicProblem = dynamic_cast< DIIDynamic * >( emodel.get() );
        if ( dynamicProblem ) {
            this->giveCurrentStep()->setTimeDiscretization( dynamicProblem->giveInitialTimeDiscretization() );
        }
        NonLinearDynamic *nonlinearDynamicProblem = dynamic_cast< NonLinearDynamic * >( emodel.get() );
        if ( nonlinearDynamicProblem ) {
            this->giveCurrentStep()->setTimeDiscretization( nonlinearDynamicProblem->giveInitialTimeDiscretization() );
        }
    }
}

void
FluidStructureProblem :: solveYourselfAt(TimeStep *stepN)
{
    PFEM *pfemProblem = dynamic_cast< PFEM * >( this->giveSlaveProblem(1) );
    Domain *pfemDomain = pfemProblem->giveDomain(1);

#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %5d, time %e]\n", stepN->giveNumber(), stepN->giveTargetTime() );
#endif

    FloatArray currentInteractionParticlesVelocities;
    FloatArray currentInteractionParticlesPressures;
    FloatArray previousInteractionParticlesVelocities;
    FloatArray previousInteractionParticlesPressures;

    if ( stepN->isTheFirstStep() ) {
        int ndman = pfemDomain->giveNumberOfDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            if ( dynamic_cast< InteractionPFEMParticle * >( pfemDomain->giveDofManager(i) ) ) {
                interactionParticles.followedBy(i, 10);
            }
        }
    }

    currentInteractionParticlesVelocities.resize( 2 * interactionParticles.giveSize() );
    previousInteractionParticlesVelocities.resize( 2 * interactionParticles.giveSize() );
    currentInteractionParticlesPressures.resize( interactionParticles.giveSize() );
    previousInteractionParticlesPressures.resize( interactionParticles.giveSize() );

    double velocityDifference = 1.0;
    double pressureDifference = 1.0;

    iterationNumber = 0;
    while ( ( ( velocityDifference > rtolv ) || ( pressureDifference > rtolp ) ) && iterationNumber < maxiter ) {
        previousInteractionParticlesPressures = currentInteractionParticlesPressures;
        previousInteractionParticlesVelocities = currentInteractionParticlesVelocities;

        for ( auto &emodel: emodelList ) {
            emodel->solveYourselfAt(stepN);
        }

        for ( int i = 1; i <= interactionParticles.giveSize(); i++ ) {
            currentInteractionParticlesPressures.at(i) = pfemProblem->giveUnknownComponent( VM_Total, stepN, pfemDomain, pfemDomain->giveDofManager( interactionParticles.at(i) )->giveDofWithID(P_f) );
            InteractionPFEMParticle *interactionParticle = dynamic_cast< InteractionPFEMParticle * >( pfemDomain->giveDofManager( interactionParticles.at(i) ) );
            FloatArray velocities;
            interactionParticle->giveCoupledVelocities(velocities, stepN);
            currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 1) = velocities.at(1);
            currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 2) = velocities.at(2);
        }

        pressureDifference = 0.0;
        velocityDifference = 0.0;
        for ( int i = 1; i <= currentInteractionParticlesPressures.giveSize(); i++ ) {
            pressureDifference += ( currentInteractionParticlesPressures.at(i) - previousInteractionParticlesPressures.at(i) ) * ( currentInteractionParticlesPressures.at(i) - previousInteractionParticlesPressures.at(i) );
            velocityDifference += ( currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 1) - previousInteractionParticlesVelocities.at(2 * ( i - 1 ) + 1) ) * ( currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 1) - previousInteractionParticlesVelocities.at(2 * ( i - 1 ) + 1) );
            velocityDifference += ( currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 2) - previousInteractionParticlesVelocities.at(2 * ( i - 1 ) + 2) ) * ( currentInteractionParticlesVelocities.at(2 * ( i - 1 ) + 2) - previousInteractionParticlesVelocities.at(2 * ( i - 1 ) + 2) );
        }
        pressureDifference = sqrt(pressureDifference);
        velocityDifference = sqrt(velocityDifference);

        if ( iterationNumber > 0 ) {
            pressureDifference /= previousInteractionParticlesPressures.computeNorm();
            velocityDifference /= previousInteractionParticlesVelocities.computeNorm();
        }

        OOFEM_LOG_RELEVANT("%3d %le %le\n", iterationNumber++, pressureDifference, velocityDifference);
    }
    if ( iterationNumber > maxiter ) {
        OOFEM_ERROR("Maximal fluid-structure interface iteration count exceded");
    }
    stepN->incrementStateCounter();
}

void
FluidStructureProblem :: preInitializeNextStep()
{
    for ( auto &emodel: emodelList ) {
        emodel->preInitializeNextStep();
    }
}

#if 0
void
FluidStructureProblem :: postInitializeCurrentStep()
{
    for ( int i = 1; i <= nModels; i++ ) {
        DIIDynamic *dynamicProblem = dynamic_cast< DIIDynamic * >( this->giveSlaveProblem(i) );
        if ( dynamicProblem ) {
            this->giveCurrentStep()->setTimeDiscretization( dynamicProblem->giveInitialTimeDiscretization() );
        }
    }
}
#endif
} // end namespace oofem

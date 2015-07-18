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

#include "structuralmaterialevaluator.h"
#include "inputrecord.h"
#include "timestep.h"
#include "domain.h"
#include "gausspoint.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "function.h"
#include "classfactory.h"

#include <fstream>

namespace oofem {
REGISTER_EngngModel(StructuralMaterialEvaluator);

StructuralMaterialEvaluator :: StructuralMaterialEvaluator(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->ndomains = 1;
}

StructuralMaterialEvaluator :: ~StructuralMaterialEvaluator()
{ }

IRResultType StructuralMaterialEvaluator :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->deltaT, _IFT_StructuralMaterialEvaluator_deltat);
    IR_GIVE_FIELD(ir, this->numberOfSteps, _IFT_StructuralMaterialEvaluator_numberOfTimeSteps);

    //IR_GIVE_FIELD(ir, this->ndim, _IFT_StructuralMaterialEvaluator_nDimensions);

    IR_GIVE_FIELD(ir, this->cmpntFunctions, _IFT_StructuralMaterialEvaluator_componentFunctions);
    IR_GIVE_FIELD(ir, this->sControl, _IFT_StructuralMaterialEvaluator_stressControl);

    tolerance = 1.0;
    if ( this->sControl.giveSize() > 0 ) {
        IR_GIVE_FIELD(ir, this->tolerance, _IFT_StructuralMaterialEvaluator_tolerance);
    }

    IR_GIVE_FIELD(ir, this->vars, _IFT_StructuralMaterialEvaluator_outputVariables);

    // Compute the strain control (everything not controlled by stress)
    for ( int i = 1; i <= 6; ++i ) {
        if ( !sControl.contains(i) ) {
            eControl.followedBy(i);
        }
    }

    return IRRT_OK;
}


void StructuralMaterialEvaluator :: solveYourself()
{
    Domain *d = this->giveDomain(1);

    MaterialMode mode = _3dMat;
    FloatArray initialStrain(6);
    gps.clear();
    gps.reserve(d->giveNumberOfMaterialModels());
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        std :: unique_ptr< GaussPoint > gp(new GaussPoint(nullptr, i, FloatArray(0), 1, mode));
        gps.emplace_back( std :: move(gp) );
        // Initialize the strain vector;
        StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( d->giveMaterial(i)->giveStatus( gps[i-1].get() ) );
        status->letStrainVectorBe(initialStrain);
    }

    std :: string outname = this->giveOutputBaseFileName() + ".matdata";
    this->outfile.open( outname.c_str() );

    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    TimeStep *tStep = giveNextStep();

    // Note, strain == strain-rate (kept as strain for brevity)
    int maxiter = 100; // User input?
    FloatArray stressC, deltaStrain, strain, stress, res;
    stressC.resize( sControl.giveSize() );
    res.resize( sControl.giveSize() );

    FloatMatrix tangent, reducedTangent;
    for ( int istep = 1; istep <= this->numberOfSteps; ++istep ) {
        this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        for ( int imat = 1; imat <= d->giveNumberOfMaterialModels(); ++imat ) {
            GaussPoint *gp = gps[imat-1].get();
            StructuralMaterial *mat = static_cast< StructuralMaterial * >( d->giveMaterial(imat) );
            StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( mat->giveStatus(gp) );

            strain = status->giveStrainVector();
            // Update the controlled parts
            for ( int j = 1; j <= eControl.giveSize(); ++j ) {
                int p = eControl.at(j);
                strain.at(p) = d->giveFunction( cmpntFunctions.at(p) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }

            for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                int p = sControl.at(j);
                stressC.at(j) = d->giveFunction( cmpntFunctions.at(p) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }

            for ( int iter = 1; iter < maxiter; iter++ ) {
                mat->giveRealStressVector_3d(stress, gp, strain, tStep);
                for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                    res.at(j) = stressC.at(j) - stress.at( sControl.at(j) );
                }

                OOFEM_LOG_RELEVANT( "Time step: %d, Material %d, Iteration: %d,  Residual = %e\n", istep, imat, iter, res.computeNorm() );
                if ( res.computeNorm() <= tolerance ) {
                    break;
                }

                mat->give3dMaterialStiffnessMatrix(tangent, TangentStiffness, gp, tStep);
                if ( res.giveSize() > 0 ) {
                    // Pick out the stress-controlled part;
                    reducedTangent.beSubMatrixOf(tangent, sControl, sControl);

                    // Update stress-controlled part of the strain
                    reducedTangent.solveForRhs(res, deltaStrain);
                    for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                        strain.at( sControl.at(j) ) += deltaStrain.at(j);
                    }
                }
            }

            if ( res.computeNorm() > tolerance ) {
                OOFEM_WARNING("Residual did not converge!");
            }

            // This material model has converged, so we update it and go on to the next.
            gp->updateYourself(tStep);
        }

        this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        this->doStepOutput(tStep);
        tStep = giveNextStep();
    }

    this->timer.stopTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    this->outfile.close();
}

int StructuralMaterialEvaluator :: checkConsistency()
{
    Domain *d =  this->giveDomain(1);
    for ( auto &mat : d->giveMaterials() ) {
        if ( !dynamic_cast< StructuralMaterial * >( mat.get() ) ) {
            OOFEM_LOG_ERROR("Material %d is not a StructuralMaterial", mat->giveNumber());
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}

void StructuralMaterialEvaluator :: doStepOutput(TimeStep *tStep)
{
    FloatArray outputValue;
    Domain *d = this->giveDomain(1);
    if ( tStep->isTheFirstStep() ) {
        this->outfile << "# Time";
        for ( int var : this->vars ) {
            this->outfile << ", " << __InternalStateTypeToString( ( InternalStateType ) var );
        }

        this->outfile << '\n';
    }

    outfile << tStep->giveIntrinsicTime();
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        Material *mat = d->giveMaterial(i);
        for ( int var : this->vars ) {
            mat->giveIPValue(outputValue, gps[i-1].get(), ( InternalStateType ) var, tStep);
            outfile << " " << outputValue;
        }
    }

    outfile << std :: endl;
}

TimeStep *StructuralMaterialEvaluator :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        //currentStep.reset( new TimeStep(*giveSolutionStepWhenIcApply()) );
        currentStep.reset( new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 1, 0., this->deltaT, 0) );
    }
    previousStep = std :: move(currentStep);
    currentStep.reset( new TimeStep(*previousStep, this->deltaT) );

    return currentStep.get();
}
} // end namespace oofem

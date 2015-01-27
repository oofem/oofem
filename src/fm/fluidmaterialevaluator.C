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

#include "fluidmaterialevaluator.h"
#include "inputrecord.h"
#include "timestep.h"
#include "domain.h"
#include "gausspoint.h"
#include "fluiddynamicmaterial.h"
#include "function.h"
#include "classfactory.h"

#include <fstream>

namespace oofem {
REGISTER_EngngModel(FluidMaterialEvaluator);

FluidMaterialEvaluator :: FluidMaterialEvaluator(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->ndomains = 1;
}

FluidMaterialEvaluator :: ~FluidMaterialEvaluator()
{ }

IRResultType FluidMaterialEvaluator :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    this->deltaT = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->deltaT, _IFT_FluidMaterialEvaluator_deltat);
    IR_GIVE_FIELD(ir, this->numberOfSteps, _IFT_FluidMaterialEvaluator_numberOfTimeSteps);

    IR_GIVE_FIELD(ir, this->ndim, _IFT_FluidMaterialEvaluator_nDimensions);

    IR_GIVE_FIELD(ir, this->cmpntFunctions, _IFT_FluidMaterialEvaluator_componentFunctions);
    IR_GIVE_FIELD(ir, this->volFunction, _IFT_FluidMaterialEvaluator_volFunction);
    IR_GIVE_FIELD(ir, this->sControl, _IFT_FluidMaterialEvaluator_stressControl);
    IR_GIVE_FIELD(ir, this->pressureControl, _IFT_FluidMaterialEvaluator_pressureControl);

    IR_GIVE_FIELD(ir, this->vars, _IFT_FluidMaterialEvaluator_outputVariables);

    // Compute the strain control (everything not controlled by stress)
    int components = ( ndim * ( ndim + 1 ) ) / 2;
    for ( int i = 1; i <= components; ++i ) {
        if ( !sControl.contains(i) ) {
            eControl.followedBy(i);
        }
    }

    return IRRT_OK;
}


void FluidMaterialEvaluator :: solveYourself()
{
    Domain *d = this->giveDomain(1);


    MaterialMode mode;
    if ( ndim == 1 ) {
        OOFEM_ERROR("1d flow not supported (should be added)")
        //mode = _1dFlow;
        mode = _Unknown;
    } else if ( ndim == 2 ) {
        mode = _2dFlow;
    } else {
        mode = _3dFlow;
    }

    int components = ( ndim * ( ndim + 1 ) ) / 2;
    FloatArray initialStrain(components);
    initialStrain.zero();
    gps.clear();
    gps.reserve(d->giveNumberOfMaterialModels());
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        std :: unique_ptr< GaussPoint > gp(new GaussPoint(nullptr, i, FloatArray(), 1, mode));
        gps.emplace_back( std :: move(gp));
        // Initialize the strain vector;
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( d->giveMaterial(i) );
        FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus * >( mat->giveStatus( &*gps[i-1] ) );
        status->letDeviatoricStrainRateVectorBe(initialStrain);
    }

    std :: string outname = this->giveOutputBaseFileName() + ".matdata";
    this->outfile.open( outname.c_str() );

    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    TimeStep *tStep = giveNextStep();

    // Note, strain == strain-rate (kept as strain for brevity)
    int maxiter = 100; // User input?
    double tolerance = 1.e-6; // Needs to be normalized somehow, or user input
    double strainVol, pressure, strainVolC = 0.;
    FloatArray stressDevC, deltaStrain, strainDev, stressDev, res;
    stressDevC.resize( sControl.giveSize() );
    res.resize( sControl.giveSize() );

    FloatMatrix tangent, reducedTangent;
    FloatArray dsdp, dedd;
    double dedp;
    for ( int istep = 1; istep <= this->numberOfSteps; ++istep ) {
        this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        for ( int imat = 1; imat <= d->giveNumberOfMaterialModels(); ++imat ) {
            GaussPoint *gp = &*gps[imat-1];
            FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( d->giveMaterial(imat) );
            FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus * >( mat->giveStatus(gp) );

            strainDev = status->giveDeviatoricStrainRateVector();
            pressure = 0.; ///@todo We should ask the material model for this initial guess.
            // Update the controlled parts
            for ( int j = 1; j <= eControl.giveSize(); ++j ) {
                int p = eControl.at(j);
                strainDev.at(p) = d->giveFunction( cmpntFunctions.at(p) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }

            for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                int p = sControl.at(j);
                stressDevC.at(j) = d->giveFunction( cmpntFunctions.at(p) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }

            if ( pressureControl ) {
                pressure = d->giveFunction(volFunction)->evaluateAtTime( tStep->giveIntrinsicTime() );
            } else {
                strainVolC = d->giveFunction(volFunction)->evaluateAtTime( tStep->giveIntrinsicTime() );
            }

            for ( int iter = 1; iter < maxiter; iter++ ) {
                mat->computeDeviatoricStressVector(stressDev, strainVol, gp, strainDev, pressure, tStep);
                for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                    res.at(j) = stressDevC.at(j) - stressDev.at( sControl.at(j) );
                }
                double resVol = 0.;
                if ( !pressureControl ) {
                    resVol = strainVolC - strainVol;
                }

                OOFEM_LOG_RELEVANT( "Time step: %d, Material %d, Iteration: %d,  Residual = %e, Resvol = %e\n", istep, imat, iter, res.computeNorm(), resVol );
                if ( res.computeNorm() <= tolerance && resVol <= tolerance ) { ///@todo More flexible control of convergence needed. Something relative?
                    break;
                }

                mat->giveStiffnessMatrices(tangent, dsdp, dedd, dedp, TangentStiffness, gp, tStep);
                if ( res.giveSize() > 0 ) {
                    // Add mean part to make it invertible
                    double norm = tangent.computeFrobeniusNorm();
                    for ( int i = 1; i <= this->ndim; ++i ) {
                        for ( int j = 1; j <= this->ndim; ++j ) {
                            tangent.at(i, j) += norm;
                        }
                    }

                    // Pick out the stress-controlled part;
                    reducedTangent.beSubMatrixOf(tangent, sControl, sControl);

                    // Update stress-controlled part of the strain
                    reducedTangent.solveForRhs(res, deltaStrain);
                    for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                        strainDev.at( sControl.at(j) ) += deltaStrain.at(j);
                    }
                }
                if ( !pressureControl ) {
                    pressure -= resVol/dedp;
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

int FluidMaterialEvaluator :: checkConsistency()
{
    Domain *d =  this->giveDomain(1);
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        if ( !dynamic_cast< FluidDynamicMaterial * >( d->giveMaterial(i) ) ) {
            return 0;
        }
    }

    return EngngModel :: checkConsistency();
}

void FluidMaterialEvaluator :: doStepOutput(TimeStep *tStep)
{
    FloatArray outputValue;
    Domain *d = this->giveDomain(1);
    if ( tStep->isTheFirstStep() ) {
        this->outfile << "# Time";
        for ( int var: this->vars ) {
            this->outfile << ", " << __InternalStateTypeToString( ( InternalStateType ) var );
        }

        this->outfile << '\n';
    }

    outfile << tStep->giveIntrinsicTime();
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        GaussPoint *gp = &*gps[i-1];
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( d->giveMaterial(i) );
        for ( int j = 1; j <= this->vars.giveSize(); ++j ) {
            mat->giveIPValue(outputValue, gp, ( InternalStateType ) this->vars.at(j), tStep);
            outfile << " " << outputValue;
        }
    }

    outfile << std :: endl;
}

TimeStep *FluidMaterialEvaluator :: giveNextStep()
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

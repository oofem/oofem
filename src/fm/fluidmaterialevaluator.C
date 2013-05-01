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

#include "fluidmaterialevaluator.h"
#include "inputrecord.h"
#include "timestep.h"
#include "domain.h"
#include "gausspoint.h"
#include "fluiddynamicmaterial.h"
#include "loadtimefunction.h"

#include <fstream>

namespace oofem {
FluidMaterialEvaluator :: FluidMaterialEvaluator(int i, EngngModel *_master) : EngngModel(i, _master)
{
    this->ndomains = 1;
}

FluidMaterialEvaluator :: ~FluidMaterialEvaluator()
{
}

IRResultType FluidMaterialEvaluator :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
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
    int components = ( ndim * (ndim + 1) ) / 2;
    for ( int i = 1; i <= components; ++i ) {
        if (!sControl.contains(i))
            eControl.followedBy(i);
    }
    
    return IRRT_OK;
}

void FluidMaterialEvaluator :: solveYourself()
{
    Domain *d = this->giveDomain(1);
    
    gps.growTo(d->giveNumberOfMaterialModels());

    MaterialMode mode;
    if ( ndim == 1 ) {
        OOFEM_ERROR("FluidMaterialEvaluator :: solveYourself - 1d flow not supported (should be added)")
        //mode = _1dFlow;
        mode = _Unknown;
    } else if ( ndim == 2 ) {
        mode = _2dFlow;
    } else {
        mode = _3dFlow;
    }

    int components = ( ndim * (ndim + 1) ) / 2;
    FloatArray initialStrain(components);
    initialStrain.zero();
    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        gps.put(i, new GaussPoint(NULL, i, NULL, 1, mode));
        // Initialize the strain vector;        
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial* >(d->giveMaterial(i));
        FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus* >( mat->giveStatus(gps.at(i)) );
        status->letTempDeviatoricStrainRateVectorBe(initialStrain);
    }

    std::string outname = this->giveOutputBaseFileName() + ".matdata";
    this->outfile.open(outname.c_str());

    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    TimeStep *tStep = giveNextStep();

    // Note, strain == strain-rate (kept as strain for brevity)
    int maxiter = 100; // User input?
    double tolerance = 1.e-6; // Needs to be normalized somehow, or user input
    double strainVol, pressure;
    FloatArray stressDevC, deltaStrain, strainDev, stressDev, res;
    stressDevC.resize(sControl.giveSize());
    res.resize(sControl.giveSize());

    FloatMatrix tangent, reducedTangent;
    for ( int istep = 1; istep <= this->numberOfSteps; ++istep ) {
        this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
        for ( int imat = 1; imat <= d->giveNumberOfMaterialModels(); ++imat ) {
            GaussPoint *gp = gps.at(imat);
            FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial* >(d->giveMaterial(imat));
            FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus* >( mat->giveStatus(gp) );

            strainDev = status->giveDeviatoricStrainRateVector();
            // Update the controlled parts
            for ( int j = 1; j <= eControl.giveSize(); ++j ) {
                int p = eControl.at(j);
                strainDev.at(p) = d->giveLoadTimeFunction(cmpntFunctions.at(p))->evaluate(tStep, VM_Total);
            }
            for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                int p = sControl.at(j);
                stressDevC.at(p) = d->giveLoadTimeFunction(cmpntFunctions.at(p))->evaluate(tStep, VM_Total);
            }
            if ( pressureControl ) {
                pressure = d->giveLoadTimeFunction(volFunction)->evaluate(tStep, VM_Total);
            } else {
                ///@todo Support volumetric strain control (which is actually quite tricky)
                OOFEM_ERROR("Volumetric strain rate control not yet implemented");
                pressure = 0.;
            }

            for ( int iter = 1; iter < maxiter; iter++ ) {
                mat->computeDeviatoricStressVector(stressDev, strainVol, gp, strainDev, pressure, tStep);
                for ( int j = 1; j <= sControl.giveSize(); ++j ) {
                    res.at(j) = stressDevC.at(j) - stressDev.at(sControl.at(j));
                }
                OOFEM_LOG_RELEVANT("Time step: %d, Material %d, Iteration: %d,  Residual = %e\n", istep, imat, iter, res.computeNorm());
                if ( res.computeNorm() <= tolerance ) { ///@todo More flexible control of convergence needed. Something relative?
                    break;
                }
                
                mat->giveDeviatoricStiffnessMatrix(tangent, TangentStiffness, gp, tStep);
                // Add mean part to make it invertible
                double norm = tangent.computeFrobeniusNorm();
                for ( int i = 1; i <= this->ndim; ++i ) {
                    for ( int j = 1; j <= this->ndim; ++j ) {
                        tangent.at(i,j) += norm;
                    }
                }
                // Pick out the stress-controlled part;
                reducedTangent.beSubMatrixOf(tangent, sControl, sControl);

                // Update stress-controlled part of the strain
                reducedTangent.solveForRhs(res, deltaStrain);
                for ( int j = 1; j <= sControl.giveSize(); ++j) {
                    strainDev.at(sControl.at(j)) += deltaStrain.at(j);
                }
            }
            if ( res.computeNorm() > tolerance ) {
                OOFEM_WARNING("Residual did not converge!");
            }
            // This material model has converged, so we update it and go on to the next.
            mat->updateYourself(gp, tStep);
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
        if ( !dynamic_cast< FluidDynamicMaterial* >(d->giveMaterial(i)) ) {
            return 0;
        }
    }
    return EngngModel :: checkConsistency();
}

void FluidMaterialEvaluator :: doStepOutput(TimeStep *tStep)
{
    FloatArray outputValue;
    Domain *d = this->giveDomain(1);
    if (tStep->isTheFirstStep()) {
        this->outfile << "# Time";
        for ( int j = 1; j <= this->vars.giveSize(); ++j ) {
            this->outfile << ", " << __InternalStateTypeToString((InternalStateType)this->vars.at(j));
        }
        this->outfile << '\n';
    }

    for ( int i = 1; i <= d->giveNumberOfMaterialModels(); i++ ) {
        GaussPoint *gp = gps.at(i);
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial* >(d->giveMaterial(i));
        outfile << tStep->giveIntrinsicTime();
        for ( int j = 1; j <= this->vars.giveSize(); ++j ) {
            mat->giveIPValue(outputValue, gp, (InternalStateType)this->vars.at(j), tStep);
            outfile << " " << outputValue;
        }
    }
    outfile << std::endl;
}

TimeStep *FluidMaterialEvaluator :: giveNextStep()
{
    if ( previousStep ) {
        delete previousStep;
    }

    if ( currentStep == NULL ) {
        int istep = this->giveNumberOfFirstStep();
        // first step -> generate initial step
        previousStep = new TimeStep(giveNumberOfTimeStepWhenIcApply(), this, 0, -this->deltaT, this->deltaT, 0);
        currentStep = new TimeStep(istep, this, 1, 0.0, this->deltaT, 1);
    } else {
        int istep =  currentStep->giveNumber() + 1;
        StateCounterType counter = currentStep->giveSolutionStateCounter() + 1;
        previousStep = currentStep;
        double dt = currentStep->giveTimeIncrement();
        double totalTime = currentStep->giveTargetTime() + dt;
        currentStep = new TimeStep(istep, this, 1, totalTime, dt, counter);
    }

    return currentStep;
}

} // end namespace oofem

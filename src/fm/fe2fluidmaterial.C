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

#include "fe2fluidmaterial.h"
#include "stokesflow.h"
#include "oofemtxtdatareader.h"
#include "domain.h"
#include "gausspoint.h"
#include "contextioerr.h"
#include "util.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "mathfem.h"

#include <sstream>

//#define DEBUG_TANGENT
#define DEBUG_ERR ( 1e-6 )

namespace oofem {
REGISTER_Material(FE2FluidMaterial);

int FE2FluidMaterial :: n = 1;

void FE2FluidMaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber( tStep->giveNumber() );
    rveTStep->setTime( tStep->giveTargetTime() );
    rveTStep->setTimeIncrement( tStep->giveTimeIncrement() );
}

void FE2FluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    double r_vol;
    this->computeDeviatoricStressVector(answer, r_vol, gp, eps, 0.0, tStep);
    if ( gp->giveMaterialMode() == _2dFlow ) {
        r_vol += eps.at(1) + eps.at(2);
    } else {
        r_vol += eps.at(1) + eps.at(2) + eps.at(3);
    }

    if ( r_vol > 1e-9 ) {
        OOFEM_ERROR("RVE seems to be compressible;"
                    " extended macro-formulation which doesn't assume incompressibility is required");
    }
}

void FE2FluidMaterial :: computeDeviatoricStressVector(FloatArray &stress_dev, double &r_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast< FE2FluidMaterialStatus * >( this->giveStatus(gp) );

    ms->setTimeStep(tStep);

    MixedGradientPressureBC *bc = ms->giveBC();

    // Set input
    bc->setPrescribedDeviatoricGradientFromVoigt(eps);
    bc->setPrescribedPressure(pressure);
    // Solve subscale problem
    ms->giveRVE()->solveYourselfAt(ms->giveRVE()->giveCurrentStep());

    bc->computeFields(stress_dev, r_vol, tStep);

    ms->letDeviatoricStressVectorBe(stress_dev);
    ms->letDeviatoricStrainRateVectorBe(eps);
    ms->letPressureBe(pressure);
    ms->markOldTangents(); // Mark this so that tangents are reevaluated if they are needed.
    // One could also just compute them here, but you don't actually need them if the problem has converged, so this method saves on that iteration.
    // Computing the tangents are often *more* expensive than computeFields, so this is well worth the time it saves
    // All the tangents are computed in one go, because they are almost always all needed, and doing so saves time.
}

void FE2FluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast< FE2FluidMaterialStatus * >( this->giveStatus(gp) );
    ms->computeTangents(tStep);
    if ( mode == TangentStiffness ) {
        answer = ms->giveDeviatoricTangent();
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

void FE2FluidMaterial :: giveStiffnessMatrices(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast< FE2FluidMaterialStatus * >( this->giveStatus(gp) );
    ms->computeTangents(tStep);
    if ( mode == TangentStiffness ) {
        dsdd = ms->giveDeviatoricTangent();
        dsdp = ms->giveDeviatoricPressureTangent();
        dedd = ms->giveVolumetricDeviatoricTangent();
        dedp = ms->giveVolumetricPressureTangent();
#if 0
        // Numerical ATS for debugging
        FloatMatrix numericalATS(6, 6);
        FloatArray dsig;
        FloatArray tempStrain(6);

        tempStrain.zero();
        FloatArray sig, strain, sigPert;
        double epspvol;
        computeDeviatoricStressVector(sig, epspvol, gp, tempStrain, 0., tStep);
        double h = 0.001; // Linear problem, size of this doesn't matter.
        for ( int k = 1; k <= 6; ++k ) {
            strain = tempStrain;
            strain.at(k) += h;
            double tmp = strain.at(1) + strain.at(2) + strain.at(3);
            strain.at(1) -= tmp/3.0;
            strain.at(2) -= tmp/3.0;
            strain.at(3) -= tmp/3.0;
            strain.printYourself();
            computeDeviatoricStressVector(sigPert, epspvol, gp, strain, 0., tStep);
            sigPert.printYourself();
            dsig.beDifferenceOf(sigPert, sig);
            numericalATS.setColumn(dsig, k);
        }
        numericalATS.times(1. / h);

        printf("Analytical deviatoric tangent = ");
        dsdd.printYourself();
        printf("Numerical deviatoric tangent = ");
        numericalATS.printYourself();
        numericalATS.subtract(dsdd);
        double norm = numericalATS.computeFrobeniusNorm();
        if ( norm > dsdd.computeFrobeniusNorm() * DEBUG_ERR && norm > 0.0 ) {
            OOFEM_ERROR("Error in deviatoric tangent");
        }
#endif
#if 0
        // Numerical ATS for debugging
        FloatArray strain(3);
        strain.zero();
        FloatArray sig, sigh;
        double epspvol, pressure = 0.0;
        double h = 1.00; // Linear problem, size of this doesn't matter.
        computeDeviatoricStressVector(sig, epspvol, gp, strain, pressure, tStep);
        computeDeviatoricStressVector(sigh, epspvol, gp, strain, pressure + h, tStep);

        FloatArray dsigh;
        dsigh.beDifferenceOf(sigh, sig);
        dsigh.times(1 / h);

        printf("Analytical deviatoric pressure tangent = ");
        dsdp.printYourself();
        printf("Numerical deviatoric pressure tangent = ");
        dsigh.printYourself();
        dsigh.subtract(dsdp);
        double norm = dsigh.computeNorm();
        if ( norm > dsdp.computeNorm() * DEBUG_ERR && norm > 0.0 ) {
            OOFEM_ERROR("Error in deviatoric pressure tangent");
        }
#endif
#if 0
        // Numerical ATS for debugging
        FloatArray tempStrain(3);
        tempStrain.zero();
        FloatArray sig, strain;
        double epspvol, epspvol11, epspvol22, epspvol12, pressure = 0.0;
        double h = 1.0; // Linear problem, size of this doesn't matter.

        computeDeviatoricStressVector(sig, epspvol, gp, tempStrain, pressure, tStep);
        strain = tempStrain;
        strain.at(1) += h;
        computeDeviatoricStressVector(sig, epspvol11, gp, strain, pressure, tStep);
        strain = tempStrain;
        strain.at(2) += h;
        computeDeviatoricStressVector(sig, epspvol22, gp, strain, pressure, tStep);
        strain = tempStrain;
        strain.at(3) += h;
        computeDeviatoricStressVector(sig, epspvol12, gp, strain, pressure, tStep);

        FloatArray dvol(3);
        dvol.at(1) = ( epspvol11 - epspvol ) / h;
        dvol.at(2) = ( epspvol22 - epspvol ) / h;
        dvol.at(3) = ( epspvol12 - epspvol ) / h;
        dvol.at(1) += 1.0;
        dvol.at(2) += 1.0;

        printf("Analytical volumetric deviatoric tangent = ");
        dedd.printYourself();
        printf("Numerical volumetric deviatoric tangent = ");
        dvol.printYourself();
        dvol.subtract(dedd);
        double norm = dvol.computeNorm();
        if ( norm > dedd.computeNorm() * DEBUG_ERR && norm > 0.0 ) {
            OOFEM_ERROR("Error in volumetric deviatoric tangent");
        }
#endif
#if 0
        // Numerical ATS for debugging
        FloatArray strain(3);
        strain.zero();
        FloatArray sig;
        double epspvol, epspvolh, pressure = 0.0;
        double h = 1.0; // Linear problem, size of this doesn't matter.

        computeDeviatoricStressVector(sig, epspvol, gp, strain, pressure, tStep);
        computeDeviatoricStressVector(sig, epspvolh, gp, strain, pressure + h, tStep);

        double dvol = -( epspvolh - epspvol ) / h;

        printf("Analytical volumetric pressure tangent = %e\n", dedp);
        printf("Numerical volumetric pressure tangent = %e\n", dvol);

        double norm = fabs(dvol - dedp);
        if ( norm > fabs(dedp) * DEBUG_ERR && norm > 0.0 ) {
            OOFEM_ERROR("Error in volumetric pressure tangent");
        }
#endif
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

IRResultType FE2FluidMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    IR_GIVE_FIELD(ir, this->inputfile, _IFT_FE2FluidMaterial_fileName);
    return FluidDynamicMaterial :: initializeFrom(ir);
}

void FE2FluidMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    FluidDynamicMaterial :: giveInputRecord(input);
    input.setField(this->inputfile, _IFT_FE2FluidMaterial_fileName);
}


int FE2FluidMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FE2FluidMaterialStatus *status = static_cast< FE2FluidMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_VOFFraction ) {
        answer = FloatArray{status->giveVOFFraction()};
        return true;
    } else if ( type == IST_Pressure ) {
        answer = FloatArray{status->givePressure()};
        return true;
    } else if ( type == IST_Undefined ) { ///@todo What should one call this value? Relation between pressure and volumetric strain-rate.
#if 0
        // Numerical ATS for debugging
        FloatArray strain(3);
        strain.zero();
        FloatArray sig;
        double epspvol, epspvolh, pressure = 0.0;
        double h = 1.0; // Linear problem, size of this doesn't matter.

        computeDeviatoricStressVector(sig, epspvol, gp, strain, pressure, tStep);
        computeDeviatoricStressVector(sig, epspvolh, gp, strain, pressure + h, tStep);

        double dvol = - ( epspvolh - epspvol ) / h;

        
        printf("Analytical volumetric pressure tangent = %f\n", status->giveVolumetricPressureTangent());
        printf("Numerical volumetric pressure tangent = %f\n", dvol);

        double norm = fabs(dvol - status->giveVolumetricPressureTangent());
        if ( norm > fabs(status->giveVolumetricPressureTangent()) * DEBUG_ERR && norm > 0.0 ) {
            OOFEM_ERROR("Error in volumetric pressure tangent");
        }
#endif
        answer = FloatArray{status->giveVolumetricPressureTangent()};
        return true;
    } else {
        return FluidDynamicMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


MaterialStatus *FE2FluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new FE2FluidMaterialStatus(n++, this->giveDomain(), gp, this->inputfile);
}

int FE2FluidMaterial :: checkConsistency()
{
    return true;
}

FE2FluidMaterialStatus :: FE2FluidMaterialStatus(int n, Domain *d, GaussPoint *gp, const std :: string &inputfile) :
    FluidDynamicMaterialStatus(n, d, gp)
{
    //this->strainVector.resize(size);
    //this->strainVector.zero();
    //this->tempStrainVector = this->strainVector;
    this->voffraction = 0.0;
    this->oldTangents = true;

    if ( !this->createRVE(n, gp, inputfile) ) {
        OOFEM_ERROR("Couldn't create RVE");
    }
}

FE2FluidMaterialStatus :: ~FE2FluidMaterialStatus()
{
}

// Uses an input file for now, should eventually create the RVE itself.
bool FE2FluidMaterialStatus :: createRVE(int n, GaussPoint *gp, const std :: string &inputfile)
{
    OOFEMTXTDataReader dr( inputfile.c_str() );
    EngngModel *em = InstanciateProblem(& dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->setProblemScale(microScale);
    em->checkProblemConsistency();
    em->initMetaStepAttributes( em->giveMetaStep(1) );
    em->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    em->init();

    this->rve.reset( em );

    std :: ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    if ( this->domain->giveEngngModel()->isParallel() && this->domain->giveEngngModel()->giveNumberOfProcesses() > 1 ) {
        name << "." << this->domain->giveEngngModel()->giveRank();
    }

    this->rve->letOutputBaseFileNameBe( name.str() );

    this->bc = dynamic_cast< MixedGradientPressureBC * >( this->rve->giveDomain(1)->giveBc(1) );
    if ( !this->bc ) {
        OOFEM_ERROR("RVE doesn't have necessary boundary condition; should have MixedGradientPressure as first b.c. (in first domain)");
    }

    return true;
}

void FE2FluidMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: printOutputAt(file, tStep);
}

void FE2FluidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    double fluid_area = this->rve->giveDomain(1)->giveSize();
    double total_area = this->bc->domainSize();
    this->voffraction = fluid_area / total_area;
    FluidDynamicMaterialStatus :: updateYourself(tStep);

    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);
}

void FE2FluidMaterialStatus :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();
}

contextIOResultType FE2FluidMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = FluidDynamicMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return this->rve->saveContext(&stream, mode, obj);
}

contextIOResultType FE2FluidMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = FluidDynamicMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    this->markOldTangents();

    return this->rve->restoreContext(&stream, mode, obj);
}

void FE2FluidMaterialStatus :: markOldTangents() { this->oldTangents = true; }

void FE2FluidMaterialStatus :: computeTangents(TimeStep *tStep)
{
    if ( !tStep->isTheCurrentTimeStep() ) {
        OOFEM_ERROR("Only current timestep supported.");
    }

    if ( this->oldTangents ) {
        bc->computeTangents(this->giveDeviatoricTangent(),
                            this->giveDeviatoricPressureTangent(),
                            this->giveVolumetricDeviatoricTangent(),
                            this->giveVolumetricPressureTangent(),
                            tStep);
    }

    this->oldTangents = false;
}

double FE2FluidMaterial :: giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *status = static_cast< FE2FluidMaterialStatus * >( this->giveStatus(gp) );
    status->computeTangents(tStep);
    const FloatMatrix &t = status->giveDeviatoricTangent();
    // Project against the normalized I_dev
    double v;
    if ( gp->giveMaterialMode() == _3dFlow ) {
        v = ( t.at(1, 1) * 2. - t.at(1, 2)    - t.at(1, 3) ) / 3. +
            ( -t.at(2, 1)    + t.at(2, 2) * 2. - t.at(2, 3) ) / 3. +
            ( -t.at(3, 1)    - t.at(3, 2)    + t.at(3, 3) * 2. ) / 3. +
            4. * ( t.at(4, 4) + t.at(5, 5) + t.at(6, 6) );
        v /= 16.;
    } else {
        v = ( t.at(1, 1) * 2. - t.at(1, 2) ) / 3. +
            ( -t.at(2, 1)    + t.at(2, 2) * 2. ) / 3. +
            4. * t.at(3, 3);
        v *= 9. / 56.;
    }

    return v;
}
} // end namespace oofem

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

#include "fe2fluidmaterial.h"
#include "structuralmaterial.h"
#include "stokesflow.h"
#include "oofemtxtdatareader.h"
#include "domain.h"
#include "gausspnt.h"
#include "engngm.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "util.h"

// Used for computing
#include "line2boundaryelement.h"

#include <cstdlib>
#include <sstream>

namespace oofem {

int FE2FluidMaterial :: n = 1;

void FE2FluidMaterialStatus :: setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber(tStep->giveNumber());
    rveTStep->setTime(tStep->giveTargetTime());
    rveTStep->setTimeIncrement(tStep->giveTimeIncrement());
}

void FE2FluidMaterial :: computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps_dev, TimeStep *tStep)
{
    // This assumes that the RVE doesn't have any pores.
    double epsp_vol;
    this->computeDeviatoricStressVector(answer, epsp_vol, gp, eps_dev, 0.0, tStep);
    if ( epsp_vol > 1e-9 ) {
        OOFEM_ERROR("FE2FluidMaterialStatus :: computeDeviatoricStressVector - RVE seems to be compressible;"
            " extended macro-formulation which doesn't assume incompressibility is required");
    }
}

void FE2FluidMaterial :: computeDeviatoricStressVector(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps_dev, double pressure, TimeStep *tStep)
{
    FloatMatrix d, tangent;
    FE2FluidMaterialStatus *ms = static_cast<FE2FluidMaterialStatus*> (this->giveStatus(gp));

    ms->setTimeStep(tStep);

    MixedGradientPressureBC *bc = ms->giveBC();
    StokesFlow *rve = ms->giveRVE();

    // Set input
    bc->setPrescribedDeviatoricGradientFromVoigt(eps_dev);
    bc->setPrescribedPressure(pressure);
    // Solve subscale problem
    rve->solveYourselfAt(tStep);

    bc->computeFields(stress_dev, epsp_vol, EID_MomentumBalance_ConservationEquation, tStep);
    ms->letTempDeviatoricStressVectorBe(stress_dev);

    // Stores temporary tangents
    // Take the solver that the RVE-problem used for linear systems (even in non-linear materials, tangent problem is still linear)
    bc->computeTangents(ms->giveDeviatoricTangent(),
                        ms->giveDeviatoricPressureTangent(),
                        ms->giveVolumetricDeviatoricTangent(),
                        ms->giveVolumetricPressureTangent(),
                        EID_MomentumBalance_ConservationEquation, tStep);
}

void FE2FluidMaterial :: giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast<FE2FluidMaterialStatus*> (this->giveStatus(gp));
    if ( mode == TangentStiffness ) {
        answer = ms->giveDeviatoricTangent();
#if 0
        // Numerical ATS for debugging
        FloatArray tempStrain(3); tempStrain.zero();
        FloatArray sig, strain, sig11, sig22, sig12;
        double epspvol;
        computeDeviatoricStressVector (sig, epspvol, gp, tempStrain, 0.0, tStep);
        printf("strain = "); tempStrain.printYourself(); printf("stress = "); sig.printYourself();
        double h = 0.001; // Linear problem, size of this doesn't matter.
        strain.resize(3);
        strain = tempStrain; strain.at(1) += h*0.5; strain.at(2) -= h*0.5;
        computeDeviatoricStressVector (sig11, epspvol, gp, strain, 0.0, tStep);
        printf("strain = "); strain.printYourself(); printf("stress = "); sig11.printYourself();
        strain = tempStrain; strain.at(1) -= h*0.5; strain.at(2) += h*0.5;
        computeDeviatoricStressVector (sig22, epspvol, gp, strain, 0.0, tStep);
        strain = tempStrain; strain.at(3) += h;
        computeDeviatoricStressVector (sig12, epspvol, gp, strain, 0.0, tStep);
        printf("strain = "); strain.printYourself(); printf("stress = "); sig12.printYourself();

        FloatArray dsig11; dsig11.beDifferenceOf(sig11,sig); dsig11.times(1/h);
        FloatArray dsig22; dsig22.beDifferenceOf(sig22,sig); dsig22.times(1/h);
        FloatArray dsig12; dsig12.beDifferenceOf(sig12,sig); dsig12.times(1/h);

        FloatMatrix numericalATS;
        numericalATS.resize(3,3);
        numericalATS.zero();
        numericalATS.setColumn(dsig11,1);
        numericalATS.setColumn(dsig22,2);
        numericalATS.setColumn(dsig12,3);
        printf("Analytical tangent = "); answer.printYourself();
        printf("Numerical tangent = "); numericalATS.printYourself();
        OOFEM_ERROR("QUIT");
#endif
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

void FE2FluidMaterial :: giveDeviatoricPressureStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast<FE2FluidMaterialStatus*> (this->giveStatus(gp));
    if ( mode == TangentStiffness ) {
        answer = ms->giveDeviatoricPressureTangent();
#if 0
        // Numerical ATS for debugging
        FloatArray strain(3); strain.zero();
        FloatArray sig, sigh;
        double epspvol, pressure = 0.0;
        double h = 1.00; // Linear problem, size of this doesn't matter.
        computeDeviatoricStressVector (sig, epspvol, gp, strain, pressure, tStep);
        computeDeviatoricStressVector (sigh, epspvol, gp, strain, pressure+h, tStep);

        FloatArray dsigh; dsigh.beDifferenceOf(sigh,sig); dsigh.times(1/h);

        printf("Analytical tangent = "); answer.printYourself();
        printf("Numerical tangent = "); dsigh.printYourself();
        OOFEM_ERROR("QUIT");
#endif
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

void FE2FluidMaterial :: giveVolumetricDeviatoricStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast<FE2FluidMaterialStatus*> (this->giveStatus(gp));
    if ( mode == TangentStiffness ) {
        answer = ms->giveVolumetricDeviatoricTangent();
#if 0
        // Numerical ATS for debugging
        FloatArray tempStrain(3); tempStrain.zero();
        FloatArray sig, strain;
        double epspvol, epspvol11, epspvol22, epspvol12, pressure = 0.0;
        double h = 1.0; // Linear problem, size of this doesn't matter.

        computeDeviatoricStressVector (sig, epspvol, gp, tempStrain, pressure, tStep);
        strain = tempStrain; strain.at(1) += h*0.5; strain.at(2) -= h*0.5;
        computeDeviatoricStressVector(sig, epspvol11, gp, strain, pressure, tStep);
        strain = tempStrain; strain.at(2) -= h*0.5; strain.at(2) += h*0.5;
        computeDeviatoricStressVector(sig, epspvol22, gp, strain, pressure, tStep);
        strain = tempStrain; strain.at(3) += h;
        computeDeviatoricStressVector(sig, epspvol12, gp, strain, pressure, tStep);

        FloatArray dvol(3);
        dvol.at(1) = (epspvol11 - epspvol)/h;
        dvol.at(2) = (epspvol22 - epspvol)/h;
        dvol.at(3) = (epspvol12 - epspvol)/h;

        printf("Analytical tangent = "); answer.printYourself();
        printf("Numerical tangent = "); dvol.printYourself();
        OOFEM_ERROR("QUIT");
#endif
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

void FE2FluidMaterial :: giveVolumetricPressureStiffness(double &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FE2FluidMaterialStatus *ms = static_cast<FE2FluidMaterialStatus*> (this->giveStatus(gp));
    if ( mode == TangentStiffness ) {
        answer = ms->giveVolumetricPressureTangent();
#if 0
        // Numerical ATS for debugging
        FloatArray strain(3); strain.zero();
        FloatArray sig;
        double epspvol, epspvolh, pressure = 0.0;
        double h = 1.0; // Linear problem, size of this doesn't matter.

        computeDeviatoricStressVector(sig, epspvol, gp, strain, pressure, tStep);
        computeDeviatoricStressVector(sig, epspvolh, gp, strain, pressure+h, tStep);

        double dvol = (epspvolh - epspvol)/h;

        printf("Analytical tangent = %e\n", answer);
        printf("Numerical tangent = %e\n", dvol);
        OOFEM_ERROR("QUIT");
#endif
    } else {
        OOFEM_ERROR("Mode not implemented");
    }
}

int FE2FluidMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _PlaneStress ) { // _2dFlow
        return 1;
    }
    return 0;
}

IRResultType FE2FluidMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    IR_GIVE_FIELD(ir, this->inputfile, IFT_MicroMaterialFileName, "inputfile");
    return this->FluidDynamicMaterial::initializeFrom(ir);
}

int FE2FluidMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FE2FluidMaterialStatus *status = static_cast<FE2FluidMaterialStatus *>(this->giveStatus(gp));
    if (type == IST_VOFFraction) {
        answer.setValues(1, status->giveVOFFraction());
        return true;
    } else {
        return FluidDynamicMaterial::giveIPValue(answer, gp, type, tStep);
    }
}

InternalStateValueType FE2FluidMaterial :: giveIPValueType(InternalStateType type)
{
    if (type == IST_VOFFraction ) {
        return ISVT_SCALAR;
    } else {
        return FluidDynamicMaterial :: giveIPValueType(type);
    }
}

int FE2FluidMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if (type == IST_VOFFraction ) {
        return 1;
    } else {
        return FluidDynamicMaterial :: giveIPValueSize(type, gp);
    }
}

int FE2FluidMaterial :: giveInputRecordString(std::string &str, bool keyword)
{
    return FluidDynamicMaterial :: giveInputRecordString(str, keyword);
}

MaterialStatus * FE2FluidMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new FE2FluidMaterialStatus(n++, this->giveDomain(), gp, this->inputfile);
}

int FE2FluidMaterial :: checkConsistency()
{
    return true;
}

FE2FluidMaterialStatus :: FE2FluidMaterialStatus(int n, Domain *d, GaussPoint *gp, const std::string &inputfile) :
            FluidDynamicMaterialStatus(n, d, gp)
{
//     MaterialMode mmode = gp->giveMaterialMode();
//     int size = 0;
//
//     if (mmode == _2dFlow) {
//         size = 3;
//     } else {
//         OOFEM_ERROR("FE2FluidMaterialStatus: unsupported material mode");
//     }

    //this->strainVector.resize(size);
    //this->strainVector.zero();
    //this->tempStrainVector = this->strainVector;
    this->voffraction = 0.0;

    if (!this->createRVE(n, gp, inputfile)) {
        OOFEM_ERROR("FE2FluidMaterialStatus :: Constructor - Couldn't create RVE");
    }
}

FE2FluidMaterialStatus :: ~FE2FluidMaterialStatus()
{
    delete this->rve;
}

// Uses an input file for now, should eventually create the RVE itself.
bool FE2FluidMaterialStatus :: createRVE(int n, GaussPoint *gp, const std::string &inputfile)
{
    OOFEMTXTDataReader dr(inputfile.c_str());
    EngngModel *em = InstanciateProblem(&dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->setProblemScale(microScale);
    em->checkProblemConsistency();
    em->initMetaStepAttributes( em->giveMetaStep( 1 ) );
    em->giveNextStep(); // Makes sure there is a timestep (which we will modify before solving a step)
    em->init();

    this->rve = dynamic_cast<StokesFlow*> (em);
    if (!this->rve) {
        return false;
    }
    std::ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    this->rve->letOutputBaseFileNameBe(name.str());

    this->bc = dynamic_cast< MixedGradientPressureBC* >(this->rve->giveDomain(1)->giveBc(1));
    if (!this->bc) {
        OOFEM_ERROR("FE2FluidMaterialStatus :: createRVE - RVE doesn't have necessary boundary condition; should have MixedGradientPressure as first b.c. (in first domain)");
    }

    return true;
}

void FE2FluidMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    FluidDynamicMaterialStatus :: printOutputAt(file, tStep);
}

double FE2FluidMaterialStatus :: computeSize()
{
    Domain *d = this->rve->giveDomain(1);
    int nsd = d->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;

    // This requires the boundary to be consistent and ordered correctly.
    for (int i = 1; i <= d->giveNumberOfElements(); ++i) {
        //BoundaryElement *e = dynamic_cast< BoundaryElement* >(d->giveElement(i));
        Line2BoundaryElement *e = dynamic_cast< Line2BoundaryElement* >(d->giveElement(i));
        if (e) {
            domain_size += e->computeNXIntegral();
        }
    }
    return domain_size/nsd;
}

void FE2FluidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    double fluid_area = this->rve->giveDomain(1)->giveArea();
    double total_area = this->computeSize();
    this->voffraction = fluid_area/total_area;
    FluidDynamicMaterialStatus::updateYourself(tStep);

    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);
}

void FE2FluidMaterialStatus :: initTempStatus()
{
    FluidDynamicMaterialStatus :: initTempStatus();
}

contextIOResultType FE2FluidMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ((iores = FluidDynamicMaterialStatus :: saveContext(stream, mode, obj)) != CIO_OK) {
        THROW_CIOERR(iores);
    }
    return CIO_OK;
}

contextIOResultType FE2FluidMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ((iores = FluidDynamicMaterialStatus :: restoreContext(stream, mode, obj)) != CIO_OK) {
        THROW_CIOERR(iores);
    }
    return CIO_OK;
}

} // end namespace oofem

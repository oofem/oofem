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

#include "fe2sinteringmaterial.h"
#include "structuralmaterial.h"
#include "oofemtxtdatareader.h"
#include "domain.h"
#include "gausspnt.h"
#include "engngm.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "util.h"

#include <cstdlib>
#include <sstream>

namespace oofem {

int FE2SinteringMaterial::n = 1;

void FE2SinteringMaterialStatus::setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber(tStep->giveNumber());
    rveTStep->setTime(tStep->giveTargetTime());
    rveTStep->setTimeIncrement(tStep->giveTimeIncrement());
}

void FE2SinteringMaterial::giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
        const FloatArray &strain, TimeStep *tStep)
{
    FloatMatrix d, tangent;
    FloatArray strainRate;
    FE2SinteringMaterialStatus *ms = static_cast<FE2SinteringMaterialStatus*> (this->giveStatus(gp));

    ms->setTimeStep(tStep);

    // First we construct the prescribed tensor
    double dt = tStep->giveTimeIncrement();
    strainRate.beDifferenceOf(strain, ms->giveStrainVector());
    strainRate.times(1/dt);

    bool success = ms->giveRVE()->computeMacroStress(answer, strainRate, tStep);

    //ms->giveBC()->setPrescribedDeviatoricGradientFromVoigt(strainRate);
    //bool success = ms->rve()->solveYourselfAt(tStep);
    //ms->giveBC()->

    if (success) {
        ms->letTempStressVectorBe(answer);
        ms->letTempStrainVectorBe(strain);
    } else {
        OOFEM_ERROR("FE2SinteringMaterial::giveRealStressVector - Failed to compute stress");
    }
}

void FE2SinteringMaterial::givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode,
        GaussPoint *gp, TimeStep *tStep)
{
    FE2SinteringMaterialStatus *ms = static_cast<FE2SinteringMaterialStatus*> (this->giveStatus(gp));
    if (mode == TangentStiffness && form == ReducedForm) {
        StokesFlowStressHomogenization *rve = ms->giveRVE();

        // Since nonlinear static calls giveRealStressVector after this function, i have to repeat it.
        // This will change eventually. giveRealStressVector should always be called first for all materials.
        //ms->setTimeStep(tStep);
        //FloatArray temp, strainRate = ms->giveTempStrainVector();
        //strainRate.times(1 / tStep->giveTimeIncrement());
        //rve->computeMacroStress(temp, strainRate, tStep);

        // This is all I should need to do.
        rve->computeMacroTangent(answer, tStep);
        answer.times(1 / tStep->giveTimeIncrement());

#if 0
        // Numerical ATS for debugging
        FloatArray tempStrain = ms->giveTempStrainVector();
        FloatArray stress, strain, sig11, sig22, sig12, sig;
        giveRealStressVector (stress, form, gp, tempStrain, tStep);
        double h = 0.001; // Linear problem, size of this doesn't matter.
        strain.resize(3);
        strain = tempStrain; strain.at(1) += h;
        giveRealStressVector (sig11, ReducedForm, gp, strain, tStep);
        strain = tempStrain; strain.at(2) += h;
        giveRealStressVector (sig22, ReducedForm, gp, strain, tStep);
        strain = tempStrain; strain.at(3) += h;
        giveRealStressVector (sig12, ReducedForm, gp, strain, tStep);

        FloatArray dsig11 = (sig11 - stress)*(1/h);
        FloatArray dsig22 = (sig22 - stress)*(1/h);
        FloatArray dsig12 = (sig12 - stress)*(1/h);

        FloatMatrix numericalATS;
        numericalATS.resize(3,3);
        numericalATS.zero();
        numericalATS.addSubVectorRow(dsig11,1,1);
        numericalATS.addSubVectorRow(dsig22,2,1);
        numericalATS.addSubVectorRow(dsig12,3,1);
        answer.printYourself();
        numericalATS.printYourself();
        OOFEM_ERROR("QUIT");
#endif
    } else {
        OOFEM_ERROR("Mode or form not implemented");
    }
}

int FE2SinteringMaterial::hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _PlaneStress ) { // _2dFlow
        return 1;
    }
    return 0;
}

IRResultType FE2SinteringMaterial::initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    IR_GIVE_FIELD(ir, this->inputfile, IFT_MicroMaterialFileName, "inputfile");
    return this->StructuralMaterial::initializeFrom(ir);
}

int FE2SinteringMaterial::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FE2SinteringMaterialStatus *status = static_cast<FE2SinteringMaterialStatus *>(this->giveStatus(gp));
    if (type == IST_VOFFraction) {
        answer.resize(1);
        answer(0) = status->giveVOFFraction();
        return true;
    } else {
        return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
    }
}

int FE2SinteringMaterial::giveInputRecordString(std::string &str, bool keyword)
{
    return StructuralMaterial::giveInputRecordString(str, keyword);
}

MaterialStatus * FE2SinteringMaterial::CreateStatus(GaussPoint *gp) const
{
    return new FE2SinteringMaterialStatus(n++, this->giveDomain(), gp, this->inputfile);
}

int FE2SinteringMaterial::checkConsistency()
{
    return true;
}

FE2SinteringMaterialStatus::FE2SinteringMaterialStatus(int n, Domain *d, GaussPoint *gp, const std::string &inputfile) :
            StructuralMaterialStatus(n, d, gp)
{
    MaterialMode mmode = gp->giveMaterialMode();
    int size = 0;

    if (mmode == _PlaneStress) { // _2dFlow
        size = 3;
    } else {
        OOFEM_ERROR("FE2SinteringMaterialStatus: unsupported material mode");
    }

    this->strainVector.resize(size);
    this->strainVector.zero();
    this->tempStrainVector = this->strainVector;
    this->voffraction = 0.0;

    if (!this->createRVE(n, gp, inputfile)) {
        OOFEM_ERROR("Couldn't create RVE");
    }
}

// Uses an input file for now, should eventually create the RVE itself.
bool FE2SinteringMaterialStatus::createRVE(int n, GaussPoint *gp, const std::string &inputfile)
{
    OOFEMTXTDataReader dr(inputfile.c_str());
    EngngModel *em = InstanciateProblem(&dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->initMetaStepAttributes(em->giveNextStep());

    this->rve = dynamic_cast<StokesFlowStressHomogenization*> (em);
    if (!this->rve) {
        return false;
    }
    std::ostringstream name;
    name << this->rve->giveOutputBaseFileName() << "-gp" << n;
    this->rve->letOutputBaseFileNameBe(name.str());
    return true;
}

void FE2SinteringMaterialStatus::printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
}

void FE2SinteringMaterialStatus::updateYourself(TimeStep *tStep)
{
    // Must be done before sms updates;
    FloatArray strainIncrement = this->tempStrainVector;
    strainIncrement.subtract(this->strainVector);

    double fluid_area = this->rve->giveDomain(1)->giveArea();
    double total_area = this->rve->computeSize(1);
    this->voffraction = fluid_area/total_area;
    StructuralMaterialStatus::updateYourself(tStep);

    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep);
}

void FE2SinteringMaterialStatus::initTempStatus()
{
    StructuralMaterialStatus::initTempStatus();
}

contextIOResultType FE2SinteringMaterialStatus::saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ((iores = StructuralMaterialStatus::saveContext(stream, mode, obj)) != CIO_OK) {
        THROW_CIOERR(iores);
    }
    return CIO_OK;
}

contextIOResultType FE2SinteringMaterialStatus::restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ((iores = StructuralMaterialStatus::restoreContext(stream, mode, obj)) != CIO_OK) {
        THROW_CIOERR(iores);
    }
    return CIO_OK;
}

} // end namespace oofem

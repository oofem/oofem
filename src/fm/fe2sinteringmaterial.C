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

#include <stdlib.h>

namespace oofem {

int FE2SinteringMaterial::n = 1;

void FE2SinteringMaterialStatus::setTimeStep(TimeStep *tStep)
{
    TimeStep *rveTStep = this->rve->giveCurrentStep(); // Should i create a new one if it is empty?
    rveTStep->setNumber(tStep->giveNumber());
    rveTStep->setTime(tStep->giveTime());
    rveTStep->setTimeIncrement(tStep->giveTimeIncrement());
}

void FE2SinteringMaterial::giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
        const FloatArray &strain, TimeStep *tStep)
{
    FE2SinteringMaterialStatus *ms = static_cast<FE2SinteringMaterialStatus*> (this->giveStatus(gp));
    FloatMatrix d, tangent;

    ms->setTimeStep(tStep);

    // First we construct the prescribed tensor
    double dt = tStep->giveTimeIncrement();
    FloatArray oldStrain = ms->giveOldStrainVector();
    FloatArray strainRate = (strain - oldStrain) * (1 / dt); // Since we update the geometry, we need to increment only, and the rate.

    bool success = ms->giveRVE()->computeMacroStress(answer, strainRate, tStep);

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

        // Numerical ATS is just silly.
        /*
         printf("*************************************************************************** FOO 1\n");
         FloatArray tempStrain = ms->giveTempStrainVector();
         FloatArray stress, strain, sig11, sig22, sig12, sig;
         giveRealStressVector (stress, form, gp, tempStrain, tStep);
         double h = 0.001; // Linear problem, size of this doesn't matter.
         strain.resize(3);
         strain.zero(); strain.at(1) += h;
         giveRealStressVector (sig11, ReducedForm, gp, strain, tStep);
         strain.zero(); strain.at(2) += h;
         giveRealStressVector (sig22, ReducedForm, gp, strain, tStep);
         strain.zero(); strain.at(3) += h;
         giveRealStressVector (sig12, ReducedForm, gp, strain, tStep);

         FloatArray dsig11 = (sig11 - stress)*(1/h);
         FloatArray dsig22 = (sig22 - stress)*(1/h);
         FloatArray dsig12 = (sig12 - stress)*(1/h);

         FloatMatrix numericalATS;
         numericalATS.resize(3,3);
         numericalATS.addSubVectorRow(dsig11,1,1);
         numericalATS.addSubVectorRow(dsig22,2,1);
         numericalATS.addSubVectorRow(dsig12,3,1);
         //tempStrain.printYourself();
         //stress.printYourself();
         answer.printYourself();
         numericalATS.printYourself();
         OOFEM_ERROR("QUIT");
         */
    } else {
        OOFEM_ERROR("Mode or form not implemented");
    }
}

int FE2SinteringMaterial::hasMaterialModeCapability(MaterialMode mode)
{
    if ((mode == _PlaneStress)) { // _2dFlow
        return 1;
    }
    return 0;
}

IRResultType FE2SinteringMaterial::initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = this->StructuralMaterial::initializeFrom(ir);
    return result;
}

int FE2SinteringMaterial::giveInputRecordString(std::string &str, bool keyword)
{
    char buff[1024];
    StructuralMaterial::giveInputRecordString(str, keyword);
    sprintf(buff, " mu %e gamma_s %e", this->mu, this->gamma_s);
    str += buff;
    return 1;
}

MaterialStatus * FE2SinteringMaterial::CreateStatus(GaussPoint *gp) const
{
    return new FE2SinteringMaterialStatus(n++, this->giveDomain(), gp);
}

int FE2SinteringMaterial::checkConsistency()
{
    return 1;
}

FE2SinteringMaterialStatus::FE2SinteringMaterialStatus(int n, Domain *d, GaussPoint *gp) :
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
    this->oldStrainVector = this->tempStrainVector = this->strainVector;
    this->volume = 2.0;

    if (!this->createRVE(n, gp)) {
        OOFEM_ERROR("Couldn't create RVE");
    }
}

// Uses an input file for now, should eventually create the RVE itself.
bool FE2SinteringMaterialStatus::createRVE(int n, GaussPoint *gp)
{
    //FloatArray *gc = gp->giveCoordinates();
    //double cx = gc->at(1);
    //double cy = gc->at(2);
    double density = 0.84;

    char filename[1024];
    //sprintf(filename,"/home/mikael/workspace/OOFEM/4particle_rve.in");
    //sprintf(filename,"/home/mikael/workspace/OOFEM/4particle_super_coarse_mesh.in");
    //sprintf(filename,"/home/mikael/workspace/OOFEM/4particle_mega_coarse_mesh.in");
    //sprintf(filename,"/home/mikael/workspace/OOFEM/4particle_coarse_mesh.in");
    //sprintf(filename,"/home/mikael/workspace/OOFEM/4particle_0.88_mesh.in");
    sprintf(filename, "/home/mikael/workspace/OOFEM/4particle_%.2f_mesh.in", density);

    OOFEMTXTDataReader dr(filename);
    EngngModel *em = InstanciateProblem(&dr, _processor, 0); // Everything but nrsolver is updated.
    dr.finish();
    em->initMetaStepAttributes(em->giveNextStep());


    this->rve = dynamic_cast<StokesFlowStressHomogenization*> (em);
    if (!this->rve) {
        printf("em = %p\n",em);
        OOFEM_ERROR("WTF");
        return false;
    }

    char tempstring[1024];
    char newName[1024];
    this->rve->giveOutputBaseFileName(tempstring, 1024);
    sprintf(newName, "%s-%d", tempstring, n);
    this->rve->letOutputBaseFileNameBe(newName);

    this->rve->activateHomogenizationMode(HT_StressDirichlet);

    return true;
}

void FE2SinteringMaterialStatus::printOutputAt(FILE *file, TimeStep *tNow)
// Prints the strains and stresses on the data file.
// Not sure i need this. I already get alot of data from every gausspoint. But the homogenized values might be of interest too.
{
}

void FE2SinteringMaterialStatus::updateYourself(TimeStep *tStep)
{
    this->oldStrainVector = this->strainVector;
    StructuralMaterialStatus::updateYourself(tStep);
    this->volume += strainVector.at(1) + strainVector.at(2) + strainVector.at(3); // must be after strain is updated.
    this->rve->updateYourself(tStep);
    this->rve->terminate(tStep); // I think updateYourself(tStep) should be part of this.
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

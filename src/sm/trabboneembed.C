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

#include "trabboneembed.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "contextioerr.h"

namespace oofem {

TrabBoneEmbed :: TrabBoneEmbed(int n, Domain *d) : StructuralMaterial(n, d)
{}


int
TrabBoneEmbed :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( mode == _3dMat ) {
        return 1;
    }

    return 0;
}


void TrabBoneEmbed :: computeCumPlastStrain(double &tempAlpha, GaussPoint *gp, TimeStep *atTime)
{
    tempAlpha = 0.;
}


void
TrabBoneEmbed :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *atTime)
{
    TrabBoneEmbedStatus *status = ( TrabBoneEmbedStatus * ) this->giveStatus(gp);

    FloatMatrix elasticity, compliance;

    this->constructIsoComplTensor(compliance, eps0, nu0);
    elasticity.beInverseOf(compliance);

    answer.resize(6, 6);
    answer = elasticity;

    status->setSmtrx(answer);
}


void
TrabBoneEmbed :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain)
{
    double tempAlpha;
    FloatArray tempPlasDef;

    TrabBoneEmbedStatus *status = ( TrabBoneEmbedStatus * ) this->giveStatus(gp);

    tempPlasDef.resize(6);
    tempAlpha = 0.;

    status->setTempPlasDef(tempPlasDef);
    status->setTempAlpha(tempAlpha);
}


double
TrabBoneEmbed :: computeDamageParam(double alpha, GaussPoint *gp)
{
    double tempDam = 0.0;

    return tempDam;
}


double
TrabBoneEmbed :: computeDamage(GaussPoint *gp,  TimeStep *atTime)
{
    double tempAlpha;

    computeCumPlastStrain(tempAlpha, gp, atTime);

    double tempDam = computeDamageParam(tempAlpha, gp);

    //  double dam=0.0;

    return tempDam;
}


void
TrabBoneEmbed :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &totalStrain,
                                      TimeStep *atTime)
{
    double tempDam, tempTSED;
    FloatArray newTotalDef, plasDef;
    FloatArray totalStress;
    double dt;
    FloatMatrix compliance, elasticity;

    this->constructIsoComplTensor(compliance, eps0, nu0);
    elasticity.beInverseOf(compliance);

    TrabBoneEmbedStatus *status = ( TrabBoneEmbedStatus * ) this->giveStatus(gp);

    this->initGpForNewStep(gp);

    performPlasticityReturn(gp, totalStrain);

    dt = atTime->giveTimeIncrement();

    tempDam = computeDamage(gp, atTime);

    plasDef.resize(6);

    totalStress.beProductOf(elasticity, totalStrain);

    tempTSED = 0.5 * totalStrain.dotProduct(totalStress);

    answer.resize(6);
    answer = totalStress;
    status->setTempDam(tempDam);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
    status->setTempTSED(tempTSED);
}


void
TrabBoneEmbed :: constructIsoComplTensor(FloatMatrix &answer, const double eps0, const double nu0)
{
    double mu0 = eps0 / ( 2 * ( 1 + nu0 ) );

    answer.resize(6, 6);
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 1 / eps0;
    answer.at(1, 2) = answer.at(2, 1) = answer.at(1, 3) = -nu0 / eps0;
    answer.at(3, 1) = answer.at(2, 3) = answer.at(3, 2) = -nu0 / eps0;
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = 1 / mu0;
}


IRResultType
TrabBoneEmbed :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Read material properties here

    IR_GIVE_FIELD(ir, eps0, IFT_TrabBoneEmbed_eps0, "eps0"); // Macro
    IR_GIVE_FIELD(ir, nu0, IFT_TrabBoneEmbed_nu0, "nu0"); // Macro

    return StructuralMaterial :: initializeFrom(ir);
}


int
TrabBoneEmbed :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    TrabBoneEmbedStatus *status = ( TrabBoneEmbedStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    } else if ( type == IST_BoneVolumeFraction ) {
        answer.resize(1);
        answer.at(1) = 1.;
        return 1;
    } else if ( type == IST_PlasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    } else if ( type == IST_ElasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = status->giveTempTSED();
        return 1;
    } else if ( type == IST_TotalStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = status->giveTempTSED();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
TrabBoneEmbed :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_BoneVolumeFraction ) || ( type == IST_PlasStrainEnerDens ) || ( type == IST_ElasStrainEnerDens ) || ( type == IST_TotalStrainEnerDens ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
TrabBoneEmbed :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_BoneVolumeFraction ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_PlasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_ElasStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( type == IST_TotalStrainEnerDens ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}

int
TrabBoneEmbed :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_DamageScalar ) || ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_BoneVolumeFraction ) || ( type == IST_ElasStrainEnerDens ) || ( type == IST_PlasStrainEnerDens ) || ( type == IST_TotalStrainEnerDens ) ) {
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        return 6;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}

/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////

TrabBoneEmbedStatus :: TrabBoneEmbedStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    alpha = 0.0;
    dam = 0.0;
    tsed = 0.0;
    tempAlpha = 0.0;
    tempDam = 0.0;
    tempTSED = 0.0;
    smtrx.resize(6, 6);
}


TrabBoneEmbedStatus :: ~TrabBoneEmbedStatus()
{}


double
TrabBoneEmbedStatus :: giveTempTSED()
{
    return tempTSED;
}


void
TrabBoneEmbedStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    fprintf( file, "plastrains: %f  %f  %f  %f  %f  %f", this->tempPlasDef.at(1), this->tempPlasDef.at(2), this->tempPlasDef.at(3), this->tempPlasDef.at(4), this->tempPlasDef.at(5), this->tempPlasDef.at(6) );
    fprintf(file, " , alpha 0. , dam 0. , esed %f , psed 0. , tsed %f ", this->tempTSED, this->tempTSED);
    fprintf(file, "}\n");
}


void
TrabBoneEmbedStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempAlpha = this->alpha;
    this->tempDam = this->dam;
    this->tempTSED = this->tsed;
    this->tempPlasDef = this->plasDef;
}


void
TrabBoneEmbedStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->alpha = this->tempAlpha;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
}


contextIOResultType
TrabBoneEmbedStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    //if (fwrite(&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fwrite(&damage,sizeof(double),1,stream)!= 1) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


contextIOResultType
TrabBoneEmbedStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    //if (fread (&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fread (&damage,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

    return CIO_OK;
}


MaterialStatus *TrabBoneEmbed :: CreateStatus(GaussPoint *gp) const
{
    TrabBoneEmbedStatus *status =
        new  TrabBoneEmbedStatus(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}

} // end namespace oofem

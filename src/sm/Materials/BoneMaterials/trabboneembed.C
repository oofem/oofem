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

#include "trabboneembed.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(TrabBoneEmbed);

TrabBoneEmbed :: TrabBoneEmbed(int n, Domain *d) : StructuralMaterial(n, d)
{ }


double TrabBoneEmbed :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    return 0.;
}


FloatMatrixF<6,6>
TrabBoneEmbed :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    // 'auto status = static_cast< TrabBoneEmbedStatus * >( this->giveStatus(gp) );

    auto compliance = this->constructIsoComplTensor(eps0, nu0);
    auto elasticity = inv(compliance);

    return elasticity;
}


void
TrabBoneEmbed :: performPlasticityReturn(GaussPoint *gp, const FloatArrayF<6> &totalStrain) const
{
    auto status = static_cast< TrabBoneEmbedStatus * >( this->giveStatus(gp) );

    status->setTempPlasDef(zeros<6>());
    status->setTempAlpha(0.);
}


double
TrabBoneEmbed :: computeDamageParam(double alpha, GaussPoint *gp) const
{
    double tempDam = 0.0;

    return tempDam;
}


double
TrabBoneEmbed :: computeDamage(GaussPoint *gp,  TimeStep *tStep) const
{
    double tempAlpha = computeCumPlastStrain(gp, tStep);
    double tempDam = computeDamageParam(tempAlpha, gp);

    //  double dam=0.0;

    return tempDam;
}


FloatArrayF<6>
TrabBoneEmbed :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                         TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneEmbedStatus * >( this->giveStatus(gp) );

    auto compliance = this->constructIsoComplTensor(eps0, nu0);
    auto elasticity = inv(compliance);

    this->initTempStatus(gp);

    performPlasticityReturn(gp, strain);

    double tempDam = computeDamage(gp, tStep);
    //FloatArrayF<6> plasDef;

    auto stress = dot(elasticity, strain);

    double tempTSED = 0.5 * dot(strain, stress);

    status->setTempDam(tempDam);
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(stress);
    status->setTempTSED(tempTSED);
    return stress;
}


FloatMatrixF<6,6>
TrabBoneEmbed :: constructIsoComplTensor(double eps0, double nu0)
{
    double mu0 = eps0 / ( 2 * ( 1 + nu0 ) );

    FloatMatrixF<6,6> c;
    c.at(1, 1) = c.at(2, 2) = c.at(3, 3) = 1 / eps0;
    c.at(1, 2) = c.at(2, 1) = c.at(1, 3) = -nu0 / eps0;
    c.at(3, 1) = c.at(2, 3) = c.at(3, 2) = -nu0 / eps0;
    c.at(4, 4) = c.at(5, 5) = c.at(6, 6) = 1 / mu0;
    return c;
}


void
TrabBoneEmbed :: initializeFrom(InputRecord &ir)
{
    StructuralMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, eps0, _IFT_TrabBoneEmbed_eps0);
    IR_GIVE_FIELD(ir, nu0, _IFT_TrabBoneEmbed_nu0);
}


int
TrabBoneEmbed :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto status = static_cast< TrabBoneEmbedStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = 0.;
        OOFEM_WARNING("No damage is exported (why?!)");
        return 1;
    } else if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasDef();
        OOFEM_WARNING("Unsure what components are stored in the plastic strain tensor");
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
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}


/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////

TrabBoneEmbedStatus :: TrabBoneEmbedStatus(GaussPoint *g) : StructuralMaterialStatus(g)
{
}


void
TrabBoneEmbedStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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
TrabBoneEmbedStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->alpha = this->tempAlpha;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
}


void
TrabBoneEmbedStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    // save parent class status
    StructuralMaterialStatus :: saveContext(stream, mode);

    // write a raw data
    //if (fwrite(&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fwrite(&damage,sizeof(double),1,stream)!= 1) THROW_CIOERR(CIO_IOERR);
}


void
TrabBoneEmbedStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    // read parent class status
    StructuralMaterialStatus :: restoreContext(stream, mode);

    // read raw data
    //if (fread (&kappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
    //if (fread (&damage,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
}


std::unique_ptr<MaterialStatus> TrabBoneEmbed :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<TrabBoneEmbedStatus>(gp);
}
} // end namespace oofem

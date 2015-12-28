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

#include "trabbonenlembed.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "nonlocalmaterialext.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////

REGISTER_Material(TrabBoneNLEmbed);

TrabBoneNLEmbed :: TrabBoneNLEmbed(int n, Domain *d) : TrabBoneEmbed(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{
    R = 0.;
}

TrabBoneNLEmbed :: ~TrabBoneNLEmbed()
{ }

void
TrabBoneNLEmbed :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray SDstrainVector;
    double cumPlastStrain;
    TrabBoneNLEmbedStatus *nlstatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    nlstatus->letTempStrainVectorBe(strainVector);

    this->performPlasticityReturn(gp, strainVector);
    this->computeLocalCumPlastStrain(cumPlastStrain, strainVector, gp, tStep);

    nlstatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}

void
TrabBoneNLEmbed :: giveRealStressVector_3d(FloatArray &answer,
                                           GaussPoint *gp,
                                           const FloatArray &strainVector,
                                           TimeStep *tStep)
{
    TrabBoneNLEmbedStatus *nlStatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );

    double tempDam, tempTSED;
    FloatArray plasDef, totalStress;
    FloatMatrix compliance, elasticity;

    compliance.resize(6, 6);
    this->constructIsoComplTensor(compliance, eps0, nu0);
    elasticity.beInverseOf(compliance);

    tempDam = 0;
    plasDef.resize(6);

    totalStress.beProductOf(elasticity, strainVector);

    tempTSED = 0.5 * strainVector.dotProduct(totalStress);

    answer.resize(6);
    answer = totalStress;

    nlStatus->setTempDam(tempDam);
    nlStatus->letTempStrainVectorBe(strainVector);
    nlStatus->letTempStressVectorBe(answer);
    nlStatus->setTempTSED(tempTSED);
}

void
TrabBoneNLEmbed :: computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep)
{
    double nonlocalContribution, nonlocalCumPlastStrain = 0.0;
    TrabBoneNLEmbedStatus *nonlocStatus, *status = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    std :: list< localIntegrationRecord > *list = status->giveIntegrationDomainList();

    for ( auto &lir: *list ) {
        nonlocStatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(lir.nearGp) );
        nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= lir.weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain /= status->giveIntegrationScale(); ///@todo This is never used.

    //  double localCumPlastStrain = status->giveLocalCumPlastStrainForAverage();
    //  alpha = mParam*nonlocalCumPlastStrain +(1-mParam)*localCumPlastStrain ;
    alpha = 0.;
}

Interface *
TrabBoneNLEmbed :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

IRResultType
TrabBoneNLEmbed :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    result = TrabBoneEmbed :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, R, _IFT_TrabBoneNLEmbed_r);
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 1.5;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneNLEmbed_m);

    return IRRT_OK;
}


void
TrabBoneNLEmbed :: giveInputRecord(DynamicInputRecord &input)
{
    TrabBoneEmbed :: giveInputRecord(input);
    input.setField(this->R, _IFT_TrabBoneNLEmbed_r);
    input.setField(this->mParam, _IFT_TrabBoneNLEmbed_m);
}


double
TrabBoneNLEmbed :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    double dist = src.distance(coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}

TrabBoneNLEmbedStatus :: TrabBoneNLEmbedStatus(int n, Domain *d, GaussPoint *g) :
    TrabBoneEmbedStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlastStrainForAverage = 0.0;
}

TrabBoneNLEmbedStatus :: ~TrabBoneNLEmbedStatus()
{ }

void
TrabBoneNLEmbedStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    for ( auto &val : plasDef ) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, " ,");
    fprintf(file, " alpha %.4e ,", alpha);
    fprintf(file, " dam  %.4e ,", dam);
    fprintf(file, " esed  %.4e ,", this->tempTSED);
    fprintf(file, " psed  0. ,");
    fprintf(file, " tsed  %.4e ,", this->tempTSED);
    fprintf(file, "}\n");
}

void
TrabBoneNLEmbedStatus :: initTempStatus()
{
    TrabBoneEmbedStatus :: initTempStatus();
}

void
TrabBoneNLEmbedStatus :: updateYourself(TimeStep *tStep)
{
    TrabBoneEmbedStatus :: updateYourself(tStep);
}

contextIOResultType
TrabBoneNLEmbedStatus :: saveContext(DataStream &stream, ContextMode mode,  void *obj)
{
    return TrabBoneEmbedStatus :: saveContext(stream, mode, obj);
}

contextIOResultType
TrabBoneNLEmbedStatus :: restoreContext(DataStream &stream, ContextMode mode,  void *obj)
{
    return TrabBoneEmbedStatus :: restoreContext(stream, mode, obj);
}


Interface *
TrabBoneNLEmbedStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return NULL;
    }
}
} // end namespace oofem

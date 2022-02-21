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
}

void
TrabBoneNLEmbed :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const
{
    auto nlstatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    FloatArray SDstrainVector;
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    nlstatus->letTempStrainVectorBe(strainVector);

    this->performPlasticityReturn(gp, strainVector);
    double cumPlastStrain = this->computeLocalCumPlastStrain(strainVector, gp, tStep);

    nlstatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}

FloatArrayF<6>
TrabBoneNLEmbed :: giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                           TimeStep *tStep) const
{
    auto nlStatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );
    
    auto compliance = this->constructIsoComplTensor(eps0, nu0);
    auto elasticity = inv(compliance);

    // Unused?
    double tempDam = 0;
    FloatArrayF<6> plasDef;

    auto totalStress = dot(elasticity, strain);

    double tempTSED = 0.5 * dot(strain, totalStress);

    nlStatus->setTempDam(tempDam);
    nlStatus->letTempStrainVectorBe(strain);
    nlStatus->letTempStressVectorBe(totalStress);
    nlStatus->setTempTSED(tempTSED);
    return totalStress;
}

double
TrabBoneNLEmbed :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    auto list = status->giveIntegrationDomainList();

    double nonlocalCumPlastStrain = 0.0;
    for ( auto &lir: *list ) {
        auto nonlocStatus = static_cast< TrabBoneNLEmbedStatus * >( this->giveStatus(lir.nearGp) );
        double nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= lir.weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain /= status->giveIntegrationScale(); ///@todo This is never used.

    //  double localCumPlastStrain = status->giveLocalCumPlastStrainForAverage();
    //  alpha = mParam*nonlocalCumPlastStrain +(1-mParam)*localCumPlastStrain ;
    return 0.;
}

Interface *
TrabBoneNLEmbed :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return nullptr;
    }
}

void
TrabBoneNLEmbed :: initializeFrom(InputRecord &ir)
{
    TrabBoneEmbed :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, _IFT_TrabBoneNLEmbed_r);
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 1.5;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneNLEmbed_m);
}


void
TrabBoneNLEmbed :: giveInputRecord(DynamicInputRecord &input)
{
    TrabBoneEmbed :: giveInputRecord(input);
    input.setField(this->R, _IFT_TrabBoneNLEmbed_r);
    input.setField(this->mParam, _IFT_TrabBoneNLEmbed_m);
}


double
TrabBoneNLEmbed :: computeWeightFunction(const double R, const FloatArray &src, const FloatArray &coord) const
{
    double dist = distance(src, coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}


TrabBoneNLEmbedStatus :: TrabBoneNLEmbedStatus(GaussPoint *g) :
    TrabBoneEmbedStatus(g), StructuralNonlocalMaterialStatusExtensionInterface()
{
}


void
TrabBoneNLEmbedStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
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

void
TrabBoneNLEmbedStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    TrabBoneEmbedStatus :: saveContext(stream, mode);
}

void
TrabBoneNLEmbedStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    TrabBoneEmbedStatus :: restoreContext(stream, mode);
}


Interface *
TrabBoneNLEmbedStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return nullptr;
    }
}
} // end namespace oofem

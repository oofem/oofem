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

#include "trabbonenlembed.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "dynalist.h"
#include "nonlocalmaterialext.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////

TrabBoneNLEmbed :: TrabBoneNLEmbed(int n, Domain *d) : TrabBoneEmbed(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{
    R = 0.;
}

TrabBoneNLEmbed :: ~TrabBoneNLEmbed()
{}

void
TrabBoneNLEmbed :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    FloatArray SDstrainVector, fullSDStrainVector;
    double cumPlastStrain;
    TrabBoneNLEmbedStatus *nlstatus = ( TrabBoneNLEmbedStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, atTime, VM_Total);

    nlstatus->letTempStrainVectorBe(strainVector);

    StrainVector strain( strainVector, gp->giveMaterialMode() );

    this->performPlasticityReturn(gp, strain);
    this->computeLocalCumPlastStrain(cumPlastStrain, strain, gp, atTime);

    nlstatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}

void
TrabBoneNLEmbed :: giveRealStressVector(FloatArray &answer,
                                        MatResponseForm form,
                                        GaussPoint *gp,
                                        const FloatArray &strainVector,
                                        TimeStep *atTime)
{
    TrabBoneNLEmbedStatus *nlStatus = ( TrabBoneNLEmbedStatus * ) this->giveStatus(gp);

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
TrabBoneNLEmbed :: computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *atTime)
{
    double nonlocalContribution, nonlocalCumPlastStrain = 0.0;
    TrabBoneNLEmbedStatus *nonlocStatus, *status = ( TrabBoneNLEmbedStatus * ) this->giveStatus(gp);

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);

    dynaList< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        nonlocStatus = ( TrabBoneNLEmbedStatus * ) this->giveStatus( ( * pos ).nearGp );
        nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= ( * pos ).weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain *= 1. / status->giveIntegrationScale();

    //  double localCumPlastStrain = status->giveLocalCumPlastStrainForAverage();
    //  alpha = mParam*nonlocalCumPlastStrain +(1-mParam)*localCumPlastStrain ;
    alpha = 0.;
}

Interface *
TrabBoneNLEmbed :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return ( StructuralNonlocalMaterialExtensionInterface * ) this;
    } else if ( type == NonlocalMaterialStiffnessInterfaceType ) {
        return ( NonlocalMaterialStiffnessInterface * ) this;
    }
    //
    else {
        return NULL;
    }
}

IRResultType
TrabBoneNLEmbed :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    TrabBoneEmbed :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, IFT_TrabBoneNLEmbed_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 1.5;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_TrabBoneNLEmbed_m, "m"); // Macro

    return IRRT_OK;
}

int
TrabBoneNLEmbed :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    TrabBoneEmbed :: giveInputRecordString(str, keyword);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);
    sprintf(buff, " r %e", this->R);
    str += buff;

    return 1;
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
{}

void
TrabBoneNLEmbedStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    int n = plasDef.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasDef.at(i) );
    }

    fprintf(file, " ,");
    fprintf(file, " alpha % .4e ,", alpha);
    fprintf(file, " dam  % .4e ,", dam);
    fprintf(file, " esed  % .4e ,", this->tempTSED);
    fprintf(file, " psed  0. ,");
    fprintf(file, " tsed  % .4e ,", this->tempTSED);
    fprintf(file, "}\n");
}

void
TrabBoneNLEmbedStatus :: initTempStatus()
{
    TrabBoneEmbedStatus :: initTempStatus();
}

void
TrabBoneNLEmbedStatus :: updateYourself(TimeStep *atTime)
{
    TrabBoneEmbedStatus :: updateYourself(atTime);
}

contextIOResultType
TrabBoneNLEmbedStatus :: saveContext(DataStream *stream, ContextMode mode,  void *obj)
{
    return TrabBoneEmbedStatus :: saveContext(stream, mode, obj);
}

contextIOResultType
TrabBoneNLEmbedStatus :: restoreContext(DataStream *stream, ContextMode mode,  void *obj)
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

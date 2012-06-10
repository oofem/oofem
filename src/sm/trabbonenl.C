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

#include "trabbonenl.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "mathfem.h"
#include "dynalist.h"
#include "nonlocalmaterialext.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
/////////////////////////////////////////////////////////////////
////////////TRABECULAR BONE NONLOCAL MATERIAL////////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
//

TrabBoneNL :: TrabBoneNL(int n, Domain *d) : TrabBoneMaterial(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{
    R = 0.;
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneNL :: ~TrabBoneNL()
{}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
// update local values of accumulated pastic strain

void
TrabBoneNL :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime)
{
    FloatArray SDstrainVector, fullSDStrainVector;
    double cumPlastStrain;
    TrabBoneNLStatus *nlstatus = ( TrabBoneNLStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, atTime, VM_Total);

    nlstatus->letTempStrainVectorBe(strainVector);

    StrainVector strain( strainVector, gp->giveMaterialMode() );

    this->performPlasticityReturn(gp, strain);
    this->computeLocalCumPlastStrain(cumPlastStrain, strain, gp, atTime);
    nlstatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}

//
// END: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
//

void
TrabBoneNL :: giveRealStressVector(FloatArray &answer,
                                   MatResponseForm form,
                                   GaussPoint *gp,
                                   const FloatArray &strainVector,
                                   TimeStep *atTime)
{
    TrabBoneNLStatus *nlStatus = ( TrabBoneNLStatus * ) this->giveStatus(gp);

    double tempDamage = computeDamage(gp, atTime);
    double plasStrain = nlStatus->giveTempPlasStrainVector().at(1);
    double elasStrain = strainVector.at(1) - plasStrain;
    double effStress = E0 * elasStrain;
    double sigc = nlStatus->giveSigC();

    answer.resize(1);
    answer.at(1) = ( 1. - tempDamage ) * effStress + sigc;

    nlStatus->setTempDam(tempDamage);
    nlStatus->letTempStrainVectorBe(strainVector);
    nlStatus->letTempStressVectorBe(answer);
}

//
// END: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF NONLOCAL ALPHA EVALUATION
//

void
TrabBoneNL :: computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *atTime)
{
    double nonlocalContribution, nonlocalCumPlastStrain = 0.0;
    TrabBoneNLStatus *nonlocStatus, *status = ( TrabBoneNLStatus * ) this->giveStatus(gp);

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(atTime);

    dynaList< localIntegrationRecord > *list = status->giveIntegrationDomainList();
    dynaList< localIntegrationRecord > :: iterator pos;

    for ( pos = list->begin(); pos != list->end(); ++pos ) {
        nonlocStatus = ( TrabBoneNLStatus * ) this->giveStatus( ( * pos ).nearGp );
        nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= ( * pos ).weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain *= 1. / status->giveIntegrationScale();

    double localCumPlastStrain = status->giveLocalCumPlastStrainForAverage();
    alpha = mParam * nonlocalCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}

//
// END: SUBROUTINE OF NONLOCAL ALPHA EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERFACE ????
//

Interface *
TrabBoneNL :: giveInterface(InterfaceType type)
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

//
// END: INTERFACE ????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PARAMETERS OF INPUT FILE
//

IRResultType
TrabBoneNL :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    TrabBoneMaterial :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, IFT_TrabBoneNL_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_TrabBoneNL_m, "m"); // Macro

    return IRRT_OK;
}

//
// END: PARAMETERS OF INPUT FILE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: ??????????????
//

int
TrabBoneNL :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    TrabBoneMaterial :: giveInputRecordString(str, keyword);
    StructuralNonlocalMaterialExtensionInterface :: giveInputRecordString(str, false);
    sprintf(buff, " r %e", this->R);
    str += buff;

    return 1;
}

//
// END: ????????????????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: WEIGHT FUNCTION FOR NONLOCAL CONTRIBUTION
//

double
TrabBoneNL :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    double dist = src.distance(coord);

    if ( ( dist >= 0. ) && ( dist <= this->R ) ) {
        double help = ( 1. - dist * dist / ( R * R ) );
        return help * help;
    }

    return 0.0;
}

//
// END: WEIGHT FUNCTION FOR NONLOCAL CONTRIBUTION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
/////////TRABECULAR BONE NONLOCAL MATERIAL STATUS////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
// init state variables

TrabBoneNLStatus :: TrabBoneNLStatus(int n, Domain *d, GaussPoint *g) :
    TrabBoneMaterialStatus(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localCumPlastStrainForAverage = 0.0;
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneNLStatus :: ~TrabBoneNLStatus()
{}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PRINTOUT
//

void
TrabBoneNLStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    int n = epsp.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", epsp.at(i) );
    }

    fprintf(file, ",");
    fprintf(file, " alpha % .4e,", alpha);
    fprintf(file, " dam  % .4e", dam);
    fprintf(file, "}\n");
}

//
// END: PRINTOUT
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
// initialize temporary state variables according to equilibriated state vars

void
TrabBoneNLStatus :: initTempStatus()
{
    TrabBoneMaterialStatus :: initTempStatus();
}

//
// END: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STEP
// Called when equlibrium reached, set equilibriated vars according to temporary (working) ones.

void
TrabBoneNLStatus :: updateYourself(TimeStep *atTime)
{
    TrabBoneMaterialStatus :: updateYourself(atTime);
}

//
// END: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERFACE ???????????
//

Interface *
TrabBoneNLStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return this;
    } else {
        return NULL;
    }
}

//
// END: INTERFACE ???????????????
/////////////////////////////////////////////////////////////////
} // end namespace oofem

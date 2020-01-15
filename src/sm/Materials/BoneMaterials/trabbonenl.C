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

#include "trabbonenl.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "mathfem.h"
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

REGISTER_Material(TrabBoneNL);

/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
//

TrabBoneNL :: TrabBoneNL(int n, Domain *d) : TrabBoneMaterial(n, d), StructuralNonlocalMaterialExtensionInterface(d)
{
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
// update local values of accumulated pastic strain

void
TrabBoneNL :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep) const
{
    auto nlstatus = static_cast< TrabBoneNLStatus * >( this->giveStatus(gp) );

    FloatArray SDstrainVector;

    this->initTempStatus(gp);
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    nlstatus->letTempStrainVectorBe(strainVector);

    this->performPlasticityReturn(gp, strainVector);
    double cumPlastStrain = this->computeLocalCumPlastStrain(strainVector, gp, tStep);
    nlstatus->setLocalCumPlastStrainForAverage(cumPlastStrain);
}

//
// END: SUBROUTINE FOR UPDATE BEFORE NON LOCAL AVERAGE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
//

FloatArrayF<1>
TrabBoneNL :: giveRealStressVector_1d(const FloatArrayF<1> &strainVector,
                                      GaussPoint *gp, TimeStep *tStep) const
{
    auto nlStatus = static_cast< TrabBoneNLStatus * >( this->giveStatus(gp) );

    double tempDamage = computeDamage(gp, tStep);
    double plasStrain = nlStatus->giveTempPlasStrainVector().at(1);
    double elasStrain = strainVector.at(1) - plasStrain;
    double effStress = E0 * elasStrain;
    double sigc = nlStatus->giveSigC();

    FloatArrayF<1> stress = {( 1. - tempDamage ) * effStress + sigc};

    nlStatus->setTempDam(tempDamage);
    nlStatus->letTempStrainVectorBe(strainVector);
    nlStatus->letTempStressVectorBe(stress);
    return stress;
}

//
// END: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF NONLOCAL ALPHA EVALUATION
//

double
TrabBoneNL :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneNLStatus * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    auto list = status->giveIntegrationDomainList();

    double nonlocalCumPlastStrain = 0.0;
    for ( auto &lir: *list ) {
        auto nonlocStatus = static_cast< TrabBoneNLStatus * >( this->giveStatus(lir.nearGp) );
        double nonlocalContribution = nonlocStatus->giveLocalCumPlastStrainForAverage();
        nonlocalContribution *= lir.weight;
        nonlocalCumPlastStrain += nonlocalContribution;
    }

    nonlocalCumPlastStrain *= 1. / status->giveIntegrationScale();

    double localCumPlastStrain = status->giveLocalCumPlastStrainForAverage();
    return mParam * nonlocalCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
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
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

//
// END: INTERFACE ????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PARAMETERS OF INPUT FILE
//

void
TrabBoneNL :: initializeFrom(InputRecord &ir)
{
    TrabBoneMaterial :: initializeFrom(ir);
    StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, _IFT_TrabBoneNL_r);
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneNL_m);
}

//
// END: PARAMETERS OF INPUT FILE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: ??????????????
//

void
TrabBoneNL :: giveInputRecord(DynamicInputRecord &input)
{
    TrabBoneMaterial :: giveInputRecord(input);
    input.setField(this->R, _IFT_TrabBoneNL_r);
    input.setField(this->mParam, _IFT_TrabBoneNL_m);
}


//
// END: ????????????????
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: WEIGHT FUNCTION FOR NONLOCAL CONTRIBUTION
//

double
TrabBoneNL :: computeWeightFunction(const FloatArray &src, const FloatArray &coord) const
{
    double dist = distance(src, coord);

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


TrabBoneNLStatus :: TrabBoneNLStatus(GaussPoint *g) :
    TrabBoneMaterialStatus(g), StructuralNonlocalMaterialStatusExtensionInterface()
{
}


void
TrabBoneNLStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, " plastrains ");
    for ( auto &val : epsp) {
        fprintf( file, " %.4e", val );
    }

    fprintf(file, ",");
    fprintf(file, " alpha %.4e,", alpha);
    fprintf(file, " dam  %.4e", dam);
    fprintf(file, "}\n");
}

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
TrabBoneNLStatus :: updateYourself(TimeStep *tStep)
{
    TrabBoneMaterialStatus :: updateYourself(tStep);
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
        return nullptr;
    }
}

//
// END: INTERFACE ???????????????
/////////////////////////////////////////////////////////////////
} // end namespace oofem

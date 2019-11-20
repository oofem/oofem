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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "latticeslip.h"
#include "latticelinearelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "floatarrayf.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(LatticeSlip);

LatticeSlip :: LatticeSlip(int n, Domain *d) : LatticeLinearElastic(n, d)
{}


bool
LatticeSlip :: hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeSlip :: initializeFrom(InputRecord &ir)
{
    LatticeLinearElastic :: initializeFrom(ir);

    //Parameter which relates the shear stiffness to the normal stiffness. Default is 1000.
    alphaOne = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeSlip_a1); // Macro

    //Parameter which is used for the definition of bending stiffness. Default is 1000.
    alphaTwo = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeSlip_a2); // Macro

    //Parameter which limits the stress in slip direction.
    IR_GIVE_FIELD(ir, tauZero, _IFT_LatticeSlip_t0); // Macro
}


MaterialStatus *
LatticeSlip :: CreateStatus(GaussPoint *gp) const
{
    LatticeSlipStatus *answer = new LatticeSlipStatus(gp);

    return answer;
}


FloatArrayF< 6 >
LatticeSlip :: giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain, GaussPoint *gp, TimeStep *atTime)
{
    FloatArray answer;
    answer.resize(6);
    answer.zero();


    LatticeSlipStatus *status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );
    status->initTempStatus();

    FloatArray strainVector;

    FloatArray tempPlasticStrain = status->givePlasticStrain();

    FloatMatrix stiffnessMatrix = this->give3dLatticeStiffnessMatrix(ElasticStiffness, gp, atTime);

    /*First component is the slip one for which the stress should be limited using plasiticity (frictional slip between fibre and matrix). The other components are kept elastic. */
    answer.at(1) = ( totalStrain.at(1) - tempPlasticStrain.at(1) ) * stiffnessMatrix.at(1, 1);

    double f = fabs(answer.at(1) ) - tauZero;

    if ( f > 0 ) {//plastic response.
        //Reduced stress by increasing plastic strain.
        tempPlasticStrain.at(1) = tempPlasticStrain.at(1) + sgn(answer.at(1) ) * f / stiffnessMatrix.at(1, 1);
        answer.at(1) = ( totalStrain.at(1) - tempPlasticStrain.at(1) ) * stiffnessMatrix.at(1, 1);

        status->setTempCrackFlag(1);
    }

    //Compute the final stress components
    for ( int i = 2; i <= 6; i++ ) {
        answer.at(i) =  stiffnessMatrix.at(i, i) * totalStrain.at(i);
    }

    //Set temp values in status needed for dissipation
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation;

    tempDeltaDissipation = computeDeltaDissipation(gp, atTime);

    tempDissipation += tempDeltaDissipation;

    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    return answer;
}


Interface *
LatticeSlip :: giveInterface(InterfaceType type)
{
    return NULL;
}


LatticeSlipStatus :: LatticeSlipStatus(GaussPoint *g) :  LatticeMaterialStatus(g)
{}


void
LatticeSlipStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();
    //Only first 3 components are used for plastic strain??
    if ( this->plasticStrain.giveSize() == 0 ) {
        this->plasticStrain.resize(6);
        this->plasticStrain.zero();
    }

    this->tempPlasticStrain = this->plasticStrain;
}


FloatArrayF< 6 >
LatticeSlip :: giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);

    //Option to add a eigendisplacement instead of strain
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    alpha += this->cAlpha / length;

    return {
               alpha, 0., 0., 0., 0., 0.
    };
}

void
LatticeSlipStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "plasticStrain %.8e, dissipation %f, deltaDissipation %f, crackFlag %d\n", this->plasticStrain.at(1), this->dissipation, this->deltaDissipation, this->crackFlag);
}


double
LatticeSlip :: computeDeltaDissipation(GaussPoint *gp,
                                       TimeStep *atTime)
{
    LatticeSlipStatus *status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );

    FloatArray plasticStrain = status->givePlasticStrain();
    FloatArray tempPlasticStrain = status->giveTempPlasticStrain();
    FloatArray tempStress = status->giveTempStressVector();

    double tempDeltaDissipation =  tempStress.at(1) * ( tempPlasticStrain.at(1) - plasticStrain.at(1) );
    return tempDeltaDissipation;
}


void
LatticeSlipStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    LatticeMaterialStatus :: saveContext(stream, mode);

    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return;
}


void
LatticeSlipStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    LatticeMaterialStatus :: restoreContext(stream, mode);

    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return;
}

void
LatticeSlipStatus :: updateYourself(TimeStep *atTime)
{
    LatticeMaterialStatus :: updateYourself(atTime);
    this->plasticStrain = this->tempPlasticStrain;
}



int
LatticeSlip :: giveIPValue(FloatArray &answer,
                           GaussPoint *gp,
                           InternalStateType type,
                           TimeStep *atTime)
{
    LatticeSlipStatus *status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );

    if ( type == IST_CrackStatuses ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveCrackFlag();
        return 1;
    } else if ( type == IST_DissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    } else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
        answer.zero();
        answer.at(1) = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        return 1;
    } else {
        return LatticeLinearElastic :: giveIPValue(answer, gp, type, atTime);
    }
}
}

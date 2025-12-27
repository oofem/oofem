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

    //Parameter which relates the shear stiffness to the normal stiffness. Default is 1000. Reread this again, because the
    alphaOne = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LatticeSlip_a1); // Macro

    //Parameter which is used for the definition of bending stiffness. Default is 1000.
    alphaTwo = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaTwo, _IFT_LatticeSlip_a2); // Macro

    //Parameter which limits the stress in slip direction.
    IR_GIVE_FIELD(ir, tauZero, _IFT_LatticeSlip_t0); // Macro
}


std::unique_ptr<MaterialStatus> 
LatticeSlip :: CreateStatus(GaussPoint *gp) const
{
    return std::make_unique<LatticeSlipStatus>(gp);
}


FloatArrayF< 6 >
LatticeSlip :: giveLatticeStress3d(const FloatArrayF< 6 > &totalStrain, GaussPoint *gp, TimeStep *atTime)
{
    auto status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );
    status->initTempStatus();

    auto tempPlasticStrain = status->givePlasticLatticeStrain();

    auto stiffnessMatrix = this->give3dLatticeStiffnessMatrix(ElasticStiffness, gp, atTime);

    /*First component is the slip one for which the stress should be limited using plasiticity (frictional slip between fibre and matrix). The other components are kept elastic. */
    FloatArrayF<6> stress;

    stress.at(1) = ( totalStrain.at(1) - tempPlasticStrain.at(1) ) * stiffnessMatrix.at(1, 1);
    double f = fabs(stress.at(1) ) - tauZero;

    if ( f > 0 ) {//plastic response.
        //Reduced stress by increasing plastic strain.
        tempPlasticStrain.at(1) = tempPlasticStrain.at(1) + sgn(stress.at(1) ) * f / stiffnessMatrix.at(1, 1);
        stress.at(1) = ( totalStrain.at(1) - tempPlasticStrain.at(1) ) * stiffnessMatrix.at(1, 1);

        status->setTempCrackFlag(1);
    }

    //Compute the final stress components
    for ( int i = 2; i <= 6; i++ ) {
        stress.at(i) =  stiffnessMatrix.at(i, i) * totalStrain.at(i);
    }

    //Set temp values in status needed for dissipation
    status->letTempPlasticLatticeStrainBe(tempPlasticStrain);
    status->letTempLatticeStrainBe(totalStrain);
    status->letTempLatticeStressBe(stress);

    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation;

    tempDeltaDissipation = computeDeltaDissipation(gp, atTime);

    tempDissipation += tempDeltaDissipation;

    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    return stress;
}


Interface *
LatticeSlip :: giveInterface(InterfaceType type)
{
    return nullptr;
}


LatticeSlipStatus :: LatticeSlipStatus(GaussPoint *g) :  LatticeMaterialStatus(g)
{}


void
LatticeSlipStatus :: initTempStatus()
{
    LatticeMaterialStatus :: initTempStatus();
    this->tempPlasticLatticeStrain = this->plasticLatticeStrain;
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
    double length = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
    alpha += this->cAlpha / length;

    return {
               alpha, 0., 0., 0., 0., 0.
    };
}

void
LatticeSlipStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "plasticStrain %.8e, dissipation %f, deltaDissipation %f, crackFlag %d\n", this->plasticLatticeStrain.at(1), this->dissipation, this->deltaDissipation, this->crackFlag);
}


double
LatticeSlip :: computeDeltaDissipation(GaussPoint *gp, TimeStep *atTime) const
{
    auto status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );

    const auto &plasticStrain = status->givePlasticLatticeStrain();
    const auto &tempPlasticStrain = status->giveTempPlasticLatticeStrain();
    const auto &tempStress = status->giveTempLatticeStress();

    return tempStress.at(1) * ( tempPlasticStrain.at(1) - plasticStrain.at(1) );
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

    if ( ( iores = this->plasticLatticeStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
LatticeSlipStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;

    LatticeMaterialStatus :: restoreContext(stream, mode);

    if ( ( iores = this->plasticLatticeStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
LatticeSlipStatus :: updateYourself(TimeStep *atTime)
{
    LatticeMaterialStatus :: updateYourself(atTime);
    this->plasticLatticeStrain = this->tempPlasticLatticeStrain;
}



int
LatticeSlip :: giveIPValue(FloatArray &answer,
                           GaussPoint *gp,
                           InternalStateType type,
                           TimeStep *atTime)
{
    auto status = static_cast< LatticeSlipStatus * >( this->giveStatus(gp) );

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

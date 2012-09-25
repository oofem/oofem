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

#include "isoheatmat.h"
#include "flotmtrx.h"
#include "gausspnt.h"

namespace oofem {
IRResultType
IsotropicHeatTransferMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // double value ;

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, conductivity, IFT_IsotropicHeatTransferMaterial_k, "k"); // Macro// conductivity
    IR_GIVE_FIELD(ir, capacity, IFT_IsotropicHeatTransferMaterial_c, "c"); // Macro// specific heat capacity

    return IRRT_OK;
}

double
IsotropicHeatTransferMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. 'k' the conductivity of the receiver).
//
{
    if ( aProperty == 'k' ) { //thermal conductivity [J/m/K]
        return conductivity;
    } else if ( aProperty == 'c' ) { //mass-specific heat capacity [J/kg]
        return capacity;
    } else if ( aProperty == HeatCapaCoeff ) { //volume-specific heat capacity [J/m3]
        return ( capacity * this->give('d', gp) );
    }

    return this->Material :: give(aProperty, gp);
}

void
IsotropicHeatTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                          MatResponseForm form,
                                                          MatResponseMode mode,
                                                          GaussPoint *gp,
                                                          TimeStep *atTime)
{
    /*
     * returns constitutive (conductivity) matrix of receiver
     */
    MaterialMode mMode = gp->giveMaterialMode();
    double cond = this->giveIsotropicConductivity(gp);
    
    /*if ( !isActivated(atTime) ) //element, which is inactive (activityLTF==0), will never go into this function
         cond = 0.;
    }
    */
    
    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
        answer.at(1, 1) = cond;
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        return;

    case _3dHeat:
        answer.resize(3, 3);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        answer.at(3, 3) = cond;
        return;
    
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


double
IsotropicHeatTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
{
    if ( mode == Capacity ) {
        return ( capacity * this->give('d', gp) );
    } else {
        _error2( "giveCharacteristicValue : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return 0.;
}

int
IsotropicHeatTransferMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_HydrationDegree ) {
        return 1;
    } else {
        return TransportMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


int
IsotropicHeatTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if (  type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1)=0.;
        return 1;
    }

    return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
}



} // end namespace oofem

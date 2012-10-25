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

#include "transportmaterial.h"
#include "gausspnt.h"
#include "contextioerr.h"

namespace oofem {
void
TransportMaterial :: updateInternalState(const FloatArray &stateVec, GaussPoint *gp, TimeStep *)
{
    TransportMaterialStatus *ms = ( TransportMaterialStatus * ) this->giveStatus(gp);
    if ( ms ) {
        ms->letTempStateVectorBe(stateVec);
    }
}


TransportMaterialStatus :: TransportMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), stateVector(), tempStateVector()
{ }

void TransportMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Print the state variable and the flow vector on the data file.
{
    int i;
    FloatArray flowVec;
    TransportElement *transpElem = ( TransportElement * ) gp->giveElement();

    MaterialStatus :: printOutputAt(File, tNow);

    fprintf(File, "  state");

    for ( i = 1; i <= stateVector.giveSize(); i++ ) {
        fprintf( File, " % .4e", stateVector.at(i) );
    }

    transpElem->computeFlow(flowVec, gp, tNow);

    fprintf(File, "   flow");
    for ( i = 1; i <= flowVec.giveSize(); i++ ) {
        fprintf( File, " % .4e", flowVec.at(i) );
    }

    fprintf(File, "\n");
}


void TransportMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);
    stateVector = tempStateVector;
}


void
TransportMaterialStatus :: initTempStatus()
//
// initialize record at the beginning of new load step
//
{
    MaterialStatus :: initTempStatus();
    tempStateVector = stateVector;
}


contextIOResultType
TransportMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stateVector.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
TransportMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stateVector.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


int
TransportMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
// IST_Humidity must be overriden!


{
    if ( ( type == IST_Temperature ) || ( type == IST_MassConcentration_1 ) || ( type == IST_Humidity ) ) {
        FloatArray vec = ( ( TransportMaterialStatus * ) this->giveStatus(aGaussPoint) )->giveStateVector();
        answer.resize(1);
        answer.at(1) = vec.at( ( type == IST_Temperature ) ? 1 : 2 );
        return 1;
    } else if ( type == IST_TemperatureFlow ) {
        TransportElement *transpElem = ( TransportElement * ) aGaussPoint->giveElement();
        transpElem->computeFlow(answer, aGaussPoint, atTime);
        return 1;
    } else {
        return Material :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
TransportMaterial :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_Temperature || type == IST_MassConcentration_1 || type == IST_Humidity ) {
        return ISVT_SCALAR;
    } else {
        return Material :: giveIPValueType(type);
    }
}


int
TransportMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_Temperature || type == IST_MassConcentration_1 || type == IST_Humidity || type == IST_HydrationDegree || type == IST_Density || type == IST_ThermalConductivityIsotropic || type == IST_HeatCapacity || type == IST_AverageTemperature  || type == IST_YoungModulusVirginPaste || type == IST_PoissonRatioVirginPaste || type == IST_YoungModulusConcrete || type == IST_PoissonRatioConcrete ) {
        answer.setValues(1, 1);
        return 1;
    } else if ( type == IST_TemperatureFlow || type == IST_MassConcentrationFlow_1 || type == IST_HumidityFlow  ) {
        answer.setValues(3, 1, 2, 3);
        return 1;
    } else {
        return Material :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
TransportMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    int size = 0;
    MaterialMode mMode = aGaussPoint->giveMaterialMode();
    switch  ( mMode ) {
    case _2dHeat:
    case _2dHeMo:
        size = 2;
        break;
    case _3dHeat:
    case _3dHeMo:
        size = 3;
        break;
    default:
        _error2( "Unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    if ( type == IST_Temperature || type == IST_MassConcentration_1 || type == IST_HydrationDegree || type == IST_Humidity || type == IST_Density || type == IST_MaterialNumber || type == IST_ElementNumber || type == IST_ThermalConductivityIsotropic || type == IST_HeatCapacity || type == IST_AverageTemperature || type == IST_YoungModulusVirginPaste || type == IST_PoissonRatioVirginPaste || type == IST_YoungModulusConcrete || type == IST_PoissonRatioConcrete ) {
        return 1;
    } else if ( type == IST_TemperatureFlow || type == IST_MassConcentrationFlow_1 || type == IST_HumidityFlow ) {
        return size;
    } else {
        return Material :: giveIPValueSize(type, aGaussPoint);
    }
}

} // end namespace oofem

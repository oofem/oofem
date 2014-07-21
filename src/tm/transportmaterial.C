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

#include "transportmaterial.h"
#include "gausspoint.h"
#include "contextioerr.h"

namespace oofem {
void
TransportMaterialStatus :: setTempGradient(FloatArray grad)
{
    this->temp_gradient = std :: move(grad);
}

void
TransportMaterialStatus :: setTempField(FloatArray field)
{
    this->temp_field = std :: move(field);
}

void
TransportMaterialStatus :: setTempFlux(FloatArray w)
{
    this->temp_flux = std :: move(w);
}

void
TransportMaterial :: updateInternalState(const FloatArray &stateVec, GaussPoint *gp, TimeStep *)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    if ( ms ) {
        ms->letTempStateVectorBe(stateVec);
    }
}


TransportMaterialStatus :: TransportMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), temp_field(), temp_gradient(), temp_flux(), field(), gradient(), flux(), maturity(0.)
{ }

void TransportMaterialStatus :: printOutputAt(FILE *File, TimeStep *tNow)
// Print the state variable and the flow vector on the data file.
{
    FloatArray flowVec;
    TransportElement *transpElem = static_cast< TransportElement * >( gp->giveElement() );

    MaterialStatus :: printOutputAt(File, tNow);

    fprintf(File, "  state");

    for ( int i = 1; i <= field.giveSize(); i++ ) {
        fprintf( File, " % .4e", field.at(i) );
    }

    transpElem->computeFlow(flowVec, gp, tNow);

    fprintf(File, "   flow");
    for ( int i = 1; i <= flowVec.giveSize(); i++ ) {
        fprintf( File, " % .4e", flowVec.at(i) );
    }

    fprintf(File, "\n");
}


void
TransportMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);
    gradient = temp_gradient;
    field = temp_field;
    flux = temp_flux;
}


void
TransportMaterialStatus :: initTempStatus()
//
// initialize record at the beginning of new load step
//
{
    MaterialStatus :: initTempStatus();
    temp_gradient = gradient;
    temp_field = field;
    temp_flux = flux;
}


contextIOResultType
TransportMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        OOFEM_ERROR("can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = gradient.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = field.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = flux.storeYourself(stream, mode) ) != CIO_OK ) {
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
        OOFEM_ERROR("can't write into NULL stream");
    }

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = gradient.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( ( iores = field.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( ( iores = flux.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


int
TransportMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
// IST_Humidity must be overriden!
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_Temperature || type == IST_MassConcentration_1 || type == IST_Humidity ) {
        FloatArray vec = ms->giveField();
        answer = FloatArray({vec.at( ( type == IST_Temperature ) ? 1 : 2 ) });
        return 1;
    } else if ( type == IST_TemperatureFlow ) {
        TransportElement *transpElem = static_cast< TransportElement * >( gp->giveElement() );
        transpElem->computeFlow(answer, gp, tStep);
        return 1;
    } else if ( type == IST_Velocity ) { ///@todo Shouldn't be named velocity.. instead, "MassFlow" or something suitable like that.
        answer = ms->giveFlux();
        answer.resizeWithValues(3);
        return 1;
    } else if ( type == IST_PressureGradient ) {
        answer = ms->giveGradient();
        answer.resizeWithValues(3);
        return 1;
    } else if ( type == IST_Density ) {
        answer = FloatArray({ this->give('d', gp) });
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer = FloatArray({ this->give('c', gp) });
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer = FloatArray({ this->give('k', gp) });
        return 1;
    } else if ( type == IST_Maturity ) {
        answer = FloatArray({ ms->giveMaturity() });
        return 1;
    }
    return Material :: giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem

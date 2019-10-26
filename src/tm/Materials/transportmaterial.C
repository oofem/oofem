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

#include "tm/Materials/transportmaterial.h"
#include "gausspoint.h"
#include "contextioerr.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "datastream.h"

namespace oofem {

TransportMaterialStatus :: TransportMaterialStatus(GaussPoint *g) :
    MaterialStatus(g)
{ }

void TransportMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep) const
{
    MaterialStatus :: printOutputAt(File, tStep);

    fprintf(File, "  state %.4e", field);

    fprintf(File, "   flow");
    for ( auto &flow : flux ) {
        fprintf( File, " %.4e", flow );
    }

    fprintf(File, "\n");
}


void
TransportMaterialStatus :: updateYourself(TimeStep *tStep)
{
    MaterialStatus :: updateYourself(tStep);
    gradient = temp_gradient;
    field = temp_field;
    flux = temp_flux;
}


void
TransportMaterialStatus :: initTempStatus()
{
    MaterialStatus :: initTempStatus();
    temp_gradient = gradient;
    temp_field = field;
    temp_flux = flux;
}


void
TransportMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = gradient.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( stream.read(field) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = flux.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
TransportMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = gradient.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( !stream.read(field) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( ( iores = flux.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}



HeMoTransportMaterialStatus :: HeMoTransportMaterialStatus(GaussPoint *g) :
    MaterialStatus(g)
{ }


void
HeMoTransportMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep) const
{
    MaterialStatus :: printOutputAt(File, tStep);

    fprintf(File, "  temperature %.4e", temperature);

    fprintf(File, "   t_flux");
    for ( auto &flow : t_flux ) {
        fprintf( File, " %.4e", flow );
    }

    fprintf(File, "  humidity %.4e", humidity);

    fprintf(File, "   h_flux");
    for ( auto &flow : h_flux ) {
        fprintf( File, " %.4e", flow );
    }

    fprintf(File, "\n");
}


void
HeMoTransportMaterialStatus :: updateYourself(TimeStep *tStep)
{
    MaterialStatus :: updateYourself(tStep);

    t_gradient = temp_t_gradient;
    temperature = temp_temperature;
    t_flux = temp_t_flux;

    h_gradient = temp_h_gradient;
    humidity = temp_humidity;
    h_flux = temp_h_flux;
}


void
HeMoTransportMaterialStatus :: initTempStatus()
{
    MaterialStatus :: initTempStatus();

    temp_t_gradient = t_gradient;
    temp_temperature = temperature;
    temp_t_flux = t_flux;

    temp_h_gradient = h_gradient;
    temp_humidity = humidity;
    temp_h_flux = h_flux;
}


void
HeMoTransportMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = t_gradient.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( stream.read(temperature) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = t_flux.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = h_gradient.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( stream.read(humidity) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = h_flux.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
HeMoTransportMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = t_gradient.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( !stream.read(temperature) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( ( iores = t_flux.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = h_gradient.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    if ( !stream.read(humidity) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( ( iores = h_flux.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}



void
TransportMaterial :: updateInternalState(const FloatArray &stateVec, GaussPoint *gp, TimeStep *)
{
    if ( stateVec.giveSize() == 1 ) {
        auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
        if ( ms ) {
            ms->setTempField(stateVec[0]);
        }
    } else if ( stateVec.giveSize() == 2 ) {
        auto ms = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );
        if ( ms ) {
            ms->setTempTemperature(stateVec[0]);
            ms->setTempHumidity(stateVec[1]);
        }
    }
}


void
TransportMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep) const
{
    if ( field.giveSize() == 1 ) {
        if ( grad.giveSize() == 3 ) {
            answer = computeFlux3D(grad, field[0], gp, tStep);
        } else if ( grad.giveSize() == 2 ) {
            answer = computeFlux2D(grad, field[0], gp, tStep);
        } else {
            answer = computeFlux1D(grad, field[0], gp, tStep);
        }
    } else {
        // Has to be a HeMo-model:
        double t = field.at(1);
        double h = field.at(2);

        int size = grad.giveSize() / 2;
        FloatArrayF<3> grad_w, grad_t;
        for (int i = 0; i < size; ++i) {
            grad_t[i] = grad[i];
        }
        for (int i = 0; i < size; ++i) {
            grad_w[i] = grad[i+size];
        }

        auto grads = computeHeMoFlux3D(grad_t, grad_w, t, h, gp, tStep);

        answer.resize(size * 2);
        for (int i = 0; i < size; ++i) {
            answer[i] = grads.first[i];
        }

        for (int i = 0; i < size; ++i) {
            answer[i+size] = grads.second[i];
        }
    }
}

void
TransportMaterial :: giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    MaterialMode mmode = gp->giveMaterialMode();
    if ( mmode == _3dHeat || mmode == _3dMTLattice ) {
        answer = computeTangent3D(mode, gp, tStep);
    } else if ( mmode == _2dHeat || mmode == _2dMTLattice ) {
        answer = computeTangent2D(mode, gp, tStep);
    } else if ( mmode == _1dHeat ) {
        answer = computeTangent1D(mode, gp, tStep);
    } else if ( mmode == _3dHeMo ) {
        answer = computeTangent3D(mode, gp, tStep);
    } else if ( mmode == _2dHeMo ) {
        answer = computeTangent2D(mode, gp, tStep);
    } else if ( mmode == _1dHeMo ) {
        answer = computeTangent1D(mode, gp, tStep);
    }
}


FloatArrayF<3>
TransportMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Method not overloaded");
}


FloatArrayF<2>
TransportMaterial :: computeFlux2D(const FloatArrayF<2> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ans = this->computeFlux3D({grad[0], grad[1], 0.}, field, gp, tStep);
    return {ans[0], ans[1]};
}


FloatArrayF<1>
TransportMaterial :: computeFlux1D(const FloatArrayF<1> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ans = this->computeFlux3D({grad[0], 0., 0.}, field, gp, tStep);
    return {ans[0]};
}


std::pair<FloatArrayF<3>, FloatArrayF<3>>
TransportMaterial :: computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Method is only validfor heat+moisture coupled materials.");
}

std::pair<FloatArrayF<2>, FloatArrayF<2>>
TransportMaterial :: computeHeMoFlux2D(const FloatArrayF<2> &grad_t, const FloatArrayF<2> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const
{
    auto grads = computeHeMoFlux3D({grad_t[0], grad_t[1], 0.}, {grad_w[0], grad_w[1], 0.}, t, h, gp, tStep);
    return {{grads.first[0], grads.first[1]}, {grads.second[0], grads.second[1]}};
}


std::pair<FloatArrayF<1>, FloatArrayF<1>>
TransportMaterial :: computeHeMoFlux1D(const FloatArrayF<1> &grad_t, const FloatArrayF<1> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const
{
    auto grads = computeHeMoFlux3D({grad_t[0], 0., 0.}, {grad_w[0], 0., 0.}, t, h, gp, tStep);
    return {{grads.first[0]}, {grads.second[0]}};
}


FloatMatrixF<2,2>
TransportMaterial :: computeTangent2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto x = this->computeTangent3D(mode, gp, tStep);
    return {
        x(0, 0), x(1, 0),
        x(0, 1), x(1, 1),
    };
}


FloatMatrixF<1,1>
TransportMaterial :: computeTangent1D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto x = this->computeTangent3D(mode, gp, tStep);
    return {x(0, 0)};
}


int
TransportMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
// IST_Humidity must be overriden!
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_Temperature || type == IST_MassConcentration_1 || type == IST_Humidity ) {
        answer = FloatArray{ms->giveField()};
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
        answer = FloatArray{ this->give('d', gp) };
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer = FloatArray{ this->give('c', gp) };
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer = FloatArray{ this->give('k', gp) };
        return 1;
    } else if ( type == IST_Maturity ) {
        answer = FloatArray{ ms->giveMaturity() };
        return 1;
    }
    return Material :: giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem

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

#include "fm/Materials/fluiddynamicmaterial.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "contextioerr.h"

namespace oofem {

std::pair<FloatArrayF<6>, double>
FluidDynamicMaterial :: computeDeviatoricStress3D(const FloatArrayF<6> &eps, double pressure, GaussPoint *gp, TimeStep *tStep) const
{
    double epsp_vol = -( eps[0] + eps[1] + eps[2] );
    auto stress_dev = this->computeDeviatoricStress3D(eps, gp, tStep);
    return {stress_dev, epsp_vol};
}


std::pair<FloatArrayF<3>, double>
FluidDynamicMaterial :: computeDeviatoricStress2D(const FloatArrayF<3> &eps, double pressure, GaussPoint *gp, TimeStep *tStep) const
{
    //auto [stress, epsv] = this->computeDeviatoricStress3D(assemble<6>(eps, {0, 1, 5}), pressure, gp, tStep);
    auto val = this->computeDeviatoricStress3D(assemble<6>(eps, {0, 1, 5}), pressure, gp, tStep);
    auto stress3 = val.first;
    auto epsv = val.second;
    return {stress3[{0, 1, 5}], epsv};
}


FloatArrayF<3>
FluidDynamicMaterial :: computeDeviatoricStress2D(const FloatArrayF<3> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    auto stress3 = this->computeDeviatoricStress3D(assemble<6>(eps, {0, 1, 5}), gp, tStep);
    return stress3[{0, 1, 5}];
}


FloatArrayF<4>
FluidDynamicMaterial :: computeDeviatoricStressAxi(const FloatArrayF<4> &eps, GaussPoint *gp, TimeStep *tStep) const
{
    auto stress3 = this->computeDeviatoricStress3D(assemble<6>(eps, {0, 1, 2, 5}), gp, tStep);
    return stress3[{0, 1, 2, 5}];
}


FloatMatrixF<3,3>
FluidDynamicMaterial :: computeTangent2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto dsdd3 = this->computeTangent3D(mode, gp, tStep);
    return dsdd3({0, 1, 5}, {0, 1, 5});
}


FloatMatrixF<4,4>
FluidDynamicMaterial :: computeTangentAxi(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto dsdd3 = this->computeTangent3D(mode, gp, tStep);
    return dsdd3({0, 1, 2, 5}, {0, 1, 2, 5});
}


FluidDynamicMaterial::Tangents<6>
FluidDynamicMaterial :: computeTangents3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    return {
        this->computeTangent3D(mode, gp, tStep),
        zeros<6>(),
        zeros<6>(),
        0.
    };
}


FluidDynamicMaterial::Tangents<3>
FluidDynamicMaterial :: computeTangents2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    Tangents<6> t = this->computeTangents3D(mode, gp, tStep);

    FloatMatrixF<3,3> dsdd = t.dsdd({0, 1, 5}, {0, 1, 5});
    FloatArrayF<3> dsdp = t.dsdp[{0, 1, 5}];
    FloatArrayF<3> dedd = t.dedd[{0, 1, 5}];
    return {dsdd, dsdp, dedd, t.dedp};
}


FluidDynamicMaterialStatus :: FluidDynamicMaterialStatus(GaussPoint *g) :
    MaterialStatus(g)
{ }

void
FluidDynamicMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    fprintf(file, "\n deviatoric stresses");
    for ( double e: deviatoricStressVector ) {
        fprintf(file, " %.4e", e );
    }

    fprintf(file, "\n");
}

int
FluidDynamicMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto status = static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DeviatoricStress ) {
        answer = status->giveDeviatoricStressVector();
        return 1;
    } else if ( type == IST_DeviatoricStrain ) {
        answer = status->giveDeviatoricStrainRateVector();
        return 1;
    } else if ( type == IST_Viscosity ) {
        answer.resize(1);
        answer.at(1) = this->giveEffectiveViscosity(gp, tStep);
        return 1;
    } else if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->give('d', gp);
        return 1;
    } else {
        return Material :: giveIPValue(answer, gp, type, tStep);
    }
}


void
FluidDynamicMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = deviatoricStressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStrainRateVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
FluidDynamicMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = deviatoricStressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStrainRateVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}
} // end namespace oofem

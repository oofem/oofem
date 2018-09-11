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

#include "fm/Materials/fluiddynamicmaterial.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "contextioerr.h"

namespace oofem {

void
FluidDynamicMaterial :: computeDeviatoricStress3D(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep)
{
    epsp_vol = -( eps[0] + eps[1] + eps[2] );
    this->computeDeviatoricStress3D(stress_dev, gp, eps, tStep);
}


void
FluidDynamicMaterial :: computeDeviatoricStress2D(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep)
{
    FloatArray stress3, eps3 = {eps[0], eps[1], 0., 0., 0., eps[2]};
    this->computeDeviatoricStress3D(stress3, epsp_vol, gp, eps3, pressure, tStep);
    stress_dev = {stress3[0], stress3[1], stress3[5]};
}


void
FluidDynamicMaterial :: computeDeviatoricStress2D(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    FloatArray stress3, eps3 = {eps[0], eps[1], 0., 0., 0., eps[2]};
    this->computeDeviatoricStress3D(stress3, gp, eps3, tStep);
    answer = {stress3[0], stress3[1], stress3[5]};
}


void
FluidDynamicMaterial :: computeDeviatoricStressAxi(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    FloatArray stress3, eps3 = {eps[0], eps[1], eps[2], 0., 0., eps[3]};
    this->computeDeviatoricStress3D(stress3, gp, eps3, tStep);
    answer = {stress3[0], stress3[1], stress3[2], stress3[5]};
}


void
FluidDynamicMaterial :: computeTangent2D(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dsdd3;
    this->computeTangent3D(dsdd3, mode, gp, tStep);
    answer = {
        {dsdd3(0, 0), dsdd3(0, 1), dsdd3(0, 5)},
        {dsdd3(1, 0), dsdd3(1, 1), dsdd3(0, 5)},
        {dsdd3(5, 0), dsdd3(5, 1), dsdd3(5, 5)}
    };
}


void
FluidDynamicMaterial :: computeTangentAxi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dsdd3;
    this->computeTangent3D(dsdd3, mode, gp, tStep);
    answer = {
        {dsdd3(0, 0), dsdd3(0, 1), dsdd3(0, 2), dsdd3(0, 5)},
        {dsdd3(1, 0), dsdd3(1, 1), dsdd3(1, 2), dsdd3(1, 5)},
        {dsdd3(2, 0), dsdd3(2, 1), dsdd3(2, 2), dsdd3(2, 5)},
        {dsdd3(5, 0), dsdd3(5, 1), dsdd3(5, 2), dsdd3(5, 5)}
    };
}


void
FluidDynamicMaterial :: computeTangents3D(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                          MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    this->computeTangent3D(dsdd, mode, gp, tStep);
    int size = dsdd.giveNumberOfRows();
    dsdp.resize(size);
    dsdp.zero();
    dedd.resize(size);
    dedd.zero();
    dedp = 0;
}


void
FluidDynamicMaterial :: computeTangents2D(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                          MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dsdd3;
    FloatArray dsdp3, dedd3;
    this->computeTangents3D(dsdd3, dsdp3, dedd3, dedp, mode, gp, tStep);

    dsdd = {
        {dsdd3(0, 0), dsdd3(0, 1), dsdd3(0, 5)},
        {dsdd3(1, 0), dsdd3(1, 1), dsdd3(1, 5)},
        {dsdd3(5, 0), dsdd3(5, 1), dsdd3(5, 5)}
    };
    dsdp = {dsdp3[0], dsdp3[1], dsdp3[5]};
    dedd = {dedd3[0], dedd3[1], dedd3[5]};
}


FluidDynamicMaterialStatus :: FluidDynamicMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), deviatoricStressVector(6), deviatoricStrainRateVector(6)
{ }

void
FluidDynamicMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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
    FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) );
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

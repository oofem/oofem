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

#include "fluiddynamicmaterial.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "contextioerr.h"

namespace oofem {
void
FluidDynamicMaterial :: updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep)
{ }


void
FluidDynamicMaterial :: computeDeviatoricStressVector(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep)
{
    if ( gp->giveMaterialMode() == _2dFlow ) {
        epsp_vol = -( eps.at(1) + eps.at(2) );
    } else {
        epsp_vol = -( eps.at(1) + eps.at(2) + eps.at(3) );
    }

    this->computeDeviatoricStressVector(stress_dev, gp, eps, tStep);
}

void
FluidDynamicMaterial :: giveStiffnessMatrices(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                              MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveDeviatoricStiffnessMatrix(dsdd, mode, gp, tStep);
    int size = dsdd.giveNumberOfRows();
    dsdp.resize(size);
    dsdp.zero();
    dedd.resize(size);
    dedd.zero();
    dedp = 0;
}


FluidDynamicMaterialStatus :: FluidDynamicMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), deviatoricStressVector(), deviatoricStrainRateVector()
{ }

void
FluidDynamicMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
{
    fprintf(File, "\n deviatoric stresses");
    for ( double e: deviatoricStressVector ) {
        fprintf( File, " %.4e", e );
    }

    fprintf(File, "\n");
}

void
FluidDynamicMaterialStatus :: initTempStatus()
{
    MaterialStatus :: initTempStatus();
}


int
FluidDynamicMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FluidDynamicMaterialStatus *status = static_cast< FluidDynamicMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DeviatoricStress ) {
        MaterialMode mmode = gp->giveMaterialMode();
        const FloatArray &vec = status->giveDeviatoricStressVector();
        if ( mmode == _2dFlow ) {
            answer.resize(6);
            answer.at(1) = vec.at(1);
            answer.at(2) = vec.at(2);
            answer.at(3) = -( vec.at(1) + vec.at(2) ); ///@todo Verify that this is correct for for all models.
            answer.at(4) = 0.;
            answer.at(5) = 0.;
            answer.at(6) = vec.at(3);
            return 1;
        } else if ( mmode == _2dAxiFlow ) {
            answer.resize(6);
            answer.at(1) = vec.at(1);
            answer.at(2) = vec.at(2);
            answer.at(3) = vec.at(3);
            answer.at(4) = 0;
            answer.at(5) = 0;
            answer.at(6) = vec.at(6);
            return 1;
        } else if ( mmode == _3dFlow ) {
            answer = vec;
            return 1;
        } else {
            OOFEM_ERROR("material mode not supported");
            return 0;
        }

    } else if ( type == IST_DeviatoricStrain ) {
        MaterialMode mmode = gp->giveMaterialMode();
        const FloatArray &vec = status->giveDeviatoricStrainRateVector();
        if ( mmode == _2dFlow ) {
            answer.resize(6);
            answer.at(1) = vec.at(1);
            answer.at(2) = vec.at(2);
            answer.at(3) = 0.;
            answer.at(4) = vec.at(3);
            answer.at(5) = 0.;
            answer.at(6) = 0.;
            answer.copySubVector(vec, 1);
            return 1;
        } else if ( mmode == _2dAxiFlow ) {
            answer.resize(6);
            answer.at(1) = vec.at(1);
            answer.at(2) = vec.at(2);
            answer.at(3) = vec.at(3);
            answer.at(4) = 0;
            answer.at(5) = 0;
            answer.at(6) = vec.at(6);
            return 1;
        } else if ( mmode == _3dFlow ) {
            answer = vec;
            return 1;
        } else {
            OOFEM_ERROR("material mode not supported");
            return 0;
        }

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


contextIOResultType
FluidDynamicMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStrainRateVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
FluidDynamicMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = deviatoricStrainRateVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem

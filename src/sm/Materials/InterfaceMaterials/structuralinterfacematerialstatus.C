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


#include "structuralinterfacematerial.h"
#include "contextioerr.h"

#include "sm/Materials/structuralmaterial.h"
#include "structuralinterfacematerialstatus.h"
#include "gausspoint.h"

namespace oofem {

StructuralInterfaceMaterialStatus :: StructuralInterfaceMaterialStatus(GaussPoint *g) :
    MaterialStatus(g)
{
    F = eye<3>();
    tempF = F;
}


void StructuralInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    MaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "  jump ");
    for ( auto &val : this->jump ) {
        fprintf(file, " %.4e", val );
    }

    fprintf(file, "\n              traction ");
    for ( auto &val : this->traction ) {
        fprintf(file, " %.4e", val );
    }
    fprintf(file, "\n");
}


void StructuralInterfaceMaterialStatus :: updateYourself(TimeStep *tStep)
{
    MaterialStatus :: updateYourself(tStep);

    this->jump            = this->tempJump;
    this->traction        = this->tempTraction;
    this->firstPKTraction = this->tempFirstPKTraction;
    this->F               = this->tempF;
}


void StructuralInterfaceMaterialStatus :: initTempStatus()
{
    MaterialStatus :: initTempStatus();

    this->tempJump            = this->jump;
    this->tempTraction        = this->traction;
    this->tempFirstPKTraction = this->firstPKTraction;
    this->tempF               = this->F;
}


void
StructuralInterfaceMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
#if 0
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
#endif
}


void
StructuralInterfaceMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
#if 0
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
#endif
}

void StructuralInterfaceMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    const StructuralInterfaceMaterialStatus &structStatus = static_cast< const StructuralInterfaceMaterialStatus & >(iStatus);

    jump                    = structStatus.giveJump();
    traction                = structStatus.giveTraction();
    tempTraction            = structStatus.giveTempTraction();
    tempJump                = structStatus.giveTempJump();
    firstPKTraction         = structStatus.giveFirstPKTraction();
    tempFirstPKTraction     = structStatus.giveTempFirstPKTraction();
    F                       = structStatus.giveF();
    tempF                   = structStatus.giveTempF();
    mNormalDir              = structStatus.giveNormal();
}


void StructuralInterfaceMaterialStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    OOFEM_ERROR("not implemented");
}
} // end namespace oofem

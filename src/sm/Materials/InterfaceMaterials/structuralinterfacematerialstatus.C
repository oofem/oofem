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


#include "structuralinterfacematerial.h"
#include "contextioerr.h"

#include "../sm/Materials/structuralmaterial.h"
#include "structuralinterfacematerialstatus.h"
#include "gausspoint.h"
namespace oofem {
StructuralInterfaceMaterialStatus :: StructuralInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g) :
    MaterialStatus(n, d, g), jump(3), traction(3), tempTraction(3), tempJump(3), firstPKTraction(3), tempFirstPKTraction(3), F(3, 3), tempF(3, 3),
    mNewlyInserted(true)
{
    this->F.beUnitMatrix();
    this->tempF.beUnitMatrix();
}


StructuralInterfaceMaterialStatus :: ~StructuralInterfaceMaterialStatus() { }


void StructuralInterfaceMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
#if 0
    MaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "  jump ");
    for ( auto &val : this->jump ) {
        fprintf(file, " %.4e", val );
    }

    fprintf(File, "\n              traction ");
    for ( auto &val : this->traction ) {
        fprintf( File, " %.4e", val );
    }
    fprintf(file, "\n");
#endif
}


void StructuralInterfaceMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);

    this->jump            = this->tempJump;
    this->traction        = this->tempTraction;
    this->firstPKTraction = this->tempFirstPKTraction;
    this->F               = this->tempF;
}


void StructuralInterfaceMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    MaterialStatus :: initTempStatus();

    // see if vectors describing reached equilibrium are defined
    if ( this->giveJump().giveSize() == 0 ) {
        
        this->jump.resize( this->giveDomain()->giveNumberOfSpatialDimensions() );
        //this->jump.resize(3);
        this->jump.zero();
    }

    if ( this->giveTraction().giveSize() == 0 ) {
        this->traction.resize( this->giveDomain()->giveNumberOfSpatialDimensions() );
        //this->traction.resize(3);
        this->traction.zero();
    }

    // reset temp vars.
    this->tempJump            = this->jump;
    this->tempTraction        = this->traction;
    this->tempFirstPKTraction = this->firstPKTraction;
    this->tempF               = this->F;
}


contextIOResultType
StructuralInterfaceMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
#if 0
    contextIOResultType iores;

    if ( ( iores = MaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
#endif
    return CIO_OK;
}


contextIOResultType
StructuralInterfaceMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
#if 0
    contextIOResultType iores;

    if ( ( iores = MaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = strainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
#endif
    return CIO_OK;
}

void StructuralInterfaceMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const StructuralInterfaceMaterialStatus &structStatus = dynamic_cast< StructuralInterfaceMaterialStatus & >(tmpStat);

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

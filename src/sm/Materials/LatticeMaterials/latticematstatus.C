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
 *               Copyright (C) 1993 - 2019   Borek Patzak
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

#include "latticematstatus.h"
#include "../structuralmaterial.h"
#include "contextioerr.h"
#include "datastream.h"
#include "gausspoint.h"

namespace oofem {
LatticeMaterialStatus :: LatticeMaterialStatus(GaussPoint *g) : MaterialStatus(g), RandomMaterialStatusExtensionInterface()
{ }


void
LatticeMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    MaterialStatus :: initTempStatus();

    this->tempLatticeStrain = this->latticeStrain;

    this->tempLatticeStress = this->latticeStress;

    this->tempReducedLatticeStrain = this->reducedLatticeStrain;

    this->tempNormalLatticeStress = this->normalLatticeStress;

    this->tempPlasticLatticeStrain = this->plasticLatticeStrain;

    this->tempDissipation = this->dissipation;
    this->tempDeltaDissipation = this->deltaDissipation;
    this->tempCrackFlag = this->crackFlag;
    this->tempCrackWidth = this->crackWidth;

    this->tempDamageLatticeStrain = this->damageLatticeStrain;

    this->updateFlag = 0;
}


void
LatticeMaterialStatus :: updateYourself(TimeStep *atTime)
{
    MaterialStatus :: updateYourself(atTime);

    this->latticeStress = this->tempLatticeStress;

    this->latticeStrain = this->tempLatticeStrain;

    this->plasticLatticeStrain = this->tempPlasticLatticeStrain;

    this->reducedLatticeStrain = this->tempReducedLatticeStrain;

    this->damageLatticeStrain = this->tempDamageLatticeStrain;

    this->dissipation = this->tempDissipation;

    this->deltaDissipation = this->tempDeltaDissipation;

    this->crackFlag = this->tempCrackFlag;

    this->crackWidth = this->tempCrackWidth;

    this->updateFlag = 1;
}


void
LatticeMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    MaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, " latticestrain ");
    for ( double s : this->latticeStrain ) {
        fprintf(file, "% .4e ", s);
    }
    fprintf(file, " latticestress ");
    for ( double s : this->latticeStress ) {
        fprintf(file, "% .4e ", s);
    }
    fprintf(file, " reducedlatticestrain ");
    for ( double s : this->reducedLatticeStrain ) {
        fprintf(file, "% .4e ", s);
    }
}


Interface *
LatticeMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return ( RandomMaterialStatusExtensionInterface * ) this;
    } else {
        return nullptr;
    }
}


void
LatticeMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = latticeStress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = latticeStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = reducedLatticeStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticLatticeStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = damageLatticeStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaDissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(crackFlag) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}


void
LatticeMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = damageLatticeStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = latticeStress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = latticeStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = reducedLatticeStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticLatticeStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaDissipation) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(crackFlag) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem

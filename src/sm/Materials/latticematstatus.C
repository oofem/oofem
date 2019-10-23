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
#include "structuralmaterial.h"
#include "contextioerr.h"
#include "datastream.h"
#include "gausspoint.h"

namespace oofem {
LatticeMaterialStatus :: LatticeMaterialStatus(GaussPoint *g) : StructuralMaterialStatus(g), RandomMaterialStatusExtensionInterface(), reducedStrain(), tempReducedStrain(), plasticStrain(), tempPlasticStrain(), oldPlasticStrain()
{
    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector(gp->giveMaterialMode() );
    plasticStrain.resize(rsize);
    plasticStrain.zero();
    tempPlasticStrain.resize(rsize);
    tempPlasticStrain.zero();
    reducedStrain.resize(rsize);
    reducedStrain.zero();
    oldPlasticStrain.resize(rsize);
    oldPlasticStrain.zero();
    tempReducedStrain.resize(rsize);
    tempReducedStrain.zero();

    normalStress = tempNormalStress = 0;

    le = 0.0;
    crackFlag = tempCrackFlag = 0;
    crackWidth = tempCrackWidth = 0;
    dissipation = tempDissipation = 0.;
    deltaDissipation = tempDeltaDissipation = 0.;
}


void
LatticeMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    StructuralMaterialStatus :: initTempStatus();

    this->oldPlasticStrain = this->plasticStrain;

    this->tempPlasticStrain = this->plasticStrain;

    this->tempReducedStrain = this->reducedStrain;

    this->tempNormalStress = this->normalStress;

    this->tempDissipation = this->dissipation;
    this->tempDeltaDissipation = this->deltaDissipation;
    this->tempCrackFlag = this->crackFlag;
    this->tempCrackWidth = this->crackWidth;

    this->updateFlag = 0;
}


void
LatticeMaterialStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    this->plasticStrain = this->tempPlasticStrain;

    this->reducedStrain = this->tempReducedStrain;

    this->dissipation = this->tempDissipation;

    this->deltaDissipation = this->tempDeltaDissipation;


    this->crackFlag = this->tempCrackFlag;

    this->crackWidth = this->tempCrackWidth;

    this->updateFlag = 1;
}

void
LatticeMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    // print only if reducedStrains are not empty
    if ( !reducedStrain.containsOnlyZeroes() ) {
        fprintf(file, "reduced strains ");
        int rSize = reducedStrain.giveSize();
        for ( int k = 1; k <= rSize; k++ ) {
            fprintf(file, "% .4e ", reducedStrain.at(k) );
        }
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
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;

    // write a raw data
    if ( ( iores = reducedStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write a raw data
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
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
    StructuralMaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;

    if ( ( iores = reducedStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
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

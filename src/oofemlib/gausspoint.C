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

#include "gausspoint.h"
#include "matstatus.h"
#include "material.h"

#include <memory>

namespace oofem {
//GaussPoint :: GaussPoint(IntegrationRule *ir, int n, FloatArray *a, double w, MaterialMode mode) : statusDict()
GaussPoint :: GaussPoint(IntegrationRule *ir, int n, FloatArray iNaturalCoord, double w, MaterialMode mode) :
    number(n),
    irule(ir),
    naturalCoordinates(std :: move(iNaturalCoord)),
    weight(w),
    materialMode(mode)
{
    subPatchCoordinates = NULL;
    globalCoordinates = NULL;
    materialStatus = NULL;
}

GaussPoint :: GaussPoint(IntegrationRule * ir, int n, double w, MaterialMode mode) :
    number(n),
    irule(ir),
    weight(w),
    materialMode(mode)
{
    subPatchCoordinates = NULL;
    globalCoordinates = NULL;
    materialStatus = NULL;
}

GaussPoint :: ~GaussPoint()
{
    for ( GaussPoint *gp: gaussPoints ) {
        delete gp;
    }

    delete subPatchCoordinates;
    delete globalCoordinates;
    delete materialStatus;
}


void GaussPoint :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
{
    int iruleNumber = 0;

    if ( irule ) {
        iruleNumber = irule->giveNumber();
    }

    fprintf(File, "  GP %2d.%-2d :", iruleNumber, number);

    // invoke printOutputAt method for all managed statuses
    IntegrationPointStatus *status = this->giveMaterialStatus();
    if ( status ) {
        status->printOutputAt(File, tStep);
    }

    if ( gaussPoints.size() != 0 ) { // layered material
        fprintf(File, "Layers report \n{\n");
        for ( GaussPoint *gp: gaussPoints ) {
            gp->printOutputAt(File, tStep);
        }

        fprintf(File, "} end layers report\n");
    }
}


GaussPoint *GaussPoint :: giveSlaveGaussPoint(int index)
// returns receivers slave gauss point
// 'slaves' are introduced in order to support various type
// of cross sections models (for example layered material, where
// each separate layer has its own slave gp.)
//
{
    if ( gaussPoints.size() == 0 ) {
        return NULL;
    }

    if ( ( index < 0 ) || ( index >= (int)gaussPoints.size() ) ) {
        OOFEM_ERROR("index out of range");
    }

    return gaussPoints [ index ];
}

bool GaussPoint :: hasSlaveGaussPoint() {
    if ( this->gaussPoints.size() != 0 ) {// layered material
        return true;
    }
    return false;
}

size_t GaussPoint :: findFirstIndexOfSlaveGaussPoint(GaussPoint *gp) {
    std::vector< GaussPoint * >::iterator it;
    it = find( gaussPoints.begin(), gaussPoints.end(), gp );
    if ( it != gaussPoints.end() ) {
        return it - gaussPoints.begin();
    } else {
        OOFEM_ERROR("Slave Gauss point not found.");
    }
}

void GaussPoint :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    IntegrationPointStatus *status = this->giveMaterialStatus();
    if ( status ) {
        status->updateYourself(tStep);
    }

    for ( GaussPoint *gp: gaussPoints ) {
        gp->updateYourself(tStep);
    }
}



/*
 * contextIOResultType
 * GaussPoint :: saveContext (FILE* stream, void *obj)
 * //
 * // saves full gp context (saves state variables, that completely describe
 * // current state)
 * // does not saves the slave - records
 * // this task is done at the layeredCrossSection level
 * {
 *
 * contextIOResultType iores;
 *
 * if ((iores = this->giveMaterial()->saveContext(stream,(void*) this)) != CIO_OK) THROW_CIOERR(iores);
 * // if (matStatusDict->saveContext(stream,obj) != 1)
 * //   error ("saveContext io error encountered");
 *
 * return CIO_OK;
 *
 * }
 *
 *
 * contextIOResultType
 * GaussPoint :: restoreContext (FILE* stream, void *obj)
 * //
 * // restores full material context (saves state variables, that completely describe
 * // current state)
 * // does not restores the slave - records
 * // this task is done at the layeredCrossSection level
 * //
 * {
 *
 * contextIOResultType iores;
 * if ((iores = this->giveMaterial()->restoreContext(stream,(void*) this)) != CIO_OK) THROW_CIOERR(iores);
 * //if (matStatusDict->restoreContext(stream,obj) != 1)
 * //  error ("restoreContext io error encountered");
 *
 * return CIO_OK;
 *
 * }
 */
} // end namespace oofem

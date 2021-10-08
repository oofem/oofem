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
}

GaussPoint :: GaussPoint(IntegrationRule * ir, int n, double w, MaterialMode mode) :
    number(n),
    irule(ir),
    weight(w),
    materialMode(mode)
{
}

GaussPoint::~GaussPoint()
{
    for ( GaussPoint *gp: gaussPoints ) {
        delete gp;
    }
}


  void GaussPoint :: printOutputAt(FILE *file, TimeStep *tStep, const char* indent)
{
    int iruleNumber = irule ? irule->giveNumber() : 0;

    fprintf(file, "%s  GP %2d.%-2d :", indent, iruleNumber, number);

    IntegrationPointStatus *status = this->giveMaterialStatus();
    if ( status ) {
        status->printOutputAt(file, tStep);
    }

    if ( gaussPoints.size() != 0 ) { // layered material
        fprintf(file, "\tLayers report {\n");
        for ( auto &gp: gaussPoints ) {
          gp->printOutputAt(file, tStep, "\t\t");
        }

        fprintf(file, "\t} end layers report\n");
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
        return nullptr;
    }

    if ( ( index < 0 ) || ( index >= (int)gaussPoints.size() ) ) {
        OOFEM_ERROR("index out of range");
    }

    return gaussPoints [ index ];
}

bool GaussPoint :: hasSlaveGaussPoint()
{
    if ( this->gaussPoints.size() != 0 ) {// layered material
        return true;
    }
    return false;
}

size_t GaussPoint :: findFirstIndexOfSlaveGaussPoint(GaussPoint *gp)
{
    auto it = find( gaussPoints.begin(), gaussPoints.end(), gp );
    if ( it != gaussPoints.end() ) {
        return it - gaussPoints.begin();
    } else {
        OOFEM_ERROR("Slave Gauss point not found.");
    }
}

void GaussPoint :: updateYourself(TimeStep *tStep)
{
    IntegrationPointStatus *status = this->giveMaterialStatus();
    if ( status ) {
        status->updateYourself(tStep);
    }

    for ( auto &gp: gaussPoints ) {
        gp->updateYourself(tStep);
    }
}

} // end namespace oofem

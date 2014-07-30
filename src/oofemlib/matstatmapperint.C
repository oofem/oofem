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

/*
 * matstatmapperint.C
 *
 *  Created on: Nov 6, 2013
 *      Author: Erik Svenning
 */

#include "matstatmapperint.h"

#include "error.h"
#include "materialmappingalgorithm.h"
#include "mmaclosestiptransfer.h"
#include "matstatus.h"

namespace oofem {
MaterialStatusMapperInterface :: MaterialStatusMapperInterface()
{
    mpMaterialMapper = new MMAClosestIPTransfer();
}

MaterialStatusMapperInterface :: ~MaterialStatusMapperInterface()
{
    if ( mpMaterialMapper != NULL ) {
        delete mpMaterialMapper;
        mpMaterialMapper = NULL;
    }
}

int MaterialStatusMapperInterface :: MSMI_map(const GaussPoint &iGP, const Domain &iOldDom, Set &sourceSet, const TimeStep &iTStep, MaterialStatus &oStatus)
{
    // Mapping of "regular" Gauss points
    int result = 1;

    Domain *dold = const_cast< Domain * >(& iOldDom);
    GaussPoint *gp = const_cast< GaussPoint * >(& iGP);
    TimeStep *tStep = const_cast< TimeStep * >(& iTStep);
    IntArray type;

    mpMaterialMapper->init(dold, type, gp, sourceSet, tStep);

    mpMaterialMapper->mapStatus(oStatus);


    return result;
}

int MaterialStatusMapperInterface :: MSMI_map_cz(const GaussPoint &iGP, const Domain &iOldDom, Set &sourceSet, const TimeStep &iTStep, MaterialStatus &oStatus)
{
    // Mapping of cohesive zone Gauss points
    int result = 1;

    Domain *dold = const_cast< Domain * >(& iOldDom);
    GaussPoint *gp = const_cast< GaussPoint * >(& iGP);
    TimeStep *tStep = const_cast< TimeStep * >(& iTStep);
    IntArray type;
    bool gpBelongsToCohesiveZone = true;

    mpMaterialMapper->init(dold, type, gp, sourceSet, tStep, gpBelongsToCohesiveZone);

    mpMaterialMapper->mapStatus(oStatus);

    return result;
}


int MaterialStatusMapperInterface :: MSMI_update(const GaussPoint &iGP, const TimeStep &iTStep)
{
    int result = 1;



    return result;
}

int MaterialStatusMapperInterface :: MSMI_finish(const TimeStep &iTStep)
{
    int result = 1;



    return result;
}
} /* namespace oofem */

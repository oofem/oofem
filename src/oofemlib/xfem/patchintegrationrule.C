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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "patchintegrationrule.h"
#include "xfemelementinterface.h"
#include "patch.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "geometry.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "datastream.h"
#include "gausspoint.h"

namespace oofem {
PatchIntegrationRule :: PatchIntegrationRule(int n, Element *e, Patch *patch) : GaussIntegrationRule(n, e)
{
    this->patch = patch;
}

PatchIntegrationRule :: ~PatchIntegrationRule()
{
    delete patch;
}

int
PatchIntegrationRule :: SetUpPointsOnTriangle(int nPoints, MaterialMode mode)
{
    numberOfIntegrationPoints = GaussIntegrationRule :: SetUpPointsOnTriangle(nPoints, mode);
    firstLocalStrainIndx = 1;
    lastLocalStrainIndx = 3;
    // convert patch coordinates into element based, update weights accordingly
    for ( int j = 0; j <  numberOfIntegrationPoints; j++ ) {
        GaussPoint *gp = this->gaussPointArray[ j ];
        patch->convertGPIntoParental(gp); // convert coordinates into parental
        Element *elg = patch->giveParent();
        double parentArea = elg->computeArea();
        Triangle *tr = ( Triangle * ) patch; ///@todo This must be a serious bug. This will do a reinterpret_cast from Patch to Triangle, which makes no sense.
        gp->setWeight(8.0 * gp->giveWeight() * tr->getArea() / parentArea); // update integration weight
    }

    return numberOfIntegrationPoints;
}

contextIOResultType
PatchIntegrationRule :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    //
    // saves full  context (saves state variables, that completely describe
    // current state)
    //

    // save parent data
    contextIOResultType iores;

    if ( ( iores = IntegrationRule :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // save patch data
    if ( this->patch ) {
        // store patch type
        int _type = this->patch->givePatchType();
        if ( !stream->write(& _type, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        patch->saveContext(stream, mode, obj);
    } else {
        OOFEM_ERROR("saveContex : can't store NULL patch");
    }

    return CIO_OK;
}

contextIOResultType
PatchIntegrationRule :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    //
    // restores full element context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;

    if ( stream == NULL ) {
        OOFEM_ERROR("restoreContex : can't write into NULL stream");
    }

    if ( ( iores = IntegrationRule :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // restore patch data
    if ( this->patch ) {
        delete this->patch;
    }

    int _ptype;
    if ( !stream->read(& _ptype, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // create new patch
    this->patch = classFactory.createPatch( ( Patch :: PatchType ) _ptype, this->giveElement() );
    this->patch->restoreContext(stream, mode, obj);

    return CIO_OK;
}
} // end namespace oofem


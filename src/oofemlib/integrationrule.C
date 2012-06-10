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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "integrationrule.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
IntegrationRule :: IntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic)
{
    number = n;
    elem = e;
    numberOfIntegrationPoints = 0;
    gaussPointArray = NULL;
    firstLocalStrainIndx = startIndx;
    lastLocalStrainIndx  = endIndx;
    isDynamic = dynamic;
    intdomain = _Unknown_integrationDomain;
}

IntegrationRule :: IntegrationRule(int n, Element *e)
{
    number = n;
    elem = e;
    numberOfIntegrationPoints = 0;
    gaussPointArray = NULL;
    firstLocalStrainIndx = lastLocalStrainIndx = 0;
    isDynamic = false;
    intdomain = _Unknown_integrationDomain;
}


IntegrationRule :: ~IntegrationRule()
{
    this->clear();
}


void
IntegrationRule :: clear()
{
    int i;
    if ( gaussPointArray ) {
        for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
            delete gaussPointArray [ i ];
        }

        delete[] gaussPointArray;
    }

    gaussPointArray = NULL;
    numberOfIntegrationPoints = 0;
}


GaussPoint *
IntegrationRule :: getIntegrationPoint(int i)
{
#  ifdef DEBUG
    if ( ( i < 0 ) || ( i >= numberOfIntegrationPoints ) ) {
        OOFEM_ERROR2("IntegrationRule::getIntegrationPoint - request out of bounds (%d)", i);
    }

#  endif
    return gaussPointArray [ i ];
}

void
IntegrationRule :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

    for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
        gaussPointArray [ i ]->printOutputAt(file, stepN);
    }
}

void
IntegrationRule :: updateYourself(TimeStep *tStep)
{
    // Updates the receiver at end of step.
    int i;

    for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
        gaussPointArray [ i ]->updateYourself(tStep);
    }
}


void
IntegrationRule :: initForNewStep()
{
    // initializes receiver to new time step or can be used
    // if current time step must be restarted
    //
    // call material->initGpForNewStep() for all GPs.
    //

    int i;

    for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
        gaussPointArray [ i ]->giveMaterial()->initGpForNewStep(gaussPointArray [ i ]);
    }
}


contextIOResultType
IntegrationRule :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    //
    // saves full  context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;
    int i;
    GaussPoint *gp;

    if ( stream == NULL ) {
        OOFEM_ERROR("saveContex : can't write into NULL stream");
    }

    int isdyn = isDynamic;
    if ( !stream->write(& isdyn, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( isDynamic ) {
        mode |= CM_Definition;          // store definition if dynamic
    }

    if ( mode & CM_Definition ) {
        if ( !stream->write(& numberOfIntegrationPoints, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // write first and last integration indices
        if ( !stream->write(& firstLocalStrainIndx, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& lastLocalStrainIndx, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
        gp = gaussPointArray [ i ];
        if ( mode & CM_Definition ) {
            // write gp weight, coordinates, element number, and material mode
            double dval = gp->giveWeight();
            if ( !stream->write(& dval, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            if ( ( iores = gp->giveCoordinates()->storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }

            //int ival = gp->giveElement()->giveNumber();
            //if (!stream->write(&ival,1)) THROW_CIOERR(CIO_IOERR);
            int mmode = gp->giveMaterialMode();
            if ( !stream->write(& mmode, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        // write gp data
        if ( ( iores = gp->giveCrossSection()->saveIPContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

contextIOResultType
IntegrationRule :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    //
    // restores full element context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;
    int i, size;
    bool __create = false;
    GaussPoint *gp;
    //Element* __parelem = (Element*) obj;

    if ( stream == NULL ) {
        OOFEM_ERROR("restoreContex : can't write into NULL stream");
    }

    int isdyn;
    if ( !stream->read(& isdyn, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    isDynamic = ( bool ) isdyn;

    if ( isDynamic ) {
        mode |= CM_Definition;          // store definition if dynamic
    }

    if ( mode & CM_Definition ) {
        if ( !stream->read(& size, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // read first and last integration indices
        if ( !stream->read(& firstLocalStrainIndx, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& lastLocalStrainIndx, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( numberOfIntegrationPoints != size ) {
            this->clear();
            __create = true;
            gaussPointArray = new GaussPoint * [ size ];
            numberOfIntegrationPoints = size;
        }
    }

    for ( i = 0; i < numberOfIntegrationPoints; i++ ) {
        if ( mode & CM_Definition ) {
            // read weight
            double w;
            if ( !stream->read(& w, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            // read coords
            FloatArray c;
            if ( ( iores = c.restoreYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }

            // restore element and material mode
            //int n;
            //if (!stream->read(&n,1)) THROW_CIOERR(CIO_IOERR);
            MaterialMode m;
            int _m;
            if ( !stream->read(& _m, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            m = ( MaterialMode ) _m;
            // read dynamic flag

            if ( __create ) {
                gaussPointArray [ i ] = new GaussPoint(this, i + 1, (new FloatArray(c)), w, m);
            } else {
                gp = gaussPointArray [ i ];
                gp->setWeight(w);
                gp->setCoordinates(c);
                //gp->setElement (__parelem);
                gp->setMaterialMode(m);
            }
        }

        // read gp data
        gp = gaussPointArray [ i ];
        if ( ( iores = gp->giveCrossSection()->restoreIPContext(stream, mode, gp) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


int
IntegrationRule :: setUpIntegrationPoints(integrationDomain mode, int nPoints,
                                          MaterialMode matMode)
{
    intdomain = mode;

    switch ( mode ) {
    case _Line:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOnLine(nPoints, matMode, & gaussPointArray) );

    case _Triangle:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOnTriangle(nPoints, matMode, & gaussPointArray) );

    case _Square:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOnSquare(nPoints, matMode, & gaussPointArray) );

    case _Cube:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOnCube(nPoints, matMode, & gaussPointArray) );

    case _Tetrahedra:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOnTetrahedra(nPoints, matMode, & gaussPointArray) );

    default:
        OOFEM_ERROR("IntegrationRule::setUpIntegrationPoints - unknown mode");
    }

    return 0;
}

int
IntegrationRule :: setUpEmbeddedIntegrationPoints(integrationDomain mode, int nPoints, MaterialMode matMode,
                                                  const FloatArray **coords)
{
    intdomain = mode;

    switch ( mode ) {
    case _Embedded2dLine:
        return  ( numberOfIntegrationPoints = this->SetUpPointsOn2DEmbeddedLine(nPoints, matMode, & gaussPointArray, coords) );

    default:
        OOFEM_ERROR("IntegrationRule::setUpEmbeddedIntegrationPoints - unknown mode");
    }

    return 0;
}
} // end namespace oofem

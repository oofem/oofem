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

#include "integrationrule.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

IntegrationRule :: IntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic)
{
    number = n;
    elem = e;
    firstLocalStrainIndx = startIndx;
    lastLocalStrainIndx  = endIndx;
    isDynamic = dynamic;
    intdomain = _Unknown_integrationDomain;
}

IntegrationRule :: IntegrationRule(int n, Element *e)
{
    number = n;
    elem = e;
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
    for ( GaussPoint *gp: *this ) {
        delete gp;
    }

    gaussPoints.clear();
}


GaussPoint *
IntegrationRule :: getIntegrationPoint(int i)
{
#  ifdef DEBUG
    if ( ( i < 0 ) || ( i >= this->giveNumberOfIntegrationPoints() ) ) {
        OOFEM_ERROR("request out of bounds (%d)", i);
    }

#  endif
    return gaussPoints [ i ];
}


GaussPoint *
IntegrationRule :: findIntegrationPointClosestTo(const FloatArray &lcoord)
{
    double mindist = -1.;
    GaussPoint *minGp = NULL;
    for ( GaussPoint *gp: *this ) {
        double dist = lcoord.distance_square(gp->giveNaturalCoordinates());
        if ( dist <= mindist || mindist < 0. ) {
            mindist = dist;
            minGp = gp;
        }
    }
    return minGp;
}


void
IntegrationRule :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    for ( GaussPoint *gp: *this ) {
        gp->printOutputAt(file, tStep);
    }
}

void
IntegrationRule :: updateYourself(TimeStep *tStep)
{
    // Updates the receiver at end of step.
    for ( GaussPoint *gp: *this ) {
        gp->updateYourself(tStep);
    }
}


contextIOResultType
IntegrationRule :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    //
    // saves full  context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;

    int isdyn = isDynamic;
    if ( !stream.write(isdyn) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( isDynamic ) {
        mode |= CM_Definition;          // store definition if dynamic
    }

    if ( mode & CM_Definition ) {
        int numberOfIntegrationPoints = (int)this->gaussPoints.size();
        if ( !stream.write(numberOfIntegrationPoints) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // write first and last integration indices
        if ( !stream.write(firstLocalStrainIndx) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(lastLocalStrainIndx) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    for ( GaussPoint *gp: *this ) {
        if ( mode & CM_Definition ) {
            // write gp weight, coordinates, element number, and material mode
            double dval = gp->giveWeight();
            if ( !stream.write(dval) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            if ( ( iores = gp->giveNaturalCoordinates().storeYourself(stream) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }

            //int ival = gp->giveElement()->giveNumber();
            //if (!stream.write(&ival,1)) THROW_CIOERR(CIO_IOERR);
            int mmode = gp->giveMaterialMode();
            if ( !stream.write(mmode) ) {
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
IntegrationRule :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    //
    // restores full element context (saves state variables, that completely describe
    // current state)
    //

    contextIOResultType iores;
    int size;

    int isdyn;
    if ( !stream.read(isdyn) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    isDynamic = ( bool ) isdyn;

    if ( isDynamic ) {
        mode |= CM_Definition;          // store definition if dynamic
    }

    if ( mode & CM_Definition ) {
        if ( !stream.read(size) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // read first and last integration indices
        if ( !stream.read(firstLocalStrainIndx) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(lastLocalStrainIndx) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        this->clear();
        
        this->gaussPoints.resize(size);
    }

    int i = 1;
    for ( GaussPoint *&gp: *this ) {
        if ( mode & CM_Definition ) {
            // read weight
            double w;
            if ( !stream.read(w) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            // read coords
            FloatArray c;
            if ( ( iores = c.restoreYourself(stream) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }

            // restore element and material mode
            //int n;
            //if (!stream.read(n)) THROW_CIOERR(CIO_IOERR);
            MaterialMode m;
            int _m;
            if ( !stream.read(_m) ) {
                THROW_CIOERR(CIO_IOERR);
            }

            m = ( MaterialMode ) _m;
            // read dynamic flag

            gp = new GaussPoint(this, i, std :: move(c), w, m);
            i++;
        }

        // read gp data
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
        return  this->SetUpPointsOnLine(nPoints, matMode);

    case _Triangle:
        return  this->SetUpPointsOnTriangle(nPoints, matMode);

    case _Square:
        return  this->SetUpPointsOnSquare(nPoints, matMode);

    case _Cube:
        return  this->SetUpPointsOnCube(nPoints, matMode);

    case _Tetrahedra:
        return  this->SetUpPointsOnTetrahedra(nPoints, matMode);

    case _Wedge:
        // Limited wrapper for now;
        if ( nPoints == 6 ) {
            return this->SetUpPointsOnWedge(3, 2, matMode);
        } else {
            return this->SetUpPointsOnWedge(3, 3, matMode);
        }

    default:
        OOFEM_ERROR("unknown mode (%d)", mode);
    }

    return 0;
}

int
IntegrationRule :: setUpEmbeddedIntegrationPoints(integrationDomain mode, int nPoints, MaterialMode matMode,
                                                  const std :: vector< FloatArray > &coords)
{
    intdomain = mode;

    switch ( mode ) {
    case _Embedded2dLine:
        if ( coords.size() != 2 ) {
            OOFEM_ERROR("Exactly 2 coordinates are required for 2D embedded lines!");
        }
        return  this->SetUpPointsOn2DEmbeddedLine(nPoints, matMode, coords[0], coords[1]);

    default:
        OOFEM_ERROR("unknown mode");
    }

    return 0;
}


int IntegrationRule :: SetUpPoint(MaterialMode mode)
{
    this->gaussPoints.resize(1);
    this->gaussPoints [ 0 ] = new GaussPoint(this, 1, FloatArray(0), 1.0, mode);
    this->intdomain = _Point;
    return 1;
}
} // end namespace oofem

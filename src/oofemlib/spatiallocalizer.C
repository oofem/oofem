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

#include "spatiallocalizer.h"
#include "connectivitytable.h"
#include "element.h"
#include "node.h"
#include "mathfem.h"
#include "error.h"
#include "floatarray.h"
#include "intarray.h"
#include "feinterpol.h"

namespace oofem {
#define POINT_TOL 1.e-6

void
SpatialLocalizerInterface :: SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1)
{
    Element *element = this->SpatialLocalizerI_giveElement();
    FloatArray *coordinates;

    coordinates = element->giveNode(1)->giveCoordinates();
    bb1 = bb0 = *element->giveNode(1)->giveCoordinates();

    for ( int i = 2; i <= element->giveNumberOfNodes(); ++i ) {
        coordinates = element->giveNode(i)->giveCoordinates();
        bb0.beMinOf(bb0, *coordinates);
        bb1.beMaxOf(bb1, *coordinates);
    }
}


double
SpatialLocalizerInterface :: SpatialLocalizerI_giveClosestPoint(FloatArray &lcoords, FloatArray &closest, const FloatArray &gcoords)
{
    Element *e = this->SpatialLocalizerI_giveElement();
    FEInterpolation *interp = e->giveInterpolation();
    
    if ( !interp->global2local(lcoords, gcoords, FEIElementGeometryWrapper(e)) ) { // Outside element
        interp->local2global(closest, lcoords, FEIElementGeometryWrapper(e));
        return closest.distance(gcoords);
    } else {
        closest = gcoords;
        return 0.0;
    }
}


int
SpatialLocalizerInterface :: SpatialLocalizerI_BBoxContainsPoint(const FloatArray &coords)
{
  FloatArray coordMin, coordMax;
  this->SpatialLocalizerI_giveBBox(coordMin, coordMax);

   int size = min( coordMin.giveSize(), coords.giveSize() );
    for ( int j = 1; j <= size; j++ ) {
      if ( coords.at(j) < coordMin.at(j) - POINT_TOL ) {
            return 0;
        }

      if ( coords.at(j) > coordMax.at(j) + POINT_TOL ) {
            return 0;
        }
    }

    return 1;
}



void
SpatialLocalizer :: giveAllElementsWithNodesWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                      const double radius)
{
    nodeContainerType nodesWithinBox;
    nodeContainerType :: iterator it;
    const IntArray *dofmanConnectivity;

    elemSet.clear();

    ConnectivityTable *ct = domain->giveConnectivityTable();

    this->giveAllNodesWithinBox(nodesWithinBox, coords, radius);

    for ( it = nodesWithinBox.begin(); it != nodesWithinBox.end(); ++it ) {
        dofmanConnectivity = ct->giveDofManConnectivityArray(* it);
        for ( int i = 1; i <= dofmanConnectivity->giveSize(); i++ ) {
            elemSet.insert( dofmanConnectivity->at(i) );
        }
    }
}
} // end namespace oofem

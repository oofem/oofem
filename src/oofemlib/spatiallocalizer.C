/* $Header: /home/cvs/bp/oofem/oofemlib/src/spatiallocalizer.C,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "conTable.h"
#include "element.h"
#include "node.h"
#include "mathfem.h"
#include "error.h"


#define POINT_TOL 1.e-6

int
SpatialLocalizerInterface :: SpatialLocalizerI_BBoxContainsPoint(const FloatArray &coords)
{
    Element *element = this->SpatialLocalizerI_giveElement();
    FloatArray *coordinates;
    double coordMin [ 3 ], coordMax [ 3 ], val;
    int i, j, size;

    coordinates = element->giveNode(1)->giveCoordinates();
    size = coordinates->giveSize();
    for ( j = 1; j <= size; j++ ) {
        val = coordinates->at(j);
        coordMin [ j - 1 ] = coordMax [ j - 1 ] = val;
    }

    for ( i = 2; i <= element->giveNumberOfNodes(); i++ ) {
        coordinates = element->giveNode(i)->giveCoordinates();
        if ( coordinates->giveSize() != size ) {
            OOFEM_ERROR("SpatialLocalizerInterface::SpatialLocalizerI_BBoxContainsPoint: coordinates size mismatch");
        }

        for ( j = 1; j <= size; j++ ) {
            val = coordinates->at(j);
            if ( val < coordMin [ j - 1 ] ) {
                coordMin [ j - 1 ] = val;
                continue;
            }

            if ( val > coordMax [ j - 1 ] ) {
                coordMax [ j - 1 ] = val;
                continue;
            }
        }
    }

    /*
     * if(coords.giveSize() < size){
     * fprintf(stderr, "SpatialLocalizerInterface::SpatialLocalizerI_BBoxContainsPoint: coordinates size mismatch");
     * exit(1);
     * }
     */
    size = min( size, coords.giveSize() );
    for ( j = 1; j <= size; j++ ) {
        if ( coords.at(j) < coordMin [ j - 1 ] - POINT_TOL ) {
            return 0;
        }

        if ( coords.at(j) > coordMax [ j - 1 ] + POINT_TOL ) {
            return 0;
        }
    }

    return 1;
}



void 
SpatialLocalizer::giveAllElementsWithNodesWithinBox(elementContainerType &elemSet, const FloatArray &coords,
                                                    const double radius)
{
  nodeContainerType nodesWithinBox;
  nodeContainerType::iterator it;
  const IntArray* dofmanConnectivity;
  int i;

  elemSet.clear();

  ConnectivityTable* ct = domain->giveConnectivityTable();

  this->giveAllNodesWithinBox (nodesWithinBox, coords, radius);
  
  for (it=nodesWithinBox.begin(); it != nodesWithinBox.end(); ++it) {
   dofmanConnectivity = ct->giveDofManConnectivityArray (*it);
   for (i=1; i<=dofmanConnectivity->giveSize(); i++) elemSet.insert(dofmanConnectivity->at(i));
  }
}

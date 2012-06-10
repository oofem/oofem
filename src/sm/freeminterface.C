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

#include "freeminterface.h"
#include "errorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "dynalist.h"
#include "conTable.h"
#include "mathfem.h"

namespace oofem {
MesherInterface :: returnCode
FreemInterface :: createMesh(TimeStep *stepN, int domainNumber, int domainSerNum, Domain **dNew)
{
    * dNew = NULL;
    if ( this->createInput(this->domain, stepN) ) {
        return MI_NEEDS_EXTERNAL_ACTION;
    } else {
        return MI_FAILED;
    }
}

int
FreemInterface :: createInput(Domain *d, TimeStep *stepN)
{
    int i;
    int nnodes = d->giveNumberOfDofManagers(), nelem = d->giveNumberOfElements();
    double density;
    FILE *outputStrem;
    Node *inode;
    Element *ielem;

    outputStrem = fopen("freem.bmf", "w");
    // print header for 2D
    fprintf(outputStrem, "nbnodes %d nbelem %d \n", nnodes, nelem);

    /* mesh densities smoothing */
    // query nodal absolute densities
    FloatArray nodalDensities(nnodes);
    for ( i = 1; i <= nnodes; i++ ) {
        nodalDensities.at(i) = d->giveErrorEstimator()->giveRemeshingCrit()->giveRequiredDofManDensity(i, stepN);
    }

    this->smoothNodalDensities(d,  nodalDensities, stepN);
    /* end of smoothing */

    // loop over nodes
    for ( i = 1; i <= nnodes; i++ ) {
        //density = d->giveErrorEstimator ()->giveRemeshingCrit()->giveRequiredDofManDensity (i, stepN, 1);
        //density = d->giveErrorEstimator ()->giveRemeshingCrit()->giveDofManDensity (i) / nodalDensities.at(i);
        density = nodalDensities.at(i);
        inode = d->giveNode(i);
        fprintf(outputStrem, "backgroungMeshNode %d x %e y %e density %e\n", i, inode->giveCoordinate(1), inode->giveCoordinate(2), density);
    }

    for ( i = 1; i <= nelem; i++ ) {
        ielem = d->giveElement(i);
        if ( ielem->giveClassID() != PlaneStress2dClass ) {
            OOFEM_ERROR("FreemInterface::createInput : unsupported element type");
        }

        fprintf( outputStrem, "backgroundMeshElem %d  nodes 4 %d %d %d %d\n", i,
                ielem->giveNode(1)->giveNumber(), ielem->giveNode(2)->giveNumber(),
                ielem->giveNode(3)->giveNumber(), ielem->giveNode(4)->giveNumber() );
    }

    fclose(outputStrem);

    OOFEM_LOG_INFO("freem.bmf file created\n");
    return 1;
}

void
FreemInterface :: smoothNodalDensities(Domain *d,  FloatArray &nodalDensities, TimeStep *stepN)
{
    int i, j, k, neighbour, candidate, found, jelemNodes;
    int nnodes = d->giveNumberOfDofManagers();
    double dist;
    const IntArray *candidateConnectivity;
    FloatArray *neighbourCoords;
    Element *jelem;
    Node *candNode;
    dynaList< int >queue;
    dynaList< int > :: iterator pos;


    // loop over nodes
    for ( i = 1; i <= nnodes; i++ ) {
        if ( !( ( d->giveDofManager(i)->giveClassID() == NodeClass ) || ( d->giveDofManager(i)->giveClassID() == RigidArmNodeClass ) ) ) {
            continue;
        }

        queue.clear();
        queue.pushFront(i);

        while ( !queue.isEmpty() ) {
            // extract candidate
            candidate = * ( queue.begin() );
            queue.erase( queue.begin() );

            candNode  = ( Node * ) d->giveDofManager(candidate);
            // find candidate neighbours
            candidateConnectivity = d->giveConnectivityTable()->giveDofManConnectivityArray(candidate);
            for ( j = 1; j <= candidateConnectivity->giveSize(); j++ ) {
                jelem = d->giveElement( candidateConnectivity->at(j) );
                jelemNodes = jelem->giveNumberOfNodes();
                for ( k = 1; k <= jelemNodes; k++ ) {
                    neighbour = jelem->giveNode(k)->giveNumber();
                    if ( neighbour == candidate ) {
                        continue;
                    }

                    // neighbour found, check if smoothing necessary
                    neighbourCoords = ( ( Node * ) jelem->giveNode(k) )->giveCoordinates();
                    dist = candNode->giveCoordinates()->distance(neighbourCoords);
                    // overshoot criteria
                    if ( ( ( nodalDensities.at(neighbour) / nodalDensities.at(candidate) ) > 1.3 ) &&
                        ( nodalDensities.at(neighbour) > 1.0 * dist ) ) {
                        // increase candidate nodal density
                        nodalDensities.at(neighbour) = max( 1.0 * dist, nodalDensities.at(candidate) );
                        // printf ("o");
                        // put candidate into queue if not yet added present
                        for ( found = 0, pos = queue.begin(); pos != queue.end(); ++pos ) {
                            if ( * pos == neighbour ) {
                                found = 1;
                                break;
                            }
                        }

                        if ( !found ) {
                            queue.pushFront(neighbour);
                        }

                        // end overshoot criteria
                    } else if ( ( nodalDensities.at(neighbour) - nodalDensities.at(candidate) ) / dist  > 2.5 ) {
                        // max gradient criteria
                        nodalDensities.at(neighbour) = nodalDensities.at(candidate) + 2.2 * dist;
                        //printf ("g");
                        // put candidate into queue if not yet added present
                        for ( found = 0, pos = queue.begin(); pos != queue.end(); ++pos ) {
                            if ( * pos == neighbour ) {
                                found = 1;
                                break;
                            }
                        }

                        if ( !found ) {
                            queue.pushFront(neighbour);
                        }
                    }
                }
            }
        } // end loop over queue

    }
}
} // end namespace oofem

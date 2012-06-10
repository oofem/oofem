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

#include "t3dinterface.h"
#include "errorestimator.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "engngm.h"

namespace oofem {
MesherInterface :: returnCode
T3DInterface :: createMesh(TimeStep *stepN, int domainNumber, int domainSerNum, Domain **dNew)
{
    * dNew = NULL;
    if ( this->createInput(this->domain, stepN) ) {
        return MI_NEEDS_EXTERNAL_ACTION;
    } else {
        return MI_FAILED;
    }
}

int
T3DInterface :: createInput(Domain *d, TimeStep *stepN)
{
    int i, j;
    int nnodes = d->giveNumberOfDofManagers(), nelem = d->giveNumberOfElements();
    double density;
    FILE *outputStrem;
    Node *inode;
    Element *ielem;
    int edges, trias, quads, tetras, pyrams, wedges, hexas;
    IntArray edgeIdArray, triaIdArray, quadIdArray, tetraIdArray, pyramIdArray, wedgeIdArray, hexaIdArray;
    bool tri_tetra = false;
    char fileName [ 32 ];

    edges = trias = quads = tetras = pyrams = wedges = hexas = 0;
    for ( i = 1; i <= nelem; i++ ) {
        ielem = d->giveElement(i);
        switch ( ielem->giveGeometryType() ) {
        case EGT_point:
            break;
        case EGT_line_1:
        case EGT_line_2:
            edges++;
            break;
        case EGT_triangle_1:
        case EGT_triangle_2:
            trias++;
            break;
        case EGT_quad_1:
        case EGT_quad_2:
            quads++;
            break;
        case EGT_tetra_1:
        case EGT_tetra_2:
            tetras++;
            break;
        case EGT_hexa_1:
        case EGT_hexa_2:
            hexas++;
            break;
        case EGT_unknown:
        case EGT_Composite:
            OOFEM_ERROR2( "T3DInterface::createInput unknown element type (%s)",
                         __Element_Geometry_TypeToString( ielem->giveGeometryType() ) );
        }
    }

    edgeIdArray.resize(edges);
    triaIdArray.resize(trias);
    quadIdArray.resize(quads);
    tetraIdArray.resize(tetras);
    pyramIdArray.resize(pyrams);
    wedgeIdArray.resize(wedges);
    hexaIdArray.resize(hexas);

    edges = trias = quads = tetras = pyrams = wedges = hexas = 0;
    for ( i = 1; i <= nelem; i++ ) {
        ielem = d->giveElement(i);
        switch ( ielem->giveGeometryType() ) {
        case EGT_point:
            break;
        case EGT_line_1:
        case EGT_line_2:
            edgeIdArray.at(++edges) = i;
            break;
        case EGT_triangle_1:
        case EGT_triangle_2:
            triaIdArray.at(++trias) = i;
            break;
        case EGT_quad_1:
        case EGT_quad_2:
            quadIdArray.at(++quads) = i;
            break;
        case EGT_tetra_1:
        case EGT_tetra_2:
            tetraIdArray.at(++tetras) = i;
            break;
        case EGT_hexa_1:
        case EGT_hexa_2:
            hexaIdArray.at(++hexas) = i;
            break;
        case EGT_unknown:
        case EGT_Composite:
            break;
        }
    }

    if ( quads + hexas + pyrams + wedges == 0 ) {
        tri_tetra = true;
    }

#ifdef __PARALLEL_MODE
    sprintf( fileName, "%s.%d", BMF_FILENAME, d->giveEngngModel()->
            giveRank() );
#else
    sprintf(fileName, "%s", BMF_FILENAME);
#endif

    outputStrem = fopen(fileName, "w");
    if ( tri_tetra == true ) {
        fprintf(outputStrem, "3 1\n");
        fprintf(outputStrem, "%d %d %d %d\n", nnodes, edges, trias, tetras);
    } else {
        fprintf(outputStrem, "7 1\n");
        fprintf(outputStrem, "%d %d %d %d %d %d %d %d\n", nnodes, edges, trias, quads, tetras, pyrams, wedges, hexas);
    }

    // loop over nodes
    for ( i = 1; i <= nnodes; i++ ) {
        density = d->giveErrorEstimator()->giveRemeshingCrit()->giveRequiredDofManDensity(i, stepN);
        inode = d->giveNode(i);
        fprintf(outputStrem, "%d %e %e %e  %e\n", i, inode->giveCoordinate(1), inode->giveCoordinate(2), inode->giveCoordinate(3), density);
    }

    // loop separately for each type of element
    // since T3ds support only linear elements in bg mesh, only corner nodes are stored
    // alternatively quadratic elements may be splitted to linear (not supported now)

    if ( edges != 0 ) {
        for ( i = 1; i <= edges; i++ ) {
            ielem = d->giveElement( edgeIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 2; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( trias != 0 ) {
        for ( i = 1; i <= trias; i++ ) {
            ielem = d->giveElement( triaIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 3; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( quads != 0 ) {
        for ( i = 1; i <= quads; i++ ) {
            ielem = d->giveElement( quadIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 4; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( tetras != 0 ) {
        for ( i = 1; i <= tetras; i++ ) {
            ielem = d->giveElement( tetraIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 4; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( pyrams != 0 ) {
        for ( i = 1; i <= pyrams; i++ ) {
            ielem = d->giveElement( pyramIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 5; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( wedges != 0 ) {
        for ( i = 1; i <= wedges; i++ ) {
            ielem = d->giveElement( wedgeIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 6; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    if ( hexas != 0 ) {
        for ( i = 1; i <= hexas; i++ ) {
            ielem = d->giveElement( hexaIdArray.at(i) );
            fprintf(outputStrem, "%d", i);
            for ( j = 1; j <= 8; j++ ) {
                fprintf( outputStrem, " %d", ielem->giveNode(j)->giveNumber() );
            }

            fprintf(outputStrem, "\n");
        }
    }

    fclose(outputStrem);

    OOFEM_LOG_INFO("t3d.bmf file created\n");
    return 1;
}
} // end namespace oofem

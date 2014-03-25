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

#include "eleminterpunknownmapper.h"
#include "eleminterpmapperinterface.h"
#include "element.h"
#include "domain.h"
#include "engngm.h"
#include "spatiallocalizer.h"
#include "node.h"
#include "dof.h"
#include "connectivitytable.h"

namespace oofem {
EIPrimaryUnknownMapper :: EIPrimaryUnknownMapper() : PrimaryUnknownMapper()
{ }

#define OOFEM_MAPPING_CHECK_REGIONS

int
EIPrimaryUnknownMapper :: mapAndUpdate(FloatArray &answer, ValueModeType mode,
                                       Domain *oldd, Domain *newd,  TimeStep *tStep)
{
    int inode, nd_nnodes = newd->giveNumberOfDofManagers();
    int nsize = newd->giveEngngModel()->giveNumberOfDomainEquations( newd->giveNumber(), EModelDefaultEquationNumbering() );
    FloatArray unknownValues;
    IntArray dofMask, locationArray;
    IntArray reglist;
#ifdef OOFEM_MAPPING_CHECK_REGIONS
    ConnectivityTable *conTable = newd->giveConnectivityTable();
    const IntArray *nodeConnectivity;
#endif

    answer.resize(nsize);
    answer.zero();

    for ( inode = 1; inode <= nd_nnodes; inode++ ) {
        /* HUHU CHEATING */
#ifdef __PARALLEL_MODE
        if ( ( newd->giveNode(inode)->giveParallelMode() == DofManager_null ) ||
            ( newd->giveNode(inode)->giveParallelMode() == DofManager_remote ) ) {
            continue;
        }

#endif

#ifdef OOFEM_MAPPING_CHECK_REGIONS
        // build up region list for node
        nodeConnectivity = conTable->giveDofManConnectivityArray(inode);
        reglist.resize( nodeConnectivity->giveSize() );
        reglist.clear();
        for ( int indx = 1; indx <= nodeConnectivity->giveSize(); indx++ ) {
            reglist.insertSortedOnce( newd->giveElementGeometry( nodeConnectivity->at(indx) )->giveRegionNumber() );
        }

#endif

        if ( this->evaluateAt(unknownValues, dofMask, mode, oldd, * newd->giveNode(inode)->giveCoordinates(), reglist, tStep) ) {
            //
            //  WARNING !! LIMITED IMPLEMENTATION HERE !!
            //
            // possible source of error -> general service allowing to request all DOFs interpolated by element
            // should be there, but newNode can accommodate only certain dofs.
            //
            //
            newd->giveNode(inode)->giveLocationArray( dofMask, locationArray, EModelDefaultEquationNumbering() );
            if ( newd->giveNode(inode)->hasAnySlaveDofs() ) {
                for ( int ii = 1; ii <= dofMask.giveSize(); ii++ ) { ///@todo How should be deal with slave dofs and such? dofMask.size() != locationArray.size() will happen
                    // exclude slaves; they are determined from masters
                    if ( newd->giveNode(inode)->giveDof(ii)->isPrimaryDof() ) {
                        answer.at( locationArray.at(ii) ) += unknownValues.at(ii);
                    }
                }
            } else {
                // assemble the interpolated values to global vector using locationArray of new node
                answer.assemble(unknownValues, locationArray);
            }
        } else {
            OOFEM_ERROR2("EIPrimaryUnknownMapper :: mapAndUpdate - evaluateAt service failed for node %d", inode);
        }
    }

    return 1;
}


int
EIPrimaryUnknownMapper :: evaluateAt(FloatArray &answer, IntArray &dofMask, ValueModeType mode,
                                     Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep)
{
    ElementGeometry *oelemGeometry;
    EIPrimaryUnknownMapperInterface *interface;
    SpatialLocalizer *sl = oldd->giveSpatialLocalizer();

    ///@todo Change to the other version after checking that it works properly. Will render "giveElementCloseToPoint" obsolete (superseeded by giveElementClosestToPoint).
#if 1
    if ( regList.isEmpty() ) {
        oelemGeometry = sl->giveElementContainingPoint(coords);
    } else {
        oelemGeometry = sl->giveElementContainingPoint(coords, & regList);
    }
    if ( !oelemGeometry ) {
        if ( regList.isEmpty() ) {
            oelemGeometry = oldd->giveSpatialLocalizer()->giveElementCloseToPoint(coords);
        } else {
            oelemGeometry = oldd->giveSpatialLocalizer()->giveElementCloseToPoint(coords, & regList);
        }
        if ( !oelemGeometry ) {
            OOFEM_WARNING("EIPrimaryUnknownMapper :: evaluateAt - Couldn't find any element containing point.");
            return false;
        }
    }
#else
    FloatArray lcoords, closest;
    if ( regList.isEmpty() ) {
        oelem = sl->giveElementClosestToPoint(lcoords, closest, coords, 0);
    } else {
        // Take the minimum of any region
        double mindist = 0.0, distance;
        oelem = NULL;
        for ( int i = 1; i < regList.giveSize(); ++i ) {
            Element *tmpelem = sl->giveElementClosestToPoint( lcoords, closest, coords, regList.at(i) );
            distance = closest.distance_square(coords);
            if ( tmpelem != NULL ) {
                distance = closest.distance_square(coords);
                if ( distance < mindist || i == 1 ) {
                    mindist = distance;
                    oelem = tmpelem;
                    if ( distance == 0.0 ) {
                        break;
                    }
                }
            }
        }
    }
    if ( !oelemGeometry ) {
        OOFEM_WARNING("EIPrimaryUnknownMapper :: evaluateAt - Couldn't find any element containing point.");
        return false;
    }
#endif

    interface = static_cast< EIPrimaryUnknownMapperInterface * >( oelemGeometry->giveInterface(EIPrimaryUnknownMapperInterfaceType) );
    if ( interface ) {
        interface->EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(dofMask);
#if 1
        interface->EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(mode, tStep, coords, answer);
#else
        interface->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
#endif
    } else {
        OOFEM_ERROR("EIPrimaryUnknownMapper :: evaluateAt - Element does not support EIPrimaryUnknownMapperInterface");
    }

    return true;
}
} // end namespace oofem

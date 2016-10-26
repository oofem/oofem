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
#include "element.h"
#include "domain.h"
#include "engngm.h"
#include "spatiallocalizer.h"
#include "node.h"
#include "dof.h"
#include "connectivitytable.h"
#include "unknownnumberingscheme.h"

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
    IntArray dofidMask;
    IntArray reglist;
#ifdef OOFEM_MAPPING_CHECK_REGIONS
    ConnectivityTable *conTable = newd->giveConnectivityTable();
    const IntArray *nodeConnectivity;
#endif

    answer.resize(nsize);
    answer.zero();

    for ( inode = 1; inode <= nd_nnodes; inode++ ) {
        DofManager *node = newd->giveNode(inode);
        /* Process local and shared nodes only */
        if ( (node->giveParallelMode() != DofManager_local ) &&
             (node->giveParallelMode() != DofManager_shared)) {
            continue;
        }

#ifdef OOFEM_MAPPING_CHECK_REGIONS
        // build up region list for node
        nodeConnectivity = conTable->giveDofManConnectivityArray(inode);
        reglist.resize( nodeConnectivity->giveSize() );
        reglist.clear();
        for ( int indx = 1; indx <= nodeConnectivity->giveSize(); indx++ ) {
            reglist.insertSortedOnce( newd->giveElement( nodeConnectivity->at(indx) )->giveRegionNumber() );
        }

#endif
        ///@todo Shouldn't we pass a primary field or something to this function?
        if ( this->evaluateAt(unknownValues, dofidMask, mode, oldd, * node->giveCoordinates(), reglist, tStep) ) {
            ///@todo This doesn't respect local coordinate systems in nodes. Supporting that would require major reworking.
            for ( int ii = 1; ii <= dofidMask.giveSize(); ii++ ) {
                // exclude slaves; they are determined from masters
                auto it = node->findDofWithDofId((DofIDItem)dofidMask.at(ii));
                if ( it != node->end() ) {
                    Dof *dof = *it;
                    if ( dof->isPrimaryDof() ) {
                        int eq = dof->giveEquationNumber(EModelDefaultEquationNumbering());
                        if (eq)
                          answer.at( eq ) += unknownValues.at(ii);
                    }
                }
            }
        } else {
            OOFEM_ERROR("evaluateAt service failed for node %d", inode);
        }
    }

    return 1;
}


int
EIPrimaryUnknownMapper :: evaluateAt(FloatArray &answer, IntArray &dofMask, ValueModeType mode,
                                     Domain *oldd, FloatArray &coords, IntArray &regList, TimeStep *tStep)
{
    Element *oelem;
    SpatialLocalizer *sl = oldd->giveSpatialLocalizer();

    FloatArray lcoords, closest;
    if ( regList.isEmpty() ) {
        oelem = sl->giveElementClosestToPoint(lcoords, closest, coords, 0);
    } else {
        // Take the minimum of any region
        double mindist = 0.0, distance;
        oelem = NULL;
        for ( int i = 1; i <= regList.giveSize(); ++i ) {
            Element *tmpelem = sl->giveElementClosestToPoint( lcoords, closest, coords, regList.at(i) );
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
    if ( !oelem ) {
        OOFEM_WARNING("Couldn't find any element containing point.");
        return false;
    }

    oelem->giveElementDofIDMask(dofMask);
    oelem->computeField(mode, tStep, lcoords, answer);

    return true;
}
} // end namespace oofem

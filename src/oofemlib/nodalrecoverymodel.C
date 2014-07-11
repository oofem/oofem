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

#include "nodalrecoverymodel.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
#endif


namespace oofem {
NodalRecoveryModel :: NodalRecoveryModel(Domain *d) : nodalValList()
{
    stateCounter = 0;
    domain = d;
    this->valType = IST_Undefined;

#ifdef __PARALLEL_MODE
    communicator = NULL;
    commBuff = NULL;
    initCommMap = true;
#endif
}


NodalRecoveryModel :: ~NodalRecoveryModel()
{
#ifdef __PARALLEL_MODE
    delete communicator;
    delete commBuff;
#endif
}


int
NodalRecoveryModel :: clear()
{
    this->nodalValList.clear();
    return 1;
}

int
NodalRecoveryModel :: giveNodalVector(const FloatArray * &answer, int node)
{
    std :: map< int, FloatArray > :: iterator it = this->nodalValList.find(node);
    if ( it != this->nodalValList.end() ) {
        answer = & it->second;
        if ( answer->giveSize() ) {
            return 1;
        }
    } else {
        answer = NULL;
    }

    return 0;
}

int
NodalRecoveryModel :: updateRegionRecoveredValues(const IntArray &regionNodalNumbers,
                                                  int regionValSize, const FloatArray &rhs)
{
    int nnodes = domain->giveNumberOfDofManagers();

    // update recovered values
    for ( int node = 1; node <= nnodes; node++ ) {
        // find nodes in region
        if ( regionNodalNumbers.at(node) ) {
            FloatArray &nodalVal = this->nodalValList [ node ];
            nodalVal.resize(regionValSize);
            for ( int i = 1; i <= regionValSize; i++ ) {
                nodalVal.at(i) = rhs.at( ( regionNodalNumbers.at(node) - 1 ) * regionValSize + i );
            }
        }
    } // end update recovered values

    return 1;
}

int
NodalRecoveryModel :: initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans, Set &region)
{
    int nnodes = domain->giveNumberOfDofManagers();
    IntArray elementRegion = region.giveElementList();

    regionNodalNumbers.resize(nnodes);
    regionNodalNumbers.zero();
    regionDofMans = 0;

    for ( int i = 1; i <= elementRegion.giveSize(); i++ ) {
        int ielem = elementRegion.at(i);
        Element *element = domain->giveElement(ielem);

        int elemNodes = element->giveNumberOfDofManagers();

        // determine local region node numbering
        for ( int elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            int node = element->giveDofManager(elementNode)->giveNumber();
            if ( regionNodalNumbers.at(node) == 0 ) { // assign new number
                regionNodalNumbers.at(node) = ++regionDofMans;
            }
        }
    }

    return 1;
}

int
NodalRecoveryModel :: giveRegionRecordSize()
{
    if ( this->nodalValList.begin() != this->nodalValList.end() ) {
        // the container is not empty
        return this->nodalValList.begin()->second.giveSize();
    } else {
        OOFEM_WARNING("data not yet initialized");
        return 0;
    }
}
} // end namespace oofem

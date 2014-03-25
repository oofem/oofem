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
NodalRecoveryModel :: NodalRecoveryModel(Domain *d) : nodalValList(0), virtualRegionMap(0)
{
    IntArray map;

    stateCounter = 0;
    domain = d;
    this->setRecoveryMode(0, map); // set whole domain recovery by default

    this->init();

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
    int nnodes = this->nodalValList.size();
    this->nodalValList.clear();
    this->nodalValList.resize(nnodes);
    return 1;
}

int
NodalRecoveryModel :: giveNodalVector(const FloatArray * &answer, int node, int region)
{
    std::map< int, FloatArray >::iterator it = this->nodalValList[node-1].find(region);
    if ( it != this->nodalValList[node-1].end() ) {
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
NodalRecoveryModel :: init()
{
    /* Initializes receiver for given domain */
    int nnodes = domain->giveNumberOfDofManagers();

    // allocate array of nodal dictionaries, containing nodal values for each region
    this->nodalValList.clear();
    this->nodalValList.resize(nnodes);

    this->stateCounter = 0;
    this->valType = IST_Undefined;
    //
    return 1;
}

int
NodalRecoveryModel :: updateRegionRecoveredValues(const int ireg, const IntArray &regionNodalNumbers,
                                                  int regionValSize, const FloatArray &rhs)
{
    int nnodes = domain->giveNumberOfDofManagers();

    // update recovered values
    for ( int node = 1; node <= nnodes; node++ ) {
        // find nodes in region
        if ( regionNodalNumbers.at(node) ) {
            FloatArray &nodalVal = this->nodalValList[node-1][ireg];
            nodalVal.resize(regionValSize);
            for ( int i = 1; i <= regionValSize; i++ ) {
                nodalVal.at(i) = rhs.at( ( regionNodalNumbers.at(node) - 1 ) * regionValSize + i );
            }
        }
    } // end update recovered values

    return 1;
}

int
NodalRecoveryModel :: initRegionNodeNumbering(IntArray &regionNodalNumbers, int &regionDofMans, int reg)
{
    int nelem = domain->giveNumberOfElements();
    int nnodes = domain->giveNumberOfDofManagers();

    regionNodalNumbers.resize(nnodes);
    regionNodalNumbers.zero();
    regionDofMans = 0;

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        ElementGeometry *elementGeometry = domain->giveElementGeometry(ielem);
        if ( this->giveElementVirtualRegionNumber(ielem) != reg ) {
            continue;
        }

        int elemNodes = elementGeometry->giveNumberOfDofManagers();

        // determine local region node numbering
        for ( int elementNode = 1; elementNode <= elemNodes; elementNode++ ) {
            int node = elementGeometry->giveDofManager(elementNode)->giveNumber();
            if ( regionNodalNumbers.at(node) == 0 ) { // assign new number
                regionNodalNumbers.at(node) = ++regionDofMans;
            }
        }
    }

    return 1;
}

int
NodalRecoveryModel :: giveRegionRecordSize(int reg, InternalStateType type)
{
    int nelem = domain->giveNumberOfElements();

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        if ( this->giveElementVirtualRegionNumber(ielem) == reg ) {
            //return domain->giveElement(ielem)->giveIPValueType(type);
        }
    }

    OOFEM_WARNING2("NodalRecoveryModel::giveRegionRecordSize: bad region number (%d) or no element in given region\n", reg);
    return 0;
}

int
NodalRecoveryModel :: giveElementVirtualRegionNumber(int ielem)
{
    return this->virtualRegionMap.at( domain->giveElementGeometry(ielem)->giveRegionNumber() );
}

void
NodalRecoveryModel :: setRecoveryMode(int nvr, const IntArray &vrmap)
{
    if ( nvr > 0 ) { // virtual regions, use provided mapping
        if ( vrmap.giveSize() != domain->giveNumberOfRegions() ) {
            //OOFEM_ERROR ("NodalRecoveryModel::setRecoveryMode: invalid size of virtualRegionMap");
        }

        this->virtualRegionMap = vrmap;
        this->numberOfVirtualRegions = nvr;
    } else if ( nvr < 0 ) { // use real regions, set up map accordingly
        int i, nreg = domain->giveNumberOfRegions();
        this->virtualRegionMap.resize(nreg);
        for ( i = 1; i <= nreg; i++ ) {
            this->virtualRegionMap.at(i) = i;
        }

        this->numberOfVirtualRegions = nreg;
    } else { // nvr == 0 (whole domain recovery, emulated by a single virtual region)
        int i, nreg = domain->giveNumberOfRegions();
        this->virtualRegionMap.resize(nreg);
        for ( i = 1; i <= nreg; i++ ) {
            this->virtualRegionMap.at(i) = 1;
        }

        this->numberOfVirtualRegions = 1;
    }
}
} // end namespace oofem

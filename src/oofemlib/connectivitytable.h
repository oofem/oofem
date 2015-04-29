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

#ifndef contable_h
#define contable_h

#include "oofemcfg.h"

#include <vector>

namespace oofem {
class IntArray;
class Domain;

/**
 * Class representing connectivity table. Usually attribute of domain. Provides
 * selected connectivity information services for domain.
 *
 * Its tasks are
 * - Creating connectivity table - method InstanciateYourself.
 * - Returning number of Elements belonging to node.
 * - Returning j-th element belonging to node i.
 */
class OOFEM_EXPORT ConnectivityTable
{
private:
    /// Pointer to domain to which receiver belongs to.
    Domain *domain;

    /// Nodal connectivity table for domain.
    std::vector< IntArray > nodalConnectivity;
    /// Flag indicating assembled connectivity table for domain.
    int nodalConnectivityFlag;

public:
    /**
     * Constructor. Creates new Connectivity table belonging to given domain.
     */
    ConnectivityTable(Domain * d) : domain(d), nodalConnectivity(), nodalConnectivityFlag(0) { }
    /// Destructor
    ~ConnectivityTable() { }
    /// reset receiver to an initial state (will force table update, when needed next time)
    void reset();

    void setDomain(Domain *ipDomain) { domain = ipDomain; }

    /**
     * Builds connectivity table. This table contains for each dofManager the list of
     * elements sharing it.
     */
    void instanciateConnectivityTable();
    /**
     * @param dofman DofManger number.
     * @return Connectivity array for dofman.
     */
    const IntArray *giveDofManConnectivityArray(int dofman);
    /**
     * Returns list of neighboring elements to given elements (they are included too).
     * Neighbor is defined as element sharing the given element node.
     * @param answer List of neighbors + given elements, every element contained only once.
     * @param elemList List of elements, which neighborhood is searched.
     */
    void giveElementNeighbourList(IntArray &answer, IntArray &elemList);
    /**
     * Returns list of elements sharing given nodes.
     * @param answer List of elements, every element contained only once.
     * @param nodeList List of nodes, which neighborhood is searched.
     */
    void giveNodeNeighbourList(IntArray &answer, IntArray &nodeList);
};
} // end namespace oofem
#endif // conTable_h

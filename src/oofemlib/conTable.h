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

#ifndef contable_h
#define contable_h

#include "alist.h"

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
class ConnectivityTable
{
private:
    /// Pointer to domain to which receiver belongs to.
    Domain *domain;

    /// Array of connectivities for dofmanagers.
    int *dofManagersConnectivity;

    /// Nodal connectivity table for domain.
    AList< IntArray >nodalConnectivity;
    /// Flag indicating assembled connectivity table for domain.
    int nodalConnectivityFlag;

public:
    /**
     * Constructor. Creates new Connectivity table belonging to given domain.
     */
    ConnectivityTable(Domain *d) : nodalConnectivity(0)
    {
        domain = d;
        dofManagersConnectivity = NULL;
        nodalConnectivityFlag = 0;
    }
    /// Destructor
    ~ConnectivityTable();
    /// reset receiver to an initial state (will force table update, when needed next time)
    void reset();

#ifdef __OOFEG
    /*
     * Creates connectivity table.
     * Alocates space for dofManagersConnectivity and dofManagersValues attributes.
     * Availbale only in oofeg.
     */
    // void               allocateConnectivityTable ();
    /*
     * Builds the dofManagersConnectivity and dofManagersValues attributes.
     * Availbale only in oofeg.
     */
    // void               instanciateYourself(oofegGraphicContext& gc);
    // Resets dofManagersConnectivity and dofManagersValues attributes (Availbale only in oofeg).

    // void               resetYourself ();
    // Returns Number Of Elements sharing the given DofManager  (Availbale only in oofeg).
    // int                giveNumberOfElementsInDofMngr (int i);

    // Returns Value assigned to particular DofManager using dofManagersValues attribute (Availbale only in oofeg).
    // double             giveDofMngrValue (oofegGraphicContext& gc, int i);
    // Returns min and max values for values stored in dofManagersValues (Availbale only in oofeg).
    // void               giveMinMaxVal (oofegGraphicContext& gc, double*, double*);
#endif

    /*
     * void               instanciateReactionForceTable (TimeStep*);
     * int                giveNumRestrDofMngrs() {return numRestrManagers;}
     * int                giveRestrDofManager (int i) {return restrDofManager[i];}
     * int                giveNumOfSharedElements (int inode);
     * int                giveSharedElement (int inode, int ielem);
     */

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

    /// Prints receiver contents on output.
    void printYourself();
};
} // end namespace oofem
#endif // conTable_h

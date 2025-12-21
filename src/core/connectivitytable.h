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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include <vector>
#include <list>
#include <memory>

#include "oofemenv.h"
#include "intarray.h"
#include "domain.h"
#include "dofmanager.h"
#include "elementgeometrytype.h"


#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
class Domain;

/**
 * @NOTE Consider just keeping the connectivity information
 *  (list of shared elements)
 *  and allocate dofmanages in Domain and keep nodes on element standard DofManList.  
 *  
 */
class SharedBoundaryEntity {
public:
    struct elementRec {
        int elementID;
        int boundaryID;
    };
    // boundary entity internal; DOfManagers (ordering is interpolation dependent) 
    // number of dofmans determines the interpolation order on the entity, however, 
    // std::list<std::unique_ptr<DofManager>> dofMans; 
    // we also store the interpolation order explicitly
    int interpolationOrder;
    // elements sharing the boundary entity (OPTIONAL, required by DG)
    // @TODO: We also need mapping from elements to boundary entities
    std::list<elementRec> elements; 
    // @NOTE this may be redundant, 
    IntArray nodes;  // nodes defining the boundary entity and its orientation
    Element_Geometry_Type geomType; // geometry type of the boundary entity
    int spatialDimension; // entity spatial dimension (1 for edge, 2 for surface)
};


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

    /// Element colors
    IntArray elementColoring;
    /// flag indicating assembled element coloring
    bool elementColoringFlag;   

    // shared element boundary entities (edges & surfaces)
    std::vector<std::unique_ptr<SharedBoundaryEntity>> sharedBoundaryEntities;

#ifdef _OPENMP
    omp_lock_t initLock;
#endif
#ifdef _MSC_VER
    // workaround for MSVC (compile error)
    // make the class noncopyable (since it contains vector<unique_ptr<SharedBoundaryEntity>)
    ConnectivityTable(ConnectivityTable const &)=delete;
    void operator=(ConnectivityTable const &)=delete;
#endif
public:
    /**
     * Constructor. Creates new Connectivity table belonging to given domain.
     */
      ConnectivityTable(Domain * d) ;
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
    void giveElementNeighbourList(IntArray &answer, const IntArray &elemList);
    /**
     * Returns list of elements sharing any given nodes.
     * @param answer List of elements, every element contained only once.
     * @param nodeList List of nodes, which neighborhood is searched.
     */
    void giveNodeNeighbourList(IntArray &answer, IntArray &nodeList);
    /** 
     * Return list of elements sharing all given nodes
     * @param answer List of elements that have all given nodes
     * @param nodeList List of nodes
    */
    void giveElementsWithNodes(IntArray &answer, const IntArray& nodes);
    
    
    /** Builds element coloring, assigning to element color, such that no neighboring elements have the same color
     */
    void buildElementColoring();
    /** Returns element color as determined by coloring algorithm*/
    int getElementColor(int e);


    /** Build list of shared element edges and surfaces */
    void buildSharedBoundaryEntities(Domain *d); 
    /** Access SharedBoundaryEntities */
    SharedBoundaryEntity* giveBoundaryEntity(int id);
};
} // end namespace oofem
#endif // conTable_h

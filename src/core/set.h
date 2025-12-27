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

#ifndef set_h
#define set_h

#include "femcmpnn.h"
#include "intarray.h"

#include <list>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace oofem {
///@name Input fields for Set
//@{
#define _IFT_Set_Name "set"
#define _IFT_Set_nodes "nodes" ///< List of specific node indices.
#define _IFT_Set_allNodes "allnodes" ///< List of specific node indices.
#define _IFT_Set_nodeRanges "noderanges" ///< List of node index ranges.
#define _IFT_Set_elements "elements" ///< List of specific element indices.
#define _IFT_Set_allElements "allelements" ///< Will generate a list of all elements in the domain.
#define _IFT_Set_elementRanges "elementranges" ///< List of element index ranges.
#define _IFT_Set_elementBoundaries "elementboundaries" ///< Interleaved array of element index + boundary number
#define _IFT_Set_elementEdges "elementedges" ///< Interleaved array of element index + edge number
#define _IFT_Set_elementSurfaces "elementsurfaces" ///< Interleaved array of element index + surface number
#define _IFT_Set_internalElementNodes "internalelementnodes" ///< Interleaved array of element index + internal node number
#define _IFT_Set_elementprops "elemprops" ///< Element properties 
#define _IFT_Set_dofmanprops "dofmanprops" ///< Dof manager properties
//@}

class EntityRenumberingFunction;
class Range;

/**
 * Set of elements, boundaries, edges and/or nodes.
 * Describes a collection of components which are given easy access to for example boundary conditions.
 * 
 * @note bp: sets are defined using component labels (global numbers). In the initialization, the sets are resolved to local component numbers.
 * In parallel mode (MPI), the sets are resolved to local component numbers on each processor, but set definition is TYPICALLY duplicated on all processors.
 * In turn, the sets are resolved to local component numbers on each processor. As the set definition can contain global component, 
 * hat are not present on the processor, these entries are removed from the set during resolution process.
 * To support also dynamic load balancing, the set resolution has to be done after each load balancing step and thus we need to keep the original set definition.
 * This can be achieved by keeping the original input record of the set.
 *   
 * @author Mikael Öhman
 */
class OOFEM_EXPORT Set : public FEMComponent
{
protected:
    IntArray elements; ///< Element numbers.
    mutable bool mElementListIsSorted;
    mutable IntArray mElementsSorted;
    IntArray elementBoundaries; /// Element numbers + boundary numbers (interleaved).
    IntArray elementEdges; /// Element numbers + edge numbers (interleaved).
    IntArray elementSurfaces; /// Element numbers + surface numbers (interleaved).
    IntArray elementInternalNodes; /// Element numbers + internal dof manager numbers (interleaved).
    IntArray nodes; ///< Node numbers.
    IntArray totalNodes; ///< Unique set of nodes (computed).
    std::string elementprops; ///< Element properties.
    std::string dofmanprops; ///< Dof manager properties.
#if 0
    /// receiver original input record
    std::unique_ptr<InputRecord> inputRec;
#endif
    bool nodalListInitialized = false;
#ifdef _OPENMP
    omp_lock_t initLock;
#endif
public:
    /**
     * Creates a empty set with given number and belonging to given domain.
     * @param n Set number.
     * @param d Domain to which component belongs to.
     */
    Set(int n, Domain * d) ;
    virtual ~Set() { }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    /**
     * Returns list of elements within set.
     * @return List of element numbers.
     */
    const IntArray &giveElementList();
    /**
     * Returns list of element boundaries within set.
     * Boundaries are either surfaces (3d), edges (2d), or corners (1d).
     * @return List of element boundaries.
     */
    const IntArray &giveBoundaryList();
    /**
     * Returns list of element edges within set (must be edges of 3D elements).
     * @return List of element edges.
     */
    const IntArray &giveEdgeList();
    /**
     * Returns list of element surfaces within set 
     * @return List of element surfaces.
     */
    const IntArray &giveSurfaceList();
    /**
     * Returns list of internal element dof managers within set 
     * @return List of internal element dofManagers.
     */
    const IntArray &giveInternalElementDofManagerList();
    /**
     * Returns list of all nodes within set.
     * This list is computed automatically, based on all elements, boundaries, edges, and specified nodes within the set.
     * @return List of node numbers.
     */
    const IntArray &giveNodeList();
    /**
     * Returns list of all directly specified nodes (excluding those generated from elements).
     * This list is exactly the list given in the input.
     * @note This is useful in for example, remeshing code, and should rarely be used elsewhere.
     * @return List of node numbers.
     */
    const IntArray &giveSpecifiedNodeList();
    /**
     * Sets list of elements within set.
     */
    void setElementList(IntArray newElements);
    /**
     * Sets list of element boundaries within set.
     */
    void setBoundaryList(IntArray newBoundaries);
    /**
     * Sets list of element edges within set (must be edges of 3D elements).
     */
    void setEdgeList(IntArray newEdges);
    /**
     * Sets list of nodes within set.
     */
    void setNodeList(IntArray newNodes);

    /**
     * Clears the entire set.
     */
    void clear();
    /**
     * Initialize the element set to contain all elements in the receiver domain
     */
    void addAllElements();
    /// Return True if given element is contained
    bool hasElement(int elem) const;

    std::string giveElementProperties() const { return this->elementprops; }
    std::string giveDofManProperties() const { return this->dofmanprops; }


    void updateLocalNumbering(EntityRenumberingFunctor &f) override;
    /**
     * Renumbering of nodes (could change due to load balancing).
     */
    void updateLocalNodeNumbering(EntityRenumberingFunctor &f);
    /**
     * Renumbering of nodes (could change due to load balancing).
     */
    void updateLocalElementNumbering(EntityRenumberingFunctor &f);

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "Set"; }
    const char *giveInputRecordName() const override { return _IFT_Set_Name; }

protected:
    /**
     * Converts list ranges to list of individual values + individually specified values.
     */
    void computeIntArray(IntArray &answer, const IntArray &specified, std :: list< Range >ranges);
};

extern Set dummySet;

} // end namespace oofem
#endif // set_h

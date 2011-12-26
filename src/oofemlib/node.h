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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef node_h
#define node_h

#include "dofmanager.h"
#include "domain.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
class Dof;
class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Class implementing node in finite element mesh. Node possess degrees of freedom
 * (see base class DofManager).
 * Node is attribute of few elements and it is managed by domain.
 * Node manages its position in space, and if specified
 * local coordinate system in node. If local coordinate system is defined, all
 * equilibrium equations are assembled in this system and therefore all DOFs and
 * applied  boundary and initial conditions apply in this local coordinate system.
 * By default, global coordinate system is assumed in each node.
 * For description, how to prescribe local coordinate system in node, see
 * input file description section.
 *
 * Tasks include:
 * - Managing its position in space. In geometrically linear analysis, this
 *   only amounts to managing its coordinates.
 * - Managing its local coordinate system. If it is defined, all equilibrium
 *   equations are assembled within it.
 *   This system is defined by triplet of unit vectors of local coordinates
 *   expressed in terms of global axes. This tripled is stored in
 *   localCoordinateSystem, its item ij is angle between e'(i) and e(j),
 *   where e' is local axis.
 * - Managing its degrees of freedom (method giveDof) ;
 * - Calculating its nodal load vector;
 * - Printing and updating at end of step ;
 * - Managing its swapping to and from disk.
 */
class Node : public DofManager
{
protected:
    /// Array storing nodal coordinates.
    FloatArray coordinates;
    /**
     * Triplet defining the local coordinate system in node.
     * Value at position (i,j) represents angle between e'(i) and e(j),
     * where e' is base vector of local coordinate system and e is
     * base vector of global c.s.
     */
    FloatMatrix *localCoordinateSystem;

public:

    /**
     * Constructor. Creates a node belonging to domain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    Node(int n, Domain *aDomain);
    /// Destructor.
    ~Node();

    bool hasCoordinates() { return true; }
    double giveCoordinate(int i);
    FloatArray *giveCoordinates() { return & coordinates; }
    /**
     * Sets node coordinates to given array.
     * @param coords New coordinates for node.
     */
    void setCoordinates(const FloatArray &coords) { this->coordinates = coords; }
    /**
     * Returns updated ic-th coordinate of receiver. Return value is computed
     * as coordinate + scale * displacement, where corresponding displacement is obtained
     * from corresponding nodal DOF. Local coordinate system is taken into account.
     * Useful mainly for postprocessing.
     * @param ic Index of coordinate.
     * @param tStep Time step for the displacement.
     * @param type Equation ID for displacement of node.
     * @param scale Scaling of displacement.
     */
    virtual double giveUpdatedCoordinate(int ic, TimeStep *tStep,
                                               EquationID type, double scale = 1.);

    // local coordinate system
    /// Returns nonzero if node has prescribed local coordinate system.
    bool hasLocalCS() { return ( localCoordinateSystem != NULL ); }
    /** Returns pointer to local coordinate triplet in node.
     * If not defined, returns NULL.
     * @return Triplet defining the local coordinate system in node.
     * Value at position (i,j) represents angle between e'(i) and e(j),
     * where e' is base vector of local coordinate system and e is
     * base vector of global c.s.
     */
    FloatMatrix *giveLocalCoordinateTriplet() { return localCoordinateSystem; }
    /** Returns true, if the local coordinate systems of receiver and given node are the same */
    bool hasSameLCS(Node *remote);

    void computeDofTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    void computeLoadTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    int requiresTransformation() { return ( this->hasLocalCS() || hasSlaveDofs ); }
    virtual void computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    /**
     * Computes transformation matrix for given dofs from global c.s to
     * rotated coordinate system (given by localCoordinateSystem member value).
     * @param answer Result of [map.giveSize(), map.giveSize()] size.
     * @param map Array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transformation matrix is
     * assembled. If dofIDArry is NULL, then all receiver dofs are assumed.
     */
    void computeGNDofTransformation(FloatMatrix &answer, const IntArray *map);
    /**
     * Updates receiver at end of time step (i.e. after equilibrium has been reached).
     * If EngngModel formulation ( see giveFormulation() member function) returns actualized
     * Lagrange mode, receiver updates its coordinates according to solution.
     * @param tStep Time step for which to update.
     * @see EngngModel::giveFormulation
     */
    void updateYourself(TimeStep *tStep);

    // miscellaneous
    const char *giveClassName() const { return "Node"; }
    classType giveClassID() const { return NodeClass; }
    IRResultType initializeFrom(InputRecord *ir);
    /// Prints receiver state on stdout. Useful for debugging.
    void printYourself();
    virtual int checkConsistency();
    /// Returns true if dof of given type is allowed to be associated to receiver
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave || type == DT_active ); }

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &);
#endif
};
} // end namespace oofem
#endif // node_h

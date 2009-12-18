/* $Header: /home/cvs/bp/oofem/oofemlib/src/node.h,v 1.11 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ******************
//   *** CLASS NODE ***
//   ******************


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
 * Class implementing node in finite element mesh. Node posses degrees of freedom
 * (see base class DofManager).
 * Node is atribute of few elements and it is managed by domain.
 * Node manages its positon in space, and if specified
 * local coordinate system in node. If local coordinate system is defined, all
 * equilibrium equations are assembled in this system and therefore all DOFs and
 * applied  boundary and initial conditions apply in this local coordinate system.
 * By default, global coordinate system is assumed in each node.
 * For description, how to prescribe local coordinate system in node, see
 * input file description section.
 */
class Node : public DofManager
{
    /*
     * This class implements a node in a finite element mesh. A node is an attri-
     * bute of a domain. It is usually also attribute of a few elements.
     * DESCRIPTION
     * The node possesses 'numberOfDofs' degrees of freedom, stored in 'dofArray'.
     * In 'loadArray' it stores the number of every nodal load it is subjected to
     * (typically, concentrated forces and moments).
     * In 'locationArray' the node stores the equation number of each of its dofs.
     * This location array is used by the node for assembling its load vector to
     * the right-hand side of the linear system ; it is also used by elements for
     * calculating their own location arrays.
     * TASKS
     * - managing its position in space. In geometrically linear analysis, this
     *   only amounts to managing its coordinates ;
     * - managing its local coordinate system. If it is defined, all equilibrium
     * equations are assemblebled within it.
     * This system is defined by triplet of unit vectors of local coordinates
     * expressed in terms of global axes. This tripled is stored in
     * localCoordinateSystem, its item ij is angle between e'(i) and e(j),
     * where e' is local axis.
     * - managing its degrees of freedom (method 'giveDof') ;
     * - calculating its nodal load vector;
     * - printing and updating at end of step ;
     * - managing its swapping to and from disk.
     */

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
     * @param n node number in domain aDomain
     * @param aDomain domain to which node belongs
     */
    Node(int n, Domain *aDomain);                      // constructor
    /// Destructor.
    ~Node();                                           // destructor

    // coordinates
    bool        hasCoordinates() { return true; }
    /// Returns i-th coordinate of node.
    double      giveCoordinate(int i);
    /// Returns pointer to node coordinate array.
    FloatArray *giveCoordinates() { return & coordinates; }
    /// Sets node coordinates to given array
    void        setCoordinates (FloatArray& _coords) {this->coordinates = _coords;}
    /**
     * Returns updated ic-th coordinate of receiver. Return value is computed
     * as coordinate + scale * displacement, where corresponding displacement is obtained
     * from cooresponding nodal DOF. Local coodinate system is taken into account.
     * Usefull mainly for postprocessing.
     */
    virtual double       giveUpdatedCoordinate(int ic, TimeStep *tStep,
                                               EquationID type, double scale = 1.);

    // local coordinate system
    /// Returns nonzero if node has prescribed  local coordinate system.
    int          hasLocalCS() { return ( localCoordinateSystem != NULL ); }
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

    /** Computes receiver DOF transformation matrix from global cs. to dofManager specific
     * coordinate system - if mode == _toNodalCS, otherwise reverse transformation is computed.
     * (In the dofManager specifis cs the governing equations are assembled, for example the
     * local coordinate system in node). This transformation may not be ortogonal.
     * @param answer computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, squre answer is produced anly if all
     * dof transformation is required).
     * @param dofIDArry array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     */
    void computeDofTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    /** Computes receiver LOAD transformation matrix from global cs. to dofManager specific
     * coordinate system - if mode == _toNodalCS, otherwise reverse transformation is computed.
     * (In the dofManager specifis cs the governing equations are assembled, for example the
     * local coordinate system in node). This transformation may not be ortogonal.
     * @param answer computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, squre answer is produced anly if all
     * dof transformation is required).
     * @param dofIDArry array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     */
    void computeLoadTransformation(FloatMatrix &answer, const IntArray *dofIDArry, DofManTransfType mode);
    /**
     * Indicates, whether dofManager requires the transformation from global c.s. to
     * dof manager specific coordinate system.
     * @return nonzero if transformation is necessary, even for single dof.
     */
    int requiresTransformation() { return ( this->hasLocalCS() || hasSlaveDofs ); }
    /** Computes the load vector of receiver in given time.
     *  @param answer load vector.
     *  @param stepN time step when answer is computed.
     *  @param mode determines response mode.
     */
    virtual void         computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    /**
     * Computes transformation matrix for given dofs from global c.s to
     * rotated coordinate system (given by localCoordinateSystem member value).
     * @param answer result of [map.giveSize(), map.giveSize()] size.
     * @param map array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transfromation mtrx is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     */
    void computeGNDofTransformation(FloatMatrix &answer, const IntArray *map);
    // time step termination
    /** Updates receiver at end of time step (i.e. after equilibrium has been reached).
     * If EngngModel formulation ( see giveFormulation() member function) returns actualized
     * Lagrange mode, node updates its coordinates according to solution.
     * @see EngngModel::giveFormulation().
     */
    void         updateYourself(TimeStep *);

    // miscellaneous
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "Node"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType    giveClassID() const { return NodeClass; }
    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);
    //virtual IntArray* ResolveDofIDArray (char* initString);
    /// prints receiver state on stdout. Usefull for debuging.
    void         printYourself();
    /**
     * Checks internal data consistency in node.
     * Current implementation checks (when receiver has slave dofs) if receiver has the same
     * coordinate system as master dofManager of slave dof.
     * @return nonzero if receiver check is o.k.
     */
    virtual int    checkConsistency();
    /// Returns true if dof of given type is allowed to be associated to receiver
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_simpleSlave ); }

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


#ifdef __OOFEG
    void         drawYourself(oofegGraphicContext &);
#endif
};

} // end namespace oofem
#endif // node_h

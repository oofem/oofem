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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef node_h
#define node_h

#include "dofmanager.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include <memory>

///@name Input fields for Node
//@{
#define _IFT_Node_Name "node"
#define _IFT_Node_coords "coords"
#define _IFT_Node_lcs "lcs"
//@}

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
class OOFEM_EXPORT Node : public DofManager
{
protected:
    /**
     * Triplet defining the local coordinate system in node.
     * Value at position (i,j) represents angle between e'(i) and e(j),
     * where e' is base vector of local coordinate system and e is
     * base vector of global c.s.
     */
    std::unique_ptr<FloatMatrix> localCoordinateSystem;

public:
    /**
     * Constructor. Creates a node belonging to domain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    Node(int n, Domain * aDomain);

    /**
     * Sets node coordinates to given array.
     * @param coords New coordinates for node.
     */
    void setCoordinates(FloatArray coords) { this->coordinates = std :: move(coords); }
    /**
     * Returns updated ic-th coordinate of receiver. Return value is computed
     * as coordinate + scale * displacement, where corresponding displacement is obtained
     * from corresponding nodal DOF. Local coordinate system is taken into account.
     * Useful mainly for postprocessing.
     * @param ic Index of coordinate.
     * @param tStep Time step for the displacement.
     * @param scale Scaling of displacement.
     */
    virtual double giveUpdatedCoordinate(int ic, TimeStep *tStep, double scale = 1.);
    /**
     * Returns updated coordinate of receiver. Return value is computed
     * as coordinate + scale * displacement, where corresponding displacement is obtained
     * from corresponding nodal DOF.
     * Local coordinate system is taken into account and the answer is given in global coordinates.
     * @param answer Updated coordinate.
     * @param tStep Time step for the displacement.
     * @param scale Scaling of displacement.
     */
    virtual void giveUpdatedCoordinates(FloatArray &answer, TimeStep *tStep, double scale = 1.);

    // local coordinate system
    /// Returns nonzero if node has prescribed local coordinate system.
    bool hasLocalCS() { return ( localCoordinateSystem != nullptr ); }
    /** Returns pointer to local coordinate triplet in node.
     * If not defined, returns NULL.
     * @return Triplet defining the local coordinate system in node.
     * Value at position (i,j) represents angle between e'(i) and e(j),
     * where e' is base vector of local coordinate system and e is
     * base vector of global c.s.
     */
    FloatMatrix *giveLocalCoordinateTriplet() { return localCoordinateSystem.get(); }
    /** Returns true, if the local coordinate systems of receiver and given node are the same */
    bool hasSameLCS(Node *remote);

    bool computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry) override;
    bool requiresTransformation() override { return ( this->hasLocalCS() || hasSlaveDofs ); }
    void computeLoadVector(FloatArray &answer, Load *load, CharType type, TimeStep *tStep, ValueModeType mode) override;

    /**
     * Updates receiver at end of time step (i.e. after equilibrium has been reached).
     * If EngngModel formulation ( see giveFormulation() member function) returns actualized
     * Lagrange mode, receiver updates its coordinates according to solution.
     * @param tStep Time step for which to update.
     * @see EngngModel::giveFormulation
     */
    void updateYourself(TimeStep *tStep) override;

    // miscellaneous
    const char *giveClassName() const override { return "Node"; }
    const char *giveInputRecordName() const override { return _IFT_Node_Name; }
    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    void printYourself() override;
    int checkConsistency() override;
    bool isDofTypeCompatible(dofType type) const override { return ( type == DT_master || type == DT_simpleSlave || type == DT_active ); }
    virtual int giveQcNodeType() { return 0; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif
};
} // end namespace oofem
#endif // node_h

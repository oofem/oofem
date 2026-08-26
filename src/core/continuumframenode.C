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

#include "continuumframenode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "node.h"
#include "element.h"
#include "elementgeometrytype.h"
#include "feinterpol.h"
#include "classfactory.h"

namespace oofem {
REGISTER_DofManager(ContinuumFrameNode);

// Adds one master-DOF contribution to a rotational slave-DOF constraint.
// Master nodes that do not carry the required translational DOF are skipped.
static void
addRotationTerm(FloatArray &coeffs, IntArray &masterNodes, IntArray &masterDofIDs,
                double coeff, Node *masterNode, DofIDItem masterDofID)
{
    if ( !masterNode->hasDofID(masterDofID) ) {
        return;
    }
    coeffs.resizeWithValues( coeffs.giveSize() + 1 );
    coeffs.at( coeffs.giveSize() ) = coeff;
    masterNodes.followedBy( masterNode->giveNumber() );
    masterDofIDs.followedBy( (int) masterDofID );
}

ContinuumFrameNode :: ContinuumFrameNode(int n, Domain *aDomain) : HangingNode(n, aDomain)
{ }

void ContinuumFrameNode :: postInitialize()
{
    // First let HangingNode set up the master element, the local coordinates and
    // the translational slave DOFs by shape-function interpolation.
    HangingNode :: postInitialize();

    // Rotational DOFs (R_u, R_v, R_w) that are defined as SLAVE DOFs are constrained to the
    // infinitesimal rotation omega = 1/2 curl(u) of the host continuum, so a beam/frame node
    // embedded in a solid mesh inherits the local element rotation and the frame element node
    // order no longer has to match the continuum. Only slaved rotational DOFs are treated here:
    // if the rotational DOFs are free (master) or fixed - e.g. a frame node whose torsion is
    // fixed and whose bending rotations are carried by the beam element itself - they are left
    // untouched. Implemented for the linear tetrahedron, whose constant shape-function gradients
    // give a constant rotation.
    SlaveDof *rotU = this->hasDofID(R_u) ? dynamic_cast< SlaveDof * >( this->giveDofWithID(R_u) ) : nullptr;
    SlaveDof *rotV = this->hasDofID(R_v) ? dynamic_cast< SlaveDof * >( this->giveDofWithID(R_v) ) : nullptr;
    SlaveDof *rotW = this->hasDofID(R_w) ? dynamic_cast< SlaveDof * >( this->giveDofWithID(R_w) ) : nullptr;

    if ( !( rotU || rotV || rotW ) ) {
        return;
    }

    Element *e = this->giveDomain()->giveGlobalElement(this->masterElement);
    if ( !e ) {
        OOFEM_ERROR("Master element %d doesn't exist.", this->masterElement);
    }
    if ( e->giveGeometryType() != EGT_tetra_1 ) {
        OOFEM_ERROR("Continuum frame node %d has slaved rotational DOFs, but master element %d is not a "
                    "linear tetrahedron (EGT_tetra_1); the continuum rotational constraint is only "
                    "implemented for linear tetrahedra. Use a plain hangingnode (with the rotational DOFs "
                    "free or fixed) to embed a frame node without this constraint.",
                    this->giveNumber(), this->masterElement);
    }

    FEInterpolation *fei = e->giveInterpolation();
    FloatArray lcoords;
    fei->global2local(lcoords, this->coordinates, FEIElementGeometryWrapper(e));
    FloatMatrix dNdX;
    fei->evaldNdx(dNdX, lcoords, FEIElementGeometryWrapper(e));
    const int nnodes = e->giveNumberOfNodes();

    // theta_x = 1/2 (du_z/dy - du_y/dz)
    if ( rotU ) {
        FloatArray coeffs;
        IntArray masterNodeIDs, masterDofIDs;
        for ( int i = 1; i <= nnodes; ++i ) {
            Node *masterNode = e->giveNode(i);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs,  0.5 * dNdX.at(i, 2), masterNode, D_w);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs, -0.5 * dNdX.at(i, 3), masterNode, D_v);
        }
        rotU->initialize(masterNodeIDs, masterDofIDs, coeffs);
    }

    // theta_y = 1/2 (du_x/dz - du_z/dx)
    if ( rotV ) {
        FloatArray coeffs;
        IntArray masterNodeIDs, masterDofIDs;
        for ( int i = 1; i <= nnodes; ++i ) {
            Node *masterNode = e->giveNode(i);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs,  0.5 * dNdX.at(i, 3), masterNode, D_u);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs, -0.5 * dNdX.at(i, 1), masterNode, D_w);
        }
        rotV->initialize(masterNodeIDs, masterDofIDs, coeffs);
    }

    // theta_z = 1/2 (du_y/dx - du_x/dy)
    if ( rotW ) {
        FloatArray coeffs;
        IntArray masterNodeIDs, masterDofIDs;
        for ( int i = 1; i <= nnodes; ++i ) {
            Node *masterNode = e->giveNode(i);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs,  0.5 * dNdX.at(i, 1), masterNode, D_v);
            addRotationTerm(coeffs, masterNodeIDs, masterDofIDs, -0.5 * dNdX.at(i, 2), masterNode, D_u);
        }
        rotW->initialize(masterNodeIDs, masterDofIDs, coeffs);
    }
}
} // end namespace oofem

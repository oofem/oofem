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

#include "Elements/PlaneStress/basicelement.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(BasicElement);

// Interpolator describing shape functions for the approximated unknowns.
// 1 -> first spatial index, 2 -> second spatial index
FEI2dTrLin BasicElement :: interp(1, 2);    

BasicElement :: BasicElement(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    this->numberOfDofMans  = 3;
    this->numberOfGaussPoints = 1;
}



IRResultType
BasicElement :: initializeFrom(InputRecord *ir)
{
    // Initialise the element from the input record (text line from input file)  
    IRResultType result = this->NLStructuralElement :: initializeFrom(ir); // Call initializeFrom for the base class
    return result; 
}


FEInterpolation *BasicElement :: giveInterpolation() const 
{ 
    /* Returns the interpolator used for the element which provide
     * shape functions, their derivatives, area of the element etc.
     */
    return & interp; 
}


double 
BasicElement :: computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    FloatArray *lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->interp.giveTransformationJacobian( *lCoords, FEIElementGeometryWrapper(this) ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness
    
    return detJ * thickness * weight; // dV
}


void
BasicElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Return which DofId's (unknowns) that are associated with the given DofManager/node
   
    answer = {D_u, D_v}; // Displacement in u- and v-direction. See 'dofiditem.h' for a list of basic dof types.
}

void
BasicElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    /* Compute the [3x6] strain-displacement matrix {B} for the element,
     * evaluated at the given gp.
     * {B}*{a} should provide the strains {eps} in Voigt form 
     * {eps} = {eps_xx, eps_yy, gam_xy}^T with {a} being the 
     * solution vector of the element.
     */
    
    /* Evaluate the derivatives of the shape functions at the position of the gp.     
     * dNdx = [dN1/dx1 dN1/dx2
     *         dN2/dx1 dN2/dx2
     *         dN3/dx1 dN3/dx2]
    */  
    FloatMatrix dNdx; 
    this->interp.evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    // Construct the B-matrix
    answer.resize(3, 6);

    answer.at(1, 1) = dNdx.at(1, 1);
    answer.at(1, 3) = dNdx.at(2, 1);
    answer.at(1, 5) = dNdx.at(3, 1);

    answer.at(2, 2) = dNdx.at(1, 2);
    answer.at(2, 4) = dNdx.at(2, 2);
    answer.at(2, 6) = dNdx.at(3, 2);

    answer.at(3, 1) = dNdx.at(1, 2);
    answer.at(3, 2) = dNdx.at(1, 1);
    answer.at(3, 3) = dNdx.at(2, 2);
    answer.at(3, 4) = dNdx.at(2, 1);
    answer.at(3, 5) = dNdx.at(3, 2);
    answer.at(3, 6) = dNdx.at(3, 1);
}


void BasicElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points (integration points)
    // used for integration over the element volume.
    
    if ( integrationRulesArray.size() == 0 ) {
        int iRuleNumber = 1;
        integrationRulesArray = { new GaussIntegrationRule(iRuleNumber, this, 1, 3) };
        // The cross section is responsible for setting up the integration rule - not the element.
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ iRuleNumber - 1 ], this->numberOfGaussPoints, this);
    }
}


/// Optional: Needed in order to have support for nonlinear kinematics / large deformations
void
BasicElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    /* Compute the [3x6] displacement gradient matrix {BH} for the element,
     * evaluated at the given gp.
     * {BH}*{a} should provide the displacement gradient {H} in Voigt form 
     * {H} = {du/dx, dv/dy, du/dy, dv/dx }^T with {a} being the solution vector 
     * of the element.
     */    

    // Evaluate the derivatives of the shape functions at the position of the gp.
    FloatMatrix dNdx;
    this->interp.evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    // Construct the BH-matrix
    answer.resize(4, 6);
    answer.zero();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(3, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}


/// Optional: Needed in order to support edge loads
void
BasicElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /* Returns the [2x4] shape function matrix {N} of the receiver, 
     * evaluated at the given gp.
     * {u} = {N}*{a} gives the displacements at the integration point.
     */ 
          
    // Evaluate the shape functions at the position of the gp. N = [N1, N2]
    FloatArray N;
    this->interp.edgeEvalN( N, 1, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    
    answer.resize(2, 4);
    answer.at(1, 1) = answer.at(2, 2) = N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = N.at(2);
}

/// Optional: Needed in order to support edge loads
void
BasicElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /* 
     * Provides a dof mapping from local edge dofs (1-4) to the 
     * element local dofs (1-6)
     */

    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer = {1, 2, 3, 4};
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = {3, 4, 5, 6};
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer = {5, 6, 1, 2};
    } else {
        OOFEM_ERROR("Wrong edge number (%d), element has three edges.", iEdge);
    }
    
    
}

/// Optional: Needed in order to support edge loads
double
BasicElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    /* Returns the line element ds associated with the given gp on the specific edge.
     * Note: The name is misleading since there is no volume to speak of in this case. 
     * The returned value is used for integration of a line integral (external forces).
     */
    double detJ = this->interp.edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ * gp->giveWeight();
}

/// Optional: Needed in order to support edge loads
int
BasicElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    /* Compute transformation matrix T from edge local coordinate 
     * system to element local c.s (same as global c.s in this case),
     * i.e. f(element local) = T * f(edge local)
     * with  T = [ cos -sin
     *             sin  cos]
     */

    IntArray edgeNodes;
    // give numbers of the nodes on a specific edge
    this->giveInterpolation()->boundaryGiveNodes(edgeNodes, iEdge); 
    Node *nodeA   = this->giveNode(edgeNodes.at(1));
    Node *nodeB   = this->giveNode(edgeNodes.at(2));

    double dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    double dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    double length = sqrt(dx * dx + dy * dy);
    double c = dx / length; // cosine
    double s = dy / length; // sine
    
    // Syntax 1 - set each column in the matrix directly
    answer = { {c, s}, {-s, c} }; 
    
    // Syntax 2 - set each component individually (needs resizing to appropriate size)
    // answer.resize(2, 2);
    // answer.zero();    
    // answer.at(1, 1) = c;     answer.at(1, 2) = -s;
    // answer.at(2, 1) = s;     answer.at(2, 2) = c;
    
    return 1;
}









} // end namespace oofem

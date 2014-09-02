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

#include "../sm/Elements/planestresselement.h"
//#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"

namespace oofem {
PlaneStressElement :: PlaneStressElement(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}



int 
PlaneStressElement :: computeNumberOfDofs() 
{ 
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask; 
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize(); 
  
}


IRResultType
PlaneStressElement :: initializeFrom(InputRecord *ir)
{
    // Initialise the element from the input record   
    return this->NLStructuralElement :: initializeFrom(ir); 
}


void 
PlaneStressElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
//void PlaneStressStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;

    //FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    FEInterpolation *interp = this->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    
    FloatMatrix dNdx; 
    interp->evaldNdx( d, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    
    //interp->evaldNdx( d, * gp->giveNaturalCoordinates(),
    //                 FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(3, d.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d.at(i, 1);
        answer.at(2, i * 2 - 0) = d.at(i, 2);

        answer.at(3, 2 * i - 1) = d.at(i, 2);
        answer.at(3, 2 * i - 0) = d.at(i, 1);
    }
}

void
PlaneStressElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    // @todo not checked if correct
  
    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(3, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}

void 
PlaneStressElement :: giveDofManDofIDMask(int inode, IntArray &answer) const 
{    
    answer = {D_u, D_v};
}


double 
PlaneStressElement :: computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    FloatArray *lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( *lCoords, FEIElementGeometryWrapper(this) ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness
    
    return detJ * thickness * weight; // dV
}




void 
PlaneStressElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray = { new GaussIntegrationRule(1, this, 1, 3) };
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}





// Edge support

void
PlaneStressElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /* Returns the [2x4] shape function matrix {N} of the receiver, 
     * evaluated at the given gp.
     * {u} = {N}*{a} gives the displacements at the integration point.
     */ 
          
    // Evaluate the shape functions at the position of the gp. 
    FloatArray N;
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeEvalN( N, iedge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );  
    answer.beNMatrixOf(N, 2);
}


void
PlaneStressElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    IntArray eNodes;
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->computeLocalEdgeMapping(eNodes,  iEdge);

    answer.resize( eNodes.giveSize() * 2 );
    for ( int i = 1; i <= eNodes.giveSize(); i++ ) {
        answer.at(i * 2 - 1) = eNodes.at(i) * 2 - 1;
        answer.at(i * 2) = eNodes.at(i) * 2;
    }
}


double
PlaneStressElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    /* Returns the line element ds associated with the given gp on the specific edge.
     * Note: The name is misleading since there is no volume to speak of in this case. 
     * The returned value is used for integration of a line integral (external forces).
     */
    double detJ = static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ * gp->giveWeight();
}



int
PlaneStressElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatArray normal(2);

    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeEvalNormal( normal, iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 2);
    answer.zero();        
    answer.at(1, 1) = normal.at(2);
    answer.at(1, 2) = normal.at(1);
    answer.at(2, 1) = -normal.at(1);
    answer.at(2, 2) = normal.at(2);

    return 1;
}


} // end namespace oofem

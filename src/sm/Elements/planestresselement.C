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
PlaneElement :: PlaneElement(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}



int 
PlaneElement :: computeNumberOfDofs() 
{ 
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask; 
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize(); 
  
}


IRResultType
PlaneElement :: initializeFrom(InputRecord *ir)
{
    // Initialise the element from the input record   
    return this->NLStructuralElement :: initializeFrom(ir); 
}



void 
PlaneElement :: giveDofManDofIDMask(int inode, IntArray &answer) const 
{    
    answer = {D_u, D_v};
}


double 
PlaneElement :: computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    FloatArray *lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( *lCoords, FEIElementGeometryWrapper(this) ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness
    
    return detJ * thickness * weight; // dV
}




void 
PlaneElement :: computeGaussPoints()
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
PlaneElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
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
PlaneElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
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



void
PlaneElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeLocal2global( answer, iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}



double
PlaneElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
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
PlaneElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
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




double
PlaneElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}







// Plane stress


PlaneStressElement :: PlaneStressElement(int n, Domain *aDomain) :
    PlaneElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}

void 
PlaneStressElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
//void PlaneStressStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{

    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    
    answer.resize(3, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);

        answer.at(3, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(3, 2 * i - 0) = dNdx.at(i, 1);
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









// Plane strain


PlaneStrainElement :: PlaneStrainElement(int n, Domain *aDomain) :
    PlaneElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}

void 
PlaneStrainElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)

// Returns the [ 4 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    
    
    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);

        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(4, 2 * i - 0) = dNdx.at(i, 1);
    }
}

void
PlaneStrainElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Returns the [ 5 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    // @todo not checked if correct
  
    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}






// Axisymmetry


AxisymElement :: AxisymElement(int n, Domain *aDomain) :
    PlaneElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}



double
AxisymElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForAxisymmElements(normalToCrackPlane);
}



double
AxisymElement :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    
    FloatArray N;
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        evalN( N, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );  
        
    double r = 0.0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        double x  = this->giveNode(i)->giveCoordinate(1);
        r += x * N.at(i);
    }

    double determinant = fabs(  static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );

    double weight = gp->giveWeight();
    return determinant * weight * r;

}






void
AxisymElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [6x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    double r, x;
    int size, ind = 1;
    FloatMatrix dnx;

    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    int nRows = dnx.giveNumberOfRows();
    
    if ( ui == ALL_STRAINS ) {
        size = 6;
        ui = 6;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 6 ) ) {
        OOFEM_ERROR("ComputeBmatrixAt size mismatch");
    }

    answer.resize(size, nRows*2);
    answer.zero();

    if ( ( li <= 1 ) && ( ui >= 1 ) ) {
        for ( int i = 1; i <= nRows; i++ ) {
            answer.at(ind, 2 * i - 1) = dnx.at(i, 1);
        }

        ind++;
    }

    if ( ( li <= 2 ) && ( ui >= 2 ) ) {
        for ( int i = 1; i <= nRows; i++ ) {
            answer.at(ind, 2 * i - 0) = dnx.at(i, 2);
        }

        ind++;
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        FloatArray n;
        static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
            evalN( n, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        r = 0.;
        for ( int i = 1; i <= numberOfDofMans; i++ ) {
            x  = this->giveNode(i)->giveCoordinate(1);
            r += x * n.at(i);
        }

        for ( int i = 0; i < numberOfDofMans; i++ ) {
            answer.at(ind, 2*i + 1) = n.at(i+1) / r;
        }        
        
        ind++;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        ind++;
    }

    if ( ( li <= 5 ) && ( ui >= 5 ) ) {
        ind++;
    }

    if ( ( li <= 6 ) && ( ui >= 6 ) ) {
        for ( int i = 1; i <= nRows; i++ ) {
            answer.at(ind, 2 * i - 1) = dnx.at(i, 2);
            answer.at(ind, 2 * i - 0) = dnx.at(i, 1);
        }
    }
}






void
AxisymElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
// @todo not checked if correct
{
    FloatArray n;
    FloatMatrix dnx;

    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    int nRows = dnx.giveNumberOfRows();
    answer.resize(9, nRows*2);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }


    // mode is _3dMat !!!!!! answer.at(4,*), answer.at(5,*), answer.at(7,*), and answer.at(8,*) is zero
    for ( int i = 1; i <= nRows*2; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }

    
    for ( int i = 0; i < numberOfDofMans; i++ ) {
        answer.at(3, 2*i + 1) = n.at(i+1) / r;
    }            

}


double
AxisymElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = static_cast< FEInterpolation2d* > ( this->giveInterpolation() )-> 
        edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );


    return c.at(1) * result * gp->giveWeight();
}


#if 0

void
L4Axisymm :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b, A;
    FloatArray u, Epsilon, help;
    fMode mode = domain->giveEngngModel()->giveFormulation();

    answer.resize(6);
    answer.zero();
    if ( mode == TL ) { // Total Lagrange formulation
        this->computeVectorOf({D_u, D_v}, VM_Total, tStep, u);
        // linear part of strain tensor (in vector form)

        this->computeBmatrixAt(gp, b, 1, 2);
        Epsilon.beProductOf(b, u);
        answer.at(1) = Epsilon.at(1);
        answer.at(2) = Epsilon.at(2);

        if ( numberOfFiAndShGaussPoints == 1 ) {
            //
            // if reduced integration in one gp only
            // force the evaluation of eps_fi in this gauss point
            // instead of evaluating in given gp
            //
            GaussPoint *helpGaussPoint;
            helpGaussPoint = integrationRulesArray [ 1 ]->getIntegrationPoint(0);

            this->computeBmatrixAt(helpGaussPoint, b, 3, 6);
        } else {
            OOFEM_ERROR("numberOfFiAndShGaussPoints size mismatch");
        }

        Epsilon.beProductOf(b, u);
        answer.at(3) = Epsilon.at(1);
        answer.at(6) = Epsilon.at(4);

        if ( nlGeometry ) {
            OOFEM_ERROR("only supports nlGeometry = 0");
        }
    } else if ( mode == AL ) { // actualized Lagrange formulation
        OOFEM_ERROR("unsupported mode");
    }
}


#endif







} // end namespace oofem

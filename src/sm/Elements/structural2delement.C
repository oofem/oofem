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

#include "Elements/structural2delement.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"

namespace oofem {
Structural2DElement :: Structural2DElement(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain),
    matRotation(false)
{
    cellGeometryWrapper = NULL;
}

Structural2DElement :: ~Structural2DElement()
{
    if ( cellGeometryWrapper ) delete cellGeometryWrapper;
}


IRResultType
Structural2DElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result = NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    matRotation = ir->hasField(_IFT_Structural2DElement_materialCoordinateSystem); //|| this->elemLocalCS.isNotEmpty();
    return IRRT_OK;
}


void
Structural2DElement :: postInitialize()
{
    // Element must be created before giveNumberOfNodes can be called
    StructuralElement :: postInitialize();
    this->numberOfDofMans = this->giveNumberOfNodes();
}


int 
Structural2DElement :: giveNumberOfNodes() const
{
    return this->giveInterpolation()->giveNumberOfNodes();
}


FEICellGeometry*
Structural2DElement :: giveCellGeometryWrapper() 
{
    if ( !cellGeometryWrapper ) {
        cellGeometryWrapper = new FEIElementGeometryWrapper(this);
    }

    return cellGeometryWrapper;
}
 

int 
Structural2DElement :: computeNumberOfDofs() 
{ 
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask; 
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize(); 
  
}


void 
Structural2DElement :: giveDofManDofIDMask(int inode, IntArray &answer) const 
{
    answer = {D_u, D_v};
}


double 
Structural2DElement :: computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    const FloatArray &lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( lCoords, *this->giveCellGeometryWrapper() ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness
    
    return detJ * thickness * weight; // dV
}


void 
Structural2DElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray[ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}


void
Structural2DElement :: giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords)
{
    if ( this->elemLocalCS.isNotEmpty() ) { // User specified orientation
        x = {elemLocalCS.at(1, 1), elemLocalCS.at(2, 1)};
        y = {-x(1), x(0)};
    } else {
        FloatMatrix jac;
        this->giveInterpolation()->giveJacobianMatrixAt( jac, lcoords, *this->giveCellGeometryWrapper() );
        x.beColumnOf(jac, 1); // This is {dx/dxi, dy/dxi, dz/dxi}
        x.normalize();
        y = {-x(1), x(0)};
    }
}


double
Structural2DElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}




// Edge support

void
Structural2DElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /* Returns the [2xn] shape function matrix {N} of the receiver, 
     * evaluated at the given gp.
     * {u} = {N}*{a} gives the displacements at the integration point.
     */
    
    // Evaluate the shape functions at the position of the gp. 
    FloatArray N;
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeEvalN( N, iedge, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );  
    answer.beNMatrixOf(N, 2);
}


void
Structural2DElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
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
Structural2DElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
}


double
Structural2DElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    /* Returns the line element ds associated with the given gp on the specific edge.
     * Note: The name is misleading since there is no volume to speak of in this case. 
     * The returned value is used for integration of a line integral (external forces).
     */
    double detJ = static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    return detJ * gp->giveWeight();
}


int
Structural2DElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
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
        edgeEvalNormal( normal, iEdge, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = normal.at(2);
    answer.at(1, 2) = normal.at(1);
    answer.at(2, 1) = -normal.at(1);
    answer.at(2, 2) = normal.at(2);

    return 1;
}





// Plane stress

PlaneStressElement :: PlaneStressElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
}


void 
PlaneStressElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
{

    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    
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
    /// @todo not checked if correct

    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

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
PlaneStressElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(2) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(2) * y(0) * y(1) + e(1) * y(1) * y(1),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(2) * ( x(1) * y(0) + x(0) * y(1) )
        };

        this->giveStructuralCrossSection()->giveRealStress_PlaneStress(s, gp, rotStrain, tStep);

        answer = {
            s(0) * x(0) * x(0) + 2 * s(2) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(2) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(1) * y(0) * y(1) + s(0) * x(0) * x(1) + s(2) * ( x(1) * y(0) + x(0) * y(1) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_PlaneStress(answer, gp, e, tStep);
    }
}


void
PlaneStressElement :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );

        Q = {
            { x(0) * x(0), x(1) * x(1), x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), y(0) * y(1) },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), x(1) * y(0) + x(0) * y(1) }
        };
        answer.rotatedWith(Q, 't');
    }
}





// Plane strain


PlaneStrainElement :: PlaneStrainElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
}


void 
PlaneStrainElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx, int upperIndx)
// Returns the [ 4 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    
    
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
    /// @todo not checked if correct
  
    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}


void
PlaneStrainElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(3) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(3) * y(0) * y(1) + e(1) * y(1) * y(1),
            e(2),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * ( x(1) * y(0) + x(0) * y(1) )
        };
        this->giveStructuralCrossSection()->giveRealStress_PlaneStrain(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(3) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(3) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(2),
            y(1) * ( s(3) * x(0) + s(1) * y(0) ) + x(1) * ( s(0) * x(0) + s(3) * y(0) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_PlaneStrain(answer, gp, e, tStep);
    }
}


void
PlaneStrainElement :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStrain(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        Q = {
            { x(0) * x(0), x(1) * x(1), 0, x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), 0, y(0) * y(1) },
            { 0, 0, 1, 0 },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 0, x(1) * y(0) + x(0) * y(1) }
        };

        answer.rotatedWith(Q, 't');
    }
}



// Axisymmetry

AxisymElement :: AxisymElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
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
        evalN( N, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );  

    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        double x  = this->giveNode(i)->giveCoordinate(1);
        r += x * N.at(i);
    }

    double determinant = fabs(  static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        giveTransformationJacobian( gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() ) );

    double weight = gp->giveWeight();
    return determinant * weight * r;

}


void
AxisymElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FEInterpolation *interp = this->giveInterpolation();
     
    FloatArray N;
    interp->evalN( N, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        double x = this->giveNode(i)->giveCoordinate(1);
        r += x * N.at(i);
    } 
    
    FloatMatrix dNdx;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    answer.resize(6, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);
        answer.at(3, i * 2 - 1) = N.at(i) / r;
        answer.at(6, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(6, 2 * i - 0) = dNdx.at(i, 1);
    }
  
}


void
AxisymElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
///@todo not checked if correct, is dw/dz = u/r for large deformations? /JB
{
    FloatArray n;
    FloatMatrix dnx;

    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    int nRows = dnx.giveNumberOfRows();
    answer.resize(9, nRows*2);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
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

    for ( int i = 0; i < this->giveNumberOfDofManagers(); i++ ) {
        answer.at(3, 2*i + 1) = n.at(i+1) / r;
    }

}


double
AxisymElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = static_cast< FEInterpolation2d* > ( this->giveInterpolation() )-> 
        edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );


    return c.at(1) * result * gp->giveWeight();
}


void 
AxisymElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray[ 0 ].reset( new GaussIntegrationRule(1, this, 1, 6) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

void
AxisymElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1),
            e(2),
            e(4) * y(0) + e(3) * y(1),
            e(4) * x(0) + e(3) * x(1),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) )
        };
        this->giveStructuralCrossSection()->giveRealStress_3d(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(2),
            s(4) * x(1) + s(3) * y(1),
            s(4) * x(0) + s(3) * y(0),
            y(1) * ( s(5) * x(0) + s(1) * y(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
    }
}

void
AxisymElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                 MatResponseMode rMode, GaussPoint *gp,
                                                 TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        Q = {
            { x(0) * x(0), x(1) * x(1), 0, 0, 0, x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), 0, 0, 0, y(0) * y(1) },
            { 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, y(1), y(0), 0 },
            { 0, 0, 0, x(1), x(0), 0 },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 0, 0, 0, x(1) * y(0) + x(0) * y(1) }
        };

        answer.rotatedWith(Q, 't');
    }
}

} // end namespace oofem

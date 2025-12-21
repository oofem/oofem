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

#include "sm/Elements/Shells/tr_shell11.h"
#include "fei2dtrlin.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "node.h"
#include "sm/Materials/structuralms.h"
#include "load.h"

namespace oofem {
REGISTER_Element(TR_SHELL11);

FEI2dTrLin TR_SHELL11 :: interp_lin(1, 2);                                                     

  
TR_SHELL11 :: TR_SHELL11(int n, Domain *aDomain) : StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this), ZZErrorEstimatorInterface(this), SpatialLocalizerInterface(this)
{
    numberOfDofMans = 3;
    numberOfGaussPoints = 4;
    numberOfRotGaussPoints = 1;
    area = 0.0;
}

void
TR_SHELL11 :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
  if ( integrationRulesArray.size() == 0 ) {
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
  }
}

FloatArray
TR_SHELL11 :: GiveDerivativeUX(const FloatArray &lCoords)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = lCoords.at(1);
    shapeFunct.at(2) = lCoords.at(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray nx(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        nx.at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * sin( angles.at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * sin( angles.at(k) ) );
    }

    return nx;
}


FloatArray
TR_SHELL11 :: GiveDerivativeVX(const FloatArray &lCoords)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = lCoords.at(1);
    shapeFunct.at(2) = lCoords.at(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray nx(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        nx.at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * cos( angles.at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * cos( angles.at(k) ) ) * ( -1.0 );
    }

    return nx;
}


FloatArray
TR_SHELL11 :: GiveDerivativeUY(const FloatArray &lCoords)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = lCoords.at(1);
    shapeFunct.at(2) = lCoords.at(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray ny(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        ny.at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * sin( angles.at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * sin( angles.at(k) ) );
    }

    return ny;
}


FloatArray
TR_SHELL11 :: GiveDerivativeVY(const FloatArray &lCoords)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = lCoords.at(1);
    shapeFunct.at(2) = lCoords.at(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray ny(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        ny.at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * cos( angles.at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * cos( angles.at(k) ) ) * ( -1.0 );
    }

    return ny;
}

FloatArray
TR_SHELL11 :: GivePitch()
// Returns angles between each side and global x-axis
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray angles(3);

    for ( int i = 0; i < 3; i++ ) {
        int j = i + 1 - i / 2 * 3;
        int k = j + 1 - j / 2 * 3;
        if ( x(k) == x(j) ) {
            if ( y(k) > y(j) ) {
                angles.at(i + 1) = M_PI / 2.;
            } else {
                angles.at(i + 1) = M_PI * 3. / 2.;
            }
        }

        if ( x(k) > x(j) ) {
            if ( y(k) >= y(j) ) {
                angles.at(i + 1) = atan( ( y(k) - y(j) ) / ( x(k) - x(j) ) );
            } else {
                angles.at(i + 1) = 2. * M_PI - atan( ( y(j) - y(k) ) / ( x(k) - x(j) ) );
            }
        }

        if ( x(k) < x(j) ) {
            if ( y(k) >= y(j) ) {
                angles.at(i + 1) = M_PI - atan( ( y(k) - y(j) ) / ( x(j) - x(k) ) );
            } else {
                angles.at(i + 1) = M_PI + atan( ( y(j) - y(k) ) / ( x(j) - x(k) ) );
            }
        }
    }

    return angles;
}
  
void
TR_SHELL11 :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns part of strain-displacement matrix {B} of the receiver,
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// type of this part is [3,9]  r=(u1,w1,fi1,u2,w2,fi2,u3,w3,fi3)
// evaluated at gp.
// strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy, gamma_rot}
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    FloatArray l(3), b(3), c(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
	if (i<3) {
          l.at(i) = gp->giveNaturalCoordinate(i);
        }
    }
    l.at(3) = 1.-l.at(1)-l.at(2);
    
    // Derivatives of the shape functions (for this special interpolation)
    FloatArray nx = this->GiveDerivativeUX(gp->giveNaturalCoordinates());
    FloatArray ny = this->GiveDerivativeVY(gp->giveNaturalCoordinates());

    FloatArray center = Vec2(0.33333333, 0.33333333);
    FloatArray nxRed = this->GiveDerivativeVX( center );
    FloatArray nyRed = this->GiveDerivativeUY( center );
    
    // These are the regular shape functions of a linear triangle evaluated at the center
    FloatArray lc = Vec3( center.at(1), center.at(2), 1.0 - center.at(1) - center.at(2) ); // = {0, 0, 1}
    
    double area = this->giveArea();
    double detJ = 1.0 / ( 2.0 * area );
    answer.resize(9, 18);
    for ( int i = 1; i <= 3; i++ ) {
      int j = i + 1 - i / 3 * 3;
      int k = j + 1 - j / 3 * 3;
      // eps_x
        answer.at(1, 6 * i - 5) = b.at(i) * detJ;
        answer.at(1, 6 * i - 0) = nx.at(i) * detJ;

	// eps_y
        answer.at(2, 6 * i - 4) = c.at(i) * detJ;
        answer.at(2, 6 * i - 0) = ny.at(i) * detJ;

        ///@note The third strain component is evaluated at the element center just as for reduced
        ///integration (component 4 below), should it be like this? /JB
	// eps_xy
        answer.at(3, 6 * i - 5) = c.at(i) * detJ;
        answer.at(3, 6 * i - 4) = b.at(i) * detJ;
        answer.at(3, 6 * i - 0) = ( nxRed.at(i) + nyRed.at(i) ) * detJ;

	// kappa_x
	answer.at(4, 6*i - 1) = b.at(i)*detJ;
	// kappa_y
	answer.at(5, 6*i - 2) =-c.at(i)*detJ;
	// kappa_xy
	answer.at(6, 6*i - 2) =-b.at(i)*detJ;
	answer.at(6, 6*i - 1) = c.at(i)*detJ;
	// gamma_zx
	answer.at(7, 6*i - 3) = b.at(i)*detJ;
	answer.at(7, 6*i - 2) = (-b.at(i)*b.at(k)*lc.at(j) + b.at(i)*b.at(j)*lc.at(k))*0.5*detJ;
	answer.at(7, 6*i - 1) = (-b.at(i)*c.at(k)*lc.at(j)-b.at(j)*c.at(k)*lc.at(i)+b.at(k)*c.at(j)*lc.at(i)+b.at(i)*c.at(j)*lc.at(k))*0.5*detJ+lc.at(i)*2.0*area*detJ;
	// gamma_zy
	answer.at(8, 6*i - 3) = c.at(i)*detJ;
	answer.at(8, 6*i - 2) = (-b.at(k)*c.at(i)*lc.at(j)-b.at(k)*c.at(j)*lc.at(i)+b.at(j)*c.at(k)*lc.at(i)+b.at(j)*c.at(i)*lc.at(k))*0.5*detJ-lc.at(i)*2.0*area*detJ;
	answer.at(8, 6*i - 1) = (-c.at(i)*c.at(k)*lc.at(j)+c.at(i)*c.at(j)*lc.at(k))*0.5*detJ;
	
        // Reduced integration of the fourth strain component
	// gamma_rot
        answer.at(9, 6 * i - 5) = -1. * c.at(i) * 1.0 / 4.0 / area;
        answer.at(9, 6 * i - 4) = b.at(i) * 1.0 / 4.0 / area;
        answer.at(9, 6 * i - 0) = ( -4. * area * lc.at(i) + nxRed.at(i) - nyRed.at(i) ) * 1.0 / 4.0 / area;
            
    }

}




void
TR_SHELL11 :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver
// evaluated at gp.
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    double l1, l2, l3, f11, f12, f13, f21, f22, f23;

    l1 = iLocCoord.at(1);
    l2 = iLocCoord.at(2);
    l3 = 1. - l1 - l2;

    f11 = d.at(2) / 2. *l1 *l3 *sin( angles.at(2) ) - d.at(3) / 2. *l2 *l1 *sin( angles.at(3) );
    f12 = d.at(3) / 2. *l2 *l1 *sin( angles.at(3) ) - d.at(1) / 2. *l3 *l2 *sin( angles.at(1) );
    f13 = d.at(1) / 2. *l3 *l2 *sin( angles.at(1) ) - d.at(2) / 2. *l1 *l3 *sin( angles.at(2) );

    f21 = d.at(3) / 2. *l2 *l1 *cos( angles.at(3) ) - d.at(2) / 2. *l1 *l3 *cos( angles.at(2) );
    f22 = d.at(1) / 2. *l3 *l2 *cos( angles.at(1) ) - d.at(3) / 2. *l2 *l1 *cos( angles.at(3) );
    f23 = d.at(2) / 2. *l1 *l3 *cos( angles.at(2) ) - d.at(1) / 2. *l3 *l2 *cos( angles.at(1) );

    //
    answer.resize(6, 18);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 6) = f11;
    answer.at(1, 7) = l2;
    answer.at(1, 12) = f12;
    answer.at(1, 13) = l3;
    answer.at(1, 18) = f13;

    answer.at(2, 2) = l1;
    answer.at(2, 6) = f21;
    answer.at(2, 8) = l2;
    answer.at(2, 12) = f22;
    answer.at(2, 14) = l3;
    answer.at(2, 18) = f23;

    answer.at(3, 3) = l1; // d_w
    answer.at(3, 4) = l1 * ( l2 * b.at(3) - l3 * b.at(2) ) * 0.5;
    answer.at(3, 5) = l1 * ( l2 * c.at(3) - l3 * c.at(2) ) * 0.5;
    answer.at(3, 9) = l2;
    answer.at(3, 10) = l2 * ( l3 * b.at(1) - l1 * b.at(3) ) * 0.5;
    answer.at(3, 11) = l2 * ( l3 * c.at(1) - l1 * c.at(3) ) * 0.5;
    answer.at(3, 15) = l3;
    answer.at(3, 16) = l3 * ( l1 * b.at(2) - l2 * b.at(1) ) * 0.5;
    answer.at(3, 17) = l3 * ( l1 * c.at(2) - l2 * c.at(1) ) * 0.5;

    answer.at(4, 4)  = l1; // fi_x
    answer.at(4, 10) = l2;
    answer.at(4, 16) = l3;

    answer.at(5, 5)  = l1; // fi_y
    answer.at(5, 11) = l2;
    answer.at(5, 17) = l3;

    answer.at(6, 6) = l1;  // fi_y
    answer.at(12, 12) = l2;
    answer.at(18, 18) = l3;
}


double
TR_SHELL11 :: giveArea()
// returns the area occupied by the receiver
{
   ///@todo replace with call to linear triangle interpolator 
  if ( fabs(area) > 0 ) { // check if previously computed
        return area;
    }

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    double x1, x2, x3, y1, y2, y3;
    x1 = x.at(1);
    x2 = x.at(2);
    x3 = x.at(3);

    y1 = y.at(1);
    y2 = y.at(2);
    y3 = y.at(3);

    //return ( area = fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) ) );
    return ( area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );
}

void
TR_SHELL11 :: giveLocalCoordinates(FloatArray &answer, const FloatArray &global)
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test the parameter
    if ( global.giveSize() != 3 ) {
        OOFEM_ERROR("cannot transform coordinates - size mismatch");
    }

    // first ensure that receiver's GtoLRotationMatrix[3,3] is defined
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        this->computeGtoLRotationMatrix();
    }

    FloatArray offset = global;
    offset.subtract( this->giveNode(1)->giveCoordinates() );
    answer.beProductOf(GtoLRotationMatrix, offset);
}


double
TR_SHELL11 :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    FloatArray x, y;
    this->giveNodeCoordinates(x, y);
    std :: vector< FloatArray > lc = {Vec2(x[0], y[0]), Vec2(x[1], y[1]), Vec2(x[2], y[2])};

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc, this->giveGeometryType()) ) );
    return detJ * weight; // * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
}


void
TR_SHELL11 :: giveNodeCoordinates(FloatArray &x, FloatArray &y)
{
    FloatArray nc1(3), nc2(3), nc3(3);

    this->giveLocalCoordinates( nc1, this->giveNode(1)->giveCoordinates() );
    this->giveLocalCoordinates( nc2, this->giveNode(2)->giveCoordinates() );
    this->giveLocalCoordinates( nc3, this->giveNode(3)->giveCoordinates() );

    x.resize(3);
    x.at(1) = nc1.at(1);
    x.at(2) = nc2.at(1);
    x.at(3) = nc3.at(1);

    y.resize(3);
    y.at(1) = nc1.at(2);
    y.at(2) = nc2.at(2);
    y.at(3) = nc3.at(2);

    //if (z) {
    //  z[0] = nc1->at(3);
    //  z[1] = nc2->at(3);
    //  z[2] = nc3->at(3);
    //}
}


void
TR_SHELL11 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w, R_u, R_v, R_w};
}


void
TR_SHELL11 :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  answer = this->giveStructuralCrossSection()->give3dShellRotStiffMtrx(rMode, gp, tStep);
}


void
TR_SHELL11 :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
  answer = this->giveStructuralCrossSection()->giveGeneralizedStress_ShellRot(strain, gp, tStep);
}


const FloatMatrix *
TR_SHELL11 :: computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{

  if ( !GtoLRotationMatrix.isNotEmpty() ) {
        FloatArray e1, e2, e3, help;

        // compute e1' = [N2-N1]  and  help = [N3-N1]
        e1.beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
        help.beDifferenceOf( this->giveNode(3)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );

        // let us normalize e1'
        e1.normalize();

        // compute e3' : vector product of e1' x help
        e3.beVectorProductOf(e1, help);
        // let us normalize
        e3.normalize();

        // now from e3' x e1' compute e2'
        e2.beVectorProductOf(e3, e1);

        //
        GtoLRotationMatrix.resize(3, 3);

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix.at(1, i) = e1.at(i);
            GtoLRotationMatrix.at(2, i) = e2.at(i);
            GtoLRotationMatrix.at(3, i) = e3.at(i);
        }
    }

    return &GtoLRotationMatrix;
}


bool
TR_SHELL11 :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [18,18]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v,w, r1, r2, r3} = T * {u,v,w,r1,r2,r3}
{
    // test if pereviously computed
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(18, 18);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
      answer.setSubMatrix(GtoLRotationMatrix, 3*i-2, 3*i-2);
    }

    return 1;
}


void
TR_SHELL11 :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
// returns characteristic tensor of the receiver at given gp and tStep
// strain vector = (Eps_X, Eps_y, Gamma_xy, Kappa_z)
{
    FloatArray charVect;
    StructuralMaterialStatus *ms = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

    answer.resize(3, 3);
    answer.zero();

    if ( type == LocalForceTensor || type == GlobalForceTensor ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ms->giveStressVector();

        double h = this->giveStructuralCrossSection()->give(CS_Thickness, gp);
        answer.at(1, 1) = charVect.at(1) * h;
        answer.at(2, 2) = charVect.at(2) * h;
        answer.at(1, 2) = charVect.at(3) * h;
        answer.at(2, 1) = charVect.at(3) * h;

	answer.at(1, 3) = charVect.at(7);
        answer.at(3, 1) = charVect.at(7);
        answer.at(2, 3) = charVect.at(8);
        answer.at(3, 2) = charVect.at(8);

    } else if ( type == LocalMomentTensor || type == GlobalMomentTensor ) {
        //this->computeStressVector(charVect, gp, tStep);
        charVect = ms->giveStressVector();

        //answer.at(3, 3) = charVect.at(4);
	answer.at(1, 1) = charVect.at(4);
        answer.at(2, 2) = charVect.at(5);
        answer.at(1, 2) = charVect.at(6);
        answer.at(2, 1) = charVect.at(6);

    } else if ( type == LocalStrainTensor || type == GlobalStrainTensor ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        answer.at(1, 1) = charVect.at(1);
        answer.at(2, 2) = charVect.at(2);
        answer.at(1, 2) = charVect.at(3) / 2.;
        answer.at(2, 1) = charVect.at(3) / 2.;

        answer.at(1, 3) = charVect.at(7) / 2.;
        answer.at(3, 1) = charVect.at(7) / 2.;
        answer.at(2, 3) = charVect.at(8) / 2.;
        answer.at(3, 2) = charVect.at(8) / 2.;

	
    } else if ( type == LocalCurvatureTensor || type == GlobalCurvatureTensor ) {
        //this->computeStrainVector(charVect, gp, tStep);
        charVect = ms->giveStrainVector();

        //answer.at(3, 3) = charVect.at(4);

        answer.at(1, 1) = charVect.at(4);
        answer.at(2, 2) = charVect.at(5);
        answer.at(1, 2) = charVect.at(6) / 2.;
        answer.at(2, 1) = charVect.at(6) / 2.;

    } else {
        OOFEM_ERROR("unsupported tensor mode");
    }

    if ( type == GlobalForceTensor || type == GlobalMomentTensor ||
         type == GlobalStrainTensor || type == GlobalCurvatureTensor ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(GtoLRotationMatrix);
    }
}


int
TR_SHELL11 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if ( type == IST_CurvatureTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_CurvatureTensor ) {
            cht = GlobalCurvatureTensor;
        } else {
            cht = GlobalStrainTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = 2 * globTensor.at(2, 3); //yz
        answer.at(5) = 2 * globTensor.at(1, 3); //xz
        answer.at(6) = 2 * globTensor.at(1, 2); //xy

        return 1;
    } else if ( type == IST_ShellMomentTensor || type == IST_ShellForceTensor ) {
        if ( type == IST_ShellMomentTensor ) {
            cht = GlobalMomentTensor;
        } else {
            cht = GlobalForceTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = globTensor.at(2, 3); //yz
        answer.at(5) = globTensor.at(1, 3); //xz
        answer.at(6) = globTensor.at(1, 2); //xy

        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}

int
TR_SHELL11 :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
{
    // test if previously computed
    if ( !GtoLRotationMatrix.isNotEmpty() ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize(6, 6);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = answer.at(4, i + 3) = GtoLRotationMatrix.at(1, i);
        answer.at(2, i) = answer.at(5, i + 3) = GtoLRotationMatrix.at(2, i);
        answer.at(3, i) = answer.at(6, i + 3) = GtoLRotationMatrix.at(3, i);
    }

    return 1;
}


void
TR_SHELL11 :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system; 6 components
    forLoad->computeComponentArrayAt(force, tStep, mode);
    // get it in local c.s.
    this->computeLoadGToLRotationMtrx(T);
    force.rotatedWith(T, 'n');

    if ( force.giveSize() ) {
      GaussIntegrationRule iRule (1, this, 1, 1);
      iRule.SetUpPointsOnTriangle(1, _Unknown);
      GaussPoint *gp = iRule.getIntegrationPoint(0);

        dens = this->giveStructuralCrossSection()->give('d', gp);
        dV   = this->computeVolumeAround(gp);// * this->giveCrossSection()->give(CS_Thickness, gp);

        answer.resize(18);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(7) = load;
        answer.at(13) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(8) = load;
        answer.at(14) = load;

        load = force.at(3) * dens * dV / 3.0;
        answer.at(3) = load;
        answer.at(9) = load;
        answer.at(15) = load;

        // transform result from global cs to local element cs.
        //if ( this->computeGtoLRotationMatrix(T) ) {
        //    answer.rotatedWith(T, 'n');
        //}
    } else {
        answer.clear();          // nil resultant
    }
}

void
TR_SHELL11 :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    this->computeNmatrixAt(sgp->giveNaturalCoordinates(), answer);
}

void
TR_SHELL11 :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(18);
    answer.zero();
    if ( iSurf == 1 ) {
      for (int i=1; i<=18; i++)
	answer.at(i)=i;
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}


double
TR_SHELL11 :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    return this->computeVolumeAround(gp);
}


int
TR_SHELL11 :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int isurf, GaussPoint *gp)
{
    return 0;
}


void
TR_SHELL11 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v;

    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < (int)integrationRulesArray.size(); i++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ i ] ) {

            fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
            fprintf(file, "  strains    ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            this->giveIPValue(v, gp, IST_CurvatureTensor, tStep);
            fprintf(file, "\n              curvatures ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            // Forces - Moments
            this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
            fprintf(file, "\n              stresses   ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            this->giveIPValue(v, gp, IST_ShellMomentTensor, tStep);
            fprintf(file, "\n              moments    ");
            for ( auto &val : v ) fprintf(file, " %.4e", val);

            
            if ( gp->hasSlaveGaussPoint()) { // layered material
              fprintf(file, "\n  Layers report \n{\n");

              for (auto &sgp: gp->giveSlaveGaussPoints()) {
                sgp->printOutputAt(file, tStep);
              }
              fprintf(file, "} end layers report\n");

            }
            
            
            fprintf(file, "\n");
        }
    }
}

void
TR_SHELL11 :: initializeFrom(InputRecord &ir, int priority)
{
    StructuralElement :: initializeFrom(ir, priority);
}

Interface *
TR_SHELL11 :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return static_cast< LayeredCrossSectionInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    }


    return NULL;
}

void
TR_SHELL11 :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
							 InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    if ( ( type == IST_ShellForceTensor ) || ( type == IST_ShellMomentTensor ) ||
        ( type == IST_ShellStrainTensor )  || ( type == IST_CurvatureTensor ) ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        this->giveIPValue(answer, gp, type, tStep);
    } else {
        answer.clear();
    }
}


void
TR_SHELL11 :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}


void
TR_SHELL11 :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}


SPRPatchType
TR_SHELL11 :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


//
// layered cross section support functions
//
void
TR_SHELL11 :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
					 GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep)
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
// // strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy, gamma_rot}
{
    double layerZeta, layerZCoord, top, bottom;

    top    = this->giveCrossSection()->give(CS_TopZCoord, masterGp);
    bottom = this->giveCrossSection()->give(CS_BottomZCoord, masterGp);
    layerZeta = slaveGp->giveNaturalCoordinate(3);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(5); // {Exx,Eyy,GMyz,GMzx,GMxy}

    answer.at(1) = masterGpStrain.at(1)+masterGpStrain.at(4) * layerZCoord;
    answer.at(2) = masterGpStrain.at(2)+masterGpStrain.at(5) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3)+masterGpStrain.at(6) * layerZCoord;
    answer.at(3) = masterGpStrain.at(8);
    answer.at(4) = masterGpStrain.at(7);
}

double
TR_SHELL11 :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


} // end namespace oofem

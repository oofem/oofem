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

#include "../sm/Elements/PlaneStress/trplanrot.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "verbose.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(TrPlaneStrRot);

TrPlaneStrRot :: TrPlaneStrRot(int n, Domain *aDomain) :
    TrPlaneStress2d(n, aDomain)
{
    numberOfDofMans        = 3;
    numberOfGaussPoints    = 4;
    numberOfRotGaussPoints = 1;
}


void
TrPlaneStrRot :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
TrPlaneStrRot :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns part of strain-displacement matrix {B} of the receiver,
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// type of this part is [3,9]  r=(u1,w1,fi1,u2,w2,fi2,u3,w3,fi3)
// evaluated at gp.
{
#if 1 
    // New version (13-09-2014 /JB)
    // Computes the B-matrix, directly taking into account the reduced 
    // integration of the fourth strain component.
    
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    FloatArray b(3), c(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
    }
    
    // Derivatives of the shape functions (for this special interpolation)
    FloatArray nx = this->GiveDerivativeUX(gp->giveNaturalCoordinates());
    FloatArray ny = this->GiveDerivativeVY(gp->giveNaturalCoordinates());

    FloatArray center = {0.0, 0.0};
    FloatArray nxRed = this->GiveDerivativeVX( center );
    FloatArray nyRed = this->GiveDerivativeUY( center );
    
    // These are the regular shape functions of a linear triangle evaluated at the center
    FloatArray shapeFunct = { center.at(1), center.at(2), 1.0 - center.at(1) - center.at(2) }; // = {0, 0, 1}
    
    double area = this->giveArea();
    double detJ = 1.0 / ( 2.0 * area );
    answer.resize(4, 9);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, 3 * i - 2) = b.at(i) * detJ;
        answer.at(1, 3 * i - 0) = nx.at(i) * detJ;

        answer.at(2, 3 * i - 1) = c.at(i) * detJ;
        answer.at(2, 3 * i - 0) = ny.at(i) * detJ;

        ///@note The third strain component is evaluated at the element center just as for reduced
        ///integration (component 4 below), should it be like this? /JB
        answer.at(3, 3 * i - 2) = c.at(i) * detJ;
        answer.at(3, 3 * i - 1) = b.at(i) * detJ;
        answer.at(3, 3 * i - 0) = ( nxRed.at(i) + nyRed.at(i) ) * detJ;
        
        // Reduced integration of the fourth strain component
        answer.at(4, 3 * i - 2) = -1. * c.at(i) * 1.0 / 4.0 / area;
        answer.at(4, 3 * i - 1) = b.at(i) * 1.0 / 4.0 / area;
        answer.at(4, 3 * i - 0) = ( -4. * area * shapeFunct.at(i) + nxRed.at(i) - nyRed.at(i) ) * 1.0 / 4.0 / area;
            
    }

#else
    // OLD version - commented out 13-09-2014 // JB

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
    }

    //
    double area;
    area = this->giveArea();

    //
    int size;
    if ( ui == ALL_STRAINS ) {
        size = 4;
        ui = 4;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 4 ) ) {
        OOFEM_ERROR("ComputeBmatrixAt size mismatch");
    }

    //
    int ind = 1;

    answer.resize(size, 9);

    if ( ( li <= 2 ) ) {
        FloatArray nx = this->GiveDerivativeUX(gp->giveNaturalCoordinates());
        FloatArray ny = this->GiveDerivativeVY(gp->giveNaturalCoordinates());

        if ( ( li <= 1 ) && ( ui >= 1 ) ) {
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(1, 3 * i - 2) = b.at(i) * 1. / ( 2. * area );
                answer.at(1, 3 * i - 0) = nx.at(i) * 1. / ( 2. * area );
            }

            ind++;
        }

        if ( ( li <= 2 ) && ( ui >= 2 ) ) {
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(2, 3 * i - 1) = c.at(i) * 1. / ( 2. * area );
                answer.at(2, 3 * i - 0) = ny.at(i) * 1. / ( 2. * area );
            }

            ind++;
        }
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        GaussIntegrationRule ir(1, this, 1, 3);
        ir.SetUpPointsOnTriangle(1, _PlaneStress);

        FloatArray nx = this->GiveDerivativeVX( *ir.getIntegrationPoint(0)->giveNaturalCoordinates() );
        FloatArray ny = this->GiveDerivativeUY( *ir.getIntegrationPoint(0)->giveNaturalCoordinates() );

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(3, 3 * i - 2) = c.at(i) * 1. / ( 2. * area );
            answer.at(3, 3 * i - 1) = b.at(i) * 1. / ( 2. * area );
            answer.at(3, 3 * i - 0) = ( nx.at(i) + ny.at(i) ) * 1. / ( 2. * area );
        }

        ind++;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
      if ( numberOfRotGaussPoints == 1 ) {
      //reduced integration, evaluate at coordinate (0,0)
        FloatArray center = {0.0, 0.0};
        FloatArray shapeFunct(3);

        FloatArray nx = this->GiveDerivativeVX( center );
        FloatArray ny = this->GiveDerivativeUY( center );

        shapeFunct.at(1) = center.at(1);
        shapeFunct.at(2) = center.at(2);
        shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = -1. * c.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 1) = b.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 0) = ( -4. * area * shapeFunct.at(i) + nx.at(i) - ny.at(i) ) * 1.0 / 4.0 / area;
        }
      }
    }


#endif


}


void
TrPlaneStrRot :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
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
    answer.resize(3, 9);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 3) = f11;
    answer.at(1, 4) = l2;
    answer.at(1, 6) = f12;
    answer.at(1, 7) = l3;
    answer.at(1, 9) = f13;

    answer.at(2, 2) = l1;
    answer.at(2, 3) = f21;
    answer.at(2, 5) = l2;
    answer.at(2, 6) = f22;
    answer.at(2, 8) = l3;
    answer.at(2, 9) = f23;

    answer.at(3, 3) = l1;
    answer.at(3, 6) = l2;
    answer.at(3, 9) = l3;
}


double
TrPlaneStrRot :: giveArea()
// returns the area occupied by the receiver
{
   ///@todo replace with call to linear triangle interpolator 
    if ( area > 0 ) { // check if previously computed
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

    return ( area = fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) ) );
    //    return ( area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );
}


void
TrPlaneStrRot :: giveNodeCoordinates(FloatArray &x, FloatArray &y)
{
    FloatArray *nc1, *nc2, *nc3;
    nc1 = this->giveNode(1)->giveCoordinates();
    nc2 = this->giveNode(2)->giveCoordinates();
    nc3 = this->giveNode(3)->giveCoordinates();

    x.at(1) = nc1->at(1);
    x.at(2) = nc2->at(1);
    x.at(3) = nc3->at(1);

    y.at(1) = nc1->at(2);
    y.at(2) = nc2->at(2);
    y.at(3) = nc3->at(2);

    //if (z) {
    //  z[0] = nc1->at(3);
    //  z[1] = nc2->at(3);
    //  z[2] = nc3->at(3);
    //}
}


FloatArray
TrPlaneStrRot :: GivePitch()
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


FloatArray
TrPlaneStrRot :: GiveDerivativeUX(const FloatArray &lCoords)
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
TrPlaneStrRot :: GiveDerivativeVX(const FloatArray &lCoords)
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
TrPlaneStrRot :: GiveDerivativeUY(const FloatArray &lCoords)
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
TrPlaneStrRot :: GiveDerivativeVY(const FloatArray &lCoords)
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


IRResultType
TrPlaneStrRot :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    numberOfGaussPoints = 4;
    result = TrPlaneStress2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    numberOfRotGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRotGaussPoints, _IFT_TrPlaneStrRot_niprot);

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 7 ) ) ) {
        numberOfGaussPoints = 4;
    }

    // According to the implementation of the B-matrix only one gp is supported, 
    // so it shouldn't be an optional parameter //JB
    if ( numberOfRotGaussPoints != 1 ) {
        OOFEM_ERROR("numberOfRotGaussPoints size mismatch - must be equal to one");
    }
    
    return IRRT_OK;
}


void
TrPlaneStrRot :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    cs->giveMembraneRotStiffMtrx(answer, rMode, gp, tStep);
}


void
TrPlaneStrRot :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    cs->giveGeneralizedStress_MembraneRot(answer, gp, strain, tStep);
}


void
TrPlaneStrRot :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, R_w};
}


void
TrPlaneStrRot :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
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

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);

    if ( force.giveSize() ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveStructuralCrossSection()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);

        answer.resize(9);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(4) = load;
        answer.at(7) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(5) = load;
        answer.at(8) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.clear();          // nil resultant
    }
}


double
TrPlaneStrRot :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


int
TrPlaneStrRot :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_StressTensor ) {
        const FloatArray &help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        answer = {help.at(1), help.at(2), 0., 0., 0.,  help.at(3)};
        return 1;
    } else if ( type == IST_StrainTensor ) {
        const FloatArray &help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        answer = {help.at(1), help.at(2), 0., 0., 0.,  help.at(3)};
        return 1;
    } else {
        return TrPlaneStress2d :: giveIPValue(answer, gp, type, tStep);
    }
}

} // end namespace oofem

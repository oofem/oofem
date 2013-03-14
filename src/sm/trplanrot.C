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

#include "trplanrot.h"

#include "node.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "verbose.h"
#include "engngm.h"
#include "load.h"
#include "mathfem.h"

namespace oofem {
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
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 2;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _PlaneStressRot);

        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 4, 4);
        integrationRulesArray [ 1 ]->setUpIntegrationPoints(_Triangle, numberOfRotGaussPoints, _PlaneStressRot);
    }
}


void
TrPlaneStrRot :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns part of strain-displacement matrix {B} of the receiver,
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// type of this part is [3,9]  r=(u1,w1,fi1,u2,w2,fi2,u3,w3,fi3)
// evaluated at aGaussPoint.
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
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
        _error("ComputeBmatrixAt size mismatch");
    }

    //
    int ind = 1;
    FloatArray *nx, *ny;

    answer.resize(size, 9);

    if ( ( li <= 2 ) ) {
        nx = this->GiveDerivativeUX(aGaussPoint);
        ny = this->GiveDerivativeVY(aGaussPoint);

        if ( ( li <= 1 ) && ( ui >= 1 ) ) {
            for ( i = 1; i <= 3; i++ ) {
                answer.at(ind, 3 * i - 2) = b.at(i) * 1. / ( 2. * area );
                answer.at(ind, 3 * i - 0) = nx->at(i) * 1. / ( 2. * area );
            }

            ind++;
        }

        if ( ( li <= 2 ) && ( ui >= 2 ) ) {
            for ( i = 1; i <= 3; i++ ) {
                answer.at(ind, 3 * i - 1) = c.at(i) * 1. / ( 2. * area );
                answer.at(ind, 3 * i - 0) = ny->at(i) * 1. / ( 2. * area );
            }

            ind++;
        }

        delete nx;
        delete ny;
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        GaussIntegrationRule ir(1, this, 1, 3);
        ir.setUpIntegrationPoints(_Triangle, 1, _PlaneStress);

        nx = this->GiveDerivativeVX( ir.getIntegrationPoint(0) );
        ny = this->GiveDerivativeUY( ir.getIntegrationPoint(0) );

        for ( i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = c.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 1) = b.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 0) = ( nx->at(i) + ny->at(i) ) * 1. / ( 2. * area );
        }

        ind++;
        delete nx;
        delete ny;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        FloatArray shapeFunct(3);

        nx = this->GiveDerivativeVX(aGaussPoint);
        ny = this->GiveDerivativeUY(aGaussPoint);

        shapeFunct.at(1) = aGaussPoint->giveCoordinate(1);
        shapeFunct.at(2) = aGaussPoint->giveCoordinate(2);
        shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

        for ( i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = -1. * c.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 1) = b.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 0) = ( -4. * area * shapeFunct.at(i) + nx->at(i) - ny->at(i) ) * 1.0 / 4.0 / area;
        }

        delete nx;
        delete ny;
    }
}


void
TrPlaneStrRot :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver
// evaluated at aGaussPoint.
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray *angles;
    angles = this->GivePitch();

    //
    double l1, l2, l3, f11, f12, f13, f21, f22, f23;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1. - l1 - l2;

    f11 = d.at(2) / 2. *l1 *l3 *sin( angles->at(2) ) - d.at(3) / 2. *l2 *l1 *sin( angles->at(3) );
    f12 = d.at(3) / 2. *l2 *l1 *sin( angles->at(3) ) - d.at(1) / 2. *l3 *l2 *sin( angles->at(1) );
    f13 = d.at(1) / 2. *l3 *l2 *sin( angles->at(1) ) - d.at(2) / 2. *l1 *l3 *sin( angles->at(2) );

    f21 = d.at(3) / 2. *l2 *l1 *cos( angles->at(3) ) - d.at(2) / 2. *l1 *l3 *cos( angles->at(2) );
    f22 = d.at(1) / 2. *l3 *l2 *cos( angles->at(1) ) - d.at(3) / 2. *l2 *l1 *cos( angles->at(3) );
    f23 = d.at(2) / 2. *l1 *l3 *cos( angles->at(2) ) - d.at(1) / 2. *l3 *l2 *cos( angles->at(1) );

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

    delete angles;
}


double
TrPlaneStrRot :: giveArea()
// returns the area occupied by the receiver
{
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


FloatArray *
TrPlaneStrRot :: GivePitch()
// Returns angles between each side and global x-axis
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray *angles;
    angles = new FloatArray(3);

    for ( i = 0; i < 3; i++ ) {
        j = i + 1 - i / 2 * 3;
        k = j + 1 - j / 2 * 3;
        if ( x(k) == x(j) ) {
            if ( y(k) > y(j) ) {
                angles->at(i + 1) = 3.14159265358979 / 2.;
            } else {
                angles->at(i + 1) = 3.14159265358979 * 3. / 2.;
            }
        }

        if ( x(k) > x(j) ) {
            if ( y(k) >= y(j) ) {
                angles->at(i + 1) = atan( ( y(k) - y(j) ) / ( x(k) - x(j) ) );
            } else {
                angles->at(i + 1) = 2. * 3.14159265358979 - atan( ( y(j) - y(k) ) / ( x(k) - x(j) ) );
            }
        }

        if ( x(k) < x(j) ) {
            if ( y(k) >= y(j) ) {
                angles->at(i + 1) = 3.14159265358979 - atan( ( y(k) - y(j) ) / ( x(j) - x(k) ) );
            } else {
                angles->at(i + 1) = 3.14159265358979 + atan( ( y(j) - y(k) ) / ( x(j) - x(k) ) );
            }
        }
    }

    return angles;
}


FloatArray *
TrPlaneStrRot :: GiveDerivativeUX(GaussPoint *aGaussPoint)
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray *angles;
    angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct.at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray *nx;
    nx = new FloatArray(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        nx->at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * sin( angles->at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * sin( angles->at(k) ) );
    }

    delete angles;
    return nx;
}


FloatArray *
TrPlaneStrRot :: GiveDerivativeVX(GaussPoint *aGaussPoint)
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray *angles;
    angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct.at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray *nx;
    nx = new FloatArray(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        nx->at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * cos( angles->at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * cos( angles->at(k) ) ) * ( -1.0 );
    }

    delete angles;
    return nx;
}


FloatArray *
TrPlaneStrRot :: GiveDerivativeUY(GaussPoint *aGaussPoint)
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray *angles;
    angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct.at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray *ny;
    ny = new FloatArray(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        ny->at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * sin( angles->at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * sin( angles->at(k) ) );
    }

    delete angles;
    return ny;
}


FloatArray *
TrPlaneStrRot :: GiveDerivativeVY(GaussPoint *aGaussPoint)
{
    int i, j, k;

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray *angles;
    angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct.at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray *ny;
    ny = new FloatArray(3);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        ny->at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * cos( angles->at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * cos( angles->at(k) ) ) * ( -1.0 );
    }

    delete angles;
    return ny;
}


IRResultType
TrPlaneStrRot :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_Element_nip, "nip");

    numberOfRotGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRotGaussPoints, IFT_TrPlaneStrRot_niprot, "niprot");

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 7 ) ) ) {
        numberOfGaussPoints = 4;
    }

    if ( !( ( numberOfRotGaussPoints == 1 ) ||
           ( numberOfRotGaussPoints == 4 ) ||
           ( numberOfRotGaussPoints == 7 ) ) ) {
        numberOfRotGaussPoints = 1;
    }

    return IRRT_OK;
}


void
TrPlaneStrRot :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u, Epsilon;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);

    answer.resize(4);
    answer.zero();

    this->computeBmatrixAt(gp, b, 1, 3);
    Epsilon.beProductOf(b, u);
    answer.at(1) = Epsilon.at(1);
    answer.at(2) = Epsilon.at(2);
    answer.at(3) = Epsilon.at(3);

    if ( numberOfRotGaussPoints == 1 ) {
        //
        // if reduced integration in one gp only
        // force the evaluation of eps_fi in this gauss point
        // instead of evaluating in given gp
        //
        GaussPoint *helpGaussPoint;
        helpGaussPoint = integrationRulesArray [ 1 ]->getIntegrationPoint(0);

        this->computeBmatrixAt(helpGaussPoint, b, 4, 4);
    } else {
        _error("ComputeStrainVector: numberOfRotGaussPoints size mismatch");
    }

    Epsilon.beProductOf(b, u);
    answer.at(4) = Epsilon.at(1);
}


void
TrPlaneStrRot :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_u, D_v, R_w);
}


void
TrPlaneStrRot :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at stepN.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    GaussPoint *gp = NULL;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);

    if ( force.giveSize() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveMaterial()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

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
        answer.resize(0);          // nil resultant
    }
}
} // end namespace oofem

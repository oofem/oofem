/* $Header: /home/cvs/bp/oofem/sm/src/trplanrot.C,v 1.3 2003/04/06 14:08:32 bp Exp $ */
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

//   file TRPLANROT.CC
//
//  triangular element with rotational degrees of freedom
//  for plane stress
//
// 5.5.1995
//

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
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

TrPlaneStrRot :: TrPlaneStrRot(int n, Domain *aDomain) :
    TrPlaneStress2d(n, aDomain)
    // Constructor.
{
    numberOfDofMans        = 3;
}

TrPlaneStrRot :: ~TrPlaneStrRot()
// destructor
{ }

FloatArray *
TrPlaneStrRot :: GivePitch()
// Returns angles between each side and global x-axis
{
    int i, j, k;
    double x [ 3 ], y [ 3 ];
    FloatArray *angles;

    angles = new FloatArray(3);

    for ( i = 1; i <= 3; i++ ) {
        x [ i - 1 ] = this->giveNode(i)->giveCoordinate(1);
        y [ i - 1 ] = this->giveNode(i)->giveCoordinate(2);
    }

    for ( i = 0; i < 3; i++ ) {
        j = i + 1 - i / 2 * 3;
        k = j + 1 - j / 2 * 3;
        if ( x [ k ] == x [ j ] ) {
            if ( y [ k ] > y [ j ] ) {
                angles->at(i + 1) = 3.14159265358979 / 2.;
            } else {
                angles->at(i + 1) = 3.14159265358979 * 3. / 2.;
            }
        }

        if ( x [ k ] > x [ j ] ) {
            if ( y [ k ] >= y [ j ] ) {
                angles->at(i + 1) = atan( ( y [ k ] - y [ j ] ) / ( x [ k ] - x [ j ] ) );
            } else {
                angles->at(i + 1) = 2. * 3.14159265358979 - atan( ( y [ j ] - y [ k ] ) / ( x [ k ] - x [ j ] ) );
            }
        }

        if ( x [ k ] < x [ j ] ) {
            if ( y [ k ] >= y [ j ] ) {
                angles->at(i + 1) = 3.14159265358979 - atan( ( y [ k ] - y [ j ] ) / ( x [ j ] - x [ k ] ) );
            } else {
                angles->at(i + 1) = 3.14159265358979 + atan( ( y [ j ] - y [ k ] ) / ( x [ j ] - x [ k ] ) );
            }
        }
    }

    return angles;
}



FloatArray *
TrPlaneStrRot :: GiveDerivativeUX(GaussPoint *aGaussPoint)
{
    int i, j, k;
    FloatArray *angles, *shapeFunct, *nx, *b, *c, *d;

    shapeFunct = new FloatArray(3);
    nx = new FloatArray(3);
    b = new FloatArray(3);
    c = new FloatArray(3);
    d = new FloatArray(3);

    angles = this->GivePitch();

    shapeFunct->at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct->at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct->at(3) = 1.0 - shapeFunct->at(1) - shapeFunct->at(2);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b->at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c->at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
        d->at(i) = sqrt( b->at(i) * b->at(i) + c->at(i) * c->at(i) );
    }


    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        nx->at(i) = ( d->at(j) / 2. * ( b->at(k) * shapeFunct->at(i) + shapeFunct->at(k) * b->at(i) ) * sin( angles->at(j) ) -
                     d->at(k) / 2. * ( b->at(i) * shapeFunct->at(j) + shapeFunct->at(i) * b->at(j) ) * sin( angles->at(k) ) );
    }

    delete angles;
    delete shapeFunct;
    delete b;
    delete c;
    delete d;
    return nx;
}



FloatArray *
TrPlaneStrRot :: GiveDerivativeVX(GaussPoint *aGaussPoint)
{
    int i, j, k;
    FloatArray *shapeFunct, *angles, *nx, *b, *c, *d;

    shapeFunct = new FloatArray(3);
    nx = new FloatArray(3);
    b  = new FloatArray(3);
    c  = new FloatArray(3);
    d  = new FloatArray(3);

    angles = this->GivePitch();

    shapeFunct->at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct->at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct->at(3) = 1.0 - shapeFunct->at(1) - shapeFunct->at(2);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b->at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c->at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
        d->at(i) = sqrt( b->at(i) * b->at(i) + c->at(i) * c->at(i) );
    }

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        nx->at(i) = ( d->at(j) / 2. * ( b->at(k) * shapeFunct->at(i) + shapeFunct->at(k) * b->at(i) ) * cos( angles->at(j) ) -
                     d->at(k) / 2. * ( b->at(i) * shapeFunct->at(j) + shapeFunct->at(i) * b->at(j) ) * cos( angles->at(k) ) ) * ( -1.0 );
    }

    delete angles;
    delete shapeFunct;
    delete b;
    delete c;
    delete d;
    return nx;
}



FloatArray *
TrPlaneStrRot :: GiveDerivativeUY(GaussPoint *aGaussPoint)
{
    int i, j, k;
    FloatArray *ny, *b, *c, *d, *shapeFunct, *angles;

    shapeFunct = new FloatArray(3);
    ny = new FloatArray(3);
    b  = new FloatArray(3);
    c  = new FloatArray(3);
    d  = new FloatArray(3);

    angles = this->GivePitch();

    shapeFunct->at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct->at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct->at(3) = 1.0 - shapeFunct->at(1) - shapeFunct->at(2);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b->at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c->at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
        d->at(i) = sqrt( b->at(i) * b->at(i) + c->at(i) * c->at(i) );
    }


    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        ny->at(i) = ( d->at(j) / 2. * ( c->at(k) * shapeFunct->at(i) + shapeFunct->at(k) * c->at(i) ) * sin( angles->at(j) ) -
                     d->at(k) / 2. * ( c->at(i) * shapeFunct->at(j) + shapeFunct->at(i) * c->at(j) ) * sin( angles->at(k) ) );
    }

    delete angles;
    delete shapeFunct;
    delete b;
    delete c;
    delete d;
    return ny;
}



FloatArray *
TrPlaneStrRot :: GiveDerivativeVY(GaussPoint *aGaussPoint)
{
    int i, j, k;
    FloatArray *ny, *shapeFunct, *angles, *b, *c, *d;

    shapeFunct = new FloatArray(3);
    ny = new FloatArray(3);
    b  = new FloatArray(3);
    c  = new FloatArray(3);
    d  = new FloatArray(3);

    angles = this->GivePitch();

    shapeFunct->at(1) = aGaussPoint->giveCoordinate(1);
    shapeFunct->at(2) = aGaussPoint->giveCoordinate(2);
    shapeFunct->at(3) = 1.0 - shapeFunct->at(1) - shapeFunct->at(2);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b->at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c->at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
        d->at(i) = sqrt( b->at(i) * b->at(i) + c->at(i) * c->at(i) );
    }

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        ny->at(i) = ( d->at(j) / 2. * ( c->at(k) * shapeFunct->at(i) + shapeFunct->at(k) * c->at(i) ) * cos( angles->at(j) ) -
                     d->at(k) / 2. * ( c->at(i) * shapeFunct->at(j) + shapeFunct->at(i) * c->at(j) ) * cos( angles->at(k) ) ) * ( -1.0 );
    }

    delete angles;
    delete shapeFunct;
    delete b;
    delete c;
    delete d;
    return ny;
}



void
TrPlaneStrRot :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns part of strain-displacement matrix {B} of the receiver,
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// type of this part is [3,9]  r=(u1,w1,fi1,u2,w2,fi2,u3,w3,fi3)
// evaluated at aGaussPoint.
{
    int i, j, k, size, ind = 1;
    double area;
    FloatArray *nx, *ny;
    // FloatMatrix *gm;
    // GaussPoint  *helpGaussPoint;


    FloatArray b(3);
    FloatArray c(3);
    area = this->giveArea();

    if ( ui == ALL_STRAINS ) {
        size = 4;
        ui = 4;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 4 ) ) {
        _error("ComputeBmatrixAt size mismatch");
    }

    // gm = new FloatMatrix (size,9);
    answer.resize(size, 9);

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c.at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
    }

    if ( ( li <= 2 ) ) {
        nx   = this->GiveDerivativeUX(aGaussPoint);
        ny   = this->GiveDerivativeVY(aGaussPoint);


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
        //helpGaussPoint     = new GaussPoint(this,1,coord,2.0,_PlaneStress);

        nx = this->GiveDerivativeVX( ir.getIntegrationPoint(0) );
        ny = this->GiveDerivativeUY( ir.getIntegrationPoint(0) );

        for ( i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = c.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 1) = b.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 0) = ( nx->at(i) + ny->at(i) ) * 1. / ( 2. * area );
        }

        ind++;
        //delete helpGaussPoint;
        delete nx;
        delete ny;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        FloatArray *shapeFunct;
        shapeFunct = new FloatArray(3);

        nx = this->GiveDerivativeVX(aGaussPoint);
        ny = this->GiveDerivativeUY(aGaussPoint);

        shapeFunct->at(1) = aGaussPoint->giveCoordinate(1);
        shapeFunct->at(2) = aGaussPoint->giveCoordinate(2);
        shapeFunct->at(3) = 1.0 - shapeFunct->at(1) - shapeFunct->at(2);

        for ( i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = -1. * c.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 1) = b.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 0) = ( -4. * area * shapeFunct->at(i) + nx->at(i) - ny->at(i) ) * 1.0 / 4.0 / area;
        }

        delete shapeFunct;
        delete nx;
        delete ny;
    }

    // delete b;  delete c;

    return;
}



void
TrPlaneStrRot :: computeGaussPoints()
//  Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 2;
    integrationRulesArray = new IntegrationRule * [ 2 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _PlaneStressRot);

    integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 4, 4);
    integrationRulesArray [ 1 ]->setUpIntegrationPoints(_Triangle, numberOfRotGaussPoints, _PlaneStressRot);
}

void
TrPlaneStrRot :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver
// evaluated at aGaussPoint.
{
    int i, j, k;
    double l1, l2, l3, f11, f12, f13, f21, f22, f23;
    FloatArray b(3), c(3), d(3), * angles;
    // FloatMatrix  *answer;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1. - l1 - l2;

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;
        b.at(i) = this->giveNode(j)->giveCoordinate(2) - this->giveNode(k)->giveCoordinate(2);
        c.at(i) = this->giveNode(k)->giveCoordinate(1) - this->giveNode(j)->giveCoordinate(1);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    angles = this->GivePitch();

    f11 = d.at(2) / 2. *l1 *l3 *sin( angles->at(2) ) - d.at(3) / 2. *l2 *l1 *sin( angles->at(3) );
    f12 = d.at(3) / 2. *l2 *l1 *sin( angles->at(3) ) - d.at(1) / 2. *l3 *l2 *sin( angles->at(1) );
    f13 = d.at(1) / 2. *l3 *l2 *sin( angles->at(1) ) - d.at(2) / 2. *l1 *l3 *sin( angles->at(2) );

    f21 = d.at(3) / 2. *l2 *l1 *cos( angles->at(3) ) - d.at(2) / 2. *l1 *l3 *cos( angles->at(2) );
    f22 = d.at(1) / 2. *l3 *l2 *cos( angles->at(1) ) - d.at(3) / 2. *l2 *l1 *cos( angles->at(3) );
    f23 = d.at(2) / 2. *l1 *l3 *cos( angles->at(2) ) - d.at(1) / 2. *l3 *l2 *cos( angles->at(1) );

    //answer = new FloatMatrix (3,9);
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
    return;
}

IRResultType
TrPlaneStrRot :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_TrPlaneStrRot_nip, "nip"); // Macro

    numberOfRotGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRotGaussPoints, IFT_TrPlaneStrRot_niprot, "niprot"); // Macro

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

    this->computeGaussPoints();
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
    if ( this->updateRotationMatrix() ) {
        u.rotatedWith(this->rotationMatrix, 'n');
    }

    answer.resize(4);
    answer.zero();

    this->computeBmatrixAt(gp, b, 1, 3);
    Epsilon.beProductOf(b, u);
    answer.at(1) = Epsilon.at(1);
    answer.at(2) = Epsilon.at(2);
    answer.at(3) = Epsilon.at(3);
    // delete b;  delete Epsilon;

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
    // delete b;  delete Epsilon;

    // delete u;
    return;
}


void
TrPlaneStrRot ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    //IntArray* answer = new IntArray (3);
    answer.resize(3);

    answer.at(1) = D_u;
    answer.at(2) = D_v;
    answer.at(3) = R_w;

    return;
}



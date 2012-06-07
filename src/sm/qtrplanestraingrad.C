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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "qtrplanestraingrad.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <math.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

QTrPlaneStrainGrad :: QTrPlaneStrainGrad(int n, Domain *aDomain) : QTrPlaneStrain( n,aDomain),GradDpElement()
// Constructor.
{
    nPrimNodes = 6;
    nPrimVars = 2;
    nSecNodes = 3;
    nSecVars = 1;
    totalSize = nPrimVars*nPrimNodes+nSecVars*nSecNodes;
    locSize   = nPrimVars*nPrimNodes;
    nlSize    = nSecVars*nSecNodes;
}


void
QTrPlaneStrainGrad ::   giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{

    if ( inode<=nSecNodes ) {
        answer.setValues(3, D_u, D_v, G_0);
    }
    else {
        answer.setValues(2, D_u, D_v);
    }
}
IRResultType
QTrPlaneStrainGrad :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                 // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QPlaneStrain_nip, "nip"); // Macro


    this->computeGaussPoints();
    return IRRT_OK;
    this->computeGaussPoints();
    return IRRT_OK;
}

void
QTrPlaneStrainGrad :: computeGaussPoints()
{

    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule* [numberOfIntegrationRules];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _PlaneStrainGrad);
    }
}

void
QTrPlaneStrainGrad :: computeNkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double l1, l2, l3;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = l2;
    answer.at(1, 3) = l3;
}

void
QTrPlaneStrainGrad :: computeBkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the [2x6] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3, area;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(2, 3);

    answer.at(1, 1) = y2 - y3;
    answer.at(1, 2) = y3 - y1;
    answer.at(1, 3) = y1 - y2;
    answer.at(2, 1) = x3 - x2;
    answer.at(2, 2) = x1 - x3;
    answer.at(2, 3) = x2 - x1;


    answer.times( 1. / ( 2. * area ) );
}

}


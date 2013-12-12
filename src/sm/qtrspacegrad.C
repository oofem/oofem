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

#include "qtrspacegrad.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "fei3dtetlin.h"
#include "classfactory.h"

#include <cstdio>

namespace oofem {
REGISTER_Element(QTRSpaceGrad);

FEI3dTetLin QTRSpaceGrad :: interpolation;

QTRSpaceGrad :: QTRSpaceGrad(int n, Domain *aDomain) :  QTRSpace(n, aDomain), GradDpElement()
    // Constructor.
{
    nPrimNodes = 10;
    nPrimVars = 3;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
}


IRResultType
QTRSpaceGrad :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = this->NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    return IRRT_OK;
}


void
QTRSpaceGrad :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( inode <= nSecNodes ) {
        answer.resize(4);
        answer.at(1) = D_u;
        answer.at(2) = D_v;
        answer.at(3) = D_w;
        answer.at(4) = G_0;
    } else {
        answer.resize(3);
        answer.at(1) = D_u;
        answer.at(2) = D_v;
        answer.at(3) = D_w;
    }
}

void
QTRSpaceGrad :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 7);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}


void
QTRSpaceGrad :: computeNkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatArray n(4);
    this->interpolation.evalN( n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(1, 4);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
QTRSpaceGrad :: computeBkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
{
    FloatMatrix dnx;
    answer.resize(3, 4);
    answer.zero();

    this->interpolation.evaldNdx( dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
        answer.at(3, i) = dnx.at(i, 3);
    }
}



void
QTRSpaceGrad :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i)
// Returns the [45x45] nonlinear part of strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint

{
    FloatMatrix dnx;

    // compute the derivatives of shape functions
    this->interpolation.evaldNdx( dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(30, 30);
    answer.zero();

    // put the products of derivatives of shape functions into the "nonlinear B matrix",
    // depending on parameter i, which is the number of the strain component
    if ( i <= 3 ) {
        for ( int k = 0; k < 10; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 30; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(i, k + 1) * dnx.at(i, ( j - 1 ) / 3 + 1);
                }
            }
        }
    } else if ( i == 4 ) {
        for ( int k = 0; k < 10; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 30; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(2, k + 1) * dnx.at(3, ( j - 1 ) / 3 + 1) + dnx.at(3, k + 1) * dnx.at(2, ( j - 1 ) / 3 + 1);
                }
            }
        }
    } else if ( i == 5 ) {
        for ( int k = 0; k < 10; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 30; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(1, k + 1) * dnx.at(3, ( j - 1 ) / 3 + 1) + dnx.at(3, k + 1) * dnx.at(1, ( j - 1 ) / 3 + 1);
                }
            }
        }
    } else if ( i == 6 ) {
        for ( int k = 0; k < 10; k++ ) {
            for ( int l = 0; l < 3; l++ ) {
                for ( int j = 1; j <= 30; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(1, k + 1) * dnx.at(2, ( j - 1 ) / 3 + 1) + dnx.at(2, k + 1) * dnx.at(1, ( j - 1 ) / 3 + 1);
                }
            }
        }
    }
}
}

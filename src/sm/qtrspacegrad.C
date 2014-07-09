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
QTRSpaceGrad :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= nSecNodes ) {
        answer = {D_u, D_v, D_w, G_0};
    } else {
        answer = {D_u, D_v, D_w};
    }
}

void
QTRSpaceGrad :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 7);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}


void
QTRSpaceGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 1);
}

void
QTRSpaceGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    answer.resize(3, 4);
    answer.zero();

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
        answer.at(3, i) = dnx.at(i, 3);
    }
}

}

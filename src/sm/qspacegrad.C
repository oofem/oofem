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

#include "qspacegrad.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QSpaceGrad);

FEI3dHexaLin QSpaceGrad :: interpolation;

QSpaceGrad :: QSpaceGrad(int n, Domain *aDomain) :  QSpace(n, aDomain), GradDpElement()
    // Constructor.
{
    nPrimNodes = 8;
    nPrimVars = 2;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
}


IRResultType
QSpaceGrad :: initializeFrom(InputRecord *ir)
{
    IRResultType result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( ( numberOfGaussPoints != 8 ) && ( numberOfGaussPoints != 14 ) && ( numberOfGaussPoints != 27 ) && ( numberOfGaussPoints != 64 ) ) {
        numberOfGaussPoints = 27;
    }

    return IRRT_OK;
}


void
QSpaceGrad :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= nSecNodes ) {
        answer = {D_u, D_v, D_w, G_0};
    } else {
        answer = {D_u, D_v, D_w};
    }
}


void
QSpaceGrad :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 7);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}



void
QSpaceGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 1);
}

void
QSpaceGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    IntArray a(8);
    for ( int i = 1; i < 9; i++ ) {
        a.at(i) = dofManArray.at(i);
    }

    answer.resize(3, 8);
    answer.zero();

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
        answer.at(3, i) = dnx.at(i, 3);
    }
}

}

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

#include "qplanestraingrad.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "engngm.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(QPlaneStrainGrad);

FEI2dQuadLin QPlaneStrainGrad :: interpolation(1, 2);

QPlaneStrainGrad :: QPlaneStrainGrad(int n, Domain *aDomain) : QPlaneStrain(n, aDomain), GradDpElement()
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


void
QPlaneStrainGrad :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const

{
    if ( inode <= nSecNodes ) {
        answer.setValues(3, D_u, D_v, G_0);
    } else {
        answer.setValues(2, D_u, D_v);
    }
}


IRResultType
QPlaneStrainGrad :: initializeFrom(InputRecord *ir)
{
    IRResultType result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    numberOfGaussPoints = 4;

    if ( !( ( numberOfGaussPoints == 1 ) ||
            ( numberOfGaussPoints == 4 ) ||
            ( numberOfGaussPoints == 9 ) ||
            ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    return IRRT_OK;
}

void
QPlaneStrainGrad :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QPlaneStrainGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    int i;
    FloatArray n(4);

    answer.resize(1, 4);
    answer.zero();
    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
QPlaneStrainGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [1x8] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 4);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, i) = dnx.at(i, 1);
        answer.at(2, i) = dnx.at(i, 2);
    }
}
}

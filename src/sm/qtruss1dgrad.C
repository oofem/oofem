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

#include "qtruss1dgrad.h"
#include "fei1dlin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(QTruss1dGrad);

FEI1dLin QTruss1dGrad :: interpolation(1);

QTruss1dGrad :: QTruss1dGrad(int n, Domain *aDomain) : QTruss1d(n, aDomain), GradDpElement()
    // Constructor.
{
    nPrimNodes = 3;
    nPrimVars = 1;
    nSecNodes = 2;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
}


void
QTruss1dGrad :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( inode < 3 ) {
        answer.setValues(2, D_u, G_0);
    } else   {
        answer.setValues(1, D_u);
    }
}


IRResultType
QTruss1dGrad :: initializeFrom(InputRecord *ir)
{
    IRResultType result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    return IRRT_OK;
}


void
QTruss1dGrad :: computeGaussPoints()
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 1);
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}


void
QTruss1dGrad :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    GradDpElement :: computeStiffnessMatrix(answer, rMode, tStep);
}


void
QTruss1dGrad :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    GradDpElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
}


void
QTruss1dGrad :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    GradDpElement :: computeForceLoadVector(answer, tStep, mode);
}


void
QTruss1dGrad :: computeNkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatArray n;
    answer.resize(1, 2);
    answer.zero();

    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
}


void
QTruss1dGrad :: computeBkappaMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    answer.resize(1, 2);
    answer.zero();
    FloatMatrix b;
    this->interpolation.evaldNdx( b, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.at(1, 1) = b.at(1, 1);
    answer.at(1, 2) = b.at(2, 1);
}
}

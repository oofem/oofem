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

#include "../sm/Elements/GradientDamage/Bars/qtruss1dgraddamage.h"
#include "fei1dlin.h"
#include "fei1dquad.h"
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
REGISTER_Element(QTruss1dGradDamage);

FEI1dLin QTruss1dGradDamage :: interpolation_lin(1);

QTruss1dGradDamage :: QTruss1dGradDamage(int n, Domain *aDomain) : QTruss1d(n, aDomain), GradientDamageElement()
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
QTruss1dGradDamage :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode < 3 ) {
        answer = {D_u, G_0};
    } else {
        answer = {D_u};
    }
}


  void
QTruss1dGradDamage :: giveDofManDofIDMask_u(IntArray &answer) const
{
  answer = {D_u};
}


void
QTruss1dGradDamage :: giveDofManDofIDMask_d(IntArray &answer) const
{
  /*
    if ( inode <= 3 ) {
      answer = {G_0};
    } else {
      answer = {};
    }
  */
}

  


void
QTruss1dGradDamage :: initializeFrom(InputRecord &ir, int priority)
{
    StructuralElement :: initializeFrom(ir, priority);
}


void
QTruss1dGradDamage :: computeGaussPoints()
{
    integrationRulesArray.resize( 1 );
    integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 1) );
    this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
}


void
QTruss1dGradDamage :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    GradientDamageElement :: computeStiffnessMatrix(answer, rMode, tStep);
}


void
QTruss1dGradDamage :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    GradientDamageElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
}


void
QTruss1dGradDamage :: computeNdMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    this->interpolation_lin.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


void
QTruss1dGradDamage :: computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    this->interpolation_lin.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dnx);
}

void QTruss1dGradDamage :: computeField(ValueModeType mode, TimeStep* tStep, const FloatArray& lcoords, FloatArray& answer)
{
    FloatArray n, unknown;
    this->interpolation_lin.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->computeVectorOf({D_u}, mode, tStep, unknown);
    answer.at(1) = n.dotProduct(unknown);

    this->interpolation.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->computeVectorOf({G_0}, mode, tStep, unknown);
    answer.at(2) = n.dotProduct(unknown);
}

}

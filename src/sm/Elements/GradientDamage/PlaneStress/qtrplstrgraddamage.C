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

#include "Elements/GradientDamage/PlaneStress/qtrplstrgraddamage.h"
#include "fei2dtrlin.h"
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
REGISTER_Element(QTrPlaneStressGradDamage);

FEI2dTrLin QTrPlaneStressGradDamage :: interpolation_lin(1, 2);

QTrPlaneStressGradDamage :: QTrPlaneStressGradDamage(int n, Domain *aDomain) : QTrPlaneStress2d(n, aDomain), GradientDamageElement()
    // Constructor.
{
    nPrimNodes = 6;
    nPrimVars = 2;
    nSecNodes = 3;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
    numberOfGaussPoints = 4;
}


void
QTrPlaneStressGradDamage :: giveDofManDofIDMask(int inode, IntArray &answer) const

{
    if ( inode <= 3 ) {
        answer = {D_u, D_v, G_0};
    } else {
        answer = {D_u, D_v};
    }
}


  void
QTrPlaneStressGradDamage :: giveDofManDofIDMask_u(IntArray &answer) const
{
  answer = {D_u, D_v};
}


void
QTrPlaneStressGradDamage :: giveDofManDofIDMask_d(IntArray &answer) const
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
QTrPlaneStressGradDamage :: initializeFrom(InputRecord &ir, int priority)
{
    QTrPlaneStress2d :: initializeFrom(ir, priority);
}


void
QTrPlaneStressGradDamage :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QTrPlaneStressGradDamage :: computeNdMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    this->interpolation_lin.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

void
QTrPlaneStressGradDamage :: computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    this->interpolation_lin.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dnx);
}
}

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

#include "../sm/Elements/GradientDamage/PlaneStress/qplanestressgraddamage.h"
#include "fei2dquadlin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QPlaneStressGradDamage);

FEI2dQuadLin QPlaneStressGradDamage :: interpolation_lin(1, 2);
IntArray QPlaneStressGradDamage :: locationArray_u = {1,2, 4,5, 7,8, 10,11, 13,14,15,16,17,18,19,20};
IntArray QPlaneStressGradDamage :: locationArray_d = {3,6,9,12};


  
QPlaneStressGradDamage :: QPlaneStressGradDamage(int n, Domain *aDomain) : QPlaneStress2d(n, aDomain), GradientDamageElement()
    // Constructor.
{
    nPrimNodes = 8;
    nPrimVars = 2;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars * nPrimNodes + nSecVars * nSecNodes;
    locSize   = nPrimVars * nPrimNodes;
    nlSize    = nSecVars * nSecNodes;
    numberOfGaussPoints = 4;
}


void
QPlaneStressGradDamage :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 4 ) {
        answer = {D_u, D_v, G_0};
    } else {
        answer = {D_u, D_v};
    }
}



  void
QPlaneStressGradDamage :: giveDofManDofIDMask_u(IntArray &answer) const
{
  answer = {D_u, D_v};
}


void
QPlaneStressGradDamage :: giveDofManDofIDMask_d(IntArray &answer) const
{
  /*
  if ( inode <= 4 ) {
    answer = {G_0};
  } else {
    answer = {};
  }
  */
}

  


void
QPlaneStressGradDamage :: initializeFrom(InputRecord &ir, int priority)
{
    QPlaneStress2d :: initializeFrom(ir, priority);
}


void
QPlaneStressGradDamage :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QPlaneStressGradDamage :: computeNdMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    this->interpolation_lin.evalN(answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

void
QPlaneStressGradDamage :: computeBdMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;
    this->interpolation_lin.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dnx);
}

void
QPlaneStressGradDamage :: giveLocationArray_u(IntArray &answer)
{
  answer = locationArray_u;
}

  
void
QPlaneStressGradDamage :: giveLocationArray_d(IntArray &answer)
{
  answer = locationArray_d;
}


void
QPlaneStressGradDamage :: postInitialize()
{
  GradientDamageElement:: postInitialize();
  QPlaneStress2d :: postInitialize();
}

  
}

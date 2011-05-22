/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.C,v 1.6 2003/04/06 14:08:32 bp Exp $ */
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

//   file Truss1d.C

#include "qplanestraingrad.h"
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
 
  FEI2dQuadLin QPlaneStrainGrad :: interpolation(1, 2);	

  QPlaneStrainGrad :: QPlaneStrainGrad(int n, Domain *aDomain) : QPlaneStrain( n,aDomain),GradDpElement()
 // Constructor.
{
  nPrimNodes = 8; 
  nPrimVars = 2;
  nSecNodes = 4;
  nSecVars = 1;
  totalSize = nPrimVars*nPrimNodes+nSecVars*nSecNodes;
  locSize   = nPrimVars*nPrimNodes;
  nlSize    = nSecVars*nSecNodes;
  
    
}

 
void
QPlaneStrainGrad ::   giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const

{

  if(inode<=nSecNodes)
    {
      answer.resize(3);
      answer.at(1) = D_u;
      answer.at(2) = D_v;
      answer.at(3) = G_0;
    }
  else
    {
      answer.resize(2);
      answer.at(1) = D_u;
      answer.at(2) = D_v;
    }
  return;
}
IRResultType
QPlaneStrainGrad :: initializeFrom(InputRecord *ir)
{
  //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
  //IRResultType result;                 // Required by IR_GIVE_FIELD macro
  this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QPlaneStrain_nip, "nip"); // Macro

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}

void 
QPlaneStrainGrad :: computeGaussPoints()
{

  if ( !integrationRulesArray ) {
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule* [numberOfIntegrationRules];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStrainGrad);
  }
}

void
QPlaneStrainGrad :: computeNkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatArray n(4);

    answer.resize(1, 4);
    answer.zero();
    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1,i) = n.at(i);
    }

    return;
}

  void
  QPlaneStrainGrad :: computeBkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
  // Returns the [1x8] strain-displacement matrix {B} of the receiver, eva-
  // luated at aGaussPoint.
  {
    int i;
    FloatMatrix dnx;
    
    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);
    
    answer.resize(2, 4);
    answer.zero();
    
    for ( i = 1; i <= 4; i++ ) {
      answer.at(1, i) = dnx.at(i, 1);
      answer.at(2, i) = dnx.at(i, 2);
    }
    
    
    
  }
}


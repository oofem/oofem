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

#include "qtrspacegrad.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "cltypes.h"
#include "structuralms.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "fei3dtetlin.h"

#include <cstdio>

namespace oofem {

FEI3dTetLin QTRSpaceGrad :: interpolation;

QTRSpaceGrad :: QTRSpaceGrad (int n, Domain* aDomain) :  QTRSpace(n, aDomain),GradDpElement()
// Constructor.
{
    nPrimNodes = 10; 
    nPrimVars = 3;
    nSecNodes = 4;
    nSecVars = 1;
    totalSize = nPrimVars*nPrimNodes+nSecVars*nSecNodes;
    locSize   = nPrimVars*nPrimNodes;
    nlSize    = nSecVars*nSecNodes;
}


IRResultType
QTRSpaceGrad :: initializeFrom (InputRecord* ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->NLStructuralElement :: initializeFrom (ir);
    IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_Element_nip, "nip");
        numberOfGaussPoints = 4;

    return IRRT_OK;
}


void
QTRSpaceGrad :: giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const
{
    if ( inode <= nSecNodes ) {
        answer.resize (4);
        answer.at(1) = D_u;
        answer.at(2) = D_v;
        answer.at(3) = D_w; 
        answer.at(4) = G_0;
    } else {
        answer.resize (3);
        answer.at(1) = D_u;
        answer.at(2) = D_v;
        answer.at(3) = D_w; 
    }
}

void
QTRSpaceGrad :: computeGaussPoints ()
  // Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule* [numberOfIntegrationRules];
    integrationRulesArray[0] = new GaussIntegrationRule (1,this,1, 7);
    MaterialMode mode = _3dMatGrad; // material model is based on strain (standard approach)
    if ( nlGeometry > 1 )
        mode = _3dMatGrad_F; // material model is based on deformation gradient, not on strain
    integrationRulesArray[0]->setUpIntegrationPoints (_Tetrahedra, numberOfGaussPoints, mode);
}


void
QTRSpaceGrad :: computeNkappaMatrixAt (GaussPoint* aGaussPoint,FloatMatrix& answer)
  // Returns the displacement interpolation matrix {N} of the receiver, eva-
  // luated at aGaussPoint.
{
    FloatArray n(4);
    this->interpolation.evalN (n, *aGaussPoint->giveCoordinates(),FEIElementGeometryWrapper(this));
    answer.resize(1,4);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
QTRSpaceGrad :: computeBkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix& answer)
{
    FloatMatrix dnx;
    answer.resize(3,4);
    answer.zero();

    this->interpolation.evaldNdx (dnx, *aGaussPoint->giveCoordinates(),FEIElementGeometryWrapper(this));
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = dnx.at(i,1);
        answer.at(2, i) = dnx.at(i,2);
        answer.at(3, i) = dnx.at(i,3);
    }
}

}

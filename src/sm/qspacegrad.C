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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "qspacegrad.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

namespace oofem {

FEI3dHexaLin QSpaceGrad :: interpolation;

QSpaceGrad :: QSpaceGrad (int n, Domain* aDomain) :  QSpace(n, aDomain), GradDpElement()
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


IRResultType
QSpaceGrad :: initializeFrom (InputRecord* ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom (ir);
    IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_QSpaceGrad_nip, "nip"); // Macro

    if ((numberOfGaussPoints != 8) && (numberOfGaussPoints != 14) && (numberOfGaussPoints != 27) && (numberOfGaussPoints != 64)) numberOfGaussPoints = 27;

    // set - up Gaussian integration points
    this -> computeGaussPoints();

    return IRRT_OK;
}


void
QSpaceGrad :: giveDofManDofIDMask (int inode, EquationID ut, IntArray& answer) const
{
    if ( inode<=nSecNodes ) {
        answer.setValues(4, D_u, D_v, D_w, G_0);
    } else {
        answer.setValues(3, D_u, D_v, D_w);
    }
}

void
QSpaceGrad :: computeGaussPoints ()
// Sets up the array containing the four Gauss points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule* [numberOfIntegrationRules];
    integrationRulesArray[0] = new GaussIntegrationRule (1,this,1, 7);
    integrationRulesArray[0]->setUpIntegrationPoints (_Cube, numberOfGaussPoints, _3dMatGrad);
}


void
QSpaceGrad :: computeNkappaMatrixAt (GaussPoint* aGaussPoint,FloatMatrix& answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatArray n(8);
    this->interpolation.evalN (n, *aGaussPoint->giveCoordinates(),FEIElementGeometryWrapper(this));
    answer.resize(1,8);
    answer.zero();

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
QSpaceGrad :: computeBkappaMatrixAt(GaussPoint *aGaussPoint, FloatMatrix& answer)
{
    FloatMatrix dnx;
    IntArray a(8);
    for( int i = 1; i < 9; i++ ) {
        a.at(i) = dofManArray.at(i);
    }
    answer.resize(3,8);
    answer.zero();

    this->interpolation.evaldNdx (dnx, *aGaussPoint->giveCoordinates(),FEIElementGeometryWrapper(this));
    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = dnx.at(i,1);
        answer.at(2, i) = dnx.at(i,2);
        answer.at(3, i) = dnx.at(i,3);
    }

}


}

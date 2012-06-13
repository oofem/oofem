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

#include "qtruss1d.h"
#include "fei1dquad.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

FEI1dQuad QTruss1d :: interpolation(1);

QTruss1d :: QTruss1d(int n, Domain *aDomain) : StructuralElement(n, aDomain)
// Constructor.
{
    numberOfDofMans = 3;
}


IRResultType
QTruss1d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    this->StructuralElement :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QTruss1d_nip, "nip"); // Macro

    this->computeGaussPoints();
    return IRRT_OK;
}

void
QTruss1d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(1, D_u);
}


double
QTruss1d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = fabs( this->interpolation.giveTransformationJacobian(* gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    double weight  = gp->giveWeight();
    return  detJ * weight * this->giveCrossSection()->give(CS_Area);
}


void QTruss1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.XF
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, numberOfGaussPoints, _1dMat);
    }
}

void
QTruss1d :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    answer.resize(1,3);
    answer.zero();

    this->interpolation.evalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    answer.at(1,1) = n.at(1);
    answer.at(1,2) = n.at(2);
    answer.at(1,3) = n.at(3);
}


void
QTruss1d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    answer.resize(1,3);
    answer.zero();

    this->interpolation.evaldNdx(answer, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}

} // end namespace oofem

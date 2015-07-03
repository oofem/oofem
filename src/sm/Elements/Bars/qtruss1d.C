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

#include "../sm/Elements/Bars/qtruss1d.h"
#include "fei1dquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(QTruss1d);

FEI1dQuad QTruss1d :: interpolation(1);

QTruss1d :: QTruss1d(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 3;
}


IRResultType
QTruss1d :: initializeFrom(InputRecord *ir)
{
    return StructuralElement :: initializeFrom(ir);
}

void
QTruss1d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u};
}

int
QTruss1d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interpolation.local2global( answer, lcoords, FEIElementGeometryWrapper(this) );
    return 1;
}

void
QTruss1d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_1d(answer, gp, strain, tStep);
}

double
QTruss1d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    double weight  = gp->giveWeight();
    return detJ *weight *this->giveCrossSection()->give(CS_Area, gp);
}


void QTruss1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.XF
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QTruss1d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    this->interpolation.evaldNdx( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

void
QTruss1d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    this->computeBmatrixAt(gp, answer);
}
} // end namespace oofem

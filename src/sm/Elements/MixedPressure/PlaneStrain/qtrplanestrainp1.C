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

#include "../sm/Elements/MixedPressure/PlaneStrain/qtrplanestrainp1.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"



namespace oofem {
REGISTER_Element(QTrPlaneStrainP1);

FEI2dTrLin QTrPlaneStrainP1 :: interpolation_lin(1, 2);

QTrPlaneStrainP1 :: QTrPlaneStrainP1(int n, Domain *aDomain) : QTrPlaneStrain(n, aDomain), BaseMixedPressureElement()
{
    displacementDofsOrdering = {
        1, 2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15
    };
    pressureDofsOrdering = {
        3, 6, 9
    };
}


void
QTrPlaneStrainP1 :: computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &answer, NLStructuralElement *elem)
{
    answer.resize(12);
    FloatMatrix dN;
    elem->giveInterpolation()->evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
        answer(k)     = dN(j, 0);
        answer(k + 1) = dN(j, 1);
    }
}

void
QTrPlaneStrainP1 :: computePressureNMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    this->interpolation_lin.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


void
QTrPlaneStrainP1 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 3 ) {
        answer = {
            D_u, D_v, P_f
        };
    } else {
        answer = {
            D_u, D_v
        };
    }
}



void
QTrPlaneStrainP1 :: giveDofManDofIDMask_u(IntArray &answer)
{
    answer = {
        D_u, D_v
    };
}


void
QTrPlaneStrainP1 :: giveDofManDofIDMask_p(IntArray &answer)
{
    answer = {
        P_f
    };
}



void
QTrPlaneStrainP1 ::  postInitialize()
{
    BaseMixedPressureElement :: postInitialize();
    QTrPlaneStrain :: postInitialize();
}
} // end namespace oofem

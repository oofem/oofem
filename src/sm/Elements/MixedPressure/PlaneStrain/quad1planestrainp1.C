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

#include "Elements/MixedPressure/PlaneStrain/quad1planestrainp1.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(Quad1PlaneStrainP1);


Quad1PlaneStrainP1 :: Quad1PlaneStrainP1(int n, Domain *aDomain) : Quad1PlaneStrain(n, aDomain), BaseMixedPressureElement()
{
    displacementDofsOrdering = {
        1, 2, 4, 5, 7, 8, 10, 11
    };
    pressureDofsOrdering = {
        3, 6, 9, 12
    };
}



void
Quad1PlaneStrainP1 :: computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &answer, NLStructuralElement *elem)
{
    answer.resize(8);
    FloatMatrix dN;
    elem->giveInterpolation()->evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int j = 0, k = 0; j < 4; j++, k += 2 ) {
        answer(k)     = dN(j, 0);
        answer(k + 1) = dN(j, 1);
    }
}

void
Quad1PlaneStrainP1 :: computePressureNMatrixAt(GaussPoint *gp, FloatArray &answer)
{
    NLStructuralElement *elem = this->giveElement();
    elem->giveInterpolation()->evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


void
Quad1PlaneStrainP1 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v, P_f
    };
}


void
Quad1PlaneStrainP1 :: giveDofManDofIDMask_u(IntArray &answer)
{
    answer = {
        D_u, D_v
    };
}


void
Quad1PlaneStrainP1 :: giveDofManDofIDMask_p(IntArray &answer)
{
    answer = {
        P_f
    };
}

void
Quad1PlaneStrainP1 ::  postInitialize()
{
    BaseMixedPressureElement :: postInitialize();
    Quad1PlaneStrain :: postInitialize();
}
} // end namespace oofem

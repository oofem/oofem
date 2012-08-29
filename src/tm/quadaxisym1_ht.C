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

#include "quadaxisym1_ht.h"
#include "node.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
QuadAxisym1_ht :: QuadAxisym1_ht(int n, Domain *aDomain) : Quad1_ht(n, aDomain)
{ }

QuadAxisym1_hmt :: QuadAxisym1_hmt(int n, Domain *aDomain) : QuadAxisym1_ht(n, aDomain)
{
    this->emode = HeatMass1TransferEM; // This could be done in a better way.
}

QuadAxisym1_ht :: ~QuadAxisym1_ht()
// Destructor
{ }

double
QuadAxisym1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );

    weight = aGaussPoint->giveWeight();
    volume = determinant * weight * this->computeRadiusAt(aGaussPoint);

    return volume;
}

double
QuadAxisym1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double radius;
    FloatArray gcoords;
    this->interpolation.edgeLocal2global(gcoords, iEdge, *gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this));
    radius = gcoords.at(1);
    
    double detJ = fabs( this->interpolation.edgeGiveTransformationJacobian(iEdge, * gp->giveCoordinates(),
                                                                     FEIElementGeometryWrapper(this)) );
    return detJ*gp->giveWeight()*radius;
}

double
QuadAxisym1_ht :: computeRadiusAt(GaussPoint *gp)
{
    FloatArray gcoords;
    this->interpolation.local2global(gcoords, *gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this));
    return gcoords.at(1);
}
} // end namespace oofem

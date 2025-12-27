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

#include "tm/Elements/quadaxisym1_ht.h"
#include "fei2dquadlin.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(QuadAxisym1_ht);
REGISTER_Element(QuadAxisym1_hmt);
REGISTER_Element(QuadAxisym1_mt);

QuadAxisym1_ht :: QuadAxisym1_ht(int n, Domain *aDomain) : Quad1_ht(n, aDomain)
{ }

QuadAxisym1_hmt :: QuadAxisym1_hmt(int n, Domain *aDomain) : QuadAxisym1_ht(n, aDomain)
{
    this->emode = HeatMass1TransferEM; // This could be done in a better way.
}

QuadAxisym1_mt :: QuadAxisym1_mt(int n, Domain *aDomain) : QuadAxisym1_ht(n, aDomain)
{
    this->emode = Mass1TransferEM;
}

double
QuadAxisym1_ht :: computeVolumeAround(GaussPoint *gp)
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );

    double weight = gp->giveWeight();
    return determinant * weight * this->computeRadiusAt(gp);
}

double
QuadAxisym1_ht :: giveThicknessAt(const FloatArray &gcoords)
{
    return gcoords.at(1);
}

double
QuadAxisym1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray gcoords;
    this->interpolation.edgeLocal2global( gcoords, iEdge, gp->giveSubPatchCoordinates(), FEIElementGeometryWrapper(this) );
    double radius = gcoords.at(1);

    double detJ = fabs( this->interpolation.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(),
                                                                           FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight() * radius;
}

double
QuadAxisym1_ht :: computeRadiusAt(GaussPoint *gp)
{
    FloatArray gcoords;
    this->interpolation.local2global( gcoords, gp->giveSubPatchCoordinates(), FEIElementGeometryWrapper(this) );
    return gcoords.at(1);
}

std::unique_ptr<IntegrationRule>
QuadAxisym1_ht :: giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
    return this->giveInterpolation()->giveBoundaryEdgeIntegrationRule(order+1, boundary, this->giveGeometryType());
}

std::unique_ptr<IntegrationRule>
QuadAxisym1_ht :: giveBoundarySurfaceIntegrationRule(int order, int boundary)
{
    return this->giveInterpolation()->giveBoundarySurfaceIntegrationRule(order+1, boundary, this->giveGeometryType());
}
} // end namespace oofem

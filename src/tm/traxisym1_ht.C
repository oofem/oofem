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

#include "traxisym1_ht.h"
#include "fei2dtrlin.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(TrAxisym1_ht);

TrAxisym1_ht :: TrAxisym1_ht(int n, Domain *aDomain) : Tr1_ht(n, aDomain)
    // Constructor.
{ }

TrAxisym1_ht :: ~TrAxisym1_ht()
// Destructor
{ }

double
TrAxisym1_ht :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight;
    determinant = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    return determinant *weight *this->computeRadiusAt(gp); ///@todo What about 2*pi ?
}

double
TrAxisym1_ht :: giveThicknessAt(const FloatArray &gcoords)
{
    return gcoords.at(1);
}

double
TrAxisym1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray gcoords;
    double determinant, radius;

    determinant = fabs( this->interp.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    this->interp.edgeLocal2global( gcoords, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    radius = gcoords.at(1);
    return determinant *radius *gp->giveWeight();
}

double
TrAxisym1_ht :: computeRadiusAt(GaussPoint *gp)
{
    FloatArray gcoords;
    this->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
    return gcoords.at(1);
}
} // end namespace oofem

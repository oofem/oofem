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

#include "trwarp.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(Tr_Warp);

FEI2dTrLin Tr_Warp :: interp(1, 2);

Tr_Warp :: Tr_Warp(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 1;
}

Tr_Warp :: ~Tr_Warp()
// Destructor
{ }


void
Tr_Warp :: computeGaussPoints()
// Sets up the array containing the Gauss point of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
Tr_Warp :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 1;
    this->StructuralElement :: initializeFrom(ir);

    return IRRT_OK;
}


void
Tr_Warp :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                    int li, int ui)
// Returns the [3x6] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 3);

    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(2, 1);
    answer.at(1, 3) = dN.at(3, 1);

    answer.at(2, 1) = dN.at(1, 2);
    answer.at(2, 2) = dN.at(2, 2);
    answer.at(2, 3) = dN.at(3, 2);
}


double
Tr_Warp :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}


double
Tr_Warp :: giveThicknessAt(const FloatArray &gcoords)
{
    return 1.;
}


double
Tr_Warp :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant = fabs( this->interp.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    FloatArray gc;
    return determinant * gp->giveWeight();
}


void
Tr_Warp :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer = {D_w};
}

Interface *
Tr_Warp :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
	//    } else if ( interface == EIPrimaryFieldInterfaceType ) {
	//        return static_cast< EIPrimaryFieldInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}

int
Tr_Warp :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Tr_Warp :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
}
} // end namespace oofem

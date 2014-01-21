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

#include "tr1_ht.h"
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
REGISTER_Element(Tr1_ht);
REGISTER_Element(Tr1_hmt);

FEI2dTrLin Tr1_ht :: interp(1, 2);

Tr1_ht :: Tr1_ht(int n, Domain *aDomain) :
    TransportElement(n, aDomain, HeatTransferEM)
    // Constructor.
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 1;
}

Tr1_hmt :: Tr1_hmt(int n, Domain *aDomain) : Tr1_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

Tr1_ht :: ~Tr1_ht()
// Destructor
{ }


void
Tr1_ht :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
Tr1_ht :: initializeFrom(InputRecord *ir)
{
    //const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    numberOfGaussPoints = 1;
    this->TransportElement :: initializeFrom(ir);

    //IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, _IFT_Element_nip);
    //if ( numberOfGaussPoints != 1) numberOfGaussPoints = 1;

    return IRRT_OK;
}


double
Tr1_ht :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight * this->giveCrossSection()->give(CS_Thickness, gp);

    return volume;
}


double
Tr1_ht :: giveThicknessAt(const FloatArray &gcoords)
{
    return this->giveCrossSection()->give(CS_Thickness, & gcoords, this, false);
}


double
Tr1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant = fabs( this->interp.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    FloatArray gc;
    double thick = this->giveCrossSection()->give(CS_Thickness, gp->giveCoordinates(), NULL, this); // 't'
    return determinant * thick * gp->giveWeight();
}


Interface *
Tr1_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >( this );
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return static_cast< EIPrimaryFieldInterface * >( this );
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    }

    return NULL;
}

int
Tr1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Tr1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
}
} // end namespace oofem

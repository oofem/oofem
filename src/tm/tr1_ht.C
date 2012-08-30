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

#include "tr1_ht.h"
#include "node.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

namespace oofem {

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
    MaterialMode mmode;

    if ( emode == HeatTransferEM ) {
        mmode = _2dHeat;
    } else {
        mmode = _2dHeMo;
    }

    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, mmode);
    }
}


IRResultType
Tr1_ht :: initializeFrom(InputRecord *ir)
{
    //const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);

    numberOfGaussPoints = 1;
    //IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, "nip"); // Macro
    //if ( numberOfGaussPoints != 1) numberOfGaussPoints = 1;

    this->computeGaussPoints();
    return IRRT_OK;
}


double
Tr1_ht :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian(* gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    weight = gp->giveWeight();
    volume = determinant * weight * this->giveCrossSection()->give(CS_Thickness);

    return volume;
}


double
Tr1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant, thick;
    determinant = fabs( this->interp.edgeGiveTransformationJacobian(iEdge, *gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    thick = this->giveCrossSection()->give(CS_Thickness);
    return determinant *thick *gp->giveWeight();
}


Interface *
Tr1_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return ( EIPrimaryFieldInterface * ) this;
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    }

    return NULL;
}

int
Tr1_ht :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
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

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

#include "qbrick1_ht.h"
#include "node.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "load.h"

namespace oofem {
FEI3dHexaQuad QBrick1_ht :: interpolation;

QBrick1_ht :: QBrick1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM), SpatialLocalizerInterface(), ZZNodalRecoveryModelInterface(), SPRNodalRecoveryModelInterface()
{
    numberOfDofMans  = 20;
    numberOfGaussPoints = 27;
}

QBrick1_hmt :: QBrick1_hmt(int n, Domain *aDomain) : QBrick1_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

QBrick1_ht :: ~QBrick1_ht()
{ }


void
QBrick1_ht :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    MaterialMode mmode;

    if ( emode == HeatTransferEM ) {
        mmode = _3dHeat;
    } else {
        mmode = _3dHeMo;
    }

    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Cube, numberOfGaussPoints, mmode);
    }
}


IRResultType
QBrick1_ht :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->TransportElement :: initializeFrom(ir);
    numberOfGaussPoints = 27;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_Element_nip, "nip");

    if ( ( numberOfGaussPoints != 8 ) && ( numberOfGaussPoints != 27 ) && ( numberOfGaussPoints != 64 ) ) {
        numberOfGaussPoints = 27;
    }

    return IRRT_OK;
}


double
QBrick1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );

    weight = aGaussPoint->giveWeight();
    volume = determinant * weight;
    return volume;
}


double
QBrick1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * gp->giveWeight();
}


IntegrationRule *
QBrick1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->setUpIntegrationPoints(_Square, npoints, _Unknown);
    return iRule;
}


double
QBrick1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
}


Interface *
QBrick1_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >( this );
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return static_cast< EIPrimaryFieldInterface * >( this );
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    }

    return NULL;
}

int
QBrick1_ht :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

int
QBrick1_ht :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    return ZZNodalRecoveryMI_giveDofManRecordSize(type);
}

void
QBrick1_ht :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(numberOfDofMans);
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QBrick1_ht :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        _error2("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: unknown node number %d", pap);
    }
}

int
QBrick1_ht :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
QBrick1_ht :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


int
QBrick1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
QBrick1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.zero();
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
}


} // end namespace oofem

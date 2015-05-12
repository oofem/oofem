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

#include "qbrick1_ht.h"
#include "fei3dhexaquad.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "load.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QBrick1_ht);
REGISTER_Element(QBrick1_hmt);

FEI3dHexaQuad QBrick1_ht :: interpolation;

QBrick1_ht :: QBrick1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this), SPRNodalRecoveryModelInterface()
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


FEInterpolation *
QBrick1_ht :: giveInterpolation() const { return & interpolation; }


void
QBrick1_ht :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
QBrick1_ht :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 27;
    return TransportElement :: initializeFrom(ir);
}


double
QBrick1_ht :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
}


double
QBrick1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}


IntegrationRule *
QBrick1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}


double
QBrick1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
}


Interface *
QBrick1_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return static_cast< EIPrimaryFieldInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    }

    return NULL;
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
        OOFEM_ERROR("unknown node number %d", pap);
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

} // end namespace oofem

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

#include "tm/Elements/qquad1_ht.h"
#include "fei2dquadquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(QQuad1_ht);
REGISTER_Element(QQuad1_hmt);
REGISTER_Element(QQuad1_mt);

FEI2dQuadQuad QQuad1_ht :: interpolation(1, 2);

QQuad1_ht :: QQuad1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this), SPRNodalRecoveryModelInterface()
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

QQuad1_hmt :: QQuad1_hmt(int n, Domain *aDomain) : QQuad1_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

QQuad1_mt :: QQuad1_mt(int n, Domain *aDomain) : QQuad1_ht(n, aDomain)
{
    emode = Mass1TransferEM;
}


FEInterpolation *
QQuad1_ht :: giveInterpolation() const { return & interpolation; }

void
QQuad1_ht :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
QQuad1_ht :: initializeFrom(InputRecord &ir)
{
    //numberOfGaussPoints = 4;
    TransportElement :: initializeFrom(ir);
}


double
QQuad1_ht :: computeVolumeAround(GaussPoint *gp)
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // 't'
    return determinant * gp->giveWeight() * thickness;
}


double
QQuad1_ht :: giveThicknessAt(const FloatArray &gcoords)
{
    return this->giveCrossSection()->give(CS_Thickness, gcoords, this, false);
}


double
QQuad1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    FloatArray gc;
    this->interpolation.edgeLocal2global( gc, iEdge, gp->giveNaturalCoordinates(),
                                         FEIElementGeometryWrapper(this) );
    // temporary gauss point on element (not edge) to evaluate thickness
    GaussPoint _gp( NULL, 1, gc, 1.0, gp->giveMaterialMode() );
    double thick = this->giveCrossSection()->give(CS_Thickness, & _gp); // 't'
    return result *thick *gp->giveWeight();
}

Interface *
QQuad1_ht :: giveInterface(InterfaceType interface)
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

    return nullptr;
}

void
QQuad1_ht :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(8);
    for ( int i = 1; i <= 8; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QQuad1_ht :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 8; i++ ) {
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
QQuad1_ht :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
QQuad1_ht :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dquadratic;
}


void
QQuad1_ht :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}




} // end namespace oofem

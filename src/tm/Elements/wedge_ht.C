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


#include "tm/Elements/wedge_ht.h"
#include "fei3dwedgelin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "load.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(Wedge_ht);

FEI3dWedgeLin Wedge_ht :: interpolation;

Wedge_ht :: Wedge_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this), SPRNodalRecoveryModelInterface()
{
    numberOfDofMans = 6;
}

Wedge_hmt :: Wedge_hmt(int n, Domain *aDomain) : Wedge_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

Wedge_mt :: Wedge_mt(int n, Domain *aDomain) : Wedge_ht(n, aDomain)
{
    emode = Mass1TransferEM;
}


void
Wedge_ht :: initializeFrom(InputRecord &ir, int priority)
{
    numberOfGaussPoints = 6;
    TransportElement :: initializeFrom(ir, priority);

}


void
Wedge_ht :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


  double
Wedge_ht :: computeVolumeAround(GaussPoint *gp)
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );

    double weight = gp->giveWeight();
    return determinant * weight;
}



double
Wedge_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result * gp->giveWeight();
}
  
  
FEInterpolation * Wedge_ht :: giveInterpolation() const { return & interpolation; }



  

Interface *
Wedge_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }

    OOFEM_LOG_INFO("Interface on Lwedge element not supported");
    return nullptr;
}

void
Wedge_ht :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(6);
    for ( int i = 1; i <= 6; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Wedge_ht :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 6; i++ ) {
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
Wedge_ht :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
Wedge_ht :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}


void
Wedge_ht :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

} // end namespace oofem

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

FEInterpolation *Tr1_ht :: giveInterpolation(DofIDItem id)
{
    if (id == T_f) {
        return &this->interp;
    } else {
        return NULL;
    }
}

void
Tr1_ht :: computeNSubMatrixAt(FloatMatrix &answer, FloatArray *coords)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    ///@todo Deal with matrix and vector (I find that the row-wise matrices should be transposed).
    //this->interp.evalN(answer, coords, FEIElementGeometryWrapper(this), 0.0);
    double l1, l2, l3;

    l1 = coords->at(1);
    l2 = coords->at(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = l2;
    answer.at(1, 3) = l3;
}

void
Tr1_ht :: computeNmatrixAt(FloatMatrix &answer, FloatArray *coords)
{
    if ( emode == HeatTransferEM ) {
        this->computeNSubMatrixAt(answer, coords);
    } else {
        FloatMatrix n;
        int i, j;

        this->computeNSubMatrixAt(n, coords);
        answer.resize(2, 6);
        for ( i = 1; i <= 2; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, ( j - 1 ) * 2 + i) = n.at(1, j);
            }
        }
    }
}


void
Tr1_ht :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    ///@todo Deal with matrix and vector (I find that the row-wise matrices should be transposed).
    //this->interp.evaldNdx(answer, gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this), 0.0);
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3, area;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(2, 3);

    answer.at(1, 1) = y2 - y3;
    answer.at(1, 2) = y3 - y1;
    answer.at(1, 3) = y1 - y2;

    answer.at(2, 1) = x3 - x2;
    answer.at(2, 2) = x1 - x3;
    answer.at(2, 3) = x2 - x1;

    answer.times( 1. / ( 2. * area ) );
}


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


void
Tr1_ht :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    if ( emode == HeatTransferEM ) {
        answer.setValues(1, T_f);
    } else if ( emode == HeatMass1TransferEM ) {
        answer.setValues(2, T_f, C_1);
    } else {
        _error("Unknown ElementMode");
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
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian(* gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    weight = gp->giveWeight();
    volume = determinant * weight * this->giveCrossSection()->give(CS_Thickness);

    return volume;
}


void
Tr1_ht :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress triangle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */
    ///@todo Use the interpolation class for this
    //this->interp.edgeEvalN(answer_vec, gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    double ksi, n1, n2;
    answer.resize(1, 2);
    answer.zero();

    ksi = gp->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.at(1, 1) = n1;
    answer.at(1, 2) = n2;
}


double
Tr1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant, thick;
    determinant = fabs( this->interp.edgeGiveTransformationJacobian(iEdge, *gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    thick = this->giveCrossSection()->give(CS_Thickness);
    return determinant *thick *gp->giveWeight();
}


void
Tr1_ht :: giveEdgeDofMapping(IntArray &answer, int iEdge)
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(2);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 2;
        answer.at(2) = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer.at(1) = 3;
        answer.at(2) = 1;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}

void
Tr1_ht :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp.edgeLocal2global(answer, iEdge, *gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this));
}

void
Tr1_ht :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
{
    if ( emode == HeatTransferEM ) {
        this->computeInternalSourceRhsSubVectorAt(answer, atTime, mode, 1);
    } else if ( emode == HeatMass1TransferEM ) {
        FloatArray subAnswer;
        int i;

        for ( i = 1; i <= 2; i++ ) {
            this->computeInternalSourceRhsSubVectorAt(subAnswer, atTime, mode, i);
            if ( subAnswer.isNotEmpty() ) {
                if ( answer.isEmpty() ) {
                    answer.resize(2);
                    answer.zero();
                }

                this->assembleLocalContribution(answer, subAnswer, 2, i, 1.0);
            }
        }
    } else {
        _error("Unknown ElementMode");
    }
}

int
Tr1_ht :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interp.local2global(answer, lcoords, FEIElementGeometryWrapper(this));
    return 1;
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
    if ( type == IST_Temperature || type == IST_HydrationDegree || type == IST_Density || type == IST_ThermalConductivityIsotropic || type == IST_HeatCapacity || type == IST_AverageTemperature || type == IST_YoungModulusVirginPaste || type == IST_PoissonRatioVirginPaste || type == IST_YoungModulusConcrete || type == IST_PoissonRatioConcrete ) {
        return 1;
    } else if ( type == IST_TemperatureFlow || type == IST_HumidityFlow ) {
        return 2;
    }

    return 0;
}

void
Tr1_ht :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    int i;
    FloatMatrix n;
    this->computeNmatrixAt( n, aGaussPoint->giveCoordinates() );

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 3);
    } else {
        return;
    }

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, i)  = n.at(1, i);
    }
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
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);

    if ( ( size = coords.giveSize() ) < ( gsize = gcoords.giveSize() ) ) {
        _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");
    }

    if ( size == gsize ) {
        dist = coords.distance(gcoords);
    } else {
        FloatArray helpCoords = coords;

        helpCoords.resize(gsize);
        dist = helpCoords.distance(gcoords);
    }

    return dist;
}


int
Tr1_ht :: computeLocalCoordinates(FloatArray &lcoords, const FloatArray &coords)
{
    return this->interp.global2local(lcoords, coords, FEIElementGeometryWrapper(this));
}

} // end namespace oofem

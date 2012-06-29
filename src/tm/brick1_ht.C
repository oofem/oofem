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

#include "brick1_ht.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"
#include "load.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "conTable.h"
#endif

namespace oofem {
FEI3dHexaLin Brick1_ht :: interpolation;

Brick1_ht :: Brick1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM)
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 8;
}

Brick1_hmt :: Brick1_hmt(int n, Domain *aDomain) : Brick1_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

Brick1_ht :: ~Brick1_ht()
{ }

void
Brick1_ht :: computeNSubMatrixAt(FloatMatrix &answer, FloatArray *coords)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    FloatArray n;
    this->interpolation.evalN(n, * coords, FEIElementGeometryWrapper(this));
    answer.resize(1, 8);

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
Brick1_ht :: computeNmatrixAt(FloatMatrix &answer, FloatArray *coords)
{
    if ( emode == HeatTransferEM ) {
        this->computeNSubMatrixAt(answer, coords);
    } else {
        FloatMatrix n;
        int i, j;

        this->computeNSubMatrixAt(n, coords);
        answer.resize(2, 16);
        for ( i = 1; i <= 2; i++ ) {
            for ( j = 1; j <= 8; j++ ) {
                answer.at(i, ( j - 1 ) * 2 + i) = n.at(1, j);
            }
        }
    }
}

void
Brick1_ht :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    answer.beTranspositionOf(dnx);
}


void
Brick1_ht :: computeGaussPoints()
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

void
Brick1_ht :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
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
Brick1_ht :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);
    numberOfGaussPoints = 8;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_Brick1_ht_nip, "nip"); // Macro

    if ( !( ( numberOfGaussPoints == 8 ) ||
           ( numberOfGaussPoints == 27 ) ) ) {
        numberOfGaussPoints = 8;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


double
Brick1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );

    weight = aGaussPoint->giveWeight();
    volume = determinant * weight;
    return volume;
}

void
Brick1_ht :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */
    FloatArray n(2);
    this->interpolation.edgeEvalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(1, 2);
    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
}


double
Brick1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * gp->giveWeight();
}


void
Brick1_ht :: giveEdgeDofMapping(IntArray &answer, int iEdge)
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
        answer.at(2) = 4;
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer.at(1) = 4;
        answer.at(2) = 1;
    } else if ( iEdge == 5 ) { // edge between nodes 1 5
        answer.at(1) = 1;
        answer.at(2) = 5;
    } else if ( iEdge == 6 ) { // edge between nodes 2 6
        answer.at(1) = 2;
        answer.at(2) = 6;
    } else if ( iEdge == 7 ) { // edge between nodes 3 7
        answer.at(1) = 3;
        answer.at(2) = 7;
    } else if ( iEdge == 8 ) { // edge between nodes 4 8
        answer.at(1) = 4;
        answer.at(2) = 8;
    } else if ( iEdge == 9 ) { // edge between nodes 5 6
        answer.at(1) = 5;
        answer.at(2) = 6;
    } else if ( iEdge == 10 ) { // edge between nodes 6 7
        answer.at(1) = 6;
        answer.at(2) = 7;
    } else if ( iEdge == 11 ) { // edge between nodes 7 8
        answer.at(1) = 7;
        answer.at(2) = 8;
    } else if ( iEdge == 12 ) { // edge between nodes 8 5
        answer.at(1) = 8;
        answer.at(2) = 5;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}

void
Brick1_ht :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global(answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}

IntegrationRule *Brick1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->setUpIntegrationPoints(_Square, npoints, _Unknown);
    return iRule;
}

void Brick1_ht :: computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n(4);
    interpolation.surfaceEvalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(1, 4);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
    answer.at(1, 3) = n.at(3);
    answer.at(1, 4) = n.at(4);
}

double Brick1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );

    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;
}

void Brick1_ht :: giveSurfaceDofMapping(IntArray &answer, int iSurf)
{
    answer.resize(4);
    if ( iSurf == 1 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 4; // node 4
        answer.at(3) = 3; // node 3
        answer.at(4) = 2; // node 2
    } else if ( iSurf == 2 ) {
        answer.at(1) = 5; // node 5
        answer.at(2) = 6; // node 6
        answer.at(3) = 7; // node 7
        answer.at(4) = 8; // node 8
    } else if ( iSurf == 3 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 2; // node 2
        answer.at(3) = 6; // node 6
        answer.at(4) = 5; // node 5
    } else if ( iSurf == 4 ) {
        answer.at(1) = 2; // node 2
        answer.at(2) = 3; // node 3
        answer.at(3) = 7; // node 7
        answer.at(4) = 6; // node 6
    } else if ( iSurf == 5 ) {
        answer.at(1) = 3; // node 3
        answer.at(2) = 4; // node 4
        answer.at(3) = 8; // node 8
        answer.at(4) = 7; // node 7
    } else if ( iSurf == 6 ) {
        answer.at(1) = 4; // node 4
        answer.at(2) = 1; // node 1
        answer.at(3) = 5; // node 5
        answer.at(4) = 8; // node 8
    } else {
        _error("giveSurfaceDofMapping: wrong surface number");
    }
}

void
Brick1_ht :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    interpolation.surfaceLocal2global(answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


void
Brick1_ht :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
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
                    answer.resize(16);
                    answer.zero();
                }

                this->assembleLocalContribution(answer, subAnswer, 2, i, 1.0);
            }
        }
    } else {
        _error("Unknown ElementMode");
    }
}


Interface *
Brick1_ht :: giveInterface(InterfaceType interface)
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
Brick1_ht :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_Temperature || type == IST_HydrationDegree || type == IST_Density || type == IST_ThermalConductivityIsotropic || type == IST_HeatCapacity || type == IST_AverageTemperature || type == IST_YoungModulusVirginPaste || type == IST_PoissonRatioVirginPaste || type == IST_YoungModulusConcrete || type == IST_PoissonRatioConcrete ) {
        return 1;
    } else if ( type == IST_TemperatureFlow || type == IST_HumidityFlow ) {
        return 3;
    }

    return 0;
}

int
Brick1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Brick1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 0.0;
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


#ifdef __OOFEG
void Brick1_ht :: drawRawGeometry(oofegGraphicContext &gc)
{
    int i;
    WCRec p [ 8 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    for ( i = 0; i < 8; i++ ) {
        p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
        p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
        p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
    }

    go =  CreateHexahedron(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Brick1_ht :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 8 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 8 ];
    double s [ 8 ];
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 8; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }

        if ( result != 8 ) {
            return;
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( (!result) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    for ( i = 1; i <= 8; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( context.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 8; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
        }

        context.updateFringeTableMinMax(s, 8);
        tr = CreateHexahedronWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem

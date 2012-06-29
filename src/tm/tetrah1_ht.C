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

#include "tetrah1_ht.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"
#include "structuralms.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "conTable.h"
#endif

namespace oofem {
FEI3dTrLin Tetrah1_ht :: interpolation;

Tetrah1_ht :: Tetrah1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM)
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 1;
}

Tetrah1_hmt :: Tetrah1_hmt(int n, Domain *aDomain) : Tetrah1_ht(n, aDomain)
{
    this->emode = HeatMass1TransferEM; // This could be done in a better way.
}

Tetrah1_ht :: ~Tetrah1_ht()
{ }

void
Tetrah1_ht :: computeNSubMatrixAt(FloatMatrix &answer, FloatArray *coords)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    FloatArray n;
    this->interpolation.evalN(n, * coords, FEIElementGeometryWrapper(this));
    answer.resize(1, 4);

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, i) = n.at(i);
    }
}

void
Tetrah1_ht :: computeNmatrixAt(FloatMatrix &answer, FloatArray *coords)
{
    if ( emode == HeatTransferEM ) {
        this->computeNSubMatrixAt(answer, coords);
    } else {
        FloatMatrix n;
        int i, j;

        this->computeNSubMatrixAt(n, coords);
        answer.resize(2, 8);
        for ( i = 1; i <= 2; i++ ) {
            for ( j = 1; j <= 4; j++ ) {
                answer.at(i, ( j - 1 ) * 2 + i) = n.at(1, j);
            }
        }
    }
}

void
Tetrah1_ht :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    answer.beTranspositionOf(dnx);
}


void
Tetrah1_ht :: computeGaussPoints()
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
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Tetrahedra, numberOfGaussPoints, mmode);
    }
}

void
Tetrah1_ht :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
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
Tetrah1_ht :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);
    numberOfGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_Tetrah1_ht_nip, "nip"); // Macro

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ) ) {
        numberOfGaussPoints = 1;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


double
Tetrah1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
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
Tetrah1_ht :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
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
Tetrah1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * gp->giveWeight();
}


void
Tetrah1_ht :: giveEdgeDofMapping(IntArray &answer, int iEdge)
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
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 3;
        answer.at(2) = 1;
    } else if ( iEdge == 4 ) { // edge between nodes 1 4
        answer.at(1) = 1;
        answer.at(2) = 4;
    } else if ( iEdge == 5 ) { // edge between nodes 2 4
        answer.at(1) = 2;
        answer.at(2) = 4;
    } else if ( iEdge == 6 ) { // edge between nodes 3 4
        answer.at(1) = 3;
        answer.at(2) = 4;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}

void
Tetrah1_ht :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global(answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}

IntegrationRule *
Tetrah1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->setUpIntegrationPoints(_Triangle, npoints, _Unknown);
    return iRule;
}

void
Tetrah1_ht :: computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n(3);
    interpolation.surfaceEvalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
    answer.at(1, 3) = n.at(3);
}

double
Tetrah1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );

    weight      = gp->giveWeight();
    volume      = 2.0 * determinant * weight;

    return volume;
}

void
Tetrah1_ht :: giveSurfaceDofMapping(IntArray &answer, int iSurf)
{
    answer.resize(3);
    if ( iSurf == 1 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 3; // node 3
        answer.at(3) = 2; // node 2
    } else if ( iSurf == 2 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 2; // node 2
        answer.at(3) = 4; // node 4
    } else if ( iSurf == 3 ) {
        answer.at(1) = 2; // node 2
        answer.at(2) = 3; // node 3
        answer.at(3) = 4; // node 4
    } else if ( iSurf == 4 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 4; // node 4
        answer.at(3) = 3; // node 3
    } else {
        _error("giveSurfaceDofMapping: wrong surface number");
    }
}

void
Tetrah1_ht :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    interpolation.surfaceLocal2global(answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


void
Tetrah1_ht :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
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
                    answer.resize(8);
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
Tetrah1_ht :: giveInterface(InterfaceType interface)
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
Tetrah1_ht :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_TemperatureFlow ) || ( type == IST_HumidityFlow ) ) {
        return 3;
    }

    return 0;
}

int
Tetrah1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Tetrah1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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
 #define TR_LENGHT_REDUCT 0.3333

void
Tetrah1_ht :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 4 ];
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
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveCoordinate(3);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void
Tetrah1_ht :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 4 ];
    double s [ 4 ];
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( (!result) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    for ( i = 1; i <= 4; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( context.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 4; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
        }

        context.updateFringeTableMinMax(s, 4);
        tr = CreateTetraWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem

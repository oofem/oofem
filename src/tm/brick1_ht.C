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


double
Brick1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * gp->giveWeight();
}


IntegrationRule *
Brick1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->setUpIntegrationPoints(_Square, npoints, _Unknown);
    return iRule;
}


double
Brick1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );
    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
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
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

int
Brick1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Brick1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.zero();
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
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

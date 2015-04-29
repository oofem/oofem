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

#include "brick1_ht.h"
#include "fei3dhexalin.h"
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

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(Brick1_ht);
REGISTER_Element(Brick1_hmt);
REGISTER_Element(Brick1_mt);

FEI3dHexaLin Brick1_ht :: interpolation;

Brick1_ht :: Brick1_ht(int n, Domain *aDomain) : TransportElement(n, aDomain, HeatTransferEM), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this), SPRNodalRecoveryModelInterface()
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 8;
}

Brick1_hmt :: Brick1_hmt(int n, Domain *aDomain) : Brick1_ht(n, aDomain)
{
    emode = HeatMass1TransferEM;
}

Brick1_mt :: Brick1_mt(int n, Domain *aDomain) : Brick1_ht(n, aDomain)
{
    emode = Mass1TransferEM;
}

Brick1_ht :: ~Brick1_ht()
{ }

FEInterpolation *
Brick1_ht :: giveInterpolation() const { return & interpolation; }

void
Brick1_ht :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
Brick1_ht :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 8;
    return TransportElement :: initializeFrom(ir);
}


double
Brick1_ht :: computeVolumeAround(GaussPoint *gp)
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
Brick1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}


IntegrationRule *
Brick1_ht :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}


double
Brick1_ht :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
}


Interface *
Brick1_ht :: giveInterface(InterfaceType interface)
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
Brick1_ht :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(numberOfDofMans);
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Brick1_ht :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
Brick1_ht :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
Brick1_ht :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}


#ifdef __OOFEG
void Brick1_ht :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
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
    for ( int i = 0; i < 8; i++ ) {
        p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
        p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
        p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
    }

    go =  CreateHexahedron(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Brick1_ht :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx, result = 0;
    WCRec p [ 8 ];
    GraphicObj *tr;
    FloatArray v [ 8 ];
    double s [ 8 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( int i = 1; i <= 8; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 8 ) {
            return;
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    for ( int i = 1; i <= 8; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( int i = 0; i < 8; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
        }

        gc.updateFringeTableMinMax(s, 8);
        tr = CreateHexahedronWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem

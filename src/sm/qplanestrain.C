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

#include "qplanestrain.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QPlaneStrain);

FEI2dQuadQuad QPlaneStrain :: interpolation(1, 2);

QPlaneStrain :: QPlaneStrain(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

Interface *
QPlaneStrain :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}


void
QPlaneStrain :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [4x16] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, 16);
    answer.zero();

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(4, 2 * i - 1) = dnx.at(i, 2);
        answer.at(4, 2 * i - 0) = dnx.at(i, 1);
    }
}

void
QPlaneStrain :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [5x16] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dnx;
    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(5, 16);
    answer.zero();

    // 3rd row is zero -> dw/dz = 0
    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dnx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dnx.at(i, 1);     // dv/dx -9
    }
}

void
QPlaneStrain :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    FloatArray n(8);

    answer.resize(2, 16);
    answer.zero();

    this->interpolation.evalN( n, iLocCoord, FEIElementGeometryWrapper(this) );

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }
}

IRResultType
QPlaneStrain :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = this->Element :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    return IRRT_OK;
}

void
QPlaneStrain :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

double
QPlaneStrain :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveCrossSection()->give(CS_Thickness, gp);
    volume      = determinant * weight * thickness;

    return volume;
}


void
QPlaneStrain :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


double
QPlaneStrain :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints() );
}


#ifdef __OOFEG
void QPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc)
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
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QPlaneStrain :: drawScalar(oofegGraphicContext &context)
{
    int i, indx;
    WCRec p [ 4 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.giveIntVarMode() == ISM_recovered ) {
        // ============ plot the recovered values (smoothed data) ===============
        /*
         * for ( i = 1; i <= 4; i++ ) {
         *    result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
         * }
         *
         * if ( result != 4 ) {
         *    return;
         * }
         *
         * indx = context.giveIntVarIndx();
         *
         * for ( i = 1; i <= 4; i++ ) {
         *    s [ i - 1 ] = v [ i - 1 ].at(indx);
         * }
         *
         * if ( context.getScalarAlgo() == SA_ISO_SURF ) {
         *    for ( i = 0; i < 4; i++ ) {
         *        if ( context.getInternalVarsDefGeoFlag() ) {
         *            // use deformed geometry
         *            defScale = context.getDefScale();
         *            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
         *            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
         *            p [ i ].z = 0.;
         *        } else {
         *            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
         *            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
         *            p [ i ].z = 0.;
         *        }
         *    }
         *
         *    //EASValsSetColor(gc.getYieldPlotColor(ratio));
         *    context.updateFringeTableMinMax(s, 4);
         *    tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
         *    EGWithMaskChangeAttributes(LAYER_MASK, tr);
         *    EMAddGraphicsToModel(ESIModel(), tr);
         * } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
         *    double landScale = context.getLandScale();
         *
         *    for ( i = 0; i < 4; i++ ) {
         *        if ( context.getInternalVarsDefGeoFlag() ) {
         *            // use deformed geometry
         *            defScale = context.getDefScale();
         *            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
         *            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
         *            p [ i ].z = s [ i ] * landScale;
         *        } else {
         *            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
         *            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
         *            p [ i ].z = s [ i ] * landScale;
         *        }
         *
         *        // this fixes a bug in ELIXIR
         *        if ( fabs(s [ i ]) < 1.0e-6 ) {
         *            s [ i ] = 1.0e-6;
         *        }
         *    }
         *
         *    if ( context.getScalarAlgo() == SA_ZPROFILE ) {
         *        EASValsSetColor( context.getDeformedElementColor() );
         *        EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
         *        tr =  CreateQuad3D(p);
         *        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
         *    } else {
         *        context.updateFringeTableMinMax(s, 4);
         *        tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
         *        EGWithMaskChangeAttributes(LAYER_MASK, tr);
         *    }
         *
         *    EMAddGraphicsToModel(ESIModel(), tr);
         * }
         */
    } else if ( context.giveIntVarMode() == ISM_local ) {
        // ========== plot the local values (raw data) =====================
        if ( numberOfGaussPoints != 4 ) {
            return;
        }

        int ip;
        GaussPoint *gp;
        IntArray ind(4);
        FloatArray *gpCoords;
        WCRec pp [ 9 ];

        for ( i = 0; i < 8; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.;

        for ( ip = 1; ip <= integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip - 1);
            gpCoords = gp->giveCoordinates();
            if ( ( gpCoords->at(1) > 0. ) && ( gpCoords->at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords->at(1) < 0. ) && ( gpCoords->at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords->at(1) < 0. ) && ( gpCoords->at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, context.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = context.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            context.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}

#endif
} // end namespace oofem

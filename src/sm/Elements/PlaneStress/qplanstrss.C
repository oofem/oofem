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

#include "../sm/Elements/PlaneStress/qplanstrss.h"
#include "fei2dquadquad.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(QPlaneStress2d);

FEI2dQuadQuad QPlaneStress2d :: interpolation(1, 2);

QPlaneStress2d :: QPlaneStress2d(int n, Domain *aDomain) :
    PlaneStressElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

Interface *
QPlaneStress2d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}

FEInterpolation *QPlaneStress2d :: giveInterpolation() const { return & interpolation; }



#ifdef __OOFEG
void QPlaneStress2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void QPlaneStress2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
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


void QPlaneStress2d :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx,  n [ 4 ], result = 0;
    WCRec p [ 4 ], pp [ 9 ];
    GraphicObj *tr;
    FloatArray v [ 8 ];
    double s [ 9 ], ss [ 4 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.giveIntVarMode() == ISM_recovered ) {
        // ============ plot the recovered values (smoothed data) ===============
        for ( i = 1; i <= 8; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 8 ) {
            return;
        }

        indx = gc.giveIntVarIndx();

        for ( i = 1; i <= 8; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        // auxiliary value at an added central node
        // computed as average of the values at all Gauss points

        s [ 8 ] = 0.;
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            s [ 8 ] +=  v [ 0 ].at(indx);
        }

        s [ 8 ] /= integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints();
        //s[8] = (s[4]+s[5]+s[6]+s[7])/4.;

        for ( i = 0; i < 8; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                // use initial geometry
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        pp [ 8 ].x = ( pp [ 4 ].x + pp [ 5 ].x + pp [ 6 ].x + pp [ 7 ].x ) / 4.;
        pp [ 8 ].y = ( pp [ 4 ].y + pp [ 5 ].y + pp [ 6 ].y + pp [ 7 ].y ) / 4.;
        pp [ 8 ].z = 0.;


        for ( int t = 1; t <= 4; t++ ) {
            if ( t == 1 ) {
                n [ 0 ] = 0;
                n [ 1 ] = 4;
                n [ 2 ] = 8;
                n [ 3 ] = 7;
            } else if ( t == 2 ) {
                n [ 0 ] = 4;
                n [ 1 ] = 1;
                n [ 2 ] = 5;
                n [ 3 ] = 8;
            } else if ( t == 3 ) {
                n [ 0 ] = 5;
                n [ 1 ] = 2;
                n [ 2 ] = 6;
                n [ 3 ] = 8;
            } else {
                n [ 0 ] = 6;
                n [ 1 ] = 3;
                n [ 2 ] = 7;
                n [ 3 ] = 8;
            }

            ss [ 0 ] = s [ n [ 0 ] ];
            ss [ 1 ] = s [ n [ 1 ] ];
            ss [ 2 ] = s [ n [ 2 ] ];
            ss [ 3 ] = s [ n [ 3 ] ];


            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ n [ i ] ].x;
                p [ i ].y = pp [ n [ i ] ].y;
                p [ i ].z = 0.;
            }

            if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
                /*
                 * for ( i = 0; i < 4; i++ ) {
                 *    if ( gc.getInternalVarsDefGeoFlag() ) {
                 *        // use deformed geometry
                 *        defScale = gc.getDefScale();
                 *        p [ i ].x = ( FPNum ) this->giveNode(n[i] + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                 *        p [ i ].y = ( FPNum ) this->giveNode(n[i] + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                 *        p [ i ].z = 0.;
                 *    } else {
                 *        // use initial geometry
                 *        p [ i ].x = ( FPNum ) this->giveNode(n[i] + 1)->giveCoordinate(1);
                 *        p [ i ].y = ( FPNum ) this->giveNode(n[i] + 1)->giveCoordinate(2);
                 *        p [ i ].z = 0.;
                 *    }
                 * }
                 */
                //EASValsSetColor(gc.getYieldPlotColor(ratio));
                gc.updateFringeTableMinMax(ss, 4);
                tr =  CreateQuadWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ], ss [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);
            } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
                //double landScale = gc.getLandScale();

                for ( i = 0; i < 4; i++ ) {
                    /*
                     * if ( gc.getInternalVarsDefGeoFlag() ) {
                     *    // use deformed geometry
                     *    defScale = gc.getDefScale();
                     *    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                     *    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                     *    p [ i ].z = ss [ i ] * landScale;
                     * } else {
                     *    // use initial geometry
                     *    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                     *    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                     *    p [ i ].z = ss [ i ] * landScale;
                     * }
                     */

                    // this fixes a bug in ELIXIR
                    if ( fabs(ss [ i ]) < 1.0e-6 ) {
                        ss [ i ] = 1.0e-6;
                    }
                }

                if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
                    EASValsSetColor( gc.getDeformedElementColor() );
                    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                    tr =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                } else {
                    gc.updateFringeTableMinMax(s, 4);
                    tr =  CreateQuadWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ], ss [ 3 ]);
                    EGWithMaskChangeAttributes(LAYER_MASK, tr);
                }

                EMAddGraphicsToModel(ESIModel(), tr);
            }
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        // ========== plot the local values (raw data) =====================
        if ( numberOfGaussPoints != 4 ) {
            return;
        }

        IntArray ind(4);
        WCRec pp [ 9 ];

        for ( i = 0; i < 8; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
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

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray& gpCoords = gp->giveNaturalCoordinates();
            if ( ( gpCoords.at(1) > 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) < 0. ) ) {
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

            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = gc.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}
#endif

void
QPlaneStress2d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                             InternalStateType type, TimeStep *tStep)
{
    if ( numberOfGaussPoints != 4 ) {
        return;
    }

    GaussPoint *gp;

    if ( node < 5 ) {
        int i = 0;
        switch ( node ) {
        case 1: i = 4;
            break;
        case 2: i = 2;
            break;
        case 3: i = 1;
            break;
        case 4: i = 3;
            break;
        }

        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i - 1);
        this->giveIPValue(answer, gp, type, tStep);
    } else {
        int i1 = 0, i2 = 0;
        switch ( node ) {
        case 5: i1 = 4;
            i2 = 2;
            break;
        case 6: i1 = 2;
            i2 = 1;
            break;
        case 7: i1 = 1;
            i2 = 3;
            break;
        case 8: i1 = 3;
            i2 = 4;
            break;
        }

        FloatArray contrib;
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i1 - 1);
        this->giveIPValue(contrib, gp, type, tStep);
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i2 - 1);
        this->giveIPValue(answer, gp, type, tStep);
        answer.add(contrib);
        answer.times(0.5);
    }
}


} // end namespace oofem

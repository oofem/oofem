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

#include "qplanstrss.h"
#include "node.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdio.h>
#endif

namespace oofem {
FEI2dQuadQuad QPlaneStress2d :: interpolation(1, 2);

QPlaneStress2d :: QPlaneStress2d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

Interface *
QPlaneStress2d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    return NULL;
}


void
QPlaneStress2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [3x16] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(3, 16);
    answer.zero();

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}

void
QPlaneStress2d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    int i;
    FloatArray n(8);

    answer.resize(2, 16);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }
}

IRResultType
QPlaneStress2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QPlaneStress2d_nip, "nip"); // Macro

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}

void
QPlaneStress2d :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStress);
    }
}

double
QPlaneStress2d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );
    weight      = aGaussPoint->giveWeight();
    thickness   = this->giveCrossSection()->give(CS_Thickness);
    volume      = determinant * weight * thickness;

    return volume;
}


void
QPlaneStress2d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


int
QPlaneStress2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this));
    return 1;
}

double
QPlaneStress2d :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    if ( normalToCrackPlane.at(3) < 0.999999 ){//ensure that characteristic length is in the plane of element
        return this->giveLenghtInDir(normalToCrackPlane) / sqrt( ( double ) this->numberOfGaussPoints );
    } else {//otherwise compute out-of-plane characteristic length from element area
        return sqrt ( this->computeVolumeAreaOrLength() / ( double ) this->numberOfGaussPoints );
    }
}



#ifdef __OOFEG
void QPlaneStress2d :: drawRawGeometry(oofegGraphicContext &gc)
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


void QPlaneStress2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QPlaneStress2d :: drawScalar(oofegGraphicContext &context)
{
    int i, indx,  n [ 4 ], result = 0;
    WCRec p [ 4 ], pp [ 9 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 8 ];
    double s [ 9 ], ss [ 4 ], defScale;
    IntArray map;
    int ip;
    GaussPoint *gp;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.giveIntVarMode() == ISM_recovered ) {
        // ============ plot the recovered values (smoothed data) ===============
        for ( i = 1; i <= 8; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }

        if ( result != 8 ) {
            return;
        }

        this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
        if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
            return;
        }

        for ( i = 1; i <= 8; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        // auxiliary value at an added central node
        // computed as average of the values at all Gauss points

        s [ 8 ] = 0.;
        for ( ip = 1; ip <= numberOfGaussPoints; ip++ ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip - 1);
            if ( giveIPValue(v [ 0 ], gp, context.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            s [ 8 ] +=  v [ 0 ].at(indx);
        }

        s [ 8 ] /= numberOfGaussPoints;
        //s[8] = (s[4]+s[5]+s[6]+s[7])/4.;

        for ( i = 0; i < 8; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
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

            if ( context.getScalarAlgo() == SA_ISO_SURF ) {
                /*
                 * for ( i = 0; i < 4; i++ ) {
                 *    if ( context.getInternalVarsDefGeoFlag() ) {
                 *        // use deformed geometry
                 *        defScale = context.getDefScale();
                 *        p [ i ].x = ( FPNum ) this->giveNode(n[i] + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                 *        p [ i ].y = ( FPNum ) this->giveNode(n[i] + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
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
                context.updateFringeTableMinMax(ss, 4);
                tr =  CreateQuadWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ], ss [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
                EMAddGraphicsToModel(ESIModel(), tr);
            } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
	      //double landScale = context.getLandScale();

                for ( i = 0; i < 4; i++ ) {
                    /*
                     * if ( context.getInternalVarsDefGeoFlag() ) {
                     *    // use deformed geometry
                     *    defScale = context.getDefScale();
                     *    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                     *    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
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

                if ( context.getScalarAlgo() == SA_ZPROFILE ) {
                    EASValsSetColor( context.getDeformedElementColor() );
                    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                    tr =  CreateQuad3D(p);
                    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                } else {
                    context.updateFringeTableMinMax(s, 4);
                    tr =  CreateQuadWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ], ss [ 3 ]);
                    EGWithMaskChangeAttributes(LAYER_MASK, tr);
                }

                EMAddGraphicsToModel(ESIModel(), tr);
            }
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        // ========== plot the local values (raw data) =====================
        if ( numberOfGaussPoints != 4 ) {
            return;
        }

        IntArray ind(4);
        FloatArray *gpCoords;
        WCRec pp [ 9 ];

        for ( i = 0; i < 8; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
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

        for ( ip = 1; ip <= numberOfGaussPoints; ip++ ) {
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

            this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
            if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
                return;
            }

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

int
QPlaneStress2d :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) || ( type == IST_DamageTensor ) ) {
        return 3;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

void
QPlaneStress2d :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    int i;
    FloatArray n;

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 8);
    } else {
        return;
    }

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, i)  = n.at(i);
    }
}


void
QPlaneStress2d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                             InternalStateType type, TimeStep *tStep)
{
    if ( numberOfGaussPoints != 4 ) {
        return;
    }

    int i, i1, i2;
    GaussPoint *gp;

    if ( node < 5 ) {
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
    } else   {
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

    /*
     * if (type == IST_StressTensor) {
     * gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
     * answer = ((StructuralMaterialStatus*) this->giveMaterial()->giveStatus(gp)) -> giveStressVector();
     * } else if (type == IST_StrainTensor) {
     * gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
     * answer = ((StructuralMaterialStatus*) this->giveMaterial()->giveStatus(gp)) -> giveStrainVector();
     * }else answer.resize(0);
     */
}

void
QPlaneStress2d :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                            InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

void
QPlaneStress2d :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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

    FloatArray n(3);
    this->interpolation.edgeEvalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(2, 6);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 3) = n.at(2);
    answer.at(1, 5) = n.at(3);
    answer.at(2, 2) = n.at(1);
    answer.at(2, 4) = n.at(2);
    answer.at(2, 6) = n.at(3);
}


void
QPlaneStress2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

  int i;
  IntArray eNodes(3);
  this->interpolation.computeLocalEdgeMapping(eNodes,  iEdge);

  answer.resize(6);
  for (i=1; i<=3; i++) {
    answer.at(i*2-1) = eNodes.at(i)*2-1;
    answer.at(i*2  ) = eNodes.at(i)*2;
  }
}

double
QPlaneStress2d ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * aGaussPoint->giveWeight();
}

void
QPlaneStress2d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global(answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


int
QPlaneStress2d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatArray normal(2);
    answer.resize(2, 2);
    answer.zero();

    this->interpolation.edgeEvalNormal(normal, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.at(1, 1) = (-1.0)*normal.at(2);
    answer.at(1, 2) = (-1.0)*normal.at(1);
    answer.at(2, 1) = normal.at(1);
    answer.at(2, 2) = (-1.0)*normal.at(2);

    return 1;
}

} // end namespace oofem

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

#include "../sm/Elements/PlaneStrain/qtrplanestrain.h"
#include "fei2dtrquad.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
 #include "oofegutils.h"
 #include "Materials/rcm2.h"
#endif


namespace oofem {
REGISTER_Element(QTrPlaneStrain);

FEI2dTrQuad QTrPlaneStrain :: interpolation(1, 2);

QTrPlaneStrain :: QTrPlaneStrain(int n, Domain *aDomain) :
    PlaneStrainElement(n, aDomain), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 6;
    numberOfGaussPoints = 4;
}


Interface *
QTrPlaneStrain :: giveInterface(InterfaceType interface)
{
    //if (interface == NodalAveragingRecoveryModelInterfaceType) return (NodalAveragingRecoveryModelInterface*) this;
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }

    return NULL;
}


#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void QTrPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QTrPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 3 ];
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
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QTrPlaneStrain :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int t, n [ 3 ], i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v [ 6 ];
    double s [ 6 ], ss [ 3 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 6; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        return;
    }

    if ( result != 6 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    for ( i = 1; i <= 6; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( t = 1; t <= 4; t++ ) {
            if ( t == 1 ) {
                n [ 0 ] = 1;
                n [ 1 ] = 4;
                n [ 2 ] = 6;
            } else if ( t == 2 ) {
                n [ 0 ] = 2;
                n [ 1 ] = 5;
                n [ 2 ] = 4;
            } else if ( t == 3 ) {
                n [ 0 ] = 3;
                n [ 1 ] = 6;
                n [ 2 ] = 5;
            } else {
                n [ 0 ] = 4;
                n [ 1 ] = 5;
                n [ 2 ] = 6;
            }


            for ( i = 0; i < 3; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = 0.;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(n [ i ])->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(n [ i ])->giveCoordinate(2);
                    p [ i ].z = 0.;
                }
            }

            //EASValsSetColor(gc.getYieldPlotColor(ratio));
            ss [ 0 ] = s [ n [ 0 ] - 1 ];
            ss [ 1 ] = s [ n [ 1 ] - 1 ];
            ss [ 2 ] = s [ n [ 2 ] - 1 ];
            gc.updateFringeTableMinMax(ss, 3);
            tr =  CreateTriangleWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }

        /* } else if (gc.getScalarAlgo() == SA_ISO_LINE) {
         *
         * EASValsSetColor(context.getActiveCrackColor());
         * EASValsSetLineWidth(OOFEG_ISO_LINE_WIDTH);
         *
         * for (t=1; t<=4; t++) {
         * if (t==1) {n[0] = 1; n[1]=4; n[2]=6;}
         * else if (t==2) {n[0]=2; n[1]=5; n[2]=4;}
         * else if (t==3) {n[0]=3; n[1]=6; n[2]=5;}
         * else {n[0]=4; n[1]=5; n[2]=6;}
         *
         *
         * for (i=0; i< 3; i++) {
         * if (gc.getInternalVarsDefGeoFlag()) {
         * // use deformed geometry
         * defScale = gc.getDefScale();
         * p[i].x = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(1,tStep,defScale);
         * p[i].y = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(2,tStep,defScale);
         * p[i].z = 0.;
         *
         * } else {
         * p[i].x = (FPNum) this->giveNode(n[i])->giveCoordinate(1);
         * p[i].y = (FPNum) this->giveNode(n[i])->giveCoordinate(2);
         * p[i].z = 0.;
         * }
         * }
         * sv[0]=s[n[0]-1];
         * sv[1]=s[n[1]-1];
         * sv[2]=s[n[2]-1];
         *
         * // isoline implementation
         * oofeg_drawIsoLinesOnTriangle (p, sv);
         * } */
    }
}

void
QTrPlaneStrain :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{ }

#endif


void
QTrPlaneStrain :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
QTrPlaneStrain :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(3);
    if ( pap == this->giveNode(1)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(4)->giveNumber();
        answer.at(3) = this->giveNode(6)->giveNumber();
    } else if ( pap == this->giveNode(2)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(5)->giveNumber();
        answer.at(3) = this->giveNode(4)->giveNumber();
    } else if ( pap == this->giveNode(3)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(6)->giveNumber();
        answer.at(3) = this->giveNode(5)->giveNumber();
    } else {
        OOFEM_ERROR("node unknown");
    }
}

int
QTrPlaneStrain :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
QTrPlaneStrain :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dquadratic;
}

} // end namespace oofem

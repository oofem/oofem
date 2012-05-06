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

#include "tr_shell01.h"

#ifdef __OOFEG
 #include "node.h"
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif


namespace oofem {
TR_SHELL01 :: TR_SHELL01(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    plate    = new CCTPlate3d(-1, aDomain);
    membrane = new TrPlaneStrRot3d(-1, aDomain);

    numberOfDofMans = 3;
}


IRResultType
TR_SHELL01 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    // proc tady neni return = this...   ??? termitovo
    this->StructuralElement :: initializeFrom(ir);

    int val = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_TrPlaneStrRot_nip, "nip"); // Macro
    if ( val != -1 ) {
        _error("key word NIP is not allowed for element TR_SHELL01");
    }

    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_TrPlaneStrRot_niprot, "niprot"); // Macro
    if ( val != -1 ) {
        _error("key word NIProt is not allowed for element TR_SHELL01");
    }

    //
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);

    return IRRT_OK;
}


void
TR_SHELL01 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver accordind to mtrx
//
{
    FloatArray aux;

    plate->giveCharacteristicVector(answer, mtrx, mode, tStep);
    membrane->giveCharacteristicVector(aux, mtrx, mode, tStep);

    answer.add(aux);
}

void
TR_SHELL01 :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    FloatMatrix aux;

    plate->giveCharacteristicMatrix(answer, mtrx, tStep);
    membrane->giveCharacteristicMatrix(aux, mtrx, tStep);

    answer.add(aux);
}


void
TR_SHELL01 :: updateInternalState(TimeStep *stepN)
// Updates the receiver at end of step.
{
    plate->updateInternalState(stepN);
    membrane->updateInternalState(stepN);
}

void
TR_SHELL01 :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);

    plate->updateYourself(tStep);
    membrane->updateYourself(tStep);
}


void
TR_SHELL01 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v, aux;
    GaussPoint *plateGP    = plate->giveMiddleGaussPoint();
    GaussPoint *membraneGP = membrane->giveMiddleGaussPoint();


#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    fprintf(file, "  GP %2d.%-2d :", 1, 1);

    // Strain - Curvature
    plate->giveIPValue(v, plateGP, IST_ShellStrainCurvatureTensor, tStep);
    membrane->giveIPValue(aux, membraneGP, IST_ShellStrainCurvatureTensor, tStep);
    v.add(aux);

    fprintf(file, "  strains ");
    fprintf( file,
            " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
            v.at(1), v.at(2), v.at(3),  2. * v.at(4), 2. * v.at(5), 2. * v.at(6),
            v.at(7), v.at(8), v.at(9),  2. * v.at(10), 2. * v.at(11), 2. * v.at(12) );

    // Strain - Curvature
    plate->giveIPValue(v, plateGP, IST_ShellForceMomentumTensor, tStep);
    membrane->giveIPValue(aux, membraneGP, IST_ShellForceMomentumTensor, tStep);
    v.add(aux);

    fprintf(file, "\n              stresses");
    fprintf( file,
            " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
            v.at(1), v.at(2), v.at(3),  v.at(4), v.at(5), v.at(6),
            v.at(7), v.at(8), v.at(9),  v.at(10), v.at(11), v.at(12) );

    fprintf(file, "\n");
}


//
// io routines
//
#ifdef __OOFEG

void
TR_SHELL01 :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL01 :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getDeformedElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

// void
// CCTPlate  :: drawScalar(oofegGraphicContext &context)
// {
//   int i, indx, result = 0;
//   WCRec p [ 3 ];
//   GraphicObj *tr;
//   TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
//   FloatArray v1, v2, v3;
//   double s [ 3 ], defScale;
//   IntArray map;
//
//   if ( !context.testElementGraphicActivity(this) ) {
//     return;
//   }
//
//   if ( !this->giveMaterial()->isActivated(tStep) ) {
//     return;
//   }
//
//   if ( context.giveIntVarMode() == ISM_recovered ) {
//     result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
//     result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
//     result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
//   } else if ( context.giveIntVarMode() == ISM_local ) {
//     GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
//     result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
//     v2 = v1;
//     v3 = v1;
//     result *= 3;
//   }
//
//   if ( result != 3 ) {
//     return;
//   }
//
//   this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
//
//   if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
//     return;
//   }
//
//   s [ 0 ] = v1.at(indx);
//   s [ 1 ] = v2.at(indx);
//   s [ 2 ] = v3.at(indx);
//
//   EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
//
//   if ( context.getScalarAlgo() == SA_ISO_SURF ) {
//     for ( i = 0; i < 3; i++ ) {
//       if ( context.getInternalVarsDefGeoFlag() ) {
//  // use deformed geometry
//  defScale = context.getDefScale();
//  p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
//  p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
//  p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
//       } else {
//  p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
//  p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
//  p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
//       }
//     }
//
//     //EASValsSetColor(gc.getYieldPlotColor(ratio));
//     context.updateFringeTableMinMax(s, 3);
//     tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
//     EGWithMaskChangeAttributes(LAYER_MASK, tr);
//     EMAddGraphicsToModel(ESIModel(), tr);
//   }
// }
//
/*
 * void
 * CCTPlate  :: drawInternalState (oofegGraphicContext & gc)
 * //
 * // Draws internal state graphics representation
 * //
 * {
 * WCRec p[3];
 * GraphicObj *tr;
 * double v1,v2,v3;
 * DrawMode mode = gc.getDrawMode();
 * TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 * double defScale = gc.getDefScale();
 *
 * if (!gc.testElementGraphicActivity(this)) return;
 *
 * // check for valid DrawMode
 * if (!((mode == mxForce) || (mode == myForce) || (mode == mxyForce) ||
 * (mode == szxForce) || (mode == syzForce))) return;
 *
 *
 *
 * EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
 * if (gc.getInternalVarsDefGeoFlag()) {
 * // use deformed geometry
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[2].z = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * } else {
 * p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
 * p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
 * p[0].z = (FPNum) this->giveNode(1)->giveCoordinate(3);
 * p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
 * p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
 * p[1].z = (FPNum) this->giveNode(2)->giveCoordinate(3);
 * p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
 * p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
 * p[2].z = (FPNum) this->giveNode(3)->giveCoordinate(3);
 * }
 *
 * int result = 0;
 * result+= this->giveInternalStateAtNode (gc, 1, &v1);
 * result+= this->giveInternalStateAtNode (gc, 2, &v2);
 * result+= this->giveInternalStateAtNode (gc, 3, &v3);
 *
 * if (result == 3) {
 * tr = CreateTriangleWD3D (p,v1,v2,v3);
 * EGWithMaskChangeAttributes(LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 * }
 * }
 */

#endif
} // end namespace oofem

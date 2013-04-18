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
#include "fei2dtrlin.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"

#ifdef __OOFEG
 #include "node.h"
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {
TR_SHELL01 :: TR_SHELL01(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    plate    = new CCTPlate3d(-1, aDomain);
    membrane = new TrPlaneStrRot3d(-1, aDomain);
    compositeIR = NULL;

    numberOfDofMans = 3;
}


IRResultType
TR_SHELL01 :: initializeFrom(InputRecord *ir)
{
    // proc tady neni return = this...   ??? termitovo
    this->StructuralElement :: initializeFrom(ir);

    /*
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_Element_nip);
    if ( val != -1 ) {
        _error("key word NIP is not allowed for element TR_SHELL01");
    }


    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_TrPlaneStrRot_niprot, "niprot");
    if ( val != -1 ) {
        _error("key word NIProt is not allowed for element TR_SHELL01");
    }
    */

    //
    plate->initializeFrom(ir);
    membrane->initializeFrom(ir);

    plate->computeGaussPoints();
    membrane->computeGaussPoints();

    // check the compatibility of irules of plate and membrane
    if (plate->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints() != membrane->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints()) {
      OOFEM_ERROR ("TR_SHELL01: incompatible integration rules detected");
    }
    
    return IRRT_OK;
}


void
TR_SHELL01 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver accordind to mtrx
//
{
  FloatArray aux;
  FloatMatrix R;

    plate->giveCharacteristicVector(answer, mtrx, mode, tStep);
    if ( answer.isNotEmpty() ) {
        if (plate->giveRotationMatrix(R, EID_MomentumBalance)) {
            answer.rotatedWith(R, 't');
        }
    }
    membrane->giveCharacteristicVector(aux, mtrx, mode, tStep);
    if ( aux.isNotEmpty() ) {
        if (membrane->giveRotationMatrix(R, EID_MomentumBalance)) {
            aux.rotatedWith(R, 't');
        }
    }

    answer.add(aux);
}

void
TR_SHELL01 :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
  FloatMatrix aux, R;

    plate->giveCharacteristicMatrix(answer, mtrx, tStep);
    if (plate->giveRotationMatrix(R, EID_MomentumBalance)) {
      answer.rotatedWith(R);
    }

    membrane->giveCharacteristicMatrix(aux, mtrx, tStep);
    if (membrane->giveRotationMatrix(R, EID_MomentumBalance)) {
     aux.rotatedWith(R);
    }

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


Interface *
TR_SHELL01 :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >( this );
    } else if ( interface == ZZRemeshingCriteriaInterfaceType ) {
        return static_cast< ZZRemeshingCriteriaInterface * >( this );
    }


    return NULL;
}

double
TR_SHELL01 :: computeVolumeAround(GaussPoint *gp) 
{
    return plate->computeVolumeAround(gp);
}

int
TR_SHELL01 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_ShellForceMomentumTensor ) {
        FloatArray aux;
        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber()-1);
	GaussPoint *plateGP = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber()-1);
        
        plate->giveIPValue(answer, plateGP, IST_ShellForceMomentumTensor, atTime);
        membrane->giveIPValue(aux, membraneGP, IST_ShellForceMomentumTensor, atTime);
        answer.add(aux);
        return 1;
    } else if ( type == IST_ShellStrainCurvatureTensor ) {
        FloatArray aux;
        GaussPoint *membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber()-1);
	GaussPoint *plateGP = plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber()-1);

        plate->giveIPValue(answer, plateGP, IST_ShellStrainCurvatureTensor, atTime);
        membrane->giveIPValue(aux, membraneGP, IST_ShellStrainCurvatureTensor, atTime);
        answer.add(aux);

        return 1;
    } else {
        return StructuralElement::giveIPValue(answer, gp, type, atTime);
    }
}


int
TR_SHELL01 :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_ShellForceMomentumTensor || type == IST_ShellStrainCurvatureTensor ) ) {
        return 12;
    } else if ((type == IST_ErrorIndicatorLevel) || (type == IST_InternalStressError)) {
        return 1;
    } else {
        return StructuralElement::giveIPValueSize(type, gp);
    }
}


int
TR_SHELL01 :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
  if ( ( type == IST_ShellForceMomentumTensor || type == IST_ShellStrainCurvatureTensor ) ) {
    answer.resize(12);
    for (int i=1; i<=12; i++) answer.at(i)=i;
    return 1;
  } else {
    return StructuralElement::giveIntVarCompFullIndx(answer, type);
  }
}


//
// The element interface required by ZZNodalRecoveryModel
//
int
TR_SHELL01 :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{

  return giveIPValueSize (type, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0));
}


void
TR_SHELL01 :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type)
// evaluates N matrix (interpolation estimated stress matrix)
// according to Zienkiewicz & Zhu paper
// N(nsigma, nsigma*nnodes)
// Definition : sigmaVector = N * nodalSigmaVector
{
    FloatArray n;
    FEI2dTrLin interp_lin(1, 2);
    interp_lin.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3);
    answer.zero();
    answer.at(1) = n.at(1);
    answer.at(2) = n.at(2);
    answer.at(3) = n.at(3);
}

double 
TR_SHELL01::ZZRemeshingCriteriaI_giveCharacteristicSize() 
{
    return sqrt(plate->computeArea() * 2.0);
}



//
// The element interface required by NodalAveragingRecoveryModel
//
void
TR_SHELL01 :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    this->giveIPValue(answer, NULL, type, tStep);
}


void
TR_SHELL01 :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                      InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}





void
TR_SHELL01 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    FloatArray v, aux;
    GaussPoint *gp, *membraneGP;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
      gp  = iRule->getIntegrationPoint(i);
      fprintf(file, "  GP %2d.%-2d :", iRule->giveNumber(), gp->giveNumber());
      membraneGP = membrane->giveDefaultIntegrationRulePtr()->getIntegrationPoint(gp->giveNumber()-1);
      // Strain - Curvature
      plate->giveIPValue(v, gp, IST_ShellStrainCurvatureTensor, tStep);
      membrane->giveIPValue(aux, membraneGP, IST_ShellStrainCurvatureTensor, tStep);
      v.add(aux);

      fprintf(file, "  strains ");
      fprintf( file,
            " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
            v.at(1), v.at(2), v.at(3),  2. * v.at(4), 2. * v.at(5), 2. * v.at(6),
            v.at(7), v.at(8), v.at(9),  2. * v.at(10), 2. * v.at(11), 2. * v.at(12) );
      
      // Strain - Curvature
      plate->giveIPValue(v, gp, IST_ShellForceMomentumTensor, tStep);
      membrane->giveIPValue(aux, membraneGP, IST_ShellForceMomentumTensor, tStep);
      v.add(aux);
      
      fprintf(file, "\n              stresses");
      fprintf( file,
            " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
            v.at(1), v.at(2), v.at(3),  v.at(4), v.at(5), v.at(6),
            v.at(7), v.at(8), v.at(9),  v.at(10), v.at(11), v.at(12) );
      
      fprintf(file, "\n");
    }
}



contextIOResultType
TR_SHELL01 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
  contextIOResultType iores;
  if ( ( iores =  StructuralElement::saveContext (stream, mode, obj) ) != CIO_OK ) {
    THROW_CIOERR(iores);
  }
  if ( ( iores =  this->plate->saveContext(stream, mode, obj) ) != CIO_OK ) { 
    THROW_CIOERR(iores);
  }
    if ( ( iores = this->membrane->saveContext(stream, mode, obj) ) != CIO_OK ) {;
    THROW_CIOERR(iores);
  }
  return iores;
}

contextIOResultType
TR_SHELL01 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores =  StructuralElement::restoreContext (stream, mode, obj) ) != CIO_OK ) {;
        THROW_CIOERR(iores);
    }
        if ( ( iores =   this->plate->restoreContext(stream, mode, obj) ) != CIO_OK ) {;
        THROW_CIOERR(iores);
    }
    if ( ( iores =  this->membrane->restoreContext(stream, mode, obj) ) != CIO_OK ) {;
        THROW_CIOERR(iores);
    }
    return iores;
}

IntegrationRule* 
TR_SHELL01 :: ZZErrorEstimatorI_giveIntegrationRule()
{
  if (this->compositeIR) {
    return this->compositeIR;
  } else {
    this->compositeIR = new GaussIntegrationRule(1, this, 1, 12);
    this->compositeIR->setUpIntegrationPoints(_Triangle, plate->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints(), _3dShell);
    return this->compositeIR;
  }
}

void 
TR_SHELL01 :: ZZErrorEstimatorI_computeLocalStress(FloatArray& answer, FloatArray& sig)  
{
  // sig is global ShellForceMomentumTensor
  FloatMatrix globTensor(3,3);
  const FloatMatrix* GtoLRotationMatrix = plate->computeGtoLRotationMatrix();
  FloatMatrix LtoGRotationMatrix ;

  answer.resize(8); // reduced, local form
  LtoGRotationMatrix.beTranspositionOf(*GtoLRotationMatrix);

  // Forces
  globTensor.at(1, 1) = sig.at(1) ; //sxForce
  globTensor.at(1, 2) = sig.at(6) ; //qxyForce
  globTensor.at(1, 3) = sig.at(5) ; //qxzForce

  globTensor.at(2, 1) = sig.at(6) ; //qxyForce
  globTensor.at(2, 2) = sig.at(2) ; //syForce
  globTensor.at(2, 3) = sig.at(4) ; //syzForce

  globTensor.at(3, 1) = sig.at(5) ; //qxzForce
  globTensor.at(3, 2) = sig.at(4) ; //syzForce
  globTensor.at(3, 3) = sig.at(3) ; //szForce

  globTensor.rotatedWith(LtoGRotationMatrix);
  // Forces: now globTensoris transformed into local c.s

  // answer should be in reduced, local  form 
  answer.at(1) = globTensor.at(1, 1); //sxForce
  answer.at(2) = globTensor.at(2, 2); //syForce
  answer.at(3) = globTensor.at(1, 2); //qxyForce
  answer.at(7) = globTensor.at(2, 3); //syzForce
  answer.at(8) = globTensor.at(1, 3); //qxzForce


  // Moments:
  globTensor.at(1, 1) = sig.at(7) ; //mxForce
  globTensor.at(1, 2) = sig.at(12); //mxyForce
  globTensor.at(1, 3) = sig.at(11); //mxzForce

  globTensor.at(2, 1) = sig.at(12); //mxyForce
  globTensor.at(2, 2) = sig.at(8) ; //myForce
  globTensor.at(2, 3) = sig.at(10); //myzForce

  globTensor.at(3, 1) = sig.at(11); //mxzForce
  globTensor.at(3, 2) = sig.at(10); //myzForce
  globTensor.at(3, 3) = sig.at(9) ; //mzForce

  globTensor.rotatedWith (LtoGRotationMatrix); 
  // now globTensoris transformed into local c.s
 
  answer.at(4)  = globTensor.at(1, 1); //mxForce
  answer.at(5)  = globTensor.at(2, 2); //myForce
  answer.at(6) = globTensor.at(1, 2); //mxyForce

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
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}

void
TR_SHELL01  :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;
    IntArray map;
    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }
    
    if ( !this->giveMaterial()->isActivated(tStep) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    
    if ( ( !result ) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }
        //     //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

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
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,defScale);
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,defScale);
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[2].z = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(3,tStep,defScale);
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

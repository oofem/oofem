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
#include "classfactory.h"

#ifdef __OOFEG
 #include "node.h"
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {

REGISTER_Element( TR_SHELL01 );

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
    if (plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints() != membrane->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints()) {
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
  IntArray loc(9);
  FloatArray aux;

  answer.resize(18);
  answer.zero();

  plate->giveCharacteristicVector(aux, mtrx, mode, tStep);
  loc.setValues(9, 3,4,5, 9,10,11, 15,16,17);
  answer.assemble(aux, loc);

  membrane->giveCharacteristicVector(aux, mtrx, mode, tStep);
  loc.setValues(9, 1,2,6, 7,8,12, 13,14,18);
  answer.assemble(aux, loc);
}

void
TR_SHELL01 :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver accordind to mtrx
//
{
    IntArray loc(9);
    FloatMatrix aux;

    answer.resize(18,18);
    answer.zero();

    plate->giveCharacteristicMatrix(aux, mtrx, tStep);
    loc.setValues(9, 3,4,5, 9,10,11, 15,16,17);
    answer.assemble(aux, loc);

    membrane->giveCharacteristicMatrix(aux, mtrx, tStep);
    loc.setValues(9, 1,2,6, 7,8,12, 13,14,18);
    answer.assemble(aux, loc);
}

bool
TR_SHELL01 :: giveRotationMatrix(FloatMatrix &answer, EquationID eid)
{
    IntArray loc(9);
    FloatMatrix aux1, aux2;
    int i,j, ncol;

    bool t1 = plate->giveRotationMatrix(aux1, eid);
    bool t2 =  membrane->giveRotationMatrix(aux2, eid);

    if (t1 != t2) OOFEM_ERROR ("Transformation demand mismatch");
    
    if (t1) {
        ncol = aux1.giveNumberOfColumns();
        answer.resize(18,ncol);
        
        loc.setValues(9, 3,4,5, 9,10,11, 15,16,17);
        for (i=1; i<=9; i++) { // row index
            for (j=1; j<=ncol; j++) {
                answer.at(loc.at(i),j)=aux1.at(i,j);
            }
        }

        loc.setValues(9, 1,2,6, 7,8,12, 13,14,18);
        for (i=1; i<=9; i++) { // row index
            for (j=1; j<=ncol; j++) {
                answer.at(loc.at(i),j)=aux2.at(i,j);
            }
        }
    }

    return t1;
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
    } else if ( interface == SpatialLocalizerInterfaceType ) {
      return static_cast< SpatialLocalizerInterface * >( this );
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


//
// The element interface required by ZZNodalRecoveryModel
//


double 
TR_SHELL01 :: ZZRemeshingCriteriaI_giveCharacteristicSize() 
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

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
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
    if ( ( iores = this->membrane->saveContext(stream, mode, obj) ) != CIO_OK ) {
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
        this->compositeIR->setUpIntegrationPoints(_Triangle, plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(), _3dShell);
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

int 
TR_SHELL01 :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return plate->computeLocalCoordinates(lcoords, coords);
}

double
TR_SHELL01 :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray c(3); 
    c.zero();
    // evaluate element center
    for (int i=1; i<=3; i++) {
        c.add(*plate->giveNode(i)->giveCoordinates());
    }
    c.times(1./3.);
    return c.distance(coords);
}

void
TR_SHELL01 :: SpatialLocalizerI_giveBBox(FloatArray &bb0, FloatArray &bb1)
{

    FloatArray lt3(3), gt3(3); // global vector in the element thickness direction of lenght thickeness/2
    const FloatMatrix* GtoLRotationMatrix = plate->computeGtoLRotationMatrix();

    // setup vector in the element local cs. perpendicular to element plane of thickness/2 length
    lt3.at(1) = 0.0;
    lt3.at(2) = 0.0;
    lt3.at(3) = 1.0; //this->giveCrossSection()->give(CS_Thickness)/2.0; // HUHU
    // transform it to globa cs
    gt3.beTProductOf (*GtoLRotationMatrix, lt3);
    
    // use gt3 to construct element bounding box respecting true element volume
    
    FloatArray *coordinates, _c(3);

    for ( int i = 1; i <= this->giveNumberOfNodes(); ++i ) {
        coordinates = this->giveNode(i)->giveCoordinates();
        
        _c = *coordinates; _c.add(gt3);
        if (i == 1) {
            bb0 = bb1 = _c;
        } else {
            bb0.beMinOf(bb0, _c);
            bb1.beMaxOf(bb1, _c);
        }

        _c = *coordinates; _c.subtract(gt3);
        bb0.beMinOf(bb0, _c);
        bb1.beMaxOf(bb1, _c);
    }
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
        int nip = plate->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
        FloatArray a, v(12);
        v.zero();
        for (int _i=1; _i<= nip; _i++) {
            this->giveIPValue(a, plate->giveDefaultIntegrationRulePtr()->getIntegrationPoint(_i-1), IST_ShellForceMomentumTensor, tStep);
            v += a;
        }
        v.times(1./nip);
        v1 = v;
        v2 =v;
        v3 =v;
    }

    indx = context.giveIntVarIndx();

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

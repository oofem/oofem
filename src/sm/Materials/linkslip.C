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

#include "linkslip.h"
#include "latticelinearelastic.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"
// #ifdef __TM_MODULE
//  #include "latticetransportelement.h"
//  #include "latticetransmat.h"
//  #include "latticelammat.h"
//  #include "pore.h"
//  #include "discretetransportproblem.h"
// #endif

namespace oofem {
REGISTER_Material(LinkSlip);

/// constructor which creates a dummy material without a status and without random extension interface
LinkSlip :: LinkSlip(int n, Domain *d, double e0, double a1, double a2) : LatticeLinearElastic(n, d)
{
    eNormalMean = e0;
    alphaOne = a1;
    alphaTwo = a2;
}



LinkSlip :: ~LinkSlip()
//
// destructor
//
{}

int
LinkSlip :: hasMaterialModeCapability(MaterialMode mode)
{
  if ( mode == _3dMat ) {
        return 1;
    }

    return 0;
}


IRResultType
LinkSlip :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro


    LatticeLinearElastic :: initializeFrom(ir);

    //axial stiffness
    IR_GIVE_FIELD(ir, eNormalMean, _IFT_LinkSlip_eNormal); // Macro

    //Ratio of lateral to axial stiffness
    alphaOne = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaOne, _IFT_LinkSlip_alphaOne); // Macro

    //Parameter which limits the stress in slip direction.
    IR_GIVE_FIELD(ir, tauZero, _IFT_LinkSlip_tauZero); // Macro
    
    return IRRT_OK;
}


MaterialStatus *
LinkSlip :: CreateStatus(GaussPoint *gp) const
{
    LinkSlipStatus *answer = new LinkSlipStatus(gp);

    return answer;
}
  

void
LinkSlip :: giveRealStressVector(FloatArray &answer,
                                             GaussPoint *gp,
                                             const FloatArray &totalStrain,
                                             TimeStep *atTime)
{
    LinkSlipStatus *status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );

    //strain has the meanig of slip. Stress is force
    
    FloatArray strainVector;

    FloatArray tempPlasticStrain = status->givePlasticStrain();
    
    //    this->initGpForNewStep(gp);

    FloatMatrix stiffnessMatrix;
    this->giveStiffnessMatrix(stiffnessMatrix, ElasticStiffness, gp, atTime);

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );

    answer.resize(rsize);
    answer.zero();

    /*First component is the slip one for which the stress should be limited using plasiticity (frictional slip between fibre and matrix). The other components are kept elastic. */
    answer.at(1) = (totalStrain.at(1)-tempPlasticStrain.at(1))*stiffnessMatrix.at(1, 1);

    double f = fabs(answer.at(1))-tauZero;

    if(f>0){//plastic response.
      //Reduced stress by increasing plastic strain.
      tempPlasticStrain.at(1) = tempPlasticStrain.at(1) + sgn(answer.at(1))*f/stiffnessMatrix.at(1, 1);
      answer.at(1) = (totalStrain.at(1)-tempPlasticStrain.at(1))*stiffnessMatrix.at(1, 1);

      status->setTempCrackFlag(1);
    }

    //Compute the final stress components
    for ( int i = 2; i <= rsize; i++ ) { // only diagonal terms matter
        answer.at(i) =  stiffnessMatrix.at(i, i) * totalStrain.at(i);
    }

    //Set temp values in status needed for dissipation
    status->letTempPlasticStrainBe(tempPlasticStrain);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);

    double tempDissipation = status->giveDissipation();
    double tempDeltaDissipation;

    tempDeltaDissipation = computeDeltaDissipation(gp, atTime);
    
    tempDissipation += tempDeltaDissipation;

    //Set all temp values
    status->setTempDissipation(tempDissipation);
    status->setTempDeltaDissipation(tempDeltaDissipation);

    return;
}

void LinkSlip :: giveRandomParameters(FloatArray &param)
{
    param.resize(3);
    param.zero();
    param.at(1) = localRandomType;

    if ( localRandomType == 1 ) { //Gaussian
        param.at(2) = coefficientOfVariation;
    } else {
        OOFEM_ERROR("Error: Unknown local random type:\n randomtype 1 = Gaussian\n");
    }

    return;
}


Interface *
LinkSlip :: giveInterface(InterfaceType type)
{
    return NULL;
}



void
LinkSlip :: give3dStiffMtrx(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime)
{
 
  /* Returns elastic moduli in reduced stress-strain space*/
    answer.resize(6, 6);
    answer.zero();

    LatticeMaterialStatus *status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );
    FloatArray tempStrain(6);
    tempStrain = status->giveTempStrainVector();
    
    answer.at(1, 1) = 1.;
    answer.at(2, 2) = this->alphaOne; // shear
    answer.at(3, 3) = this->alphaOne; // shear
    
    answer.times(this->give(eNormal_ID, gp) * this->eNormalMean);

    
}

LinkSlipStatus :: LinkSlipStatus(GaussPoint *g) :  LatticeMaterialStatus(g)
{
}
  

void
LinkSlipStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
// builds new crackMap
//
{
    LatticeMaterialStatus :: initTempStatus();
    this->tempPlasticStrain = this->plasticStrain;
}
  

void
LinkSlip :: giveThermalDilatationVector(FloatArray &answer,
                                                    GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);

    //Option to add a eigendisplacement instead of strain
    double length = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveLength();
    alpha += this->cAlpha / length;

    answer.resize(6);
    answer.zero();

    answer.at(1) = alpha;
}

void
LinkSlipStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    LatticeMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "plasticStrain %.8e, dissipation %f, deltaDissipation %f, crackFlag %d\n", this->plasticStrain.at(1), this->dissipation, this->deltaDissipation, this->crackFlag);
    return;
}
  
double
LinkSlip :: give(int aProperty, GaussPoint *gp)
{
    double answer;
    if ( RandomMaterialExtensionInterface :: give(aProperty, gp, answer) ) {
        if ( answer < 0.1 ) { //Introduce cut off to avoid numerical problems
            answer = 0.1;
        } else if ( answer > 10 ) {
            answer = 10;
        }
        return answer;
    } else if ( aProperty == eNormal_ID ) {
        return 1.;
    } else {
        return LatticeLinearElastic :: give(aProperty, gp);
    }
}


double
LinkSlip :: computeDeltaDissipation(GaussPoint *gp,
				       TimeStep *atTime)
{
    LinkSlipStatus *status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );

    FloatArray plasticStrain = status->givePlasticStrain();
    FloatArray tempPlasticStrain = status->giveTempPlasticStrain();
    FloatArray tempStress = status->giveTempStressVector();

    double tempDeltaDissipation =  tempStress.at(1)*(tempPlasticStrain.at(1)-plasticStrain.at(1));
    return tempDeltaDissipation;
}

  
contextIOResultType
LinkSlipStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    contextIOResultType iores;
    // save parent class status
    
    LatticeMaterialStatus :: saveContext(stream, mode);
      
    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}


contextIOResultType
LinkSlipStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;
    
    LatticeMaterialStatus :: restoreContext(stream, mode);
    
    //FloatArrays
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}
 
void
LinkSlipStatus :: updateYourself(TimeStep *atTime)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium.
//
{
    LatticeMaterialStatus :: updateYourself(atTime);
    this->plasticStrain = this->tempPlasticStrain;
 }

  
  
int
LinkSlip :: giveIPValue(FloatArray &answer,
                                    GaussPoint *gp,
                                    InternalStateType type,
                                    TimeStep *atTime)
{
  LinkSlipStatus *status = static_cast< LinkSlipStatus * >( this->giveStatus(gp) );
    
  if ( type == IST_CrackStatuses ) {
        answer.resize(1);
	answer.zero();
        answer.at(1) = status->giveCrackFlag();
        return 1;
  } else if ( type == IST_DissWork ) {
        answer.resize(1);
	answer.zero();
        answer.at(1) = status->giveDissipation();
        return 1;
    } else if ( type == IST_DeltaDissWork ) {
    answer.resize(1);
	answer.zero();
        answer.at(1) = status->giveDeltaDissipation();
        return 1;
    }
    else if ( type == IST_CharacteristicLength ) {
        answer.resize(1);
	answer.zero();
        answer.at(1) = static_cast< LatticeStructuralElement * >( gp->giveElement() )->giveLength();
        return 1;
    } else {
      return LatticeLinearElastic :: giveIPValue(answer, gp, type, atTime);
    }
}
    
}

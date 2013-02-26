/*
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2000   Borek Patzak
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


#include "trabbonegrad3d.h"
#include "stressvector.h"
#include "strainvector.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

#include "sparsemtrx.h"
#include "isolinearelasticmaterial.h"
#include "dynalist.h"
#include "error.h"
#ifndef __MAKEDEPEND
 #include <math.h>
#endif


namespace oofem {

TrabBoneGrad3D :: TrabBoneGrad3D(int n, Domain *d) : TrabBone3D(n, d)
{
    l = 0.;
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneGrad3D :: ~TrabBoneGrad3D()
{}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////
int
TrabBoneGrad3D :: hasMaterialModeCapability(MaterialMode mode)
{
  if ( mode == _3dMatGrad )
      {
        return 1;
      }

    return 0;
}
void
TrabBoneGrad3D :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMatGrad_F:
    case _3dMatGrad:
      if ( form == PDGrad_uu ) 
	{
	  give3dMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
	  break;
        }
      else if ( form == PDGrad_ku ) 
	{
	  give3dKappaMatrix(answer, form, rMode, gp, atTime);
	  break;
        } 
      else if ( form == PDGrad_uk ) 
	{
	  give3dGprime(answer, form, rMode, gp, atTime);
	  break;
        } 
      else if ( form == PDGrad_kk ) 
	{
	  giveInternalLength(answer, form, rMode, gp, atTime);
	  break;
        }

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
        return;
    }

    return;
}

void
TrabBoneGrad3D :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
 
  double tempDam, beta, tempKappa, kappa;
  FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
  FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;
  TrabBoneGrad3DStatus *status = ( TrabBoneGrad3DStatus * ) this->giveStatus(gp);
  
  if ( mode == ElasticStiffness ) 
    {
      this->constructAnisoComplTensor(compliance);
      elasticity.beInverseOf(compliance);
      answer = elasticity;
    } 
  else if ( mode == SecantStiffness || ( (mode == TangentStiffness) && (status->giveNsubsteps()>1))) 
    {
      if (printflag)
	printf("secant\n");
      this->constructAnisoComplTensor(compliance);
      elasticity.beInverseOf(compliance);
      tempDam = status->giveTempDam();
      answer = elasticity;
      answer.times(1.0 - tempDam);
    } 
  else if ( mode == TangentStiffness ) 
    {
      kappa = status->giveKappa();
      tempKappa = status->giveTempKappa();
      tempDam = status->giveTempDam();  
      if ( tempKappa > kappa ) 
	{ // plastic loading
	  // Imports
	  tempEffectiveStress = * status->giveTempEffectiveStress();
	  plasFlowDirec = * status->givePlasFlowDirec();
	  SSaTensor = * status->giveSSaTensor();
	  beta = status->giveBeta();
	  // Construction of the dyadic product tensor
	  prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
	  // Construction of the tangent stiffness third term
	  double dam = status->giveDam();
	  // Construction of the tangent stiffness second term
	  tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
	  secondTerm.beDyadicProductOf(tempTensor2, prodTensor);
	  secondTerm.times(-( 1.0 - tempDam ) / beta);
	  // Construction of the tangent stiffness
	  tangentMatrix = SSaTensor;
	  tangentMatrix.times(1.0 - tempDam);
	  tangentMatrix.add(secondTerm);
	  if(tempDam > status->giveDam())
	    {
	      //double nlKappa =  status->giveTempStrainVector().at(7);
	      double nlKappa =  status->giveNlKappa();
	      kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
	      FloatArray tempEffStress =  *status->giveTempEffectiveStress();   
	      double gPrime = TrabBone3D::computeDamageParamPrime(kappa);     
	      // Construction of the tangent stiffness third term
	      thirdTerm.beDyadicProductOf(tempEffectiveStress, prodTensor);
	      thirdTerm.times( -expDam * critDam * exp(-expDam * kappa) );
	      thirdTerm.times((1.-mParam) / beta);
	      tangentMatrix.add(thirdTerm);
	    }
	  answer = tangentMatrix;
	  //	    answer.beTranspositionOf(tangentMatrix);	  
        }
      else
	{ // elastic behavior with damage
	  // Construction of the secant stiffness
	  this->constructAnisoComplTensor(compliance);
	  elasticity.beInverseOf(compliance);
	  answer = elasticity;
	  answer.times(1.0 - tempDam);
        }
    }

  double g = status ->giveDensG();
  if (g <= 0)
    {
 
      double factor = gammaL0 * __OOFEM_POW(rho,rL) + gammaP0 * __OOFEM_POW(rho,rP)*(tDens-1) * __OOFEM_POW(g,tDens-2);
      tangentMatrix.resize(6,6);
      tangentMatrix.zero();
      tangentMatrix.at(1,1) = tangentMatrix.at(1,2)=tangentMatrix.at(1,3) = 1;
      tangentMatrix.at(2,1) = tangentMatrix.at(2,2)=tangentMatrix.at(2,3) = 1;
      tangentMatrix.at(3,1) = tangentMatrix.at(3,2)=tangentMatrix.at(3,3) = 1;
      tangentMatrix.times(factor);
      answer.add(tangentMatrix);		
      }
  
  status->setSmtrx(answer);
}




void
TrabBoneGrad3D :: give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{

  answer.resize(1,6);
  answer.zero();
  if(mode == TangentStiffness)
    {
      TrabBoneGrad3DStatus *status = ( TrabBoneGrad3DStatus * ) this->giveStatus(gp); 
      double beta;
      FloatArray plasFlowDirec,prodTensor;
      FloatMatrix SSaTensor;
      
      double kappa = status->giveKappa();
      double tempKappa = status->giveTempKappa();
      double dKappa = tempKappa - kappa;
      
      if ( dKappa > 0.0 )
	{

	  plasFlowDirec = * status->givePlasFlowDirec();
	  SSaTensor = * status->giveSSaTensor();
	  beta = status->giveBeta();
	  prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
	  for(int i = 1; i<= 6; i++)
	    answer.at(1,i) = prodTensor.at(i);
	  answer.times(1./beta);
	} 
    }
}



void
TrabBoneGrad3D :: give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{

  answer.resize(6, 1);
  answer.zero();
  if(mode == TangentStiffness)
    {
      TrabBoneGrad3DStatus *status = ( TrabBoneGrad3DStatus * ) this->giveStatus(gp);
      double damage, tempDamage;
      double nlKappa, kappa;
      FloatArray tempEffStress;
      double gPrime;
      double tempKappa = status->giveTempKappa();
      damage = status->giveDam();
      //nlKappa =  status->giveTempStrainVector().at(7);
      nlKappa =  status->giveNlKappa();
      kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
      tempDamage = TrabBone3D::computeDamageParam(kappa);

      
      if ( ( tempDamage - damage ) > 0 ) 
	{
	  tempEffStress =  *status->giveTempEffectiveStress();   
	  for ( int i = 1; i <= 6; i++ ) 
	    answer.at(i, 1) = tempEffStress.at(i);
	  kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
	  gPrime = TrabBone3D::computeDamageParamPrime(kappa);
	  answer.times(gPrime * mParam);
	} 
    }
    
}

void
TrabBoneGrad3D :: giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.at(1, 1) = l;
}

void
TrabBoneGrad3D :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain, TimeStep *atTime)
{
    TrabBoneGrad3DStatus *status = ( TrabBoneGrad3DStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    MaterialMode plReturnMode;
    if ( gp->giveMaterialMode() == _3dMatGrad ) 
      {
        plReturnMode = _3dMat;
      }
    
    
    double tempDam;
    FloatArray tempEffStress, totalStress, locTotalStrain;

    int size = totalStrain.giveSize();
    status->letTempStrainVectorBe(totalStrain);
    status -> setNlKappa(totalStrain.at(7));
    locTotalStrain = totalStrain;
    locTotalStrain.resize(size - 1);
    
    TrabBone3D ::  performPlasticityReturn(gp, locTotalStrain,atTime, plReturnMode);
    
    double localCumPlastStrain = status->giveTempKappa();
    tempDam = computeDamage(gp, atTime);
    tempEffStress = *status->giveTempEffectiveStress();
    answer = ( 1 - tempDam ) * tempEffStress;
    if(densCrit != 0)
      {
	FloatArray densStress;
	computeDensificationStress(densStress, gp, totalStrain, atTime);
	answer.add(densStress);
      }
    
    size = tempEffStress.giveSize();
    answer.resize(size + 1);
    answer.at(size + 1) = localCumPlastStrain;

    status->setTempDam(tempDam);
    status->setTempEffectiveStress(tempEffStress);
    status->letTempStressVectorBe(answer);

    return;
}


void
TrabBoneGrad3D :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nlCumPlastStrain;
    TrabBoneGrad3DStatus *status = ( TrabBoneGrad3DStatus * ) this->giveStatus(gp);
    double localCumPlastStrain = status->giveTempKappa();
    FloatArray strain;
    /*  strain = status->giveTempStrainVector();
    int size = strain.giveSize();
    nlCumPlastStrain = strain.at(size);*/
    nlCumPlastStrain = status->giveNlKappa();
    kappa = mParam * nlCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


IRResultType
TrabBoneGrad3D :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    TrabBone3D :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, l, IFT_TrabBoneGrad3D_l, "l"); // Macro
    if ( l < 0.0 ) {
        l = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_TrabBoneGrad3D_m, "m"); // Macro

    return IRRT_OK;
}


TrabBoneGrad3DStatus :: TrabBoneGrad3DStatus(int n, Domain *d, GaussPoint *g) :
    TrabBone3DStatus(n, d, g)
{
  nlKappa = 0;
  

}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

TrabBoneGrad3DStatus :: ~TrabBoneGrad3DStatus()
{}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PRINTOUT
//

void
TrabBoneGrad3DStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
  TrabBone3DStatus:: printOutputAt(file, tStep);
}


void
TrabBoneGrad3DStatus :: initTempStatus()
{
   TrabBone3DStatus :: initTempStatus();
 }


void
TrabBoneGrad3DStatus :: updateYourself(TimeStep *atTime)
{
  StructuralMaterialStatus :: updateYourself(atTime);
  this->kappa = this->tempKappa;
  this->dam = this->tempDam;
  this->tsed = this->tempTSED;
  this->plasDef = this->tempPlasDef;
  nss = 1;
  //TrabBone3DStatus :: updateYourself(atTime);
}
} // end namespace oofem

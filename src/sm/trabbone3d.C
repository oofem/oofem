/* $Header: /home/cvs/bp/oofem/sm/src/isodamagemodel.C,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#include "trabbone3d.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "internalstatetype.h"
#include "contextioerr.h"

namespace oofem {

/////////////////////////////////////////////////////////////////
////////////////TRABECULAR BONE MATERIAL/////////////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
//

TrabBone3D :: TrabBone3D (int n, Domain *d) : StructuralMaterial (n,d)
{
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
//

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: MODE COMPATIBILITY
// returns whether receiver supports given mode

int
TrabBone3D :: hasMaterialModeCapability (MaterialMode mode)
{
  if ((mode == _3dMat)) return 1;
  return 0;
}

//
// END: MODE COMPATIBILITY
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF PLASTIC ENERGY DENSITY EVALUATION
//

void TrabBone3D::computePlasStrainEnerDensity (GaussPoint* gp, const FloatArray& totalStrain, const FloatArray& totalStress)
{
  TrabBone3DStatus *status = (TrabBone3DStatus*) this -> giveStatus (gp);

  double tsed, tempPSED, tempTSED, tempESED;
  FloatArray newTotalDef, oldTotalDef, tempPlasDef, oldStress;

  newTotalDef = totalStrain;
  oldTotalDef = status -> giveStrainVector();

  tsed = status -> giveTSED();
  tempPlasDef = *status -> giveTempPlasDef();
  oldStress = status -> giveTempStressVector();

  tempESED = dotProduct(0.5*(totalStrain-tempPlasDef),totalStress,6);
  tempTSED = tsed + dotProduct(0.5*(totalStrain-oldTotalDef),(totalStress+oldStress),6);
  tempPSED = tempTSED - tempESED;

  if (sqrt(tempPSED*tempPSED) < pow(10.,-16.))
    tempPSED = 0.;

  status -> setTempTSED(tempTSED);
  status -> setTempPSED(tempPSED);
}

//
// END:  SUBROUTINE OF PLASTIC ENERGY DENSITY EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: EVALUATION OF STIFFNESS MATRIX
//

void
TrabBone3D::give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                                   MatResponseForm form,MatResponseMode mode,GaussPoint* gp,
                                   TimeStep* atTime)
{
  TrabBone3DStatus *status = (TrabBone3DStatus*) this -> giveStatus (gp);

  double tempDam, beta, deltaKappa, tempKappa;
  FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
  FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;

  if (mode == ElasticStiffness)
  {
    this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
    elasticity.beInverseOf(compliance);

    answer = elasticity;
  }
  else if (mode == SecantStiffness)
  {
    this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
    elasticity.beInverseOf(compliance);
    tempDam = status -> giveTempDam();

    answer = elasticity;
    answer.times(1.0-tempDam);
  }
  else if (mode == TangentStiffness)
  {
    deltaKappa = status -> giveDeltaKappa();
    if (deltaKappa > 0.0)
    {
// Imports
      tempEffectiveStress = *status -> giveTempEffectiveStress();
      tempKappa = status -> giveTempKappa();
      tempDam = status -> giveTempDam();
      plasFlowDirec = *status -> givePlasFlowDirec();
      SSaTensor = *status -> giveSSaTensor();
      beta = status -> giveBeta();
// Construction of the dyadic product tensor
      prodTensor.beTProductOf(SSaTensor,plasFlowDirec);
// Construction of the tangent stiffness third term
      thirdTerm.beDyadicProductOf(tempEffectiveStress,prodTensor);
      thirdTerm.times(-expDam*critDam*exp(-expDam*tempKappa)/beta);
// Construction of the tangent stiffness second term
      tempTensor2.beProductOf(SSaTensor,plasFlowDirec);
      secondTerm.beDyadicProductOf(tempTensor2,prodTensor);
      secondTerm.times(-(1.0-tempDam)/beta);
// Construction of the tangent stiffness
      tangentMatrix = SSaTensor;
      tangentMatrix.times(1.0-tempDam);
      tangentMatrix.plus(secondTerm);
      tangentMatrix.plus(thirdTerm);

      answer = tangentMatrix;

    }
    else
    {
// Import of state variables
      tempDam = status -> giveTempDam();
// Construction of the tangent stiffness
      this->constructAnisoComplTensor(compliance, m1, m2, (3.0-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
      elasticity.beInverseOf(compliance);

      answer = elasticity;
      answer.times(1.0-tempDam);
    }
  }

  status -> setSmtrx(answer);

}

//
// END: EVALUATION OF STIFFNESS MATRIX
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF PLASTIC RETURN ALGORITHM
//

void
TrabBone3D :: performPlasticityReturn(GaussPoint* gp, const FloatArray& totalStrain)
{

  int flagLoop = 0, flag0 = 0;

  double tempKappa, deltaKappa, incKappa, halfSpaceCriterion, beta, tempScalar, norm;
  double plasCriterion, plasModulus, toSolveScalar, err, errTolerance = pow(10.,-15.);

  FloatArray tempPlasDef, tempEffectiveStress, incTempEffectiveStress, trialEffectiveStress, errLoop;
  FloatArray toSolveTensor, plasFlowDirec, yieldDerivative, tempTensor2, tensorFF_S;

  FloatMatrix elasticity, compliance, fabric, identity, SSaTensor, tempTensor4, normAdjust, normedFFFF, derivPlasFlowDirec;

  TrabBone3DStatus *status = (TrabBone3DStatus*) this -> giveStatus (gp);
  tempPlasDef = *status -> givePlasDef();
  tempKappa = status -> giveKappa();

  // Construction of the fourth order tensors

  this->constructAnisoComplTensor(compliance, m1, m2, (3-(m1+m2)), rho, eps0, nu0, mu0, expk, expl);
  this->constructIdentityTensor(identity);
  this->constructNormAdjustTensor(normAdjust);
  elasticity.beInverseOf(compliance);

  // Selection of the plastic criterion half-space

  trialEffectiveStress.beProductOf(elasticity,(totalStrain-tempPlasDef));
  tempEffectiveStress = trialEffectiveStress;

  halfSpaceCriterion=(pow(m1,-2.0*expq)*trialEffectiveStress.at(1)+pow(m2,-2.0*expq)*trialEffectiveStress.at(2)+pow((3.0-(m1+m2)),-2.0*expq)*trialEffectiveStress.at(3))/sqrt((pow(m1,-4.0*expq)+pow(m2,-4.0*expq)+pow((3.0-(m1+m2)),-4.0*expq)));

  if(halfSpaceCriterion < 0.0)
    this->constructAnisoFabricTensor(fabric, m1, m2, (3.0-(m1+m2)), rho, sig0Neg, chi0Neg, tau0, expk, expq);
  else 
    this->constructAnisoFabricTensor(fabric, m1, m2, (3.0-(m1+m2)), rho, sig0Pos, chi0Pos, tau0, expk, expq);

  // Evaluation of the plastic criterion

  tempTensor2.beProductOf(fabric,trialEffectiveStress);
  plasCriterion = sqrt(dotProduct(trialEffectiveStress,tempTensor2,6))-(1.0 + plasHardFactor*(1.0-exp(-tempKappa*expPlasHard)));

  if (plasCriterion <= errTolerance)
  {
    tempPlasDef = tempPlasDef;
    tempKappa = tempKappa;
    deltaKappa = 0.0;
  }
  else
  {

    // Initial values

    toSolveTensor.resize(6);
    toSolveScalar = plasCriterion;
    err = plasCriterion;
    deltaKappa = 0.0;
    SSaTensor = elasticity;

    // Construction of the norm 

    tempTensor4.beProductOf(normAdjust,fabric);
    normedFFFF.beProductOf(fabric,tempTensor4);
    tempTensor2.beProductOf(normedFFFF,tempEffectiveStress);
    norm = sqrt(dotProduct(tempEffectiveStress,tempTensor2,6));

    // Construction of the product FF.S

    tensorFF_S.beProductOf(fabric,tempEffectiveStress);

    // Construction of the direction of plastic flow

    plasFlowDirec = tensorFF_S;
    plasFlowDirec.times(1.0/norm);


      while (err > errTolerance)
      {

        plasModulus = -(plasHardFactor*expPlasHard*exp(-(tempKappa + deltaKappa)*expPlasHard));

	//*************************************
	//Evaluation of the Recursive Equations
	//*************************************

	tempTensor2.beProductOf(SSaTensor,plasFlowDirec);
	beta = dotProduct(plasFlowDirec,tempTensor2,6);
	beta -= sqrt(dotProduct(tempEffectiveStress,tensorFF_S,6))/norm*plasModulus;

	// Construction of the equation of Delta Kappa

	tempTensor2.beProductOf(SSaTensor,toSolveTensor);
	tempScalar = dotProduct(plasFlowDirec,tempTensor2,6);
	incKappa = (sqrt(dotProduct(tempEffectiveStress,tensorFF_S,6))/norm*toSolveScalar - tempScalar)/beta;
	deltaKappa += incKappa;

	// Construction of the equation of stress

	tempTensor2 = plasFlowDirec;
	tempTensor2.times(incKappa);
	tempTensor2 += toSolveTensor;
	incTempEffectiveStress.beProductOf(SSaTensor,tempTensor2);	
	tempEffectiveStress -= incTempEffectiveStress;

	//*************************************
	// Evaluation of the f and R
	//*************************************

	// Evaluation of the norm

	tempTensor2.beProductOf(normedFFFF,tempEffectiveStress);
	norm = sqrt(dotProduct(tempEffectiveStress,tempTensor2,6));

	// Construction of the direction of plastic flow

	tensorFF_S.beProductOf(fabric,tempEffectiveStress);
	plasFlowDirec = tensorFF_S;
	plasFlowDirec.times(1.0/norm);

	// Construction of the derivative of the plastic flow

	tempTensor2.beProductOf(normedFFFF,tempEffectiveStress);
	tempTensor4.beDyadicProductOf(tensorFF_S,tempTensor2);
	tempTensor4.times(-1.0/(norm*norm*norm));
	derivPlasFlowDirec = fabric;
	derivPlasFlowDirec.times(1.0/norm);
	derivPlasFlowDirec.plus(tempTensor4);

	// Construction of the gradient Nabla_S of R and SSa tensor

	tempTensor4 = derivPlasFlowDirec;
	tempTensor4.times(deltaKappa);
	tempTensor4.plus(compliance);
	SSaTensor.beInverseOf(tempTensor4);

	// Evaluation of R

        tempTensor2.beProductOf(compliance,(tempEffectiveStress-trialEffectiveStress));
	toSolveTensor = tempTensor2 + deltaKappa*plasFlowDirec;

	// Evaluation of f

	tempTensor2.beProductOf(fabric,tempEffectiveStress);
	toSolveScalar = sqrt(dotProduct(tempEffectiveStress,tempTensor2,6))-(1.0 + plasHardFactor*(1.0-exp(-(tempKappa + deltaKappa)*expPlasHard)));

	//*************************************
	// Evaluation of the error
	//*************************************

	errLoop = toSolveTensor;
	errLoop.resize(7);
	errLoop.at(7) = toSolveScalar;
        err = sqrt(dotProduct(errLoop,errLoop,7));
	flagLoop+=1;

        if (flagLoop > 10) flag0 =1;
        if (flagLoop > 1000) _error("No convergence of the stress return algorithm");

      }

    if (flag0 == 1)
      printf ("\nElement Number: %d (GP %d ), iterationNumber %d, deltaKappa %f ", gp->giveElement()->giveNumber(), gp->giveNumber(), flagLoop, deltaKappa);

    tempPlasDef -= -deltaKappa*plasFlowDirec;
    tempKappa += deltaKappa;

  status -> setBeta(beta);
  status -> setPlasFlowDirec(plasFlowDirec);
  status -> setSSaTensor(SSaTensor);

  }

  status -> setTempPlasDef(tempPlasDef);
  status -> setTempKappa(tempKappa);
  status -> setDeltaKappa(deltaKappa);
  status -> setTempEffectiveStress(tempEffectiveStress);

////  if (gp->giveNumber()==1)
////  {
////    printf ("\nParameter Number: %d ", gp->giveElement()->giveNumber());
////    printf ("\nFabric Tensor\n");
////    fabric.printYourself();
////    printf ("\nCompliance Tensor\n");
////    compliance.printYourself();
////    printf ("\nElasticity Tensor\n");
////    elasticity.printYourself();
////    printf ("\nInput Parameter (D_C, expDam, chi, expPlas)\n");
////    printf ("%f %f %f %f\n", critDam, expDam, plasHardFactor, expPlasHard);
////  }

}

//
// END: SUBROUTINE OF PLASTIC RETURN ALGORITHM
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: FUNCTION FOR DAMAGE PARAMETER
//

double  
TrabBone3D::computeDamageParam (double tempKappa, GaussPoint* gp)  
{
  double tempDam;
  if (tempKappa > 0.)
  {
    tempDam = critDam*(1.0-exp(-expDam*tempKappa));
  }
  else
  { 
    tempDam = 0.0;
  }
  return tempDam;
}

//
// END: FUNCTION FOR DAMAGE PARAMETER
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: FUNCTION FOR DAMAGE EVALUATION
//

double
TrabBone3D::computeDamage (GaussPoint* gp,  TimeStep* atTime)
{
  double tempKappa;

  computeCumPlastStrain(tempKappa, gp, atTime);

  double tempDam = computeDamageParam (tempKappa, gp);

  return tempDam;

}

//
// END: FUNCTION FOR DAMAGE EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF ALPHA EVALUATION
//

void TrabBone3D::computeCumPlastStrain (double& tempKappa, GaussPoint* gp, TimeStep* atTime)
{
  TrabBone3DStatus *status = (TrabBone3DStatus*) this -> giveStatus (gp);
  tempKappa = status -> giveTempKappa();
}

//
// END: SUBROUTINE OF ALPHA EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE OF DENSIFICATOR EVALUATION
//

void TrabBone3D::computeDensificationStress (FloatArray& answer, GaussPoint* gp, const FloatArray& totalStrain, TimeStep* atTime)
{
  double J;
  FloatArray Id(6);
  FloatMatrix C(3,3), invC(3,3), U(3,3), invU(3,3);

// The following development is valid only with the nonlinear geometries (Green-Lagrange: E=0.5(C-I))
  C.at(1,1) = 1. + 2.*totalStrain.at(1); 
  C.at(2,2) = 1. + 2.*totalStrain.at(2); 
  C.at(3,3) = 1. + 2.*totalStrain.at(3); 
  C.at(1,2) = C.at(2,1) = totalStrain.at(6);
  C.at(1,3) = C.at(3,1) = totalStrain.at(5);
  C.at(2,3) = C.at(3,2) = totalStrain.at(4);
  invC.beInverseOf(C);
  J = sqrt(C.giveDeterminant());

  if (J < JCrit)
  {
    Id.at(1) = invC.at(1,1);
    Id.at(2) = invC.at(2,2);
    Id.at(3) = invC.at(3,3);
    Id.at(4) = 2*invC.at(3,2);
    Id.at(5) = 2*invC.at(3,1);
    Id.at(6) = 2*invC.at(2,1);
    answer = gamDens*pow(J/JCrit-JCrit/J,tDens)*J*Id;
  }
  else answer.resize(6);

}

//
// END: SUBROUTINE OF ALPHA EVALUATION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
// returns real stress vector in 3d stress space of receiver according to 
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//

void
TrabBone3D :: giveRealStressVector (FloatArray& answer, MatResponseForm form, GaussPoint* gp, 
                                    const FloatArray& totalStrain, TimeStep* atTime)
{
  TrabBone3DStatus *status = (TrabBone3DStatus*) this -> giveStatus (gp);
  this->initGpForNewStep(gp);

  int i;
  double tempDam;
  FloatArray totalStress, effStress, densStress;

  performPlasticityReturn(gp, totalStrain);
  tempDam = computeDamage(gp, atTime);
  effStress = *status -> giveTempEffectiveStress();

  totalStress = (1-tempDam)*effStress;

  for(i=1;i<=6;i++){
     if (sqrt(totalStress.at(i)*totalStress.at(i))< pow(10.,-16.))
	totalStress.at(i)=0.;
    }

  computePlasStrainEnerDensity(gp, totalStrain, totalStress);

  if (JCrit != 0.) computeDensificationStress(densStress, gp, totalStrain, atTime);
  else densStress.resize(6);

  answer = totalStress + densStress;

  status -> setTempDam(tempDam);
  status -> letTempStrainVectorBe(totalStrain);
  status -> letTempStressVectorBe(answer);

////   double tempKappa = status -> giveTempKappa();
////   double deltaKappa = status -> giveDeltaKappa();
////   printf("\nInternal Variable, Element %d (GP %d )\n", gp->giveElement()->giveNumber(), gp->giveNumber());
////   printf("deltaKappa %10.16e , kappa %10.16e , damage %10.16e , nlKappa %10.16e \n", deltaKappa, tempKappa, tempDam, tempKappa);
////   printf("Plastic Stress");
////   effStress.printYourself();

}

//
// END: SUBROUTINE FOR EVALUATION OF TOTAL STRESS
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: MATRIX DEFINITION
void
TrabBone3D :: constructAnisoComplTensor (FloatMatrix& answer, const double m1, const  double m2, const  double m3, const  double rho, const  double eps0, const  double nu0, const  double mu0, const  double expk, const  double expl)
{
  answer.resize (6,6);

  answer.at(1,1)=1/(eps0*pow(rho,expk)*pow(m1,2*expl));
  answer.at(2,2)=1/(eps0*pow(rho,expk)*pow(m2,2*expl));
  answer.at(3,3)=1/(eps0*pow(rho,expk)*pow(m3,2*expl));
  answer.at(1,2)=-nu0/(eps0*pow(rho,expk)*pow(m1,expl)*pow(m2,expl));
  answer.at(2,1)=-nu0/(eps0*pow(rho,expk)*pow(m2,expl)*pow(m1,expl));
  answer.at(1,3)=-nu0/(eps0*pow(rho,expk)*pow(m1,expl)*pow(m3,expl));
  answer.at(3,1)=-nu0/(eps0*pow(rho,expk)*pow(m3,expl)*pow(m1,expl));
  answer.at(2,3)=-nu0/(eps0*pow(rho,expk)*pow(m2,expl)*pow(m3,expl));
  answer.at(3,2)=-nu0/(eps0*pow(rho,expk)*pow(m3,expl)*pow(m2,expl));
  answer.at(4,4)=1/(mu0*pow(rho,expk)*pow(m2,expl)*pow(m3,expl));
  answer.at(5,5)=1/(mu0*pow(rho,expk)*pow(m3,expl)*pow(m1,expl));
  answer.at(6,6)=1/(mu0*pow(rho,expk)*pow(m1,expl)*pow(m2,expl));

   return  ;
 }

void
TrabBone3D :: constructAnisoFabricTensor (FloatMatrix& answer, const double m1, const  double m2, const  double m3, const  double rho, const  double sig0, const  double chi0, const  double tau0, const  double expp, const  double expq)
{
  answer.resize (6,6);

  answer.at(1,1)=1/pow(sig0*pow(rho,expp)*pow(m1,2*expq),2);
  answer.at(2,2)=1/pow(sig0*pow(rho,expp)*pow(m2,2*expq),2);
  answer.at(3,3)=1/pow(sig0*pow(rho,expp)*pow(m3,2*expq),2);
  answer.at(1,2)=-chi0/pow(sig0*pow(rho,expp)*pow(m1,expq)*pow(m2,expq),2);
  answer.at(2,1)=-chi0/pow(sig0*pow(rho,expp)*pow(m2,expq)*pow(m1,expq),2);
  answer.at(1,3)=-chi0/pow(sig0*pow(rho,expp)*pow(m1,expq)*pow(m3,expq),2);
  answer.at(3,1)=-chi0/pow(sig0*pow(rho,expp)*pow(m3,expq)*pow(m1,expq),2);
  answer.at(2,3)=-chi0/pow(sig0*pow(rho,expp)*pow(m2,expq)*pow(m3,expq),2);
  answer.at(3,2)=-chi0/pow(sig0*pow(rho,expp)*pow(m3,expq)*pow(m2,expq),2);
  answer.at(4,4)=1/(pow(tau0*pow(rho,expp)*pow(m2,expq)*pow(m3,expq),2));
  answer.at(5,5)=1/(pow(tau0*pow(rho,expp)*pow(m3,expq)*pow(m1,expq),2));
  answer.at(6,6)=1/(pow(tau0*pow(rho,expp)*pow(m1,expq)*pow(m2,expq),2));

   return  ;
 }

void
TrabBone3D :: constructIdentityTensor (FloatMatrix& answer)
{
   int i, j;

  answer.resize (6,6);

  for (i=1 ; i<=6 ; i++)
  for (j=1 ; j<=6 ; j++)
   {
    if (i==j)
    answer.at(i,j)=1.;
    else
    answer.at(i,j)=0.;
   }

   return  ;
 }

void
TrabBone3D :: constructNormAdjustTensor (FloatMatrix& answer)
{
   int i, j;

  answer.resize (6,6);

  for (i=1 ; i<=3 ; i++)
  for (j=1 ; j<=3 ; j++)
   {
    if (i==j)
    answer.at(i,j)=1.;
    else
    answer.at(i,j)=0.;
   }
  for (i=4 ; i<=6 ; i++)
  for (j=4 ; j<=6 ; j++)
   {
    if (i==j)
    answer.at(i,j)=0.5;
    else
    answer.at(i,j)=0.;
   }

   return  ;
 }
//
// END: MATRIX DEFINITION
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: PARAMETERS OF INPUT FILE
//

IRResultType
TrabBone3D :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 /// Read material properties here

 IR_GIVE_FIELD (ir, eps0, IFT_TrabBone3D_eps0, "eps0"); // Macro
 IR_GIVE_FIELD (ir, nu0, IFT_TrabBone3D_nu0, "nu0"); // Macro
 IR_GIVE_FIELD (ir, mu0, IFT_TrabBone3D_mu0, "mu0"); // Macro
 IR_GIVE_FIELD (ir, expk, IFT_TrabBone3D_expk, "expk"); // Macro
 IR_GIVE_FIELD (ir, expl, IFT_TrabBone3D_expk, "expl"); // Macro

 IR_GIVE_FIELD (ir, m1, IFT_TrabBone3D_m1, "m1"); // Macro
 IR_GIVE_FIELD (ir, m2, IFT_TrabBone3D_m2, "m2"); // Macro
 IR_GIVE_FIELD (ir, rho, IFT_TrabBone3D_rho, "rho"); // Macro

 IR_GIVE_FIELD (ir, sig0Pos, IFT_TrabBone3D_sig0Pos, "sig0pos"); // Macro
 IR_GIVE_FIELD (ir, sig0Neg, IFT_TrabBone3D_sig0Neg, "sig0neg"); // Macro
 IR_GIVE_FIELD (ir, chi0Pos, IFT_TrabBone3D_chi0Pos, "chi0pos"); // Macro
 IR_GIVE_FIELD (ir, chi0Neg, IFT_TrabBone3D_chi0Neg, "chi0neg"); // Macro
 IR_GIVE_FIELD (ir, tau0, IFT_TrabBone3D_tau0, "tau0"); // Macro
 IR_GIVE_FIELD (ir, expq, IFT_TrabBone3D_expq, "expq"); // Macro

 IR_GIVE_FIELD (ir, plasHardFactor, IFT_TrabBone3D_plasHardFactor, "plashardfactor"); // Macro
 IR_GIVE_FIELD (ir, expPlasHard, IFT_TrabBone3D_expPlasHard, "expplashard"); // Macro

 IR_GIVE_FIELD (ir, expDam, IFT_TrabBone3D_expDam, "expdam"); // Macro
 IR_GIVE_FIELD (ir, critDam, IFT_TrabBone3D_critDam, "critdam"); // Macro

 gamDens = 0.;
 tDens = 1.;
 JCrit = 0.;
 IR_GIVE_OPTIONAL_FIELD (ir, gamDens, IFT_TrabBone3D_gamDens, "gamdens"); // Macro
 IR_GIVE_OPTIONAL_FIELD (ir, tDens, IFT_TrabBone3D_tDens, "tdens"); // Macro
 IR_GIVE_OPTIONAL_FIELD (ir, JCrit, IFT_TrabBone3D_JCrit, "jcrit"); // Macro

 return StructuralMaterial::initializeFrom (ir);
} 

//
// END: PARAMETERS OF INPUT FILE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: VTK Output
//

int 
TrabBone3D::giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime)
{
 TrabBone3DStatus* status = (TrabBone3DStatus*) this -> giveStatus (aGaussPoint);
 if (type == IST_DamageScalar) {
  answer.resize(1);
  answer.at(1) = status->giveTempDam();
  return 1;
 } else if (type == IST_PlasticStrainTensor) {
  answer.resize(6);
  answer = *status->giveTempPlasDef();
  return 1;
 } else if (type == IST_MaxEquivalentStrainLevel) {
  answer.resize(1);
  answer.at(1) = status->giveTempKappa ();
  return 1;
 } else if (type == IST_BoneVolumeFraction) {
  answer.resize(1);
  answer.at(1) = rho;
  return 1;
 } else if (type == IST_PlasStrainEnerDens) {
  answer.resize(1);
  answer.at(1) = status ->giveTempPSED();
  return 1;
 } else if (type == IST_ElasStrainEnerDens) {
  answer.resize(1);
  answer.at(1) = status ->giveTempTSED() - status ->giveTempPSED();
  return 1;
 } else if (type == IST_TotalStrainEnerDens) {
  answer.resize(1);
  answer.at(1) = status ->giveTempTSED();
  return 1;
 } else return StructuralMaterial::giveIPValue (answer, aGaussPoint, type, atTime);
}

InternalStateValueType 
TrabBone3D::giveIPValueType (InternalStateType type)
{
 if (type == IST_PlasticStrainTensor) return ISVT_TENSOR_S3;
 else if ((type == IST_DamageScalar)||(type == IST_MaxEquivalentStrainLevel)||(type == IST_BoneVolumeFraction)||(type == IST_PlasStrainEnerDens)||(type == IST_ElasStrainEnerDens)||(type == IST_TotalStrainEnerDens)) return ISVT_SCALAR;
 else return StructuralMaterial::giveIPValueType (type);
}

int 
TrabBone3D::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
{
 if (type == IST_DamageScalar) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_PlasticStrainTensor) {
  answer.resize (6);
  answer.at(1) = 1;
  answer.at(2) = 2;
  answer.at(3) = 3;
  answer.at(4) = 4;
  answer.at(5) = 5;
  answer.at(6) = 6;
  return 1;
 } else if (type == IST_MaxEquivalentStrainLevel) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_BoneVolumeFraction) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_PlasStrainEnerDens) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_ElasStrainEnerDens) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_TotalStrainEnerDens) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else 
  return StructuralMaterial::giveIntVarCompFullIndx (answer, type, mmode);
}

int
TrabBone3D::giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint)
{
 if ((type == IST_DamageScalar)||(type == IST_MaxEquivalentStrainLevel)||(type == IST_BoneVolumeFraction)||(type == IST_ElasStrainEnerDens)||(type == IST_PlasStrainEnerDens)||(type == IST_TotalStrainEnerDens)) return 1;
 else if (type == IST_PlasticStrainTensor) return 6;
 else return StructuralMaterial::giveIPValueSize (type, aGaussPoint);
}

//
// END: PARAMETERS OF INPUT FILE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
//////////////////TRABECULAR BONE STATUS/////////////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CONSTRUCTOR
// init state variables

TrabBone3DStatus::TrabBone3DStatus (int n, Domain*d, GaussPoint* g):StructuralMaterialStatus(n,d,g)
{

	kappa = 0.0;
	tempKappa = 0.0;
	tempDam = 0.0;
	tempPSED = 0.0;
	tsed = 0.0;
	tempTSED = 0.0;
	deltaKappa = 0.0;
	beta = 0.0;
	tempEffectiveStress.resize(6);
	plasDef.resize(6);
	tempPlasDef.resize(6);
	plasFlowDirec.resize(6);
	smtrx.resize(6,6);
	tangentMatrix.resize(6,6);
	SSaTensor.resize(6,6);
}

//
// END: CONSTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DESTRUCTOR
// 

TrabBone3DStatus::~TrabBone3DStatus () 
{
}

//
// END: DESTRUCTOR
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: DEFINITION OF "GIVE PARAMETER"
//

double
TrabBone3DStatus :: giveKappa()
{
  return kappa;
}

double
TrabBone3DStatus :: giveTempKappa()
{
  return tempKappa;
}

double
TrabBone3DStatus :: giveDam()
{
  return dam;
}

double
TrabBone3DStatus :: giveTempDam()
{
  return tempDam;
}

double
TrabBone3DStatus :: giveDeltaKappa()
{
  return deltaKappa;
}

double
TrabBone3DStatus :: giveTempPSED()
{
  return tempPSED;
}

double
TrabBone3DStatus :: giveTSED()
{
  return tsed;
}

double
TrabBone3DStatus :: giveTempTSED()
{
  return tempTSED;
}

double
TrabBone3DStatus :: giveBeta()
{
  return beta;
}

const FloatArray*
TrabBone3DStatus :: giveTempEffectiveStress()
{
  return &tempEffectiveStress;
}

const FloatArray*
TrabBone3DStatus :: givePlasFlowDirec()
{
  return &plasFlowDirec;
}

const FloatArray*
TrabBone3DStatus :: givePlasDef()
{
  return &plasDef;
}

const FloatArray*
TrabBone3DStatus :: giveTempPlasDef()
{
  return &tempPlasDef;
}

const FloatMatrix*
TrabBone3DStatus :: giveSSaTensor()
{
  return &SSaTensor;
}

//
// END: DEFINITION OF "GIVE PARAMETER"
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: OUTPUT
// print state to output stream

void 
TrabBone3DStatus :: printOutputAt  (FILE *file, TimeStep* tStep)
{
  StructuralMaterialStatus :: printOutputAt (file, tStep);
  fprintf (file,"status { ");
  fprintf (file,"plastrains: %f  %f  %f  %f  %f  %f",this->plasDef.at(1),this->plasDef.at(2),this->plasDef.at(3),this->plasDef.at(4),this->plasDef.at(5),this->plasDef.at(6));
  fprintf (file," , kappa %f , dam %f , esed %f , psed %f , tsed %f ",this->tempKappa, this->tempDam, this->tempTSED - this->tempPSED, this->tempPSED, this->tempTSED);
  fprintf (file,"}\n");
}

//
// END: OUTPUT
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
// initialize temporary state variables according to equilibriated state vars

void 
TrabBone3DStatus::initTempStatus ()
{
  StructuralMaterialStatus :: initTempStatus();
}

//
// END: INITIALIZE TEMP VARIABLE (UPDATED DURING ITERATIONS)
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STEP
// Called when equlibrium reached, set equilibriated vars according to temporary (working) ones.

void 
TrabBone3DStatus::updateYourself(TimeStep* atTime)
{
  StructuralMaterialStatus::updateYourself(atTime);
  this->kappa = this->tempKappa;
  this->dam = this->tempKappa;
  this->tsed = this->tempTSED;
  this->plasDef= this->tempPlasDef;
}

//
// END: SETS VARIABLE EQUAL TO TEMP VARIABLE AT THE END OF THE STEP
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERRUPT RESTART UTILITY - SAVE
//

contextIOResultType
TrabBone3DStatus::saveContext (DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // save parent class status
 if ((iores = StructuralMaterialStatus :: saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 // write a raw data
 //if (fwrite(&deltaKappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 //if (fwrite(&damage,sizeof(double),1,stream)!= 1) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

//
// END: INTERRUPT RESTART UTILITY - SAVE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: INTERRUPT RESTART UTILITY - RESTORE
//

contextIOResultType
TrabBone3DStatus::restoreContext(DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // read parent class status
 if ((iores = StructuralMaterialStatus :: restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 // read raw data 
 //if (fread (&deltaKappa,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 //if (fread (&damage,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

//
// END: INTERRUPT RESTART UTILITY - RESTORE
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
// BEGIN: CREATE STATUS
//

MaterialStatus *TrabBone3D::CreateStatus (GaussPoint* gp) const
{
  TrabBone3DStatus *status =
    new  TrabBone3DStatus (1, StructuralMaterial :: giveDomain(), gp) ;
  return status ;
}

//
// END: CREATE STATUS
/////////////////////////////////////////////////////////////////

} // end namespace oofem

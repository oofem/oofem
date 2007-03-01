/* $Header: /home/cvs/bp/oofem/sm/src/maxwellChM.C,v 1.5.4.1 2004/04/05 15:19:47 bp Exp $ */
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

// file: MaxwelllChM.C
 
#ifndef __MAKEDEPEND
#include <math.h>
#endif
#include "mathfem.h"
#include "maxwellChM.h"
#include "material.h"
#include "isolinearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matstatus.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"


MaxwellChainMaterial :: MaxwellChainMaterial (int n,Domain* d) : StructuralMaterial (n,d),
 EmjuVal(), relaxationTimes(), discreteTimeScale()
{
 nChainUnits = 0;
 relMatAge = 0.0;
 linearElasticMaterial = NULL;

 //timeOfStoredEVal = timeOfStoredEmjuVal = 0.;
 Eval = 0.;
 EmjuValTime = -1.0;
}


MaxwellChainMaterial:: ~MaxwellChainMaterial () 
{
 if (linearElasticMaterial) delete linearElasticMaterial;
}


int
MaxwellChainMaterial :: hasMaterialModeCapability (MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
 if ((mode == _3dMat) || (mode == _PlaneStress) || 
   (mode == _PlaneStrain) || (mode == _1dMat) ||
   (mode == _2dPlateLayer) || (mode == _2dBeamLayer) ||
   (mode == _3dShellLayer)) return 1;

 if ((mode == _3dBeam) || (mode == _2dBeam) || (mode == _3dShell) || (mode == _2dPlate))
   return 1;
 return 0;
}



void
MaxwellChainMaterial :: giveRealStressVector (FloatArray& answer, MatResponseForm form, 
                       GaussPoint* gp, 
                       const  FloatArray& totalStrain,
                       TimeStep* atTime) 
//
// returns total stress vector (in full or reduced form - form parameter)
// of receiver according to 
// previous level of stress and current
// strain increment at time atTime.
//
// 
// may be good idea for nonlinear materials overload this function in order to 
// capture strain - stress history for proper modelling receiver behaviour
// and in order to capture possible failure of material
//
{
 //  FloatMatrix* Material :: (*dfunc) (GaussPoint*, FloatArray* = NULL);
 FloatArray  stressIncrement, stressVector, strainIncrement, reducedStrain, eigenStrain; // *deltaEps0, *help ;
 //FloatArray *answer;
 FloatMatrix Binv;
 double Emodulus;
 MaxwellChainMaterialStatus *status = (MaxwellChainMaterialStatus*) this -> giveStatus (gp);

 
 this->initTempStatus(gp);
 this->initGpForNewStep(gp);

 // substract stress independent part (temperature+shrinkage)
 this->giveStressDependentPartOfStrainVector(reducedStrain, gp, totalStrain,
                       atTime, VM_Incremental);

 strainIncrement = reducedStrain;
 strainIncrement.substract (status-> giveStrainVector());
 

 if (status-> giveStressVector().giveSize()) {
  stressVector      = status-> giveStressVector();
 } else {
  stressVector.resize (strainIncrement.giveSize());
 }

 Emodulus = this-> giveEModulus (gp, atTime);
 this-> giveGeneralizationBInvMatrix(Binv, ReducedForm, gp, atTime);
 // deltaEps0 = this->GiveEigenStrainVector(FullForm, gp, atTime);
 // help = strainIncrementIn3d->GiveCopy();
 // help -> substract(deltaEps0);
 stressIncrement.beProductOf (Binv, strainIncrement);
 stressIncrement.times (Emodulus);
 stressVector.add (stressIncrement);

 status -> letTempStrainVectorBe (totalStrain);
 status -> letTempStressVectorBe (stressVector);
 
 // update shrinkage strain if needed
 if (this->hasIncrementalShrinkageFromulation()) {
  FloatArray shv;
  this->giveShrinkageStrainVector(shv, ReducedForm, gp, atTime, VM_Total);
  status->setShrinkageStrainVector(shv);
 }

 // delete strainIncrement;
 // delete Binv;
 // delete deltaEps0;
 // delete help;
 
 if (form == ReducedForm) {answer = stressVector; return;}
 
 FloatArray helpAnswer;
 ((StructuralCrossSection*) gp -> giveElement()->giveCrossSection())
  ->giveFullCharacteristicVector(answer, gp, stressVector);
 //delete stressVector;
 return ;

}




void
MaxwellChainMaterial :: computeDiscreteRelaxationFunction (FloatArray& answer, 
                              GaussPoint *gp,
                              const FloatArray& atTimes,
                              double t0, double tr)
{
/*
 This functions solves numerically integral equation of the form

 \varepsilon(t) = \int{0}{t} J(t, \Tau) d\sigma(\Tau) + \varepsilon_n(t),
 (where \varepsilon_n(t) is stress independent deformation) for the case
 where \varepsilon(t) = 1  is keept at constant value in time.
 
 input variables are:
  
 t0 - age of material, when load has been applied.
 tr - age of material, when relaxation has begun.
 sig0 - stress at time t0
 atTimes - at which times the relaxation function will be evaluated.
           WARNING: atTimes should be uniformly distributed in log time scale.
              and relatively dense (100 intervals) in order to achieve reasonable
        accuracy.
*/
 
 int i,k, size,nsteps,si;
 double taui, tauk, Jtrt0, totalDeltaSigma;
 double sig0;
 double sum;
 
 size = atTimes.giveSize();
 nsteps = size ;
 FloatArray deltaSigma (size);
 answer.resize (size); answer.zero();
 
// Jtrt0 = this->computeCreepFunction(gp, tr, t0);
// sig0  = 1./Jtrt0;
// si = 1;
//
  Jtrt0 = this->computeCreepFunction(gp, t0+atTimes.at(1), t0);
  sig0  = 1./ Jtrt0;
  answer.at(1) = sig0; 
 si = 2;
// 
 totalDeltaSigma = 0.;
 for (k=si; k<= nsteps; k++) {
  for (sum=0., i=si; i<= k-1; i++) {
   if (i==1) {taui = 0.5*(t0+atTimes.at(i) + tr );}
   else {taui = t0+0.5*(atTimes.at(i) + atTimes.at(i-1));}
   sum += deltaSigma.at(i) * this->computeCreepFunction (gp, t0+atTimes.at(k), taui);
  }
  if (k==1) {tauk = 0.5*(t0+atTimes.at(k) + tr );}
  else {tauk = t0+0.5*(atTimes.at(k) + atTimes.at(k-1));}
  deltaSigma.at(k) = (sig0*(this->computeCreepFunction (gp, t0+atTimes.at(k),t0)-Jtrt0)
         -sum) / this->computeCreepFunction (gp, t0+atTimes.at(k),tauk);
  
  totalDeltaSigma += deltaSigma.at(k);
  answer.at(k) = sig0 - totalDeltaSigma;
 }
 
 // delete working arrays
 return ;
}


void
MaxwellChainMaterial :: generateLogTimeScale (FloatArray& answer, double from, double to, 
                       int nsteps,
                       int fromIncluded)
{
/* 
 function generates discrete times starting from time "from" to time "to"
   uniformely distributed in log time scale. The time interval (to-from) is
 divided to nsteps intervals. if from is wanted to be included in return
 array fromIncluded must be set to 1, otherwise we return times starting
 from ("from" + first increment);
 */
 int size = nsteps+1;
 int i,i0 = 1;
 double help, almostlogZero = -1.0;

 if (fromIncluded) {
  size ++; 
  i0 = 2;
  // ensure that really fromIncluded == 1
  fromIncluded = 1;
 }

 answer.resize (size); answer.zero();

 help = (log10(to-from) - almostlogZero) / (double) nsteps;
 for (i=0; i<= nsteps; i++)
  answer.at(i+i0) = __OOFEM_POW(10.,(almostlogZero + help * i)) + from ;
 if(fromIncluded) answer.at(1) = from;

/* 
 help = (fabs(almostlogZero) + log(to-from)) / (double) nsteps;
 for (i=1; i<= nsteps; i++) {
  answer->at(i+i0) = exp(almostlogZero + help * i) + from;
 }

 if(fromIncluded) answer->at(1) = from;
*/
 return ;
}


const FloatArray&
MaxwellChainMaterial :: giveDiscreteTimes ()
{
/*
 Function returns generated relative time scale (iuniformely distributed 
 in log time scale).
 Because of multiple use of this data, once generated, they are
 stored 
*/
   if (discreteTimeScale.isNotEmpty()) return discreteTimeScale;

 double endTime;
 endTime = this->giveEndOfTimeOfInterest() ;

 this->generateLogTimeScale (discreteTimeScale, 0., endTime, MNC_NPOINTS-1,0);
 return discreteTimeScale;
}


void
MaxwellChainMaterial :: computeCharCoeficients (FloatArray& answer, GaussPoint* gp,
                        double atTime )
{
/*
 function computes matrix of characteristic components of dirichlet series
 (computes Emju values).
 
 i-th value in returned array corresponds to coefficient for i-th Maxwell unit
 These values are obtained using least squared method by 
 minimizing following functional (atTime = t_0)
 
 $$ F=\sum^{k}_{r=1} \left[ \sum^{N}_{\mju=1} E_m(t_0) \exp^{-(t_r-t_0)/\tau_{\mju}
    - \bar{R}(t_r, t_0) \right]^2 = min $$

 
 INPUTS:
  atTime = age of material, when load has been applied.
  
 TASKS:

*/
 int i,j,r,rSize;
 //double endTime;
 double taui, tauj, sum, tti, ttj, sumRhs;
 FloatArray rhs (this->nChainUnits), discreteRelaxFunctionVal ;
 FloatMatrix A (this->nChainUnits, this->nChainUnits);


 // endTime = this->giveEndOfTimeOfInterest()+atTime ;
 // first generate discrete sampling points used to aproximate Relaxation function.
 const FloatArray& rTimes = this->giveDiscreteTimes ();
 rSize = rTimes.giveSize();
 // compute discrete values of relaxation function at times rTimes
 // from creep function.
 // (or we can use direct call to relaxation function, if it is available).
 this->computeDiscreteRelaxationFunction (discreteRelaxFunctionVal,
                      gp, rTimes,
                      atTime,
                      atTime);

 // assemble equations for computing unknown char. components of series
 for (i=1; i<= this->nChainUnits; i++) {
  taui = this->giveRelaxationTime(i);
  for (j=1; j<= this->nChainUnits; j++) {
   tauj = this->giveRelaxationTime(j);
   for (sum=0., r=1; r<= rSize; r++) {
    tti = __OOFEM_POW((atTime+rTimes.at(r))/taui,giveRelaxationTimeExponent(i)) - 
     __OOFEM_POW(atTime/taui,giveRelaxationTimeExponent(i));
    ttj = __OOFEM_POW((atTime+rTimes.at(r))/tauj,giveRelaxationTimeExponent(j)) - 
     __OOFEM_POW(atTime/tauj,giveRelaxationTimeExponent(j));
    sum += exp(-tti-ttj);
   }
   A.at(i,j) = sum; 
  }
  // assemble rhs
  for (sumRhs=0., r=1; r<= rSize; r++) {
   tti = __OOFEM_POW((atTime+rTimes.at(r))/taui,giveRelaxationTimeExponent(i)) -
    __OOFEM_POW(atTime/taui,giveRelaxationTimeExponent(i));
   sumRhs += exp(-tti) * discreteRelaxFunctionVal.at(r);
  }
  rhs.at(i)=sumRhs;
 }
 // solve assembled linear system
 A.solveForRhs (rhs, answer);
 
 return ;

}


void
MaxwellChainMaterial :: giveGeneralizationBInvMatrix (FloatMatrix& answer,
                           MatResponseForm form,
                           GaussPoint* gp,
                           TimeStep* tStep)
{
/*
 Returns generalization matrix, which based on principle of superposition
 (we have in mind Boltzmann one) this generalizes 1-d stiffness to
 stress-strain space given by gp.

 we assume, that Poisson ratio is constant (this is valid only for
 problems at constant humidity).

 This matrix is identical to li-iso stiffness matrix of stress strain space 
 defined for this gp, taking E = 1;

 We call crossSection, beacuse we need Binv matrix also for modes
 defined at the crosssection level.
 (Hidenn stresses are are stored in MatStatus level, but they are
 defined by stressStrain mode in gp, which may be of type suuported 
 at the crosssection level).
*/

 ((StructuralCrossSection*)gp->giveCrossSection())
  ->giveCharMaterialStiffnessMatrixOf (answer,
                     form, TangentStiffness,gp,
                     this->giveLinearElasticMaterial(),
                     tStep);
 return;
}

void
MaxwellChainMaterial :: giveGeneralizationBMatrix (FloatMatrix& answer,
                          MatResponseForm form,
                          GaussPoint* gp,
                          TimeStep* tStep)
{
/*
 Returns generalization matrix, which based on principle of superposition
 (we have in mind Boltzmann one) this generalizes 1-d stiffness to
 stress-strain space given by gp.

 we assume, that Poisson ratio is constant (this is valid only for
 problems at constant humidity).

 This matrix is identical to li-iso compliance matrix of stress strain space 
 defined for this gp, taking E = 1;
*/

 ((StructuralCrossSection*)gp->giveCrossSection())->
  giveCharMaterialComplianceMatrixOf (answer, form, TangentStiffness,gp,
                    this->giveLinearElasticMaterial(),
                    tStep);
 return;
}



double
MaxwellChainMaterial :: giveEModulus (GaussPoint* gp, TimeStep* atTime)
{
/*
 this function returns time dependent E'' modulus  for given time.
 more precisely, E" is evaluated at time (atTime-"1/2") 
 because we assume that change of loading is applied at the midle of
 time interval (atTime-1, atTime).
 E" may also be dependent on specimen geometry (gp - dependence).
 
 E" is stored for further expected requests from other gaussPoints.

 note: time -1 refer to previous time                           (")
*/
 int i;
 double lambdaMju, Emju, deltaYmju;
 double E = 0.0;

 // if (!eValueValid(atTime->giveTime())) {
 this->updateEmjuModuluses (gp, relMatAge + (atTime->giveTime()-0.5*atTime->giveTimeIncrement())/timeFactor);
 for (i=1; i<= nChainUnits; i++) {
  deltaYmju = atTime->giveTimeIncrement()/ timeFactor / this->giveRelaxationTime(i);
  if (deltaYmju <=0.0) deltaYmju = 1.e-3;
  deltaYmju = __OOFEM_POW(deltaYmju, this->giveRelaxationTimeExponent(i));
  
  lambdaMju = (1.0-exp(-deltaYmju)) / deltaYmju;
  Emju      = this->giveEmjuModulus(i); // previously updated by updateEmjuModuluses
  
  E+= lambdaMju * Emju;
  
 }
 Eval = E;
 return Eval;
}


double 
MaxwellChainMaterial :: giveEmjuModulus (int iChain)
{
/* returns value of iChain EMju modulus, previosly computed by updateEmjuModuluses()
 function
*/
 return EmjuVal.at(iChain);
}


void
MaxwellChainMaterial :: updateEmjuModuluses (GaussPoint* gp, double atTime)
{
/*
 Computes char. coefficients of dirichlet series - Emju values
 for relaxation function.

 INPUTS: 
 
   atTime - age of material, when load has been applied.
 
 DESCRIPTION:
   we store  computed values, because they will be used by other material points in subsequent 
   calculations. Their computation is very costly.
   for active time step, there should be requsts only to Emju values in time
   (currentTime+prevTime)*0.5. Previous values are not needed.

*/
 // compute new values and store them in temporary array for further use
  if (fabs(atTime-EmjuValTime) > TIME_DIFF) {
    //if (EmjuVal) delete EmjuVal;
    if (atTime < 0) this->computeCharCoeficients (EmjuVal, gp, 1.e-3);
    else this->computeCharCoeficients (EmjuVal, gp, atTime);
    EmjuValTime = atTime;
  } 
}


void 
MaxwellChainMaterial :: computeTrueStressIndependentStrainVector (FloatArray& answer, 
                                 GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
{
 // computes true stress independent strain vector, containing temperature and shrinkage effects
 FloatArray e0;
 
 this->giveShrinkageStrainVector(answer, ReducedForm, gp, stepN, mode);
 // take into account temperature
 StructuralMaterial::computeStressIndependentStrainVector (e0, gp, stepN, mode);
 answer.add(e0);
}


void
MaxwellChainMaterial :: computeStressIndependentStrainVector (FloatArray& answer, 
                               GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
//
// returns initial strain vector induced by stress independent effects
// like temperatue or shrinkage.
// takes into account form of load vector assumed by engngModel (Incremental or Total Load form).
// 
{
 FloatArray e0;
 this->computeTrueStressIndependentStrainVector (answer, gp, stepN, mode);
 // add eigen strain
 this->giveEigenStrainVector(e0, ReducedForm, gp, stepN, mode);
 if (e0.giveSize()) {
  answer.add(e0);
 }
}



void
MaxwellChainMaterial :: giveEigenStrainVector (FloatArray& answer, MatResponseForm form,
                        GaussPoint* gp, TimeStep* atTime, ValueModeType mode)
{
/*
 Return vector of increments of eigen strains which are induced by stress history
 at time atTime.
*/
 int i;
 double E ;
 double deltaYmju;
 // double deltaEps=0.0;
 FloatArray *sigmaMju, help, reducedAnswer;
 FloatMatrix B;
 MaxwellChainMaterialStatus *status = (MaxwellChainMaterialStatus*) this -> giveStatus (gp);
 
 if (mode==VM_Incremental) {
  this -> giveGeneralizationBMatrix (B, ReducedForm, gp, atTime);
  //reducedAnswer = new FloatArray (B.giveNumberOfRows());
  reducedAnswer.resize (B.giveNumberOfRows());
  
  for (i=1; i<= nChainUnits; i++) {
   deltaYmju = atTime->giveTimeIncrement() / timeFactor / this->giveRelaxationTime(i);
   deltaYmju = __OOFEM_POW(deltaYmju, this->giveRelaxationTimeExponent(i));
   sigmaMju  = status->giveHiddenStressVector (i);
   if (sigmaMju) {
    help.beProductOf (B, *sigmaMju);  // B can be moved before sum !!!
    help.times (1.0-exp(-deltaYmju));
    reducedAnswer.add (help);  
    // delete help;
   }
  }
  E = this -> giveEModulus (gp,atTime);
  reducedAnswer.times (1.0/E);
  
  // delete B;
  
  if (form == ReducedForm) {answer =  reducedAnswer; return ;}
  ((StructuralCrossSection*)gp->giveCrossSection())->
   giveFullCharacteristicVector(answer, gp, reducedAnswer) ;
  // delete reducedAnswer;
 } else {
  /* error - not implemented now */
  _error ("giveEigenStrainVector - mode is not supported");
 } 
 
 return ;
}


double
MaxwellChainMaterial :: giveRelaxationTime (int i)
{
/* returns relaxation time for i-th unit */

 if (!relaxationTimes.isNotEmpty()) this->computeRelaxationTimes();
 if ((i <= 0) || (i > nChainUnits)) 
  _error ("giveRelaxationTime - no such unit defined");
 return relaxationTimes.at(i);
}

void
MaxwellChainMaterial :: computeRelaxationTimes()
{
/*
 This function generates a discrete relaxation times
 accordig to rules, which garantee good aproximation
 of relaxation function by Dirichlet series

 First relaxation time Tau(1) is choosen to be equal to 0.1 day
 The biggest relaxation time Tau(n) is choosem 1.0e30
 (to approxime precisely also relaxation in time infinity.
 the second biggest time Tau(n-1) is choosen as 0.75 tmax, where tmax is 
 lifetime of structure or end of time of interest.
 times  Tau(2) .. Tau(n-2) are defined by uniform division to n-2 steps in
 log scale. It is necesary to check condition a <= 10, where Tau(k) = a Tau(k-1)


*/
# define a 10.

 int size, nsteps, i;
 double endTime, Taun1, Tau1, help;

 endTime = this->giveEndOfTimeOfInterest()+relMatAge;
 Taun1 = 0.75 * endTime; 
// Tau1  = 0.1;
 Tau1  = 0.1;
 nsteps = (int) ((log (Taun1) - log(Tau1)) / log(a) + 1.);
 if (nsteps < 8) nsteps = 8;
 this->nChainUnits = size = nsteps + 2;
 
 this->relaxationTimes.resize (size);
 
 relaxationTimes.at(1) = Tau1;
 relaxationTimes.at(size) = 1.0e10;
 help = (log (Taun1) - log(Tau1)) / (double) nsteps;
 // real used "a" is exp(help);
 for (i=1; i<= nsteps; i++) {
  relaxationTimes.at(i+1) = exp(log(Tau1) + help *i);
 }

 return ;
}


void
MaxwellChainMaterial :: giveCharacteristicMatrix (FloatMatrix& answer,
                         MatResponseForm form,
                         MatResponseMode mode,
                         GaussPoint* gp,
                         TimeStep* atTime)
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 this -> giveGeneralizationBInvMatrix(answer, form, gp, atTime);
 answer.times(this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial :: give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                            MatResponseForm form, MatResponseMode mode, 
                            GaussPoint* gp,
                            TimeStep *atTime) 
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 //FloatMatrix *answer;
 this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix(answer, form,mode,gp,
                                  atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial :: givePlaneStressStiffMtrx (FloatMatrix& answer,
                         MatResponseForm form, MatResponseMode mode,
                         GaussPoint* gp,
                         TimeStep *atTime)
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 //FloatMatrix *answer;
 this->giveLinearElasticMaterial()->givePlaneStressStiffMtrx(answer, form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}

void
MaxwellChainMaterial:: givePlaneStrainStiffMtrx (FloatMatrix& answer,
                         MatResponseForm form, MatResponseMode mode,
                         GaussPoint* gp,
                         TimeStep *atTime)
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 // FloatMatrix *answer;
 this->giveLinearElasticMaterial()->givePlaneStrainStiffMtrx(answer,form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial :: give1dStressStiffMtrx (FloatMatrix& answer,
                        MatResponseForm form, MatResponseMode mode,
                        GaussPoint* gp,
                        TimeStep *atTime)
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 // FloatMatrix *answer;
 this->giveLinearElasticMaterial()->give1dStressStiffMtrx(answer,form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial :: give2dBeamLayerStiffMtrx (FloatMatrix& answer,
                         MatResponseForm form, MatResponseMode mode,
                         GaussPoint* gp,
                         TimeStep *atTime)
{
 //
 // Returns characteristic material stiffness matrix of the receiver
 //
 // FloatMatrix *answer;
 this->giveLinearElasticMaterial()->give2dBeamLayerStiffMtrx(answer,form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial::give2dPlateLayerStiffMtrx(FloatMatrix& answer,
                        MatResponseForm form, MatResponseMode mode, 
                        GaussPoint* gp,
                        TimeStep *atTime)
{
//
// Returns characteristic material stiffness matrix of the receiver
//
 //FloatMatrix *answer;
 this->giveLinearElasticMaterial()->give2dPlateLayerStiffMtrx(answer,form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}


void
MaxwellChainMaterial::give3dShellLayerStiffMtrx(FloatMatrix& answer,
                        MatResponseForm form, MatResponseMode mode,
                        GaussPoint* gp,
                        TimeStep *atTime)
{
 //
 // Returns characteristic material stiffness matrix of the receiver
 //
 //FloatMatrix *answer;
 this->giveLinearElasticMaterial()->give3dShellLayerStiffMtrx(answer,form,mode,gp,atTime);
 answer.times (this->giveEModulus (gp,atTime));
 return ;
}



void
MaxwellChainMaterial :: updateYourself (GaussPoint* gp, TimeStep* tNow)
{
/*
 Updates hidden stresses used to effectively trace the load history
*/


 int i;
 double deltaYmju, Emju, lambdaMju;
 FloatArray  help, *ithHiddenStressVector, deltaEps0, help1;
 FloatMatrix Binv;
 MaxwellChainMaterialStatus *status = 
  (MaxwellChainMaterialStatus*) this -> giveStatus (gp);
 
 this -> giveGeneralizationBInvMatrix (Binv, ReducedForm, gp, tNow);
 //help = status->giveStrainIncrementVector();
 help = status->giveTempStrainVector();
 help.substract (status->giveStrainVector());

 // Substract stress independent part of strain
 this->computeTrueStressIndependentStrainVector(deltaEps0, gp, tNow, VM_Incremental);
 //this->giveThermalDilatationVector(deltaEps0, gp, tNow);
 if (deltaEps0.giveSize()) help.substract(deltaEps0);
 help1.beProductOf (Binv, help);
 //delete help;  // delete deltaEps0; // delete Binv;

 this->updateEmjuModuluses (gp, relMatAge + (tNow->giveTime()-0.5*tNow->giveTimeIncrement()) / timeFactor);

 for (i=1; i<= nChainUnits; i++) {

  deltaYmju = tNow->giveTimeIncrement() / timeFactor / this->giveRelaxationTime(i);
  deltaYmju = __OOFEM_POW(deltaYmju, this->giveRelaxationTimeExponent(i));
   
  lambdaMju = (1.0-exp(-deltaYmju)) / deltaYmju;
  Emju      = this->giveEmjuModulus(i);

  ithHiddenStressVector = status->giveHiddenStressVector(i);
  help = help1; help.times(lambdaMju*Emju);
  if (ithHiddenStressVector) {
    ithHiddenStressVector->times(exp(-deltaYmju));
    ithHiddenStressVector->add(&help);
  } else {
    status->letHiddenStressVectorBe (i,help.GiveCopy());
  }

  // clen up space allocated in loop 
  //delete help;
 }
 //delete help1;
 
 // We call MaxwellChainMaterialStatus->updateYourself()
 status->updateYourself(tNow);
 
}

MaterialStatus*
MaxwellChainMaterial:: CreateStatus (GaussPoint* gp) const
/* 
 creates new  material status  corresponding to this class
*/
{
 return new MaxwellChainMaterialStatus (1, this->giveDomain(), gp, nChainUnits);
}


IRResultType
MaxwellChainMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

//   if (this->readWhetherHas("d")) {
//      value = this -> read("d") ;
//      propertyDictionary -> add('d',value) ;}
 
 StructuralMaterial::initializeFrom(ir);
 IR_GIVE_FIELD (ir, nu, IFT_MaxwellChainMaterial_n, "n"); // Macro
 IR_GIVE_FIELD (ir, relMatAge, IFT_MaxwellChainMaterial_relmatage, "relmatage"); // Macro
 this->endOfTimeOfInterest = -1.0;
 IR_GIVE_OPTIONAL_FIELD (ir, endOfTimeOfInterest, IFT_MaxwellChainMaterial_endoftimeofinterest, "endoftimeofinterest"); // Macro
 IR_GIVE_FIELD (ir, timeFactor, IFT_MaxwellChainMaterial_timefactor, "timefactor"); // solution time/timeFactor should give time in reqired scale


 //nChainUnits = readInteger (initString,"nChainUnits");
 //if (nChainUnits < 4) nChainUnits = 4;
 this->computeRelaxationTimes();  // sets up nChainUnits variable

 //this -> giveLinearElasticMaterial() -> instanciateFrom(ir);
 
 return IRRT_OK;
}

LinearElasticMaterial*
MaxwellChainMaterial ::  giveLinearElasticMaterial ()
{
 if (linearElasticMaterial == NULL) {
  linearElasticMaterial = new IsotropicLinearElasticMaterial (this->giveNumber(),
                      this->giveDomain(),
                      1.0, this-> nu);
 }
 return linearElasticMaterial;
}

double
MaxwellChainMaterial ::giveEndOfTimeOfInterest() 
{
 if (this->endOfTimeOfInterest > 0.0) return this->endOfTimeOfInterest;
 else this->endOfTimeOfInterest = this->giveDomain()->giveEngngModel()->giveEndOfTimeOfInterest()/timeFactor;

 return this->endOfTimeOfInterest;
}

contextIOResultType
MaxwellChainMaterial :: saveContext (FILE* stream, void *obj)
//
// saves full status for this material, also invokes saving
// for sub-objects of this (yieldcriteria, loadingcriteria, linearElasticMaterial)
// which can have their own statuses stored in gp.
{
 contextIOResultType iores;

 if (stream == NULL) _error ("saveContex : can't write into NULL stream");
 if ((iores = Material :: saveContext (stream, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = giveLinearElasticMaterial()->saveContext(stream,obj)) != CIO_OK) THROW_CIOERR(iores);
 return CIO_OK;
}


contextIOResultType 
MaxwellChainMaterial :: restoreContext (FILE* stream, void *obj)
// 
//
// resaves full status for this material, also invokes saving
// for sub-objects of this (yieldcriteria, loadingcriteria, linearElasticMaterial)
// which can have their own statuses stored in gp.

//
{
 contextIOResultType iores;

  if ((iores = Material :: restoreContext( stream, obj)) != CIO_OK) THROW_CIOERR(iores);
  // invoke possible restoring of statuses for yield conditions and submaterials
  if ((iores = giveLinearElasticMaterial()->restoreContext(stream,obj)) != CIO_OK) THROW_CIOERR(iores);
 return CIO_OK;
}



/****************************************************************************************/

MaxwellChainMaterialStatus :: MaxwellChainMaterialStatus (int n, Domain *d,
                             GaussPoint* g, int nUnits)
: StructuralMaterialStatus(n,d,g), shrinkageStrain() {
// constructor
 int i;
 
 nMaxwellUnits = nUnits;
 
 hiddenStreses = new FloatArray* [nUnits];
 for (i=0; i < nUnits; i++)
  hiddenStreses [i] = NULL;
 
 
}


MaxwellChainMaterialStatus :: ~MaxwellChainMaterialStatus ()
{
// destructor

  if (hiddenStreses) {
    for (int i=0 ; i< nMaxwellUnits ; i++)
      delete hiddenStreses[i] ;
    delete hiddenStreses ;
  }
}

FloatArray*
MaxwellChainMaterialStatus :: letHiddenStressVectorBe (int i, FloatArray* newVector)
{
/*
 sets i-th  hidden stress vector to newVector
 deletetin previous one, if defined
*/
  if (i > nMaxwellUnits) {
  _error ("letHiddenStressVectorBe: i-th chain exceed specified limit");
  exit(1);
  }

  if (hiddenStreses[i-1]) delete hiddenStreses[i-1];
  hiddenStreses[i-1] = newVector;

  return newVector;
}
  


void 
MaxwellChainMaterialStatus :: printOutputAt (FILE *file, TimeStep* tStep)
{
 // printing of hidden stresses 
 // useful only for debugging
 int i,j; 
 FloatArray helpVec;

 StructuralMaterialStatus:: printOutputAt (file, tStep);

 fprintf(file,"{hidden stresses: ");
 for (i=1; i< nMaxwellUnits; i++) {
  ((StructuralCrossSection*)
   gp->giveCrossSection())->giveFullCharacteristicVector(helpVec, gp, *(hiddenStreses[i]));
  fprintf(file,"{ ");
  for (j=1; j<= helpVec.giveSize(); j++)
   fprintf(file,"%f ",helpVec.at(j));
  fprintf(file,"} ");
  //delete helpVec;
 }
 if (shrinkageStrain.isNotEmpty()) {
  fprintf(file,"shrinkageStrain: {");
  for (j=1; j<= shrinkageStrain.giveSize(); j++)
   fprintf(file,"%f ",shrinkageStrain.at(j));
  fprintf(file,"} ");
 }

 fprintf(file,"}\n");
}


void
MaxwellChainMaterialStatus :: updateYourself(TimeStep* tStep)
{
  StructuralMaterialStatus :: updateYourself (tStep);
}

void
MaxwellChainMaterialStatus :: initTempStatus () 
{
  StructuralMaterialStatus :: initTempStatus ();
}

contextIOResultType
MaxwellChainMaterialStatus :: saveContext (FILE* stream, void *obj)
//
// saves full information stored in this Status
// 
{
 int i;
 contextIOResultType iores;

 if (stream == NULL) _error ("saveContex : can't write into NULL stream");

 if ((iores = StructuralMaterialStatus :: saveContext (stream, obj)) != CIO_OK) THROW_CIOERR(iores);

 // write a raw data
 for (i=0; i<nMaxwellUnits; i++) {
  if (hiddenStreses[i] == NULL) hiddenStreses[i] = new FloatArray(0);
  if ((iores = hiddenStreses[i]->storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 }

 if ((iores = shrinkageStrain.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 
 
 // return result back
 return CIO_OK;
}


contextIOResultType 
MaxwellChainMaterialStatus :: restoreContext (FILE* stream, void *obj)
// 
// restore state variables from stream
//
{
 int i;
 contextIOResultType iores;

 if ((iores = StructuralMaterialStatus :: restoreContext (stream, obj))!=CIO_OK) THROW_CIOERR (iores);
 
 // read raw data 
 for (i=0; i<nMaxwellUnits; i++) {
  if (hiddenStreses[i] == NULL) hiddenStreses[i] = new FloatArray(0);
  if ((iores = hiddenStreses[i]->restoreYourself(stream))!= CIO_OK) THROW_CIOERR(iores);
 }
 
 if ((iores = shrinkageStrain.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);

 // return result back
 return CIO_OK;
} 


/*
CohesiveInterface added by Milan Jirasek on 1 Feb 2010


                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



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

#include "cohint.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {

//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterial
//---------------------------------------------------------------------------------------------------

CohesiveInterfaceMaterial :: CohesiveInterfaceMaterial (int n, Domain *d) : StructuralMaterial (n,d)
//
// constructor
//
{
  // parameters for rate-dependent version(default values lead to rate-independent response)
  damchartime = 0.0; // default value of characteristic time (damage)
  damrateexp = 1.0; // default value of rate exponent (damage)
  plchartime = 0.0; // default value of characteristic time (plasticity)
  plrateexp = 1.0; // default value of rate exponent (plasticity)
}

void
CohesiveInterfaceMaterial :: giveRealStressVector (FloatArray& answer, MatResponseForm form, GaussPoint* gp, 
                                                    const FloatArray& totalStrain, 
                                                    TimeStep* atTime)
//
// returns real stress vector in 3d stress space of receiver according to 
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
  CohesiveInterfaceMaterialStatus *status = (CohesiveInterfaceMaterialStatus*) this -> giveStatus (gp);
  double kappa, omega, tempGam1p, tempGam2p;
 
  // initialize
  this->initGpForNewStep(gp);
  answer.resize(3);

  // set history variables to values at the end of previous step
  kappa = status->giveKappa();
  omega = status->giveDamage();
 
  // compute equivalent strain
  double equivStrain = computeEquivalentStrain (totalStrain);

  // compute damage
  double om = computeDamage(equivStrain);
  if (om > omega)
    omega = om;

  // normal part of stress-strain law (elastic damage in tension, elastic in compression)
  answer.at(1) = kn * totalStrain.at(1);
  if (answer.at(1)>0.)
    answer.at(1) *= (1.0-omega);

  // shear part of stress-strain law (elastoplasticity)
  double gam1p = status->giveGam1p();
  double gam2p = status->giveGam2p();
  double tau1 = ks * (totalStrain.at(2)-gam1p);
  double tau2 = ks * (totalStrain.at(3)-gam2p);
  double tauYield = (1.0-omega) * coh - tanphi * answer.at(1); 
  if (tauYield<0.) 
    tauYield = 0.;
  double tauNorm = sqrt(tau1*tau1+tau2*tau2);
  if (tauNorm>tauYield){ // plastic flow in shear takes place
    double scale;
    if (plchartime>0.0){ // rate-dependent version
      double dt = atTime->giveTimeIncrement();      
      scale = computeViscoplasticScalingFactor(tauNorm,tauYield,dt);
    }
    else // rate-independent version
      scale = tauYield/tauNorm;
    tau1 *= scale;
    tau2 *= scale;
    tempGam1p = totalStrain.at(2) - tau1/ks;
    tempGam2p = totalStrain.at(3) - tau2/ks;
  }
  else{ // no plastic flow
    tempGam1p = gam1p;
    tempGam2p = gam2p;
  }
  answer.at(2) = tau1;
  answer.at(3) = tau2;

  // add damage overstress (rate-dependent model only)
  if (damchartime>0.0){ 
    double dt = atTime->giveTimeIncrement();
    if (dt>0.0)
      answer.at(1) += computeDamageOverstress(totalStrain.at(1),kappa,omega,dt);
  }

  // update gp
  status-> letTempStrainVectorBe (totalStrain);
  status-> letTempStressVectorBe (answer);
  status-> setTempKappa (kappa);
  status-> setTempDamage(omega);
  status-> setTempGam1p(tempGam1p);
  status-> setTempGam2p(tempGam2p);
  return ;
}

double
CohesiveInterfaceMaterial::computeViscoplasticScalingFactor
(double tauTrial, double tauYield, double dt)
{
  if (tauTrial<=tauYield)
    return 1.0;
  double c = coh * pow(plchartime/(ks*dt),plrateexp) * pow(tauTrial-tauYield,plrateexp-1.);
  double beta = solveBeta(c, plrateexp); 
  return 1. - exp(beta) * (1.-tauYield/tauTrial);
}

double
CohesiveInterfaceMaterial::solveBeta(double c, double N)
{
  const int MAXITER = 20;
  const double MAXERROR = 1.e-12;
  double f, answer = 0.0;

  for (int i=0; i<MAXITER; i++){
    double aux = c*exp(N*answer)+exp(answer);
    f = log(aux);
    if (fabs(f)<MAXERROR)
      return answer;
    double df = (c*N*exp(N*answer)+exp(answer))/aux;
    answer -= f/df;
  }
  printf("No convergence in CohesiveInterfaceMaterial::solveBeta, c=%g, answer=%g, f=%g\n",c,answer,f);
  exit(0);
}

double
CohesiveInterfaceMaterial::computeDamageOverstress(double eps, double& damstrain, double omega, double dt)
{
  if (damstrain >= eps*omega){ // unloading ... no viscous stress
    damstrain = eps*omega;
    return 0.0;
  }
  double c = e0*(1.-omega)*pow(damchartime/dt,damrateexp)*pow(eps*omega-damstrain,damrateexp-1.);
  double beta = solveBeta(c,damrateexp);
  double ddamstrain = (eps*omega-damstrain)*exp(beta);
  damstrain += ddamstrain;
  double answer = (eps*omega-damstrain)*kn;
  return answer;
}

double 
CohesiveInterfaceMaterial::computeEquivalentStrain (const FloatArray& strain)
{
  double answer = ksi*ksi*(strain.at(2)*strain.at(2)+strain.at(3)*strain.at(3));
  double epsn = strain.at(1);
  if (epsn>0.)
    answer += epsn*epsn;
  return sqrt(answer);
}

double
CohesiveInterfaceMaterial::computeDamage(double epseq)
{
  if (epseq<e0)
    return 0.;

  // simple exponential law with 2 parameters e0, ef
  double aux = (epseq-e0)/ef;
  return 1.-(e0/epseq)*exp(-aux);
}

void
CohesiveInterfaceMaterial :: giveCharacteristicMatrix (FloatMatrix& answer,
                                                        MatResponseForm form, MatResponseMode rMode,
                                                        GaussPoint* gp, TimeStep* atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
  MaterialMode mMode = gp->giveMaterialMode();
  switch(mMode) {
  case _3dInterface:
    give3dInterfaceMaterialStiffnessMatrix (answer,form,rMode,gp,atTime);
    break;
  default:
    StructuralMaterial::giveCharacteristicMatrix (answer,form,rMode,gp,atTime);
  }
}


void
CohesiveInterfaceMaterial :: give3dMaterialStiffnessMatrix (FloatMatrix& answer, 
                                                             MatResponseForm form,
                                                             MatResponseMode mode,
                                                             GaussPoint* gp,
                                                             TimeStep* atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
  _error ("give3dMaterialStiffnessMatrix: not implemented");
  return ;
}

int
CohesiveInterfaceMaterial :: giveSizeOfReducedStressStrainVector (MaterialMode mode)
  //
  // returns the size of reduced stress-strain vector
  // acording to mode given by gp.
  //
{
  
  switch(mode) {
  case _3dInterface:
    return 3;
  default:
    return StructuralMaterial::giveSizeOfReducedStressStrainVector(mode);
  }
}

int
CohesiveInterfaceMaterial :: giveStressStrainComponentIndOf (MatResponseForm form, MaterialMode mmode, int ind)
//
// this function returns index of reduced(if form == ReducedForm) 
// or Full(if form==FullForm) stressStrain component in Full or reduced 
// stressStrainVector acording to stressStrain mode of given gp.
//
{  
  if (mmode == _3dInterface) return ind;
  else return StructuralMaterial::giveStressStrainComponentIndOf(form, mmode, ind);
}

void
CohesiveInterfaceMaterial :: giveStressStrainMask(IntArray& answer, MatResponseForm form,
                                                   MaterialMode mmode) const
  //
  // this function returns mask of reduced(if form == ReducedForm) 
  // or Full(if form==FullForm) stressStrain vector in full or 
  // reduced StressStrainVector
  // acording to stressStrain mode of given gp.
  // 
  //
  // mask has size of reduced or full StressStrain Vector and  i-th component
  // is index to full or reduced StressStrainVector where corresponding 
  // stressStrain resides. 
  //
  // Reduced form is sub-vector (of stress or strain components),
  // where components corresponding to imposed zero stress (plane stress,...)
  // are not included. On the other hand, if zero strain component is imposed
  // (Plane strain, ..) this condition must be taken into account in geometrical
  // relations, and corresponding component is included in reduced vector.
  //
{
 int i;
  
 if (mmode == _3dInterface) {
   answer.resize (3);
   for (i = 1; i<=3; i++) answer.at(i) = i;
 } else StructuralMaterial::giveStressStrainMask(answer, form, mmode);
 
}

void
CohesiveInterfaceMaterial :: giveReducedCharacteristicVector (FloatArray& answer, GaussPoint* gp, 
                                              const FloatArray& charVector3d)
//
// returns reduced stressVector or strainVector from full 3d vector reduced
// to vector required by gp->giveStressStrainMode()
//
{
 MaterialMode mode = gp-> giveMaterialMode ();
 
 if (mode == _3dInterface) {
   answer = charVector3d;
   return ;
 } else StructuralMaterial::giveReducedCharacteristicVector(answer,gp,charVector3d);
}

void
CohesiveInterfaceMaterial :: giveFullCharacteristicVector (FloatArray& answer, 
                                           GaussPoint* gp, 
                                           const FloatArray& strainVector) 
//
// returns full 3d general strain vector from strainVector in reducedMode
// based on StressStrainMode in gp. Included are strains which 
// perform nonzero work.
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
// you must assigng your stress strain mode to one of the folloving modes (or add new)
// FullForm of MaterialStiffnessMatrix must have the same form.
//
{ 
 MaterialMode mode = gp-> giveMaterialMode ();
 if (mode == _3dInterface) {answer = strainVector; return; }
 else StructuralMaterial::giveFullCharacteristicVector(answer, gp, strainVector);
}
 
void 
CohesiveInterfaceMaterial::give3dInterfaceMaterialStiffnessMatrix (FloatMatrix& answer, MatResponseForm form, MatResponseMode rMode,
                                                                    GaussPoint* gp, TimeStep* atTime)
{
  double om, en;
  CohesiveInterfaceMaterialStatus *status = (CohesiveInterfaceMaterialStatus*) this -> giveStatus (gp);
  
  if ((rMode == ElasticStiffness) || (rMode == SecantStiffness) || (rMode == TangentStiffness)) {
    // assemble elastic stiffness
    answer.resize(3,3);
    answer.zero();
    answer.at(1,1) = kn;
    answer.at(2,2) = answer.at(3,3) = ks;
    
    if (rMode == ElasticStiffness) 
      return;
    if (rMode == SecantStiffness) {
      // secant stiffness
      om = status->giveTempDamage();
      en = status->giveTempStrainVector().at(1);
      // damage in tension only
      if (en >= 0) { answer.times (1.0-om);}
      return;
    } else {
      // tangent stiffness
      //      _error ("give3dInterfaceMaterialStiffnessMatrix: tangent stiffness not implemented");
      // secant stiffness
      om = status->giveTempDamage();
      en = status->giveTempStrainVector().at(1);
      // damage in tension only
      if (en >= 0) { answer.times (1.0-om);}
      return;
    }
   
  }  else _error ("give3dInterfaceMaterialStiffnessMatrix: unknown MatResponseMode");

  return ;
}

int 
CohesiveInterfaceMaterial::giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime)
{
 CohesiveInterfaceMaterialStatus* status = (CohesiveInterfaceMaterialStatus*) this -> giveStatus (aGaussPoint);
 if ((type == IST_DamageTensor)||(type == IST_PrincipalDamageTensor)) {
   answer.resize(1);
   answer.at(1) = status->giveDamage();
   return 1;
 } else if ((type == IST_DamageTensorTemp)||(type == IST_PrincipalDamageTempTensor)) {
   answer.resize(1);
   answer.at(1) = status->giveTempDamage();
   return 1;
 } else if (type == IST_MaxEquivalentStrainLevel) {
   answer.resize(1);
   answer.at(1) = status->giveKappa ();
   return 1;
 } else return StructuralMaterial::giveIPValue (answer, aGaussPoint, type, atTime);
}

InternalStateValueType 
CohesiveInterfaceMaterial::giveIPValueType (InternalStateType type)
{
 if ((type == IST_DamageTensor)||(type == IST_DamageTensorTemp) ||
   (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)) return ISVT_TENSOR_S3;
 else if (type == IST_MaxEquivalentStrainLevel) return ISVT_SCALAR;
 else return StructuralMaterial::giveIPValueType (type);
}

int 
CohesiveInterfaceMaterial::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
{
  if ((type == IST_DamageTensor)||(type ==IST_DamageTensorTemp) ||
      (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)) {
    answer.resize (9);
    answer.at(1) = 1;
    return 1;
  } else if (type == IST_MaxEquivalentStrainLevel) {
    answer.resize (1);
    answer.at(1) = 1;
    return 1;
  } else 
    return StructuralMaterial::giveIntVarCompFullIndx (answer, type, mmode);
}

int
CohesiveInterfaceMaterial::giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint)
{
 if ((type == IST_DamageTensor)||(type ==IST_DamageTensorTemp) ||
   (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)||
   (type == IST_MaxEquivalentStrainLevel)) return 1;
 else return StructuralMaterial::giveIPValueSize (type, aGaussPoint);
}

void
CohesiveInterfaceMaterial :: giveThermalDilatationVector (FloatArray& answer, 
                                                           GaussPoint * gp,  TimeStep* tStep)
   //
   // returns a FloatArray(6) of initial strain vector
   // eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
   // caused by unit temperature in direction of 
   // gp (element) local axes
   //
{
 //FloatArray *result = new FloatArray (6);
 answer.resize (6); answer.zero();
 answer.at(1) = this->tempDillatCoeff;
 answer.at(2) = this->tempDillatCoeff;
 answer.at(3) = this->tempDillatCoeff;
 
 return ;
}

IRResultType
CohesiveInterfaceMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 // elastic parameters
 IR_GIVE_FIELD (ir, kn, IFT_IsoInterfaceDamageMaterial_kn, "kn");
 IR_GIVE_FIELD (ir, ks, IFT_IsoInterfaceDamageMaterial_ks, "ks"); 

 // tensile damage parameters
 IR_GIVE_FIELD (ir, e0, IFT_IsoInterfaceDamageMaterial_ft, "e0");
 IR_GIVE_FIELD (ir, ef, IFT_IsoInterfaceDamageMaterial_gf, "ef"); 
 IR_GIVE_FIELD (ir, ksi, IFT_IsoInterfaceDamageMaterial_gf, "ksi"); 

 // plastic shear parameters
 IR_GIVE_FIELD (ir, coh, IFT_IsoInterfaceDamageMaterial_gf, "coh"); 
 IR_GIVE_FIELD (ir, tanphi, IFT_IsoInterfaceDamageMaterial_gf, "tanphi"); 

 // thermal coefficient
 tempDillatCoeff = 0.0;

 // parameters for the rate-dependent version
 IR_GIVE_OPTIONAL_FIELD (ir, damchartime, IFT_IsoInterfaceDamageMaterial_gf, "tau"); 
 IR_GIVE_OPTIONAL_FIELD (ir, damrateexp, IFT_IsoInterfaceDamageMaterial_gf, "rate"); 
 IR_GIVE_OPTIONAL_FIELD (ir, plchartime, IFT_IsoInterfaceDamageMaterial_gf, "taupl"); 
 IR_GIVE_OPTIONAL_FIELD (ir, plrateexp, IFT_IsoInterfaceDamageMaterial_gf, "ratepl"); 
 
 return StructuralMaterial::initializeFrom (ir);
} 

int
CohesiveInterfaceMaterial::giveInputRecordString(std::string &str, bool keyword)
{
 char buff[1024];

 StructuralMaterial::giveInputRecordString(str, keyword);

 sprintf(buff, " talpha %e kn %e ks %e", this -> tempDillatCoeff,kn,ks);
 str += buff;

 return 1;
}

//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterialStatus
//---------------------------------------------------------------------------------------------------
 
CohesiveInterfaceMaterialStatus::CohesiveInterfaceMaterialStatus (int n, Domain*d, GaussPoint* g) : StructuralMaterialStatus(n,d,g)
{
 kappa = tempKappa = 0.0;
 damage = tempDamage = 0.0;
 gam1p = tempGam1p = 0.0;
 gam2p = tempGam2p = 0.0;
}

CohesiveInterfaceMaterialStatus::~CohesiveInterfaceMaterialStatus () 
{}

void 
CohesiveInterfaceMaterialStatus :: printOutputAt  (FILE *file, TimeStep* tStep)
{
 
 StructuralMaterialStatus :: printOutputAt (file, tStep);
 fprintf (file,"status { ");
 if (this->damage > 0.0) {
  fprintf (file,"kappa %f, damage %f ",this->kappa, this->damage);
 }
 if (gam1p != 0.0 || gam2p != 0.0) {
   fprintf (file,"gam1p %f, gam2p %f, gamp %f ",gam1p, gam2p, sqrt(gam1p*gam1p+gam2p*gam2p));
 }
 fprintf (file,"}\n");
}
 
void 
CohesiveInterfaceMaterialStatus::initTempStatus ()
{
 StructuralMaterialStatus :: initTempStatus();
 this->tempKappa = this->kappa;
 this->tempDamage= this->damage;
 this->tempGam1p = this->gam1p;
 this->tempGam2p = this->gam2p;
}

void 
CohesiveInterfaceMaterialStatus::updateYourself(TimeStep* atTime)
{
 StructuralMaterialStatus::updateYourself(atTime);
 this->kappa = this->tempKappa;
 this->damage= this->tempDamage;
 this->gam1p = this->tempGam1p;
 this->gam2p = this->tempGam2p;
}

contextIOResultType
CohesiveInterfaceMaterialStatus::saveContext (DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // save parent class status
 if ((iores = StructuralMaterialStatus :: saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 // write raw data
 if (!stream->write(&kappa,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->write(&damage,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->write(&gam1p,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->write(&gam2p,1)) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

contextIOResultType
CohesiveInterfaceMaterialStatus::restoreContext(DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // read parent class status
 if ((iores = StructuralMaterialStatus :: restoreContext (stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);

 // read raw data 
 if (!stream->read (&kappa,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->read (&damage,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->read (&gam1p,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->read (&gam2p,1)) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

} // namespace oofem


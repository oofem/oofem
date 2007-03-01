/* $Header: /home/cvs/bp/oofem/sm/src/mazarsmodel.C,v 1.3.4.1 2004/04/05 15:19:47 bp Exp $ */
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

// file: mazarsmodel.C


#include "mazarsmodel.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"

#define _MAZAR_MODEL_ITER_TOL 1.e-15

MazarsMaterial :: MazarsMaterial (int n, Domain *d) : IsotropicDamageMaterial1 (n,d)
//
// constructor
//
{
 // force the loading / unloading controll according to max. reached damage level
 llcriteria = idm_damageLevelCR;
}


MazarsMaterial :: ~MazarsMaterial ()
//
// destructor
//
{}

IRResultType
MazarsMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro
 int ver;

 IsotropicDamageMaterial1::initializeFrom (ir);

 ver         = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, ver, IFT_MazarsMaterial_version, "version"); // Macro
 if (ver == 1) this->modelVersion = maz_modTension;
 else if (ver == 0) this->modelVersion = maz_original;
 else _error ("instanciateFrom: unknown version");
 
 IR_GIVE_FIELD (ir, this->eps_0, IFT_MazarsMaterial_e0, "e0"); // Macro
 IR_GIVE_FIELD (ir, this->Ac, IFT_MazarsMaterial_ac, "ac"); // Macro
 IR_GIVE_FIELD (ir, this->Bc, IFT_MazarsMaterial_bc, "bc"); // this->Bc = (Ac-1.0)/(Ac*eps_0)

 beta = 1.06;
 IR_GIVE_OPTIONAL_FIELD (ir, beta, IFT_MazarsMaterial_beta, "beta"); // Macro

 if (this->modelVersion == maz_original) {
  IR_GIVE_FIELD (ir, this->At, IFT_MazarsMaterial_at, "at"); // Macro
  IR_GIVE_FIELD (ir, this->Bt, IFT_MazarsMaterial_bt, "bt"); // Macro
 } else if (this->modelVersion == maz_modTension) {
  // in case of modified model store ef to At
  IR_GIVE_FIELD (ir, this->At, IFT_MazarsMaterial_ef, "ef"); // Macro
 }
 
 if (ir->hasField(IFT_MazarsMaterial_r, "r")) { // nonlocal version 
  // ask for compulsory "reference length"  
  IR_GIVE_FIELD (ir, this->hReft, IFT_MazarsMaterial_hreft, "hreft"); // Macro
  IR_GIVE_FIELD (ir, this->hRefc, IFT_MazarsMaterial_hrefc, "hrefc"); // Macro
 }
/*
 double gt, gc;

 gt   = this->readDouble(initString,"gt");
 gc   = this->readDouble(initString,"gc");

 // determine hRef (diffrenet for compression and tension?)
 hReft = gt/(eps_0*eps_0*linearElasticMaterial->give('E'));
 hRefc = gc/(eps_0*eps_0*linearElasticMaterial->give('E'));
*/ 
 return IRRT_OK;
}


void 
MazarsMaterial::computeEquivalentStrain (double& kappa, const FloatArray& strain, GaussPoint* gp, TimeStep* atTime)
{
 int i;
 LinearElasticMaterial* lmat = this->giveLinearElasticMaterial();
 double posNorm = 0.0;
 FloatArray principalStrains, strainb;
 StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
 
 if (strain.isEmpty()) {
  kappa= 0.; 
  return;
 }
 crossSection->giveFullCharacteristicVector(strainb, gp, strain);
 // if plane stress mode -> compute strain in z-direction from condition of zero stress in corresponding direction
 if (gp->giveMaterialMode() == _PlaneStress) {
  double nu = lmat->give (NYxz);
  strainb.at(3) = -nu * (strainb.at(1)+strainb.at(2))/(1.-nu);
 } else if (gp->giveMaterialMode() == _1dMat) {
  double nu = lmat->give (NYxz);
  strainb.at(2) = -nu * strainb.at(1);
  strainb.at(3) = -nu * strainb.at(1);
 }

 if (gp->giveMaterialMode() == _1dMat) {
  principalStrains.resize(3);
  for (i=1; i<=3; i++) principalStrains.at(i) = strainb.at(i);
 } else this->computePrincipalValues (principalStrains, strainb, principal_strain);

/*
 // simple check
 double i1 = strainb.at(1)+strainb.at(2)+strainb.at(3) - principalStrains.at(1)-principalStrains.at(2)-principalStrains.at(3);
 double i2 = strainb.at(1)*strainb.at(3)+strainb.at(2)*strainb.at(3)+strainb.at(1)*strainb.at(2) - 
  0.25*(strainb.at(4)*strainb.at(4)+strainb.at(5)*strainb.at(5)+strainb.at(6)*strainb.at(6)) - 
   principalStrains.at(1)*principalStrains.at(3)+principalStrains.at(2)*principalStrains.at(3)+principalStrains.at(1)*principalStrains.at(2);
 if ((fabs(i1) > 1.e-6) || (fabs(i2) > 1.e-6)) {
  printf("v");
 }
 // end simple check
*/
 for (i=1; i<=3; i++) {
  if (principalStrains.at(i) > 0.0) posNorm += principalStrains.at(i)*principalStrains.at(i);
 }

 kappa = sqrt(posNorm);

}


void 
MazarsMaterial::giveNormalElasticStiffnessMatrix (FloatMatrix& answer,
                         MatResponseMode rMode,
                         GaussPoint*gp, TimeStep* atTime)
{
//
// return Elastic Stiffness matrix for normal Stresses
 LinearElasticMaterial *lMat = this->giveLinearElasticMaterial();
 FloatMatrix de;
 int i,j;
 
 answer.resize(3,3);
 lMat -> give3dMaterialStiffnessMatrix(de, FullForm, rMode, gp, atTime);
 // fullAnswer = new FloatMatrix (3,3);
 // copy first 3x3 submatrix to answer
 for (i=1;i<=3; i++)
  for (j=1;j<=3; j++)
   answer.at(i,j)=de.at(i,j);
}



void 
MazarsMaterial::computeDamageParam (double& omega, double kappa, const FloatArray& strain, GaussPoint* gp) 
{
 int i, positive_flag=0, negat_count = 0, nite;
 FloatMatrix de, ce;
 FloatArray strainb, sigp, epsti, epsi, epst, sigppos(3);
 double gt, gc, alpha_t, alpha_c, help, kappaRefT, kappaRefC, gtCurr, gcCurr, hCurrt, hCurrc, R, dgt, dgc, equivStrain;
 LinearElasticMaterial* lmat = this->giveLinearElasticMaterial();
 StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();


 if (kappa > this->eps_0) {
  // omega = 1.0-(this->e0/kappa)*exp(-(kappa-this->e0)/(this->ef-this->e0));
  MazarsMaterialStatus *status = (MazarsMaterialStatus*) this -> giveStatus (gp);
  // do not change passed parameter
  crossSection->giveFullCharacteristicVector(strainb, gp, strain);

  // compute principal strains due to positive stresses
  if (gp->giveMaterialMode() == _PlaneStress) {
   double nu = lmat->give (NYxz);
   strainb.at(3) = -nu * (strainb.at(1)+strainb.at(2))/(1.-nu);
  } else if (gp->giveMaterialMode() == _1dMat) {
   double nu = lmat->give (NYxz);
   strainb.at(2) = -nu * strainb.at(1);
   strainb.at(3) = -nu * strainb.at(1);
  }

  this -> giveNormalElasticStiffnessMatrix (de, TangentStiffness,gp, domain->giveEngngModel()->giveCurrentStep());
  if (!strainb.isEmpty()) this->computePrincipalValues (epsi, strainb, principal_strain);
  else {epsi.resize(3); epsi.zero();}

  // compute actual (local) equivStrain
  equivStrain = 0.0;
  for (i=1; i<=3; i++) {
   if (epsi.at(i) > 0.0) equivStrain += epsi.at(i)*epsi.at(i);
  }
  equivStrain = sqrt(equivStrain);

/*

  if (gp->giveMaterialMode() == _PlaneStress) {
   FloatArray help;
   this->computePrincipalValues (help, strain, principal_strain);
   epsi.resize(3);
   epsi.at(1) = help.at(1); epsi.at(2)=help.at(2);
  }
   
*/ 

  // compute principal stresses
  sigp.beProductOf (de, epsi);
  // take positive part
  for (i=1; i<=3; i++) sigppos.at(i) = macbra(sigp.at(i));
  // compute principal strains due to positive stresses
  ce.beInverseOf (de);
  epsti.beProductOf (ce,sigppos);
 
  positive_flag = negat_count = 0;
  for (i=1; i<=3; i++) {
   if (sigp.at(i) > 1.e-6) {
    positive_flag = 1;
    break;
   } else if (sigp.at(i) < 1.e-6) 
    negat_count++;
  }
  if ((positive_flag == 0)&&(negat_count > 1)) {
   // adjust equivalent strain to improve biaxial compression
   double f = 0.0, g= 0.0 ;
   for (i=1; i<=3; i++) if (sigp.at(i)<0) {f+= sigp.at(i)*sigp.at(i); g+= fabs(sigp.at(i));}
   f = sqrt (f) / g;
   kappa *= f;
  }

  // test adjusted kappa
  if (kappa < this->eps_0) {
   omega=0.0; 
   return;
  }
  
/*
  gt = 1.0-(1.0-this->At)*this->eps_0/kappa - this->At*exp(-this->Bt*(kappa-this->eps_0));
  gc = 1.0-(1.0-this->Ac)*this->eps_0/kappa - this->Ac*exp(-this->Bc*(kappa-this->eps_0));
  if (gc > 1.0) gc = 1.0;
*/

  
  //objectivity 
  nite = 0;
  hCurrt = status->giveLe();
  hCurrc = status->giveLec();
  kappaRefT = kappaRefC = kappa; // this->eps_0;
  
  // tension objectivity
  do {
   nite++;
   if (this->modelVersion == maz_modTension) {
    gt = 1.0-(this->eps_0/kappaRefT)*exp((this->eps_0-kappaRefT)/this->At);
    dgt= this->eps_0/kappaRefT/kappaRefT+this->eps_0/kappaRefT/this->At;
    dgt*=exp((this->eps_0-kappaRefT)/this->At);
   } else { // maz_original
    gt   = 1.0-(1.0-this->At)*this->eps_0/kappaRefT - this->At*exp(-this->Bt*(kappaRefT-this->eps_0));
    dgt  = (1.0-this->At)* this->eps_0/kappaRefT/kappaRefT+this->At*exp(-this->Bt*(kappaRefT-this->eps_0))*this->Bt;
   }
   help = 1+gt*(hReft/hCurrt - 1.0);
   R    = kappaRefT*help - kappa;
   if (fabs(R) <= _MAZAR_MODEL_ITER_TOL) break; 
   if (nite > 400) 
    _error ("computeDamageParam: tension objectivity iteration internal error - no convergence");
   help+= dgt*kappaRefT*(hReft/hCurrt - 1.0);
   kappaRefT += -R/help;
  } while (1);

  nite = 0;
  // compression objectivity
  do {
   nite++;
   gc   = 1.0-(1.0-this->Ac)*this->eps_0/kappaRefC - this->Ac*exp(-this->Bc*(kappaRefC-this->eps_0));
   help = 1+gc*(hRefc/hCurrc - 1.0);
   R    = kappaRefC*help - kappa;
   if (fabs(R) <= _MAZAR_MODEL_ITER_TOL) break; 
   if (nite > 400) 
    _error ("computeDamageParam: comression objectivity iteration internal error - no convergence");
   dgc  = (1.0-this->Ac)* this->eps_0/kappaRefC/kappaRefC+this->Ac*exp(-this->Bc*(kappaRefC-this->eps_0))*this->Bc;
   help+= dgc * kappaRefC*(hRefc/hCurrc - 1.0);
   kappaRefC += -R/help;
  } while (1);


  gtCurr = gt*(hReft*kappaRefT)/(hCurrt*kappa);
  gcCurr = gc*(hRefc*kappaRefC)/(hCurrc*kappa);

  if ((gt < 0.) || (gt > 1.0)) 
   _error ("computeDamageParam: gt out of range ");
  //if (gcCurr > 1.0) 
  // _warning ("computeDamageParam: gc out of range ");
  // if (gcCurr > 1.0) gcCurr = 1.0;


  help = 0.0;
  for (i=1;i<=3;i++) if (epsi.at(i)>0.) help+= epsti.at(i)*epsi.at(i);
  help /= equivStrain*equivStrain;
  if (help>1.0) help = 1.0;
  alpha_t = __OOFEM_POW(help, this->beta);
  alpha_c = __OOFEM_POW((1.-help), this->beta);

  if ((alpha_t < 0.) || (alpha_t > 1.0)) 
   _error ("computeDamageParam: apha_t out of range");
  if ((alpha_c < 0.) || (alpha_c > 1.0)) 
   _error ("computeDamageParam: apha_t out of range");

  omega = alpha_t*gtCurr + alpha_c*gcCurr;
  if (omega > 1.0) omega = 1.0;
//  if ((omega <0.) || (omega > 1.0))
//   _error ("computeDamageParam: omega out of range");

 } else 

  omega = 0.0;
}


void 
MazarsMaterial::initDamaged (double kappa, FloatArray& totalStrainVector, GaussPoint* gp) 
{
 int i, indmin = 1, indmax = 1;
 double le;
 FloatArray principalStrains, crackPlaneNormal(3), fullstrain;
 FloatMatrix principalDir (3,3);
 MazarsMaterialStatus *status = (MazarsMaterialStatus*) this -> giveStatus (gp);
 //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();

 // crossSection->giveFullCharacteristicVector(fullstrain, gp, totalStrainVector);

 if ((kappa > this->e0) && (status->giveDamage() == 0.)) {
   //printf (":");


  if (gp->giveMaterialMode() == _1dMat) {
   crackPlaneNormal.zero(); crackPlaneNormal.at(1) = 1.0;
   le = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);
   status->setLe (le);
   status->setLec (le);
   return;
  } 

  
  if (gp->giveMaterialMode() == _PlaneStress) principalDir.resize(2,2);
  this->computePrincipalValDir (principalStrains, principalDir, totalStrainVector, principal_strain);
  if (gp->giveMaterialMode() == _PlaneStress) {
   if (principalStrains.at(1) > principalStrains.at(2)) {indmax = 1; indmin = 2;}
   else {indmax = 2; indmin = 1;}
  } else {
   // find index of max and min positive principal strains
   for (i=2; i<=3; i++)  {
    if (principalStrains.at(i) > principalStrains.at(indmax)) indmax = i;
    if (principalStrains.at(i) < principalStrains.at(indmin)) indmin = i;
   }   
  }

  for (i=1;i<=principalStrains.giveSize();i++) crackPlaneNormal.at(i) = principalDir.at(i,indmax);
  le = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);
  // remember le in cooresponding status for tension
  status->setLe (le);
  
  for (i=1;i<=principalStrains.giveSize();i++) crackPlaneNormal.at(i) = principalDir.at(i,indmin);
  le = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);
  // remember le in cooresponding status for compression
  status->setLec (le);

  // printf ("les: %e %e\n", status->giveLe(), status->giveLec());
 }

 //status->setLe(hReft); status->setLec(hRefc);

}








MazarsMaterialStatus :: MazarsMaterialStatus (int n, Domain*d, GaussPoint *g) 
: IsotropicDamageMaterial1Status(n,d,g)
{
 lec = 0.0;
}

contextIOResultType
MazarsMaterialStatus :: saveContext (FILE* stream, void *obj)
//
// saves full information stored in this Status
// no temp variables stored
//
{
 contextIOResultType iores;
 // save parent class status
 if ((iores = IsotropicDamageMaterial1Status :: saveContext (stream, obj))!= CIO_OK) THROW_CIOERR(iores);

 // write a raw data
 if (fwrite(&lec,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

contextIOResultType
MazarsMaterialStatus :: restoreContext(FILE* stream, void *obj)
//
// restores full information stored in stream to this Status
//
{
 contextIOResultType iores;
 // read parent class status
 if ((iores = IsotropicDamageMaterial1Status :: restoreContext (stream,obj)) != CIO_OK) THROW_CIOERR(iores);
 // read raw data 
 if (fread (&lec,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}


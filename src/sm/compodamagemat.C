/* v 1.8 2009/04/27 vs  */
/*initTempStatus

                   *****    *****   ******  ******  ***   ***
                 **   **  **   **  **      **      ** *** **
                **   **  **   **  ****    ****    **  *  **
               **   **  **   **  **      **      **     **
              **   **  **   **  **      **      **     **
              *****    *****   **      ******  **     **


               OOFEM : Object Oriented Finite Element Code

                 Copyright (C) 2009   Vit Smilauer



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



#include "compodamagemat.h"
#include "material.h"
#include "linearelasticmaterial.h"
#include "ortholinearelasticmaterial.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "cltypes.h"
#include "mathfem.h"
#include "contextioerr.h"

CompoDamageMat :: CompoDamageMat (int n, Domain *d) : StructuralMaterial (n,d)
{
 /// Constructor
}

CompoDamageMat :: ~CompoDamageMat()
{
 /// destructor
}

CompoDamageMat :: ~CompoDamageMat ();



IRResultType CompoDamageMat :: initializeFrom (InputRecord* ir)
{
 int i;
 double value;
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 this -> Material::initializeFrom(ir);
 //define transversely othotropic material stiffness parameters
 IR_GIVE_FIELD (ir, value, IFT_CompoDamageMat_ex, "exx"); // Macro
 propertyDictionary -> add(Ex,value);
 IR_GIVE_FIELD (ir, value, IFT_CompoDamageMat_ez, "eyyezz");
 propertyDictionary -> add(Ey,value);
 propertyDictionary -> add(Ez,value);
 IR_GIVE_FIELD (ir, value, IFT_CompoDamageMat_nyxy, "nyxynuxz");
 propertyDictionary -> add(NYxy,value);
 propertyDictionary -> add(NYxz,value);
 IR_GIVE_FIELD (ir, value, IFT_CompoDamageMat_nyyz, "nyyz");
 propertyDictionary -> add(NYyz,value);
 propertyDictionary -> add(NYzy,value);
 IR_GIVE_FIELD (ir, value, IFT_CompoDamageMat_Gxy, "gxyxz");
 propertyDictionary -> add(Gxy,value);
 propertyDictionary -> add(Gxz,value);

 //calulate remaining components
 propertyDictionary -> add(Gyz,this->give(Ey,NULL)/(1.+this->give(NYxy,NULL)));
 propertyDictionary -> add(NYyx,this->give(Ey,NULL)*this->give(NYxy,NULL)/this->give(Ex,NULL));
 //propertyDictionary -> add(Gzy,this->give(Gyz,NULL));
 propertyDictionary -> add(NYzx,this->give(NYyx,NULL));

 IR_GIVE_FIELD (ir, this->inputTension, IFT_CompoDamageMat_components, "tension_f0_gf");

 //this->inputTension.printYourself();
 IR_GIVE_FIELD (ir, this->inputCompression, IFT_CompoDamageMat_components, "compres_f0_gf");

 if(this->inputTension.giveSize()!=12)
  _error ("instanciateFrom: need 12 components for tension in pairs f0 Gf for all 6 directions");

 if(this->inputCompression.giveSize()!=12)
  _error ("instanciateFrom: need 12 components for compression in pairs f0 Gf for all 6 directions");

 for(i=1;i<=12;i+=2){
  if(this->inputTension.at(i) < 0.0)
   _error ("instanciateFrom: negative f0 detected for tension");
  if(this->inputCompression.at(i) > 0.0)
   _error ("instanciateFrom: positive f0 detected for compression");
 }
//
 //OOFEM_LOG_INFO("READ \n");
 return IRRT_OK;
}

//used in debugging only ?
int CompoDamageMat :: giveInputRecordString(std::string &str, bool keyword)
{

}

//called at the beginning of each time increment (not iteration), no influence of parameter
void  CompoDamageMat :: give3dMaterialStiffnessMatrix (FloatMatrix& answer, MatResponseForm form, MatResponseMode mode, GaussPoint* gp, TimeStep* atTime)
{
 FloatMatrix *rotationMatrix;
 StructuralElement* element = (StructuralElement*) gp-> giveElement();

 //already with reduced components
 this->giveUnrotated3dMaterialStiffnessMatrix (answer, gp);
 rotationMatrix = this->giveMatStiffRotationMatrix(gp);
 (*rotationMatrix).printYourself();
 answer.rotatedWith (*rotationMatrix);
 //answer.printYourself();
 delete rotationMatrix;
}

//called in each iteration, support for 3D and 1D material mode
void CompoDamageMat :: giveRealStressVector (FloatArray& answer,  MatResponseForm form, GaussPoint* gp, const FloatArray& totalStrain, TimeStep* atTime)
{
 int i, i_max, s, tensCompr;
 double delta, sigma, charLen, tmp;
 CompoDamageMatStatus *st = (CompoDamageMatStatus*) this -> giveStatus (gp);
 Element* element = gp->giveElement();
 FloatArray strainVectorL(6), stressVectorL(6), tempStressVectorL(6),  reducedTotalStrainVector(6), tempKappaL(12), ans, equilStressVectorL(6), equilStrainVectorL(6);
 FloatArray *stressLimit;
 FloatMatrix de, elementCs;
 MaterialMode mMode = gp->giveMaterialMode();

 //substract strain independent part - temperature, creep ..., in global c.s.
 this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
 //reducedTotalStrainVector.printYourself();

 switch(mMode) {
  case _3dMat://applies only for 3D
  {
  element -> giveMatLocalCS(elementCs);//if mlcs undefined, returns unix matrix

  //first run (can be called somewhere better ?)
  if(st->elemCharLength.at(1) == 0.)
  this->giveCharLength(st, gp, elementCs);

  //transform strain to local c.s.
  this -> transformStrainVectorTo(strainVectorL,elementCs, reducedTotalStrainVector,0);
  //strainVectorL.printYourself();

  //damage criteria based on stress, assuming same damage parameter for tension/compression
  //determine unequilibrated stress vector
  this->giveUnrotated3dMaterialStiffnessMatrix (de, gp);
  tempStressVectorL.beProductOf (de, strainVectorL);
  i_max = 6;
  break;
  }

  case _1dMat:
  {//applies only for 1D, strain vectors are already in local c.s.
  if(st->elemCharLength.at(1) == 0.){
   st->elemCharLength.at(1) = gp->giveElement()->giveCharacteristicLenght(gp, NULL);
  }

  strainVectorL.at(1) = reducedTotalStrainVector.at(1);
  tempStressVectorL.zero();
  tempStressVectorL.at(1) = this->give(Ex,NULL) * strainVectorL.at(1);
  i_max = 1;
  break;
  }
  default: _error("Compodamagemat: unsupported material mode, accept only 3D and 1D mode");
 }

 //proceed 6 components for 3D or 1 component for 1D, damage evolution is based on the evolution of strains
 //xx, yy, zz, yz, zx, xy
 for(i=1;i<=i_max;i++){

  if(tempStressVectorL.at(i)>=0.){//unequilibrated stress, tension
   stressLimit = &inputTension;//contains pairs (stress - fracture energy)
   s=0;
  }
  else{//compression
   stressLimit = &inputCompression;
   s=6;
  }
  if(fabs(tempStressVectorL.at(i))>fabs((*stressLimit).at(2*i-1)) && st->strainAtMaxStress.at(i+s)==0.){//damage starts now, can be replaced for more advanced initiation criteria, e.g. Hill's maximum combination of stresses
   //equilibrated strain and stress from the last time step, transform to local c.s.
   switch(mMode) {
    case _3dMat:{
    ans=st->giveStrainVector();
    this -> transformStrainVectorTo(equilStrainVectorL,elementCs, ans,0);
    ans=st->giveStressVector();
    this -> transformStressVectorTo(equilStressVectorL,elementCs, ans,0);
    break;
    }
    case _1dMat:{
     equilStrainVectorL = st->giveStrainVector();
     equilStressVectorL = st->giveStressVector();
     break;
    }
   }

   //subdivide last increment, interpolate, delta in range <0;1>
   delta = ((*stressLimit).at(2*i-1) - equilStressVectorL.at(i))/(tempStressVectorL.at(i)-equilStressVectorL.at(i));
   delta = min(delta, 1.); delta = max(delta, 0.);//stabilize transition from tensile to compression (denom -> 0)

   st->strainAtMaxStress.at(i+s)=equilStrainVectorL.at(i) + delta * (strainVectorL.at(i)-equilStrainVectorL.at(i));

   st->initDamageStress.at(i+s) = equilStressVectorL.at(i) + delta * (tempStressVectorL.at(i)-equilStressVectorL.at(i));

   //determine characteristic length for six stresses/strains
   switch (i){
    case 1:
    case 2:
    case 3:
     charLen=st->elemCharLength.at(i); break;
    case 4:charLen=(st->elemCharLength.at(2)+st->elemCharLength.at(3))/2.; break;//average two directions
    case 5:charLen=(st->elemCharLength.at(3)+st->elemCharLength.at(1))/2.; break;//average two directions
    case 6:charLen=(st->elemCharLength.at(1)+st->elemCharLength.at(2))/2.; break;//average two directions
   }

   st->maxStrainAtZeroStress.at(i+s) = st->strainAtMaxStress.at(i+s) + 2*fabs((*stressLimit).at(2*i))/charLen/st->initDamageStress.at(i+s);//Gf scaled for characteristic size
  }

  if(st->strainAtMaxStress.at(i+s) != 0. && fabs(strainVectorL.at(i))>fabs(st->kappa.at(i+s))){//damage started and grows
   //desired stress
   sigma = st->initDamageStress.at(i+s)*(st->maxStrainAtZeroStress.at(i+s) - strainVectorL.at(i))/(st->maxStrainAtZeroStress.at(i+s) - st->strainAtMaxStress.at(i+s));

   //check that sigma remains in tension/compression area
   if(s==0)//tension
    sigma = max(sigma,0.000001);
   else
    sigma = min(sigma,-0.000001);

   switch (i){//need to substract contributions from strains
    case 1: tmp = 1.-sigma/(this->give(Ex,NULL)*strainVectorL.at(i) + this->give(NYxy,NULL)*tempStressVectorL.at(2) + this->give(NYxz,NULL)*tempStressVectorL.at(3)); break;
    case 2: tmp = 1.-sigma/(this->give(Ey,NULL)*strainVectorL.at(i) + this->give(NYyx,NULL)*tempStressVectorL.at(1) + this->give(NYyz,NULL)*tempStressVectorL.at(3)); break;
    case 3: tmp = 1.-sigma/(this->give(Ez,NULL)*strainVectorL.at(i) + this->give(NYzx,NULL)*tempStressVectorL.at(1) + this->give(NYzy,NULL)*tempStressVectorL.at(2)); break;
    case 4: tmp = 1.-sigma/this->give(Gyz,NULL)/strainVectorL.at(i); break;
    case 5: tmp = 1.-sigma/this->give(Gxz,NULL)/strainVectorL.at(i); break;
    case 6: tmp = 1.-sigma/this->give(Gxy,NULL)/strainVectorL.at(i); break;
   }

   st->tempOmega.at(i) = max(tmp,st->omega.at(i));//damage can only grow, interval <0;1>
   st->tempOmega.at(i) = min(st->tempOmega.at(i),0.999999);
   st->tempOmega.at(i) = max(st->tempOmega.at(i),0.0);
   st->tempKappa.at(i+s) = strainVectorL.at(i);
  }
 }

 switch(mMode) {
  case _3dMat:{
   //already with reduced stiffness components in local c.s.
   this->giveUnrotated3dMaterialStiffnessMatrix (de, gp);
   //de.printYourself();
   //in local c.s.
   stressVectorL.beProductOf (de, strainVectorL);
   //stressVectorL.printYourself();
   //transform local c.s to global c.s.
   this->transformStressVectorTo(answer,elementCs,stressVectorL,1);
   break;
  }
  case _1dMat:{
   answer.resize(1);
   answer.at(1) = (1-st->tempOmega.at(1)) * this->give(Ex,NULL) * strainVectorL.at(1);//tempStress
  break;
  }
 }

 st->letTempStressVectorBe (answer);//needed in global c.s for 3D

 //not changed inside this function body
 st->letTempStrainVectorBe (totalStrain);//needed in global c.s for 3D
}

//used for output in *.hom a *.out
int CompoDamageMat :: giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime){

 CompoDamageMatStatus* status = (CompoDamageMatStatus*) this -> giveStatus (aGaussPoint);
 if (type == IST_DamageTensor){
  answer.resize(6);
  answer = status->omega;
 }
 else
  StructuralMaterial::giveIPValue (answer, aGaussPoint, type, atTime);
 return 1;
}

InternalStateValueType CompoDamageMat :: giveIPValueType (InternalStateType type)
{
 if (type == IST_DamageTensor)
  return ISVT_TENSOR_S3;
 else
  return StructuralMaterial::giveIPValueType (type);
}


int CompoDamageMat::giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint)
{
 if (type == IST_DamageTensor)
  return 1;
 else return StructuralMaterial::giveIPValueSize (type, aGaussPoint);
}

int
  CompoDamageMat::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
{
 if (type == IST_DamageTensor){
  answer.resize (9);
  answer.at(1) = 1;
  return 1;
 }
 else
 return StructuralMaterial::giveIntVarCompFullIndx (answer, type, mmode);
}


void CompoDamageMat :: giveUnrotated3dMaterialStiffnessMatrix (FloatMatrix& answer, GaussPoint* gp){
 double i;
 double denom;
 double ex, ey, ez, nxy, nxz, nyz, gyz, gzx, gxy;
 double a,b,c,d,e,f;
 FloatArray tempOmega;

 answer.resize (6,6);
 answer.zero();

 CompoDamageMatStatus *st = (CompoDamageMatStatus*) this -> giveStatus (gp);

 ex = this->give(Ex,NULL);
 ey = this->give(Ey,NULL);
 ez = this->give(Ez,NULL);
 nxy = this->give(NYxy,NULL);
 nxz = nxy;
 nyz = this->give(NYyz,NULL);
 gyz = this->give(Gyz,NULL);
 gzx = this->give(Gxz,NULL);
 gxy = this->give(Gxy,NULL);

 //xx, yy, zz, yz, zx, xy
 //assebmle stiffness matrix for transversely orthotropic material with reduced moduli, derived from compliance matrix with only reduced diagonal terms. Procedure can be used for fully orthotropic stiffness matrix as well
 a=1.-st->tempOmega.at(1);
 b=1.-st->tempOmega.at(2);
 c=1.-st->tempOmega.at(3);
 d=1.-st->tempOmega.at(4);
 e=1.-st->tempOmega.at(5);
 f=1.-st->tempOmega.at(6);

 denom = -ey*ex + ex*nyz*nyz*b*c*ez + nxy*nxy*a*b*ey*ey + 2*nxy*a*b*ey*nxz*nyz*c*ez + nxz*nxz*a*ey*c*ez;

 answer.at(1,1) = (-ey + nyz*nyz*b*c*ez)*a*ex*ex/denom;
 answer.at(1,2) = -(nxy*ey + nxz*nyz*c*ez)*ex*ey*a*b/denom;
 answer.at(1,3) = -(nxy*nyz*b + nxz)*ey*ex*a*c*ez/denom;
 answer.at(2,2) = (-ex + nxz*nxz*a*c*ez)*b*ey*ey/denom;
 answer.at(2,3) = -(nyz*ex + nxz*nxy*a*ey)*ey*b*c*ez/denom;
 answer.at(3,3) = (-ex + nxy*nxy*a*b*ey)*ey*c*ez/denom;
 answer.at(4,4) = gyz;
 answer.at(5,5) = gzx;
 answer.at(6,6) = gxy;
 answer.symmetrized ();
 //answer.printYourself();
}

//returns material rotation stiffness matrix [6x6]
FloatMatrix* CompoDamageMat :: giveMatStiffRotationMatrix(GaussPoint *gp){
 FloatMatrix t(3,3), answer, elementCs;
 int elementCsFlag;
 StructuralElement* element = (StructuralElement*) gp-> giveElement();

 element -> giveMatLocalCS(t);//if mlcs undefined, returns unix matrix
 t.printYourself();
 
 //rotate from unrotated (base) c.s. to local material c.s.
 this -> giveStrainVectorTranformationMtrx (answer, t);
 return answer.GiveCopy();
}

//determine characteristic fracture area for three orthogonal cracks, based on the size of element (crack band model). Since the orientation of cracks is aligned with the orientation of material, determination is based only on the geometry (not on the direction of principal stress etc.). Assumption that fracture localizes into all integration points on element. Material orientation in global c.s. is passed. Called in the first run
void
  CompoDamageMat::giveCharLength (CompoDamageMatStatus *status, GaussPoint* gp, FloatMatrix& elementCs){

 int i, j;
 FloatArray crackPlaneNormal(3);

 //elementCs.printYourself();

 //normal to x,y,z is the same as in elementCs

 for (i=1;i<=3;i++){
  for (j=1;j<=3;j++)
   crackPlaneNormal.at(j) = elementCs.at(j,i);

  //already corrected for the number of integration points on element
  status->elemCharLength.at(i) = gp->giveElement()->giveCharacteristicLenght (gp, crackPlaneNormal);
 }
}



// constructor
CompoDamageMatStatus :: CompoDamageMatStatus (int n, Domain*d, GaussPoint* g) : StructuralMaterialStatus (n,d,g)
{
 //largest strain ever reached [6 tension, 6 compression]
 this->kappa.resize(12); this->kappa.zero();
 this->tempKappa.resize(12); this->tempKappa.zero();

 //array of damage parameters [6] for both tension and compression
 this->omega.resize(6); this->omega.zero();
 this->tempOmega.resize(6); this->tempOmega.zero();

 this->initDamageStress.resize(12);this->initDamageStress.zero();
 this->maxStrainAtZeroStress.resize(12); this->maxStrainAtZeroStress.zero();
 this->strainAtMaxStress.resize(12); this->strainAtMaxStress.zero();

 this->elemCharLength.resize(3); this->elemCharLength.zero();
}

// destructor
CompoDamageMatStatus :: ~CompoDamageMatStatus ()
{
}


 /// Prints the receiver state to stream
void CompoDamageMatStatus :: printOutputAt (FILE *file, TimeStep* tStep)
{
 int i;
 StructuralMaterialStatus :: printOutputAt (file, tStep);
 fprintf (file,"status {");
 if(!this->omega.containsOnlyZeroes()){
  fprintf (file, " omega ");
  for(i=1;i<=6;i++)
   fprintf (file,"%.4f ", this->omega.at(i));
 }

 fprintf (file," kappa ");
 for(i=1;i<=12;i++)
  fprintf (file,"%.4f ",this->kappa.at(i));
 fprintf (file,"}\n");
}


//initializes temp variables according to variables form previous equilibrium state, resets tempStressVector, tempStrainVector
//function called at the beginning of each time increment (not iteration)
void CompoDamageMatStatus :: initTempStatus ()
{
}

// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reached equilibrium
void CompoDamageMatStatus :: updateYourself(TimeStep* atTime)
{
 //here stressVector = tempStressVector; strainVector = tempStrainVector;
 StructuralMaterialStatus :: updateYourself (atTime);//MaterialStatus::updateYourself, i.e. stressVector = tempStressVector; strainVector = tempStrainVector;
 this->kappa = this->tempKappa;
 this->omega= this->tempOmega;
}


contextIOResultType CompoDamageMatStatus :: saveContext (DataStream* stream, ContextMode mode, void *obj = NULL)
{
 contextIOResultType iores;
 // save parent class status
 if ((iores = StructuralMaterialStatus :: saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}

contextIOResultType CompoDamageMatStatus :: restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL){
 contextIOResultType iores;
 // read parent class status ?
 if ((iores = StructuralMaterialStatus :: restoreContext (stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);

 return CIO_OK;
}


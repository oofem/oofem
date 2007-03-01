/* $Header: /home/cvs/bp/oofem/sm/src/ortholinearelasticmaterial.C,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   file ortholinearelasticmaterial.cc

#include "linearelasticmaterial.h" 
#include "ortholinearelasticmaterial.h"
#include "material.h"
#include "structuralms.h"
#include "domain.h"
#include "flotmtrx.h"
#include "timestep.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <math.h>
#endif

#define ZERO_LENGTH 1.e-6

IRResultType
OrthotropicLinearElasticMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 double value ;
 double nyzx, nyzy, nyyx;
 int j,size;
 FloatArray triplets;
 
 
 this -> Material::initializeFrom(ir);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_ex, "ex"); // Macro
 propertyDictionary -> add(Ex,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_ey, "ey"); // Macro
 propertyDictionary -> add(Ey,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_ez, "ez"); // Macro
 propertyDictionary -> add(Ez,value);
 
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_nyyz, "nyyz"); // Macro
 propertyDictionary -> add(NYyz,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_nyxz, "nyxz"); // Macro
 propertyDictionary -> add(NYxz,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_nyxy, "nyxy"); // Macro
 propertyDictionary -> add(NYxy,value);
 
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_gyz, "gyz"); // Macro
 propertyDictionary -> add(Gyz,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_gxz, "gxz"); // Macro
 propertyDictionary -> add(Gxz,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_gxy, "gxy"); // Macro
 propertyDictionary -> add(Gxy,value);
 
 
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_talphax, "talphax"); // Macro
 propertyDictionary -> add(tAlphax,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_talphay, "talphay"); // Macro
 propertyDictionary -> add(tAlphay,value);
 
 IR_GIVE_FIELD (ir, value, IFT_OrthotropicLinearElasticMaterial_talphaz, "talphaz"); // Macro
 propertyDictionary -> add(tAlphaz,value);
 
 // check for suspicious parameters
 // ask for dependent parameters (symmetry conditions) and check if reasonable 
 nyzx = this->give (NYzx);
 nyzy = this->give (NYzy);
 nyyx = this->give (NYyx);
 if ((nyzx < 0.) || (nyzx > 0.5) || (nyzy < 0.) || (nyzy > 0.5) || (nyyx < 0.) || (nyyx > 0.5))
  _warning2 ("instanciateFrom: suspicious parameters", 1);
 
 // Read local coordinate system of principal axes of ortotrophy 
 // in localCoordinateSystem the unity vectors are stored
 // COLUMWISE (this is exception, but allows faster numerical
 // implementation)
 
 // try to read lcs section
 triplets.resize(0);
 IR_GIVE_OPTIONAL_FIELD (ir, triplets, IFT_OrthotropicLinearElasticMaterial_lcs, "lcs"); // Macro

 size=triplets.giveSize();
 if (!((size==0)||(size ==6))) 
   _warning2("instanciateFrom: Warning: lcs in material %d is not properly defined, will be assumed as global",
             this->giveNumber());
 
 if (size==6) {
  
  cs_type = localCS;
  double n1=0.0, n2=0.0 ;
  
  localCoordinateSystem = new FloatMatrix(3,3);
  for (j=1; j<= 3; j++) {
   localCoordinateSystem->at(j,1) = triplets.at(j);
   n1 += triplets.at(j)*triplets.at(j);
   localCoordinateSystem->at(j,2) = triplets.at(j+3);
   n2 += triplets.at(j+3)*triplets.at(j+3);
  }
  n1 = sqrt(n1); n2 = sqrt(n2);
  for (j=1; j<= 3; j++) { // normalize e1' e2'
   localCoordinateSystem->at(j,1) /= n1;
   localCoordinateSystem->at(j,2) /= n2;
  }
  // vector e3' computed from vector product of e1', e2'
  localCoordinateSystem->at(1,3) = 
   (localCoordinateSystem->at(2,1)*localCoordinateSystem->at(3,2) -
    localCoordinateSystem->at(3,1)*localCoordinateSystem->at(2,2)) ;
  localCoordinateSystem->at(2,3) = 
   (localCoordinateSystem->at(3,1)*localCoordinateSystem->at(1,2) -
    localCoordinateSystem->at(1,1)*localCoordinateSystem->at(3,2)) ;
  localCoordinateSystem->at(3,3) = 
   (localCoordinateSystem->at(1,1)*localCoordinateSystem->at(2,2) -
    localCoordinateSystem->at(2,1)*localCoordinateSystem->at(1,2)) ;
 }
 
 // try to read ElementCS section
 if (cs_type == unknownCS) {
  triplets.resize(0);
  IR_GIVE_OPTIONAL_FIELD (ir, triplets, IFT_OrthotropicLinearElasticMaterial_scs, "scs"); // cs for shells. 
  // first three numbers are direction of normal n - see orthoelasticmaterial.h for description
  // shellCS  - coordinate system of principal axes is specified in shell  coordinate system
  //            this is defined as follows: principal z-axis is perpendicular to mid-section 
  //            x-axis is perpendicular to z-axis and normal to user specified vector n.
  //            (so x-axis is parallel to plane, with n beeing normal to this plane).
  //            y-axis is then perpendicular both to x and z axes.
  //            WARNING: this definition of cs is valid only for plates and shells
  //            when vector n is paralel to z-axis an error occurs and program is terminated.
  //
  size=triplets.giveSize();
  if (!((size==0)||(size ==3))) 
    _warning2("instanciateFrom: scs in material %d is not properly defined, will be assumed as global",
              this->giveNumber());
  
  if (size==3) {
   
   cs_type = shellCS;
   triplets.normalize();
   helpPlaneNormal = new FloatArray(3);
   for (j=1; j<4;j++) helpPlaneNormal->at(j) = triplets.at(j);
   //
   // store normal defining help plane into row matrix
   // localCoordinateSystemmust be computed on demand from specific element
   //
  }
 } //
 if (cs_type == unknownCS) {  
  //
  // if no cs defined assume global one
  //
  cs_type = localCS;
  localCoordinateSystem = new FloatMatrix(3,3);
  for (j=1; j< 4; j++) localCoordinateSystem->at(j,j) = 1.0;
 }
 return IRRT_OK;
}


double
OrthotropicLinearElasticMaterial :: give (int aProperty)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
 if (aProperty == NYzx) return this->give(NYxz)*this->give(Ex) / this->give(Ez);
 if (aProperty == NYzy) return this->give(NYyz)*this->give(Ey) / this->give(Ez);
 if (aProperty == NYyx) return this->give(NYxy)*this->give(Ex) / this->give(Ey);
 
 return this -> Material :: give(aProperty);
}


void
OrthotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                                  MatResponseForm form,
                                  MatResponseMode mode,
                                  GaussPoint* gp, 
                                  TimeStep* atTime)
//
// forceElasticResponse ignored - always elastic
//
{
 FloatMatrix *rotationMatrix ;
 double eksi,nxz,nyz,nxy,nzx,nzy,nyx;
 int i,j ;
 
 nxz = this->give(NYxz);
 nyz = this->give(NYyz);
 nxy = this->give(NYxy);
 nzx = this->give(NYzx);
 nzy = this->give(NYzy);
 nyx = this->give(NYyx);
 
 eksi = 1.-(nxy*nyx+nyz*nzy+nzx*nxz)-(nxy*nyz*nzx+nyx*nzy*nxz);

 //constitutiveMatrix = new FloatMatrix(6,6) ;
 answer.resize (6,6);
 answer.zero();
 
 answer.at(1,1) =  this->give(Ex)*(1.-nzy*nyz)/eksi ;
 answer.at(1,2) =  this->give(Ey)*(nyx+nzx*nyz)/eksi;
 answer.at(1,3) =  this->give(Ez)*(nzx+nzy*nyx)/eksi; 
 answer.at(2,2) =  this->give(Ey)*(1.-nzx*nxz)/eksi; 
 answer.at(2,3) =  this->give(Ez)*(nzy+nxy*nzx)/eksi; 
 answer.at(3,3) =  this->give(Ez)*(1.-nxy*nyx)/eksi;
 
 // define the lover triangle
 for (i=1;i<4;i++)
  for (j=1; j< i; j++)
   answer.at(i,j)=answer.at(j,i);


 answer.at(4,4) =  this->give(Gyz);
 answer.at(5,5) =  this->give(Gxz);
 answer.at(6,6) =  this->give(Gxy);

 rotationMatrix = this -> GiveRotationMatrix (gp);
 answer.rotatedWith (*rotationMatrix);
 delete rotationMatrix;

 return ;
}



FloatMatrix*
OrthotropicLinearElasticMaterial :: GiveTensorRotationMatrix (GaussPoint* gp)
//
// returns [3,3] rotation matrix from local principal axes of material
// to local axes used at gp (element) level
// 
{
 int elementCsFlag;
 FloatMatrix elementCs, *t=NULL;
 StructuralElement* element = (StructuralElement*) gp-> giveElement();
 
 elementCsFlag = element -> giveLocalCoordinateSystem(elementCs);
 // 
 // in localCoordinateSystem the directional cosines are stored columwise(exception)
 // in elementCs rowwise.
 // 
 if (this->cs_type == localCS) {
  // 
  // in localCoordinateSystem are stored directional cosines
  //
  if (elementCsFlag) {
   t = elementCs.Times (this->localCoordinateSystem);
   //delete elementCs;
  } else {
   t = this->localCoordinateSystem->GiveCopy();
  }
 } else if (this->cs_type == shellCS) {
  FloatArray *elementNormal, *helpx, *helpy;
  localCoordinateSystem = new FloatMatrix(3,3);
  int i;
  
  elementNormal = element -> ComputeMidPlaneNormal(gp); 
  elementNormal -> normalize();
  helpx = this->helpPlaneNormal->VectorProduct (elementNormal);
  // test if localCoordinateSystem is uniquely 
  // defined by elementNormal and helpPlaneNormal
  if (dotProduct(helpx->givePointer(),helpx->givePointer(),3) < ZERO_LENGTH)
   _error ("GiveTensorRotationMatrix: element normal parallel to plane normal encountered");
  helpy = elementNormal->VectorProduct(helpx);
  for (i=1; i< 4; i++) {
   localCoordinateSystem->at(i,1) = helpx->at(i);
   localCoordinateSystem->at(i,2) = helpy->at(i);
   localCoordinateSystem->at(i,3) = elementNormal->at(i);
  }
  delete elementNormal;
  delete helpx;
  delete helpy;

  // 
  // possible rotation about local z-axix should be considered in future
  //
  /*
   //
   // GiveZRotationMtrx assembles rotMtrx from cs rotated from curent about rotAngle
   // to current cs
   //
   zRotMtrx = GiveZRotationMtrx (rotAngle); // rotAngle supplied by user
   rotatedLocalCoordinateSystem = localCoordinateSystem->Times (zRotMtrx);
   delete localCoordinateSystem;
   localCoordinateSystem = rotatedLocalCoordinateSystem;
   */
  if (elementCsFlag) {
   t = elementCs.Times (this->localCoordinateSystem);
   //delete elementCs;
  } else {
   t = this->localCoordinateSystem->GiveCopy();
  }
  delete localCoordinateSystem; localCoordinateSystem= NULL;

 } else {
  _error("GiveTensorRotationMatrix - internal error no cs defined");
 }
 // 
 // t at (i,j) contains cosine of angle between elementAxis(i) and localMaterialAxis(j).
 //
 return t;
}


FloatMatrix*
OrthotropicLinearElasticMaterial :: GiveRotationMatrix (GaussPoint *gp)
//
// returns [6,6] rotation matrix from local principal axes of material 
// to local axes used at the gp (element) level
// at element level is implemented next transformation to global cs.
//
//
{
 FloatMatrix *t, answer;
 t = this->GiveTensorRotationMatrix (gp);
 // 
 // t at (i,j) contains cosine of angle between elementAxis(i) and localMaterialAxis(j).
 //
 

 this -> giveStrainVectorTranformationMtrx (answer, *t);
 delete t;
 return answer.GiveCopy();
}


/*
FloatMatrix*
IsotropicLinearElasticMaterial :: GivePlaneStressStiffMtrx (MatResponseForm form,
                    GaussPoint* gp,
                    FloatArray* strainIncrement,
                    TimeStep* atTime)
// 
// strainIncrement may be used for loading/unloading criteria
// in this simple material it is not necessary
//
{
 FloatMatrix* constitutiveMatrix ;
 double e,nu,ee,shear ;
 
 e     = this -> give('E') ;
 nu    = this -> give('n') ;
 ee    = e/(1.-nu*nu);
 shear = e/(2.0*(1.+nu));

 if (form == FullForm) {
  constitutiveMatrix = new FloatMatrix(6,6);
  
  constitutiveMatrix->at(1,1) = ee;
  constitutiveMatrix->at(1,2) = nu*ee;
  constitutiveMatrix->at(2,1) = nu*ee;
  constitutiveMatrix->at(2,2) = ee;
  constitutiveMatrix->at(6,6) = shear;
 } else {
  constitutiveMatrix = new FloatMatrix(3,3);
  
  constitutiveMatrix->at(1,1) = ee;
  constitutiveMatrix->at(1,2) = nu*ee;
  constitutiveMatrix->at(2,1) = nu*ee;
  constitutiveMatrix->at(2,2) = ee;
  constitutiveMatrix->at(3,3) = shear;
 }  
   return constitutiveMatrix ;
}



FloatMatrix*
IsotropicLinearElasticMaterial :: GivePlaneStrainStiffMtrx  (MatResponseForm form,
                     GaussPoint* gp,
                     FloatArray* strainIncrement,
                     TimeStep* atTime)
// 
// strainIncrement may be used for loading/unloading criteria
// in this simple material it is not necessary
//
{
 FloatMatrix* constitutiveMatrix;
 double e,nu,ee,shear;
 
 e     = this -> give('E');
 nu    = this -> give('n');
 ee    = e/(1.0+nu)/(1.-2.0*nu);
 shear = e/(2.0*(1.+nu)) ;

 if (form == FullForm) {
  constitutiveMatrix = new FloatMatrix(6,6);
 
  constitutiveMatrix->at(1,1) = ee*(1.0-nu);
  constitutiveMatrix->at(1,2) = nu*ee;
  constitutiveMatrix->at(2,1) = nu*ee;
  constitutiveMatrix->at(2,2) = ee*(1.0-nu);
  constitutiveMatrix->at(3,1) = nu*ee;
  constitutiveMatrix->at(3,2) = nu*ee;
  constitutiveMatrix->at(1,3) = nu*ee;
  constitutiveMatrix->at(2,3) = nu*ee;
  constitutiveMatrix->at(6,6) = shear;
 } else {
  constitutiveMatrix = new FloatMatrix(3,3);
 
  constitutiveMatrix->at(1,1) = ee*(1.0-nu);
  constitutiveMatrix->at(1,2) = nu*ee;
  constitutiveMatrix->at(2,1) = nu*ee;
  constitutiveMatrix->at(2,2) = ee*(1.0-nu);
  constitutiveMatrix->at(3,3) = shear;
 }  
 
   return constitutiveMatrix ;
}


FloatMatrix*
IsotropicLinearElasticMaterial :: Give1dStressStiffMtrx (MatResponseForm form,
                   GaussPoint* gp,
                   FloatArray* strainIncrement,
                   TimeStep* atTime)
// 
// strainIncrement may be used for loading/unloading criteria
// in this simple material it is not necessary
//
{
 FloatMatrix* constitutiveMatrix ;
 double e,nu,ee,shear ;
 
 e     = this -> give('E') ;
 
 if (form == FullForm) {
  constitutiveMatrix = new FloatMatrix(6,6);
  
  constitutiveMatrix->at(1,1) = e;
 } else {
  constitutiveMatrix = new FloatMatrix(1,1);
  
  constitutiveMatrix->at(1,1) = e;
 }  
   return constitutiveMatrix ;
}
*/

void
OrthotropicLinearElasticMaterial :: giveThermalDilatationVector (FloatArray& answer, 
                                 GaussPoint* gp, TimeStep* tStep)
   //
   // returns a FloatArray(3) of coefficients of thermal dillatation in direction
   // of each (local) axisgiven by element lcs.
   //
{
 FloatMatrix *transf;
 FloatArray  help(6);
 help.at(1) = this->give(tAlpha);
 help.at(2) = this->give(tAlpha);
 help.at(3) = this->give(tAlpha);
 
 transf = this->GiveRotationMatrix (gp);
 answer.beProductOf (*transf, help);
 //delete help;
 delete transf;
 
 return ;
}
 

MaterialStatus*
OrthotropicLinearElasticMaterial:: CreateStatus (GaussPoint* gp) const
/* 
 creates new  material status  corresponding to this class
*/
{
 return new StructuralMaterialStatus (1,this->giveDomain(),gp);
}

